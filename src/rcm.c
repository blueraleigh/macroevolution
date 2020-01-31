#include "rcm.h"


static struct rcm_state *state_add(struct rcm_statelist *sl)
{
    sl->len += 1;

    if (sl->len > sl->r)
        return sl->tail;

    struct rcm_state *state = NULL;
    size_t nbytes = 3 * sl->nnode * sizeof(double);

    if (!sl->fl)
    {
        state = malloc(sizeof(struct rcm_state));
        if (!state)
            error("unable to allocate state");
        state->mem = malloc(nbytes);
        if (!state->mem)
        {
            free(state);
            error("unable to allocate state memory");
        }
        state->index = sl->len;
        state->loglk = 0;
        state->ntip = 0;
        state->dclk = (double *)(state->mem);
        state->sclk = state->dclk + sl->nnode;
        state->uclk = state->sclk + sl->nnode;
        state->next = NULL;
        state->prev = NULL;
        state->stat = sl->stat_alloc(sl->p, sl->hyperparam);

        memset(state->mem, 0, nbytes);

        state->uclk[sl->root] = 1 / (double)(sl->r);
    }
    else
    {
        state = sl->fl;
        sl->fl = state->prev;
        state->prev = NULL;
    }

    if (!sl->head)
    {
        state->prev = NULL;
        sl->head = state;
        return sl->head;
    }
    else
    {
        if (sl->tail)
        {
            sl->tail->next = state;
            state->prev = sl->tail;

            // tail state always interpreted as non-analog state
            memcpy(state->mem, sl->tail->mem, nbytes);
        }
        else
        {
            sl->head->next = state;
            state->prev = sl->head;
        }
        sl->tail = state;
    }
    // was previously this (prior to 2020-01-17):
    //     return (sl->len < sl->r) ? sl->tail->prev : sl->tail;
    // but this caused a bug when model was initialized with
    // no empty states. pretty sure this change will not now
    // cause bugs when model is initialized with empty states
    return (sl->len <= sl->r) ? sl->tail->prev : sl->tail;
}


static void state_remove(struct rcm_state *state, struct rcm_statelist *sl)
{
    struct rcm_state *prev = state->prev;
    struct rcm_state *next = state->next;

    if (sl->len <= sl->r)
    {
        if (state == sl->head)
        {
            sl->head = next;
            next->prev = NULL;

        }
        else
        {
            prev->next = next;
            next->prev = prev;
        }

        state->ntip = 0;
        state->next = NULL;
        state->prev = NULL;

        state->prev = sl->fl;
        sl->fl = state;
    }
    else
    {
        // removing a state when maximum number of states allowed by model are
        // already represented among the tips requires special handling. in this
        // case, the removed state is repositioned at the tail of the list where
        // it can be reused as the non-analog state.
        if (state == sl->head)
        {
            sl->head = next;
            next->prev = NULL;

            // move state to tail
            sl->tail->next = state;
            state->next = NULL;
            state->prev = sl->tail;
            sl->tail = state;

        }
        else if (state != sl->tail)
        {
            prev->next = next;
            next->prev = prev;

            // move state to tail
            sl->tail->next = state;
            state->next = NULL;
            state->prev = sl->tail;
            sl->tail = state;
        }
    }

    sl->len -= 1;
}


static void state_free(struct rcm_state *state)
{
    free(state->mem);
    free(state);
}


static void statelist_free(struct rcm_statelist *sl)
{
    struct rcm_state *a;
    struct rcm_state *b;

    a = b = sl->head;

    if (a)
    {
        do {
            b = a->next;
            sl->stat_free(a->stat);
            state_free(a);
            a = b;
        } while (b);
    }

    a = b = sl->fl;

    if (a)
    {
        do {
            b = a->prev;
            state_free(a);
            a = b;
        } while (b);
    }
}


void rcm_clear(struct rcm *model)
{
    if (model)
    {
        free(model->stateid);
        free(model->stateid_r);
        free(model->needsupdate_up);
        free(model->needsupdate_dn);
        free(model->lzd);
        free(model->lzu);
        statelist_free(&model->sl);
    }
}


static void rcm_branch_downpass(struct node *node, struct rcm *model)
{
    struct rcm_statelist *sl = &model->sl;
    struct rcm_state *i;
    struct rcm_state *j;
    struct rcm_state *non = sl->tail;
    int r = sl->len - 1;  // number of analog states
    int r_max = sl->r;

    if (r == r_max)
    {
        // there are zero non-analog states
        r = r_max - 1;
    }

    double D = exp(-r_max * model->rate * node->brlen);
    double pij = (1 - D) / (double)r_max;
    double pii = pij + D;

    SCLK(non, node) = pii * DCLK(non, node) + (r_max - r - 1) * pij * DCLK(non, node);

    for (i = sl->head; i != non; i = i->next)
    {
        SCLK(i, node) = pii * DCLK(i, node) + (r_max - r) * pij * DCLK(non, node);
        SCLK(non, node) += pij * DCLK(i, node);
    }

    for (i = sl->head; i != non->prev; i = i->next)
    {
        for (j = i->next; j != non; j = j->next)
        {
            SCLK(i, node) += pij * DCLK(j, node);
            SCLK(j, node) += pij * DCLK(i, node);
        }
    }
}


static void rcm_node_downpass(struct node *node, struct rcm *model)
{
    int scale;
    struct rcm_statelist *sl = &model->sl;
    struct rcm_state *i;

    struct node *lfdesc;
    struct node *rtdesc;

    lfdesc = node->lfdesc;
    rtdesc = lfdesc->next;

    model->lzd[node->index] = model->lzd[lfdesc->index]
        + model->lzd[rtdesc->index];

    for (i = sl->head; i != NULL; i = i->next)
        DCLK(i, node) = SCLK(i, lfdesc) * SCLK(i, rtdesc);

    scale = 1;
    for (i = sl->head; scale && (i != NULL); i = i->next)
        scale = (DCLK(i, node) < minlikelihood) && (DCLK(i, node) > minusminlikelihood);

    if (scale)
    {
        for (i = sl->head; i != NULL; i = i->next)
            DCLK(i, node) *= twotothe256;

        model->lzd[node->index] += 1;
    }
}


static void rcm_downpass(struct rcm *model)
{
    struct node *node;

    phy_traverse_prepare(model->phy, model->phy->root, INTERNAL_NODES_ONLY,
        POSTORDER);

    while ((node = phy_traverse_step(model->phy)) != 0)
    {
        rcm_branch_downpass(node->lfdesc, model);
        rcm_branch_downpass(node->lfdesc->next, model);
        rcm_node_downpass(node, model);
    }
}


static void rcm_node_uppass(struct node *node, struct rcm *model)
{
    int scale;
    struct rcm_statelist *sl = &model->sl;
    struct rcm_state *i;
    struct rcm_state *j;
    struct rcm_state *non = sl->tail;
    int r = sl->len - 1;  // number of analog states
    int r_max = sl->r;

    struct node *anc;
    struct node *sib;

    anc = node->anc;

    sib = anc->lfdesc;

    if (sib == node)
        sib = sib->next;

    if (r == r_max)
    {
        // there are zero non-analog states
        r = r_max - 1;
    }

    double D = exp(-r_max * model->rate * node->brlen);
    double pij = (1 - D) / (double)r_max;
    double pii = pij + D;

    UCLK(non, node) = UCLK(non, anc) * SCLK(non, sib) * pii
        + (r_max - r - 1) * UCLK(non, anc) * SCLK(non, sib) * pij;

    for (j = sl->head; j != non; j = j->next)
    {
        UCLK(j, node) = (r_max - r) * UCLK(non, anc) * SCLK(non, sib) * pij
            + UCLK(j, anc) * SCLK(j, sib) * pii;
        UCLK(non, node) += UCLK(j, anc) * SCLK(j, sib) * pij;
    }

    for (j = sl->head; j != non->prev; j = j->next)
    {
        for (i = j->next; i != non; i = i->next)
        {
            UCLK(j, node) += UCLK(i, anc) * SCLK(i, sib) * pij;
            UCLK(i, node) += UCLK(j, anc) * SCLK(j, sib) * pij;
        }
    }

    model->lzu[node->index] = model->lzu[anc->index] + model->lzd[sib->index];

    scale = 1;
    for (j = sl->head; scale && (j != NULL); j = j->next)
        scale = (UCLK(j, node) < minlikelihood) && (UCLK(j, node) > minusminlikelihood);

    if (scale)
    {
        for (j = sl->head; j != NULL; j = j->next)
            UCLK(j, node) *= twotothe256;
        model->lzu[node->index] += 1;
    }
}


static void rcm_uppass(struct rcm *model)
{
    struct node *node;

    phy_traverse_prepare(model->phy, model->phy->root, ALL_NODES, PREORDER);

    // step past root
    phy_traverse_step(model->phy);

    while ((node = phy_traverse_step(model->phy)) != 0)
    {
        rcm_node_uppass(node, model);
        model->needsupdate_up[node->index] = 0;
        model->needsupdate_dn[node->index] = 0;
    }
}


static void rcm_update_node(struct node *node, struct rcm *model)
{
    // these changes ensure that the local configuration of uppass likelihoods
    // is correct for performing updates to terminal nodes. the global
    // configuration of uppass likelihoods is NOT correct, however, and will
    // need to be reset with a full uppass from the root.

    if (model->needsupdate_up[node->anc->index])
        rcm_node_uppass(node, model);

    if (model->needsupdate_up[node->index])
    {
        struct node *sib;

        // model_update_nodes ensures node will never equal the root node so that
        // anc is always non-NULL
        sib = node->anc->lfdesc;

        // we should always be on a right descendant at this stage
        assert(sib != node);

        rcm_branch_downpass(sib, model);
        rcm_node_uppass(node, model);
    }

    if (!node->ndesc)
    {
        struct rcm_statelist *sl = &model->sl;
        struct rcm_state *z_new;
        struct rcm_state *z_old;
        struct rcm_state *i;

        double w;
        double g;
        double wmax = -HUGE_VAL;

        z_new = z_old = model->stateid[node->index];

        model->pop(node->index, model->data, z_old->stat);

        for (i = sl->head; i != NULL; i = i->next)
        {
            // we want to sample a new state from the distribution
            //
            // exp(w) / sum(exp(w))
            //
            // we do this in one pass by computing
            //
            // argmax w + g,
            //
            // where g is a Gumbel(0, 1) random variate
            //
            // for why this works see:
            // https://timvieira.github.io/blog/post/2014/07/31/gumbel-max-trick/
            // or
            // https://arxiv.org/pdf/1706.04161.pdf

            w = model->push(node->index, model->data, i->stat);
            w += log(UCLK(i, node)) + model->lzu[node->index]
                        * log(minlikelihood);

            // Gumbel(0,1) random variate
            g = -log(-log(unif_rand()));

            if ((g + w) > wmax)
            {
                wmax = g + w;
                z_new = i;
            }

            model->pop(node->index, model->data, i->stat);
        }

        if (z_new->ntip == 0 && z_old->ntip == 1)
        {
            // update transfers node from singleton state to another singleton
            // state so there is no point in changing the label. in fact, the
            // likelihoods for z_new and z_old are identical in this case.
            z_new = z_old;
        }

        model->push0(node->index, model->data, z_new->stat);

        if (z_old != z_new)
        {
            if (!z_new->ntip)
                state_add(sl);

            z_old->ntip -= 1;
            z_new->ntip += 1;

            model->stateid[node->index] = z_new;
            model->stateid_r[node->index] = z_new->index;
            DCLK(z_old, node) = 0;
            DCLK(z_new, node) = 1;

            model->dataloglk -= z_old->loglk;
            model->dataloglk -= z_new->loglk;
            z_old->loglk = sl->stat_loglk(z_old->stat, model->data);
            z_new->loglk = sl->stat_loglk(z_new->stat, model->data);
            model->dataloglk += z_old->loglk;
            model->dataloglk += z_new->loglk;

            if (!z_old->ntip)
                state_remove(z_old, sl);

            if (node->next)
                model->needsupdate_up[node->next->index] = 1;
            else
                rcm_branch_downpass(node, model); // force update now as we won't get another opportunity

            model->needsupdate_dn[node->anc->index] = 1;
        }

        while (node->anc != 0 && node->next == 0)
        {
            node = node->anc;
            if (model->needsupdate_dn[node->index])
            {
                rcm_node_downpass(node, model);
                if (node->next)
                    model->needsupdate_up[node->next->index] = 1;
                else
                    rcm_branch_downpass(node, model);
                if (node->anc)
                    model->needsupdate_dn[node->anc->index] = 1;
            }
        }
    }
}


static double rcm_root_loglk(struct rcm *model)
{
    struct rcm_statelist *sl = &model->sl;
    struct rcm_state *i;

    double lk = 0;

    for (i = sl->head; i != NULL; i = i->next)
        lk += DCLK(i, model->phy->root) / (double)(model->r);

    return log(lk) + model->lzd[model->phy->root->index] * log(minlikelihood);
}


double rcm_treeloglk(struct rcm *model)
{
    rcm_downpass(model);
    rcm_uppass(model);
    return rcm_root_loglk(model);
}


double rcm_dataloglk(struct rcm *model)
{
    struct rcm_statelist *sl = &model->sl;
    struct rcm_state *i;

    double loglk = 0;

    for (i = sl->head; i != NULL; i = i->next)
    {
        i->loglk = sl->stat_loglk(i->stat, model->data);
        loglk += i->loglk;
    }

    return loglk;
}


static void rcm_update_nodes(struct rcm *model)
{
    struct node *node;
    struct node *root = model->phy->root;

    phy_traverse_prepare(model->phy, root, ALL_NODES,
        PREORDER);

    phy_traverse_step(model->phy);
    while ((node = phy_traverse_step(model->phy)) != 0)
        rcm_update_node(node, model);

    rcm_uppass(model);

    model->treeloglk = rcm_root_loglk(model);

    #ifndef NDEBUG
        double delta = model->treeloglk - rcm_treeloglk(model);
        assert(delta < 1e-12 && delta > -1e-12);
    #endif
}


// step 1. define vertical level for slice
static double slice(struct rcm *model)
{
    return (model->treeloglk + log(1 / model->rate_max)) - rexp(1);
}


// step 2. find slice interval containing x
static void slice_locate(double x, double w, double z, double *L, double *R, struct rcm *model)
{
    double u;
    double logprior = -log(model->rate_max);

    u = unif_rand();
    *L = x - w*u;
    *R = *L + w;

    model->rate = *L;
    while (*L > 0 && (rcm_treeloglk(model) + logprior) > z) {
        *L -= w;
        model->rate = *L;
    }

    model->rate = *R;
    while (*R < model->rate_max && (rcm_treeloglk(model) + logprior) > z) {
        *R += w;
        model->rate = *R;
    }

    *L = (*L < 0) ? 1e-8 : *L;
    *R = (*R > model->rate_max) ? model->rate_max : *R;
}


// step 3. sample new point x1 from slice
static void slice_sample(double x, double z, double L, double R, struct rcm *model)
{
    double u;
    double x1;
    double Lbar = L;
    double Rbar = R;
    double logprior = -log(model->rate_max);

    do {
        u = unif_rand();
        x1 = Lbar + u * (Rbar - Lbar);

        if (x1 < x)
            Lbar = x1;
        else
            Rbar = x1;

        model->rate = x1;
        model->treeloglk = rcm_treeloglk(model);

    } while ((model->treeloglk + logprior) < z);
}


static void rcm_update_rate(double tune_rate, struct rcm *model)
{
    double x;
    double z;
    double L;
    double R;

    x = model->rate;
    z = slice(model);

    slice_locate(x, tune_rate * model->rate_max, z, &L, &R, model);
    slice_sample(x, z, L, R, model);
}


void rcm_init_states(int *stateid, struct rcm *model)
{
    int i;
    int j;
    int w;
    int inlist;
    int n = 1;
    struct rcm_state *state;
    SEXP root = PROTECT(list1(ScalarInteger(stateid[0])));
    SEXP head;
    SEXP tail;

    for (i = 1; i < model->phy->ntip; ++i)
    {
        head = root;
        tail = CDR(head);

        inlist = 0;
        for (j = 0; j < n && !inlist; ++j)
        {
            inlist = (INTEGER(CAR(head))[0] == stateid[i]);
            if (tail != R_NilValue)
            {
                head = tail;
                tail = CDR(head);
            }
        }

        if (!inlist)
        {
            SETCDR(head, list1(ScalarInteger(stateid[i])));
            ++n;
        }
    }

    // initialize the statelist
    state = state_add(&model->sl);
    state_add(&model->sl);

    head = root;
    for (i = 0; i < model->phy->ntip; ++i)
    {
        if (INTEGER(CAR(head))[0] == stateid[i])
        {
            DCLK(state, phy_getnode_with_index(model->phy, i)) = 1;
            state->ntip += 1;
            model->stateid[i] = state;
            model->stateid_r[i] = state->index;
            model->push0(i, model->data, state->stat);
        }
    }

    head = CDR(head);
    for (j = 1; j < n; ++j)
    {
        state = state_add(&model->sl);
        for (i = 0; i < model->phy->ntip; ++i)
        {
            if (INTEGER(CAR(head))[0] == stateid[i])
            {
                DCLK(state, phy_getnode_with_index(model->phy, i)) = 1;
                state->ntip += 1;
                model->stateid[i] = state;
                model->stateid_r[i] = state->index;
                model->push0(i, model->data, state->stat);
            }
        }

        head = CDR(head);
    }

    UNPROTECT(1);
}


struct rcm *rcm_init_start(
    int p,
    int r,
    double rate,
    struct phy *phy)
{
    double f;
    double tbar;
    struct rcm *model;
    struct node *node;

    if (!phy)
        error("passing a NULL pointer as phy object");

    model = malloc(sizeof(struct rcm));

    if (!model)
        return NULL;

    model->p = p;
    model->r = r;
    model->stateid         = malloc(phy->ntip * sizeof(struct rcm_state *));
    model->stateid_r       = malloc(phy->ntip * sizeof(int));
    model->needsupdate_up  = calloc(phy->nnode, sizeof(int));
    model->needsupdate_dn  = calloc(phy->nnode, sizeof(int));

    model->lzd          = calloc(phy->nnode, sizeof(int));
    model->lzu          = calloc(phy->nnode, sizeof(int));

    model->sl.len       = 0;
    model->sl.r         = r;
    model->sl.p         = p;
    model->sl.root      = phy->root->index;
    model->sl.ntip      = phy->ntip;
    model->sl.nnode     = phy->nnode;
    model->sl.head      = NULL;
    model->sl.tail      = NULL;
    model->sl.fl        = NULL;
    model->sl.hyperparam = NULL;

    model->phy = phy;

    model->rate = rate;

    if (!model->stateid || !model->needsupdate_up || !model->needsupdate_dn
        || !model->stateid_r)
    {
        free(model->stateid);
        free(model->stateid_r);
        free(model->needsupdate_up);
        free(model->needsupdate_dn);
        free(model);
        return NULL;
    }

    tbar = 0;
    f = (double)(((phy->nnode - 1) / r) + 1) / (phy->nnode - 1);

    phy_traverse_prepare(phy, phy->root, ALL_NODES, PREORDER);
    phy_traverse_step(phy);

    while ((node = phy_traverse_step(phy)) != 0)
        tbar += node->brlen / (phy->nnode - 1);

    model->rate_max = -log((f*r - 1) / (r - 1)) / (r * tbar);

    return model;
}


static void rcm_mcmc_sample(struct rcm_mcmc *mcmc)
{
    Rconnection con;
    con = R_GetConnection(mcmc->outputConn);
    con->write(&mcmc->model->dataloglk, sizeof(double), 1, con);
    con->write(&mcmc->model->treeloglk, sizeof(double), 1, con);
    con->write(&mcmc->model->rate, sizeof(double), 1, con);
    con->write(mcmc->model->stateid_r, sizeof(int), mcmc->model->phy->ntip, con);
}


static void rcm_mcmc_run(struct rcm_mcmc *mcmc)
{
    int i;
    double u;

    rcm_mcmc_sample(mcmc);

    for (i = 0; i < mcmc->niter; ++i)
    {
        u = unif_rand() * mcmc->update_weight;

        u -= mcmc->update_node;

        if (u <= 0)
            rcm_update_nodes(mcmc->model);
        else
            rcm_update_rate(mcmc->tune_rate, mcmc->model);

        if (((i+1) % mcmc->sample_freq) == 0)
            rcm_mcmc_sample(mcmc);

        if (((i+1) % 1048576) == 0)
            R_CheckUserInterrupt();
    }
}


SEXP rcm_model_mcmc_run(
    SEXP model,
    SEXP niter,
    SEXP thin,
    SEXP updatenode,
    SEXP updaterate,
    SEXP tunerate,
    SEXP outputConn)
{
    struct rcm_mcmc mcmc;

    if (R_CONNECTIONS_VERSION != 1)
        error("mkdmm package code was written for R_CONNECTIONS_VERSION=1 "
            "but this version of R is running R_CONNECTIONS_VERSION=%d", R_CONNECTIONS_VERSION);

    mcmc.niter          = INTEGER(niter)[0];
    mcmc.sample_freq    = INTEGER(thin)[0];
    mcmc.update_node    = REAL(updatenode)[0];
    mcmc.update_rate    = REAL(updaterate)[0];
    mcmc.tune_rate      = REAL(tunerate)[0];
    mcmc.model          = (struct rcm *)R_ExternalPtrAddr(model);
    mcmc.outputConn     = outputConn;

    mcmc.update_weight = mcmc.update_node + mcmc.update_rate;

    GetRNGstate();

    rcm_mcmc_run(&mcmc);

    PutRNGstate();

    return R_NilValue;
}


/******************************************************************************
** Post MCMC run functions
*******************************************************************************/

/* Compute posterior pairwise coincidence probabilities */
SEXP rcm_posterior_coincidence(SEXP stateid)
{
    int i;
    int j;
    int k;
    int item;
    double *rho;
    int *partition = INTEGER(stateid);
    int nrow = INTEGER(getAttrib(stateid, R_DimSymbol))[0];
    int ntip = INTEGER(getAttrib(stateid, R_DimSymbol))[1];
    int nelem = 0.5 * ntip * (ntip-1);

    SEXP RHO = PROTECT(allocVector(REALSXP, nelem));
    rho = REAL(RHO);
    memset(rho, 0, nelem * sizeof(double));

    for (i = 0; i < nrow; ++i)
    {
        item = 0;
        for (j = 0; j < (ntip-1); ++j)
        {
            for (k = (j+1); k < ntip; ++k)
            {
                if (partition[i + j * nrow] == partition[i + k * nrow])
                    rho[item] += 1 / (double)nrow;
                ++item;
            }
        }
    }

    UNPROTECT(1);
    return RHO;
}


/*
** a := penalty for assigning two tips to different states when they are
**      really in the same state
** b := penalty for assigning two tips to the same state when they are
**      really in different states
** Argument K to the function is then b / (a + b)
** Argument rho is the coincidence matrix computed above
** Argument stateid is the matrix of tip partitions
*/
SEXP rcm_expected_loss(SEXP K, SEXP stateid, SEXP rho)
{
    int i;
    int j;
    int k;
    int item;
    int nrow = INTEGER(getAttrib(stateid, R_DimSymbol))[0];
    int ntip = INTEGER(getAttrib(stateid, R_DimSymbol))[1];
    int *a = INTEGER(stateid);
    double f = REAL(K)[0];

    SEXP LOSS = PROTECT(allocVector(REALSXP, nrow));
    double *loss = REAL(LOSS);
    double *r = REAL(rho);

    for (i = 0; i < nrow; ++i)
    {
        item = 0;
        loss[i] = 0;
        for (j = 0; j < (ntip-1); ++j)
        {
            for (k = (j+1); k < ntip; ++k)
            {
                if (a[i + j*nrow] == a[i + k*nrow])
                    loss[i] += r[item] - f;
                ++item;
            }
        }
    }

    UNPROTECT(1);
    return LOSS;
}


/* Marginal ancestral state reconstruction */

static double asr_compute(
    struct rcm_state *j,
    struct node *node,
    struct rcm *model)
{
    struct rcm_statelist *sl = &model->sl;
    struct rcm_state *i;
    struct rcm_state *non = sl->tail;
    int r = sl->len - 1;  // number of analog states
    int r_max = sl->r;

    double lk = 0;

    struct node *anc;
    struct node *sib;

    double D = exp(-r_max * model->rate * node->brlen);
    double pij = (1 - D) / (double)r_max;
    double pii = pij + D;

    anc = node->anc;

    if (r == r_max)
    {
        // there are zero non-analog states
        r = r_max - 1;
    }

    if (anc)
    {
        sib = anc->lfdesc;

        if (sib == node)
            sib = sib->next;

        if (j != non)
        {
            lk = (r_max - r) * UCLK(non, anc) * SCLK(non, sib) * pij * DCLK(j, node);
        }
        else
        {
            lk = UCLK(non, anc) * SCLK(non, sib) * pii * DCLK(j, node)
                + (r_max - r - 1) * UCLK(non, anc) * SCLK(non, sib) * pij * DCLK(j, node);
        }

        for (i = sl->head; i != non; i = i->next)
        {
            if (i != j)
                lk += UCLK(i, anc) * SCLK(i, sib) * pij * DCLK(j, node);
            else
                lk += UCLK(i, anc) * SCLK(i, sib) * pii * DCLK(j, node);
        }

    }
    else
    {
        return DCLK(j, node);
    }

    return lk;
}


static void asr_normalize(int r, double *out)
{
    int i;
    double norm = 0;

    for (i = 0; i < r; ++i)
        norm += out[i];

    for (i = 0; i < r; ++i)
        out[i] /= norm;
}


SEXP rcm_marginal_asr(SEXP model)
{
    int i;
    int r;
    int ntip;
    struct rcm_statelist *sl;
    struct rcm_state *j;
    double *asr;

    struct rcm *m;
    struct node *node;

    SEXP ASR;

    m = (struct rcm *)R_ExternalPtrAddr(model);
    ntip = m->phy->ntip;
    sl = &m->sl;
    r = sl->len - 1;

    if (r < sl->r)
        ++r;

    ASR = PROTECT(allocMatrix(REALSXP, r, ntip - 1));
    asr = REAL(ASR);

    phy_traverse_prepare(m->phy, m->phy->root, INTERNAL_NODES_ONLY, POSTORDER);

    while ((node = phy_traverse_step(m->phy)) != 0)
    {
        for (j = sl->head, i = 0; j != NULL; j = j->next, ++i)
            asr[i + r * (node->index - ntip)] = asr_compute(j, node, m);

        asr_normalize(r, asr + r * (node->index - ntip));
    }

    UNPROTECT(1);
    return ASR;
}


/* Simulate internal node states using stochastic character mapping */

SEXP rcm_stochastic_map(SEXP nmaps, SEXP model)
{
    int i;
    int k;
    int z;
    int n;
    int r;
    int r_max;
    struct rcm_state *j;
    struct rcm_statelist *sl;
    struct node *node;
    struct rcm *m = (struct rcm *)R_ExternalPtrAddr(model);
    double g;
    double w;
    double wmax;
    double D;
    double pii;
    double pij;

    int *nodestate;
    SEXP NodeState;

    sl = &m->sl;
    r = sl->len - 1;
    r_max = sl->r;

    if (r == r_max)
    {
        // there are zero non-analog states
        r = r_max - 1;
    }

    n = INTEGER(nmaps)[0];

    NodeState = PROTECT(allocMatrix(INTSXP, m->phy->nnode, n));
    nodestate = INTEGER(NodeState);

    GetRNGstate();

    for (i = 0; i < n; ++i)
    {
        wmax = -HUGE_VAL;
        phy_traverse_prepare(m->phy, m->phy->root, ALL_NODES, PREORDER);
        node = phy_traverse_step(m->phy);

        for (j = sl->head, k = 1; j != NULL; j = j->next, ++k)
        {
            // we want to sample a new state from the distribution
            //
            // exp(w) / sum(exp(w))
            //
            // where w is the conditional likelihood
            //
            // we do this in one pass by computing
            //
            // argmax w + g,
            //
            // where g is a Gumbel(0, 1) random variate
            //
            // for why this works see:
            // https://timvieira.github.io/blog/post/2014/07/31/gumbel-max-trick/
            // or
            // https://arxiv.org/pdf/1706.04161.pdf

            w = log(DCLK(j, node));

            // Gumbel(0,1) random variate
            g = -log(-log(unif_rand()));

            if ((g + w) > wmax)
            {
                wmax = g + w;
                z = k;
            }
        }

        nodestate[node->index + i * m->phy->nnode] = z;

        while ((node = phy_traverse_step(m->phy)) != 0)
        {
            wmax = -HUGE_VAL;
            D = exp(-r_max * m->rate * node->brlen);
            pij = (1 - D) / (double)r_max;
            pii = pij + D;

            for (j = sl->head, k = 1; j != NULL; j = j->next, ++k)
            {
                w = (nodestate[node->anc->index] != k) ? log(pij * DCLK(j, node)) :
                    log(pii * DCLK(j, node));

                g = -log(-log(unif_rand()));

                if ((g + w) > wmax)
                {
                    wmax = g + w;
                    z = k;
                }
            }

            nodestate[node->index + i * m->phy->nnode] = z;
        }
    }

    PutRNGstate();
    UNPROTECT(1);
    return NodeState;
}


/* Simulation free stochastic character maps. Compute expected number of
** character change events (i.e., all state-to-state transitions) on each branch
** when rate of change is identical for all states */

static double rcm_branch_posterior(struct node *node, struct rcm *model)
{
    struct rcm_statelist *sl = &model->sl;
    struct rcm_state *i;
    struct rcm_state *j;
    struct rcm_state *non = sl->tail;
    int r = sl->len - 1;  // number of analog states
    int r_max = sl->r;

    struct node *anc;
    struct node *sib;

    anc = node->anc;

    sib = anc->lfdesc;

    if (sib == node)
        sib = sib->next;

    if (r == r_max)
    {
        // there are zero non-analog states
        r = r_max - 1;
    }

    double D = exp(-r_max * model->rate * node->brlen);
    double pij = (1 - D) / (double)r_max;
    double pii = pij + D;

    double eij = model->rate * (1 - pij) * node->brlen;
    double eii = model->rate * (1 - pii) * node->brlen;

    double n;
    double lk;

    n = UCLK(non, anc) * SCLK(non, sib) * eii * DCLK(non, node)
        + (r_max - r - 1) * UCLK(non, anc) * SCLK(non, sib) * eij * DCLK(non, node);

    lk = DCLK(non, model->phy->root) / (double)(model->r);

    for (j = sl->head; j != non; j = j->next)
    {
        n += (r_max - r) * UCLK(non, anc) * SCLK(non, sib) * eij * DCLK(j, node)
            + UCLK(j, anc) * SCLK(j, sib) * eii * DCLK(j, node);
        n += UCLK(j, anc) * SCLK(j, sib) * eij * DCLK(non, node);

        lk += DCLK(j, model->phy->root) / (double)(model->r);
    }

    for (j = sl->head; j != non->prev; j = j->next)
    {
        for (i = j->next; i != non; i = i->next)
        {
            n += UCLK(i, anc) * SCLK(i, sib) * eij * DCLK(j, node);
            n += UCLK(j, anc) * SCLK(j, sib) * eij * DCLK(i, node);
        }
    }

    return n / lk;
}


static double rcm_branch_prior(struct node *node, struct rcm *model)
{
    double r = (double)(model->r);
    double D = exp(-r * model->rate * node->brlen);
    double pij = (1 - D) / r;
    double pii = pij + D;

    double eij = model->rate * (1 - pij) * node->brlen;
    double eii = model->rate * (1 - pii) * node->brlen;

    return eii + (r-1) * eij;
}


SEXP rcm_stochastic_map_expected_counts(SEXP model)
{
    struct node *node;
    struct rcm *m = (struct rcm *)R_ExternalPtrAddr(model);

    SEXP BSTATS;
    double *bstats;

    BSTATS = PROTECT(allocMatrix(REALSXP, 2, m->phy->nnode));
    bstats = REAL(BSTATS);

    phy_traverse_prepare(m->phy, m->phy->root, ALL_NODES,
        PREORDER);

    // step past root
    phy_traverse_step(m->phy);

    while ((node = phy_traverse_step(m->phy)) != 0)
    {
        bstats[0 + node->index * 2] = rcm_branch_posterior(node, m);
        bstats[1 + node->index * 2] = rcm_branch_prior(node, m);
    }
    bstats[0 + m->phy->root->index * 2] = 0;
    bstats[1 + m->phy->root->index * 2] = 0;

    UNPROTECT(1);
    return BSTATS;
}

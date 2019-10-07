#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Connections.h>
#include "assert.h"
#include "phy.h"

/* These are the values raxml uses for scaling likelihood calculations to avoid
** underflow, taken from axml.h. The scaling technique is described in the book
** chapter by Stamatakis titled "Orchestrating the phylogenetic likelihood function
** on emerging parallel architectures" in the edited volume "Bioinformatics: High
** Performance Parallel Computer Architectures". */

#define twotothe256 \
    115792089237316195423570985008687907853269984665640564039457584007913129639936.0

#define minlikelihood  (1.0/twotothe256)
#define minusminlikelihood -minlikelihood

/* Macros for accessing the downpass, stem, and uppass
** conditional likelihood arrays */
#define DCLK(state, node) (state)->dclk[(node)->index]
#define SCLK(state, node) (state)->sclk[(node)->index]
#define UCLK(state, node) (state)->uclk[(node)->index]


/* Sparse representation of count data */
struct data {
    /* Number of species in the dataset */
    int n;

    /* Records position in the data array */
    int cursor;

    /* Records the end of data position */
    int stop;

    /* Records the current count */
    int cnt;

    /* Records the current category */
    int cat;

    /* Number of resource categories used
    ** by each species */
    int *m;

    /* Beginning of data for each species
    ** in the array of counts */
    int *ofs;

    /* An array of count data for each species.
    ** Each observation is a tuple of (cnt, cat).
    ** Observations for species are contiguous.
    ** Equivalently, a column major matrix
    ** with 2 rows and column count equal to
    ** the number of non-zero elements in a
    ** full nspecies by ncategory matrix plus
    ** one extra column for species with missing
    ** data. */
    int *data;
};


/* Returns >0 if there is data, 0 otherwise.
** So, to loop over the data for a particular
** species use the following idiom,
**
** for (cnt = data_begin(tip->index, model->data);
**          cnt > 0; cnt = data_next(model->data)) {
**      cat = data->cat;
**      // use cnt and cat for something
** }
*/
static int data_begin(int id, struct data *data)
{
    if (id < data->n) {
        data->cursor = data->ofs[id];
        data->stop = data->cursor + 2 * data->m[id];
        data->cnt = data->data[data->cursor++];
        data->cat = data->data[data->cursor++];
        return data->cnt;
    }
    data->cnt = 0;
    data->cat = 0;
    return 0;
}


/* Returns >0 if there is data, 0 otherwise */
static int data_next(struct data *data)
{
    if (data->cursor < data->stop) {
        data->cnt = data->data[data->cursor++];
        data->cat = data->data[data->cursor++];
        return data->cnt;
    }
    data->cursor = data->stop = 0;
    data->cnt = 0;
    data->cat = 0;
    return 0;
}


/* Resource states on the tree
**
** To iterate over the states on the tree we do
**
** struct statelist *sl = &model->sl;
**
** for (state = sl->head; state; state = state->next) {
**     // do something with state
** }
**
*/
struct statelist {
    /* Number of states that are currently on the tree plus 1 */
    int len;

    /* Maximum number of states allowed on the tree (= number of states in Mk model) */
    int len_max;

    /* Number of resource categories in each resource state */
    int nrcat;

    /* Number of nodes in the tree */
    int nnode;

    /* Root node index */
    int root;

    struct state *head;

    /* tail state is always non-analog state (when len < len_max) */
    struct state *tail;

    /* Pointer to the end of the list of states malloc'd but
    ** then removed from the tree that are available for reuse
    ** i.e.,
    **
    ** NULL <-- state <-- state <-- state <-- fl
    */
    struct state *fl;
};


/* Resource state */
struct state {

    /* Identifier */
    int index;

    /* Count of terminal nodes assigned to this state */
    int ntip;

    /* Counts for each resource category */
    int *count;

    /* Downpass conditional likelihood array */
    double *dclk;

    /* Stem conditional likelihood array */
    double *sclk;

    /* Uppass conditional likelihood array */
    double *uclk;

    /* Sampling weight. Dynamic variable used for Gibbs updates */
    double w;

    /* Log likelihood of count data generated from state */
    double loglk;

    /* Memory alloc'd */
    void *mem;

    /* Resource states stored as doubly linked list */
    struct state *next;

    struct state *prev;
};


static struct state *state_add(struct statelist *sl)
{
    sl->len += 1;

    if (sl->len > sl->len_max)
        return sl->tail;

    struct state *state = NULL;
    size_t nbytes = (sl->nrcat + 1) * sizeof(int) + 3 * sl->nnode * sizeof(double);

    if (!sl->fl) {

        state = malloc(sizeof(struct state));
        if (!state)
            error("unable to allocate state");
        state->mem = malloc(nbytes);
        if (!state->mem) {
            free(state);
            error("unable to allocate state memory");
        }
        state->loglk = 0;
        state->ntip = 0;
        state->w = 0;
        state->index = sl->len;
        state->count = (int *)(state->mem);
        state->dclk = (double *)((int *)(state->mem) + (sl->nrcat + 1));
        state->sclk = state->dclk + sl->nnode;
        state->uclk = state->sclk + sl->nnode;
        state->next = NULL;

        memset(state->mem, 0, nbytes);

        state->uclk[sl->root] = 1 / (double)(sl->len_max);

    } else {

        state = sl->fl;
        sl->fl = state->prev;
        state->prev = NULL;
    }

    if (!sl->head) {
        state->prev = NULL;
        sl->head = state;
        return sl->head;
    } else {
        if (sl->tail) {
            sl->tail->next = state;
            state->prev = sl->tail;
            // tail state always interpreted as non-analog state
            memcpy(state->mem, sl->tail->mem, nbytes);

            // this is necessary because model_update_node increments
            // the counts of the state before we know it is non-analog
            memset(state->count, 0, (sl->nrcat+1) * sizeof(int));
        } else {
            sl->head->next = state;
            state->prev = sl->head;
        }
        sl->tail = state;
    }

    return (sl->len < sl->len_max) ? sl->tail->prev : sl->tail;
}


static void state_free(struct state *state)
{
    free(state->mem);
    free(state);
}


static void statelist_free(struct statelist *sl)
{
    struct state *a;
    struct state *b;

    a = b = sl->head;

    if (a) {
        do {
            b = a->next;
            state_free(a);
            a = b;
        } while (b);
    }

    a = b = sl->fl;

    if (a) {
        do {
            b = a->prev;
            state_free(a);
            a = b;
        } while (b);
    }
}


static void state_remove(struct state *state, struct statelist *sl)
{
    struct state *prev = state->prev;
    struct state *next = state->next;

    if (sl->len > sl->len_max)
        error("Attempt to decrement state count after max states reached. "
            "The conditional likelihoods are now erroneous. Increase max states "
            "to avoid this error.");

    if (state == sl->head) {

        sl->head = next;
        next->prev = NULL;

    } else {

        prev->next = next;
        next->prev = prev;
    }

    state->ntip = 0;
    state->next = NULL;
    state->prev = NULL;

    state->prev = sl->fl;
    sl->fl = state;

    sl->len -= 1;
}


struct model {
    /* Number of resource-use categories in the data. */
    int nrcat;

    /* Number of resource-use classes in the model. */
    int nrcls;

    /* Each item is the current resource class assignment for a node */
    struct state **stateid;

    int *stateid_r;

    /* Encodes state information about whether the conditional likelihoods
    ** of a node need to be updated to reflect new partition of terminals */
    int *needsupdate_up;

    int *needsupdate_dn;

    /* Probability that descendant is in same resource-use
    ** class as ancestor */
    double pii;

    /* Probability that descendant is in different resource-use
    ** class than ancestor */
    double pij;

    /* Hyperparameter of the Gamma(alpha, 1) prior on branch lengths */
    double alpha;

    /* Upper bound for alpha */
    double alpha_max;

    /* Hyperparameter of the common Dirichlet prior on the
    ** multinomial distributions assumed to underlay each
    ** resource use class. */
    double *beta;

    /* log gamma function for each beta */
    double *lgamma_beta;

    /* sum of all betas */
    double beta_sum;

    /* log gamma function of beta_sum */
    double lgamma_beta_sum;

    /* downpass scaling factor exponent for each node */
    int *lzd;

    /* uppass scaling factor exponent for each node */
    int *lzu;

    double dataloglk;

    double treeloglk;

    struct phy *phy;

    /* Observed count data at the terminal nodes of phy */
    struct data data;

    struct statelist sl;
};


struct mcmc {
    int niter;

    int sample_freq;

    double update_node;

    double update_alpha;

    double update_beta;

    double update_weight;

    double tune_alpha;

    double tune_beta;

    struct model *model;

    SEXP outputConn;
};


static void model_free(struct model *model)
{
    if (model) {
        free(model->stateid);
        free(model->stateid_r);
        free(model->needsupdate_up);
        free(model->needsupdate_dn);
        free(model->beta);
        free(model->lgamma_beta);
        free(model->data.m);
        free(model->data.ofs);
        free(model->data.data);
        free(model->lzd);
        free(model->lzu);
        statelist_free(&model->sl);
        free(model);
    }
}



/* Log likelihood of count data in state z */
static double model_stateloglk(struct state *z, struct model *model)
{
    int j;
    int n;
    z->loglk = 0;
    n = z->count[model->nrcat];
    for (j = 0; j < model->nrcat; ++j)
        z->loglk += lgammafn(z->count[j] + model->beta[j]) - model->lgamma_beta[j];

    z->loglk += model->lgamma_beta_sum - lgammafn(n + model->beta_sum);

    return z->loglk;
}


static double model_dataloglk(struct model *model)
{
    struct state *i;
    double loglk = 0;
    for (i = model->sl.head; i != NULL; i = i->next)
        loglk += model_stateloglk(i, model);
    return loglk;
}



static void branch_downpass(struct node *node, struct model *model)
{
    struct statelist *sl = &model->sl;
    struct state *i;
    struct state *j;
    struct state *non = sl->tail;
    int r = sl->len - 1;  // number of analog states
    int r_max = sl->len_max;

    if (r == r_max) {
        // there are zero non-analog states
        r = r_max - 1;
    }

    SCLK(non, node) = model->pii * DCLK(non, node) + (r_max - r - 1) * model->pij * DCLK(non, node);

    for (i = sl->head; i != non; i = i->next) {
        SCLK(i, node) = model->pii * DCLK(i, node) + (r_max - r) * model->pij * DCLK(non, node);
        SCLK(non, node) += model->pij * DCLK(i, node);
    }

    for (i = sl->head; i != non->prev; i = i->next) {
        for (j = i->next; j != non; j = j->next) {
            SCLK(i, node) += model->pij * DCLK(j, node);
            SCLK(j, node) += model->pij * DCLK(i, node);
        }
    }
}


static void node_downpass(struct node *node, struct model *model)
{
    int scale;
    struct statelist *sl = &model->sl;
    struct state *i;

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

    if (scale) {
        for (i = sl->head; i != NULL; i = i->next)
            DCLK(i, node) *= twotothe256;

        model->lzd[node->index] += 1;
    }
}


static void model_downpass(struct model *model)
{
    struct node *node;

    phy_traverse_prepare(model->phy, model->phy->root, INTERNAL_NODES_ONLY,
        POSTORDER);

    while ((node = phy_traverse_step(model->phy)) != 0) {
        branch_downpass(node->lfdesc, model);
        branch_downpass(node->lfdesc->next, model);
        node_downpass(node, model);
    }
}


static void node_uppass(struct node *node, struct model *model)
{
    int scale;
    struct statelist *sl = &model->sl;
    struct state *i;
    struct state *j;
    struct state *non = sl->tail;
    int r = sl->len - 1;  // number of analog states
    int r_max = sl->len_max;

    struct node *anc;
    struct node *sib;

    anc = node->anc;

    sib = anc->lfdesc;

    if (sib == node)
        sib = sib->next;

    if (r == r_max) {
        // there are zero non-analog states
        r = r_max - 1;
    }

    UCLK(non, node) = UCLK(non, anc) * SCLK(non, sib) * model->pii
        + (r_max - r - 1) * UCLK(non, anc) * SCLK(non, sib) * model->pij;

    for (j = sl->head; j != non; j = j->next) {
        UCLK(j, node) = (r_max - r) * UCLK(non, anc) * SCLK(non, sib) * model->pij
            + UCLK(j, anc) * SCLK(j, sib) * model->pii;
        UCLK(non, node) += UCLK(j, anc) * SCLK(j, sib) * model->pij;
    }

    for (j = sl->head; j != non->prev; j = j->next) {
        for (i = j->next; i != non; i = i->next) {
            UCLK(j, node) += UCLK(i, anc) * SCLK(i, sib) * model->pij;
            UCLK(i, node) += UCLK(j, anc) * SCLK(j, sib) * model->pij;
        }
    }

    model->lzu[node->index] = model->lzu[anc->index] + model->lzd[sib->index];

    scale = 1;
    for (j = sl->head; j != NULL; j = j->next)
        scale = (UCLK(j, node) < minlikelihood) && (UCLK(j, node) > minusminlikelihood);

    if (scale) {
        for (j = sl->head; j != NULL; j = j->next)
            UCLK(j, node) *= twotothe256;
        model->lzu[node->index] += 1;
    }

}


static void model_uppass(struct model *model)
{
    struct node *node;

    phy_traverse_prepare(model->phy, model->phy->root, ALL_NODES, PREORDER);

    // step past root
    phy_traverse_step(model->phy);

    while ((node = phy_traverse_step(model->phy)) != 0) {
        node_uppass(node, model);
        model->needsupdate_up[node->index] = 0;
        model->needsupdate_dn[node->index] = 0;
    }

}


/* Given an array of conditional log likelihoods
** choose the index of one based on its weight,
** where a high log likelihood means a high
** weight */
static struct state *choose_state(struct statelist *sl)
{
    struct state *i;
    double u;
    double maxin = R_NegInf;
    double norm = 0;

    for (i = sl->head; i != NULL; i = i->next) {
        if (i->w > maxin)
            maxin = i->w;
    }

    for (i = sl->head; i != NULL; i = i->next) {
        i->w = exp(i->w - maxin);
        norm += i->w;
    }

    u = unif_rand() * norm;

    for (i = sl->head; i != NULL; i = i->next) {
        u -= i->w;
        if (u <= 0)
            return i;
    }

    return sl->tail;
}


static void model_update_node(struct node *node, struct model *model)
{
    // these changes ensure that the local configuration of uppass likelihoods
    // is correct for performing updates to terminal nodes. the global
    // configuration of uppass likelihoods is NOT correct, however, and will
    // need to be reset with a full uppass from the root.

    if (model->needsupdate_up[node->anc->index])
        node_uppass(node, model);

    if (model->needsupdate_up[node->index]) {

        struct node *sib;

        // model_update_nodes ensures node will never equal the root node so that
        // anc is always non-NULL
        sib = node->anc->lfdesc;

        // we should always be on a right descendant at this stage
        assert(sib != node);

        branch_downpass(sib, model);
        node_uppass(node, model);
    }

    if (!node->ndesc) {
        struct statelist *sl = &model->sl;
        struct state *z_new;
        struct state *z_old;
        struct state *i;

        int n;
        int m;
        int M;
        int N;

        double A;
        double B;

        z_new = z_old = model->stateid[node->index];

        for (n = data_begin(node->index, &model->data);
                n > 0; n = data_next(&model->data))
        {
            z_old->count[model->data.cat] -= n;
            z_old->count[sl->nrcat] -= n;
        }

        for (i = sl->head; i != NULL; i = i->next) {
            A = 0;
            B = 0;
            N = 0;
            M = i->count[sl->nrcat];

            for (n = data_begin(node->index, &model->data);
                    n > 0; n = data_next(&model->data))
            {
                m = i->count[model->data.cat];

                A += lgammafn(m + n + model->beta[model->data.cat]);
                A -= lgammafn(m + model->beta[model->data.cat]);

                N += n;
            }

            B = lgammafn(M + model->beta_sum)
                 - lgammafn(M + N + model->beta_sum);

            i->w = A + B + log(UCLK(i, node)) + model->lzu[node->index] * log(minlikelihood);
        }

        z_new = choose_state(sl);

        if (z_new->ntip == 0 && z_old->ntip == 1) {
            // update transfers node from singleton state to another singleton
            // state so there is no point in changing the label. in fact, the
            // likelihoods for z_new and z_old are identical in this case.
            z_new = z_old;
        }

        for (n = data_begin(node->index, &model->data);
                n > 0; n = data_next(&model->data))
        {
            z_new->count[model->data.cat] += n;
            z_new->count[sl->nrcat] += n;
        }

        if (z_old != z_new) {

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
            model->dataloglk += model_stateloglk(z_old, model);
            model->dataloglk += model_stateloglk(z_new, model);

            if (!z_old->ntip)
                state_remove(z_old, sl);

            if (node->next)
                model->needsupdate_up[node->next->index] = 1;
            else
                branch_downpass(node, model); // force update now as we won't get another opportunity

            model->needsupdate_dn[node->anc->index] = 1;

        }

        while (node->anc != 0 && node->next == 0) {
            node = node->anc;
            if (model->needsupdate_dn[node->index]) {
                node_downpass(node, model);
                if (node->next)
                    model->needsupdate_up[node->next->index] = 1;
                else
                    branch_downpass(node, model);
                if (node->anc)
                    model->needsupdate_dn[node->anc->index] = 1;
            }
        }
    }
}


static inline double root_loglk(struct model *model)
{
    struct statelist *sl = &model->sl;
    struct state *i;
    struct state *non = sl->tail;
    int r = sl->len - 1;
    int r_max = sl->len_max;

    double lk;

    if (r == r_max) {
        // there are zero non-analog states
        r = r_max - 1;
    }

    lk = (r_max - r) * DCLK(non, model->phy->root) / (double)r_max;

    for (i = sl->head; i != non; i = i->next)
        lk += DCLK(i, model->phy->root) / (double)r_max;

    return log(lk) + model->lzd[model->phy->root->index] * log(minlikelihood);
}


static double model_treeloglk(struct model *model)
{
    model_downpass(model);
    model_uppass(model);
    return root_loglk(model);
}


static void model_update_nodes(struct model *model)
{
    struct node *node;
    struct node *root = model->phy->root;

    phy_traverse_prepare(model->phy, root, ALL_NODES,
        PREORDER);

    phy_traverse_step(model->phy);
    while ((node = phy_traverse_step(model->phy)) != 0)
        model_update_node(node, model);

    model_uppass(model);

    model->treeloglk = root_loglk(model);

    #ifndef NDEBUG
        double delta = model->treeloglk - model_treeloglk(model);
        assert(delta < 1e-12 && delta > -1e-12);
    #endif
}


static void model_init_states(int *stateid, struct model *model)
{
    int i;
    int j;
    int w;
    int inlist;
    int n = 1;
    struct state *state;
    SEXP root = PROTECT(list1(ScalarInteger(stateid[0])));
    SEXP head;
    SEXP tail;

    for (i = 1; i < model->phy->ntip; ++i) {

        head = root;
        tail = CDR(head);

        inlist = 0;
        for (j = 0; j < n && !inlist; ++j) {
            inlist = (INTEGER(CAR(head))[0] == stateid[i]);
            if (tail != R_NilValue) {
                head = tail;
                tail = CDR(head);
            }
        }

        if (!inlist) {
            SETCDR(head, list1(ScalarInteger(stateid[i])));
            ++n;
        }
    }

    // initialize the statelist
    state = state_add(&model->sl);
    state_add(&model->sl);

    head = root;
    for (i = 0; i < model->phy->ntip; ++i) {
        if (INTEGER(CAR(head))[0] == stateid[i]) {
            DCLK(state, phy_getnode_with_index(model->phy, i)) = 1;
            state->ntip += 1;
            model->stateid[i] = state;
            model->stateid_r[i] = state->index;
            for (w = data_begin(i, &model->data);
                w > 0; w = data_next(&model->data))
            {
                state->count[model->data.cat] += w;
                state->count[model->nrcat] += w;
            }
        }
    }

    head = CDR(head);
    for (j = 1; j < n; ++j) {
        state = state_add(&model->sl);
        for (i = 0; i < model->phy->ntip; ++i) {
            if (INTEGER(CAR(head))[0] == stateid[i]) {
                DCLK(state, phy_getnode_with_index(model->phy, i)) = 1;
                state->ntip += 1;
                model->stateid[i] = state;
                model->stateid_r[i] = state->index;
                for (w = data_begin(i, &model->data);
                    w > 0; w = data_next(&model->data))
                {
                    state->count[model->data.cat] += w;
                    state->count[model->nrcat] += w;
                }
            }
        }
        head = CDR(head);
    }

    UNPROTECT(1);
}


static struct model *model_init(
    int nrcat,
    int nrcls,
    int nnz,        // number of non-zero values in data
    int *data,      // nnz by three matrix. row structure (species index, category index, cnt)
    int *stateid,   // initial partition for terminal nodes.
    double alpha,   // initial parameter value for gamma prior on branch lengths
    double *beta,
    struct phy *phy
){
    int i;
    int w;
    int m;
    int ofs;
    double x;
    double f;

    struct model *model;

    model = malloc(sizeof(struct model));

    if (!model)
        return NULL;

    x = pow(1/(nrcls/(double)(nrcls-1)+1), alpha);

    model->nrcat = nrcat;
    model->nrcls = nrcls;
    model->pij = (1 - x) / nrcls;
    model->pii = model->pij + x;
    model->alpha = alpha;
    model->beta_sum = 0;
    model->beta         = calloc(nrcat, sizeof(double));
    model->lgamma_beta  = calloc(nrcat, sizeof(double));
    model->stateid      = malloc(phy->ntip * sizeof(struct state *));
    model->stateid_r    = malloc(phy->ntip * sizeof(int));
    model->needsupdate_up  = calloc(phy->nnode, sizeof(int));
    model->needsupdate_dn  = calloc(phy->nnode, sizeof(int));

    model->data.n       = phy->ntip;
    model->data.m       = calloc(phy->ntip, sizeof(int));
    model->data.ofs     = calloc(phy->ntip, sizeof(int));
    model->data.data    = calloc(2 * (nnz+1), sizeof(int));

    model->lzd          = calloc(phy->nnode, sizeof(int));
    model->lzu          = calloc(phy->nnode, sizeof(int));

    model->sl.len       = 0;
    model->sl.len_max   = nrcls;
    model->sl.nrcat     = nrcat;
    model->sl.root      = phy->root->index;
    model->sl.nnode     = phy->nnode;
    model->sl.head      = NULL;
    model->sl.tail      = NULL;
    model->sl.fl        = NULL;

    model->phy = phy;

    if (!model->beta || !model->lgamma_beta || !model->stateid
        || !model->data.m || !model->data.ofs || !model->data.data
        || !model->lzd || !model->lzu || !model->needsupdate_up
        || !model->needsupdate_dn || !model->stateid_r)
    {
        model_free(model);
        return NULL;
    }

    memcpy(model->beta, beta, nrcat * sizeof(double));

    for (i = 0; i < nrcat; ++i) {
        model->lgamma_beta[i] = lgammafn(beta[i]);
        model->beta_sum += beta[i];
    }
    model->lgamma_beta_sum = lgammafn(model->beta_sum);

    // Set the initial data pointers for all species to the
    // empty element at the end of the data array. This way,
    // species in the phylogeny that are missing data all
    // point to the same element.
    for (i = 0; i < phy->ntip; ++i) {
        model->data.m[i] = 1;
        model->data.ofs[i] = 2 * nnz;
    }

    if (nnz) {
        // nnz will be 0 if matrix(0L, 0L, 3L) passed from R, in which case
        // using the data arg here is undefined behavior
        ofs = 0;
        m = 0;
        w = data[0];
        model->data.ofs[w] = 0;
        for (i = 0; i < nnz; ++i) {
            if (data[i] != w) {
                model->data.m[w] = m;
                w = data[i];
                model->data.ofs[w] = ofs;
                m = 0;
            }
            model->data.data[ofs++] = data[i + 2 * nnz];
            model->data.data[ofs++] = data[i + 1 * nnz];
            ++m;
        }
        // don't forget to record the number of resource cats
        // recorded for the last species
        model->data.m[w] = m;
    }

    model_init_states(stateid, model);

    // this is the minimum fraction of nodes in a different state than their
    // ancestor that is consistent with a finite branch length under the fully
    // symmetric model.
    f = (double)(((phy->nnode - 1) / nrcls) + 1) / (phy->nnode - 1);

    model->alpha_max = -log((f*nrcls - 1) / (nrcls - 1)) /
                            log((nrcls / (double)(nrcls-1)) + 1);

    model->dataloglk = model_dataloglk(model);
    model->treeloglk = model_treeloglk(model);

    return model;
}


static int model_update_alpha_mh(double tune_alpha, struct model *model)
{
    double x;
    double treeloglk;
    double alpha = model->alpha;
    double pij = model->pij;
    double pii = model->pii;

    model->alpha += (unif_rand() - 0.5) * tune_alpha;

    while (model->alpha > model->alpha_max || model->alpha < 0) {
        if (model->alpha > model->alpha_max)
            model->alpha = model->alpha_max - (model->alpha - model->alpha_max);

        if (model->alpha < 0)
            model->alpha = -model->alpha;
    }

    x = pow(1/(model->nrcls/(double)(model->nrcls-1)+1), model->alpha);
    model->pij = (1 - x) / model->nrcls;
    model->pii = model->pij + x;

    treeloglk = model_treeloglk(model);

    if (unif_rand() < exp(treeloglk - model->treeloglk)) {
        model->treeloglk = treeloglk;

        return 1;

    } else {
        model->alpha = alpha;
        model->pij = pij;
        model->pii = pii;

        // necessary to reset the conditional likelihood arrays as to be
        // correct for model->alpha
        model_treeloglk(model);

        return 0;
    }
}


static inline
void model_setalpha(double alpha, struct model *model)
{
    double x;

    x = pow(1 / (model->nrcls/(double)(model->nrcls-1) + 1), alpha);

    model->alpha = alpha;
    model->pij = (1 - x) / model->nrcls;
    model->pii = model->pij + x;
}

// step 1. define vertical level for slice
static double slice(struct model *model)
{
    return (model->treeloglk + log(1 / model->alpha_max)) - rexp(1);
}


// step 2. find slice interval containing x
static void slice_locate(double x, double w, double z, double *L, double *R, struct model *model)
{
    double u;
    double logprior = -log(model->alpha_max);

    u = unif_rand();
    *L = x - w*u;
    *R = *L + w;

    model_setalpha(*L, model);
    while (*L > 0 && (model_treeloglk(model) + logprior) > z) {
        *L -= w;
        model_setalpha(*L, model);
    }

    model_setalpha(*R, model);
    while (*R < model->alpha_max && (model_treeloglk(model) + logprior) > z) {
        *R += w;
        model_setalpha(*R, model);
    }

    *L = (*L < 0) ? 1e-8 : *L;
    *R = (*R > model->alpha_max) ? model->alpha_max : *R;
}


// step 3. sample new point x1 from slice
static void slice_sample(double x, double z, double L, double R, struct model *model)
{
    double u;
    double x1;
    double Lbar = L;
    double Rbar = R;
    double logprior = -log(model->alpha_max);

    do {
        u = unif_rand();
        x1 = Lbar + u * (Rbar - Lbar);

        if (x1 < x)
            Lbar = x1;
        else
            Rbar = x1;

        model_setalpha(x1, model);
        model->treeloglk = model_treeloglk(model);

    } while ((model->treeloglk + logprior) < z);
}


static void model_update_alpha(double tune_alpha, struct model *model)
{
    double x;
    double z;
    double L;
    double R;

    x = model->alpha;
    z = slice(model);

    slice_locate(x, tune_alpha * model->alpha_max, z, &L, &R, model);
    slice_sample(x, z, L, R, model);
}


static int model_update_beta_mh(double tune_beta, struct model *model)
{
    int i;
    int r;
    double dataloglk;
    double beta;
    double beta_sum = model->beta_sum;
    double lgamma_beta;
    double lgamma_beta_sum = model->lgamma_beta_sum;

    i = (int)(unif_rand() * model->nrcat);

    beta = model->beta[i];
    lgamma_beta = model->lgamma_beta[i];

    model->beta[i] += (unif_rand() - 0.5) * tune_beta;

    if (model->beta[i] < 0)
        model->beta[i] = -model->beta[i];

    model->beta_sum -= beta;
    model->beta_sum += model->beta[i];

    model->lgamma_beta[i] = lgammafn(model->beta[i]);
    model->lgamma_beta_sum = lgammafn(model->beta_sum);

    dataloglk = model_dataloglk(model);

    if (unif_rand() < exp(dataloglk - model->dataloglk)) {

        model->dataloglk = dataloglk;

        return 1;

    } else {
        model->beta[i] = beta;
        model->lgamma_beta[i] = lgamma_beta;
        model->beta_sum = beta_sum;
        model->lgamma_beta_sum = lgamma_beta_sum;

        model->dataloglk = model_dataloglk(model);

        return 0;
    }
}


static void model_setbeta(int j, double b, struct model *model)
{
    model->beta_sum = model->beta_sum - model->beta[j] + b;
    model->lgamma_beta_sum = lgammafn(model->beta_sum);
    model->lgamma_beta[j] = lgammafn(b);
    model->beta[j] = b;
}


static void model_update_beta(double tune_beta, struct model *model)
{
    int j;
    double x;
    double z;
    double u;
    double w = tune_beta * 10;
    double logprior = -log(10);
    double L;
    double R;
    double Lbar;
    double Rbar;
    double x1;

    for (j = 0; j < model->nrcat; ++j) {
        x = model->beta[j];
        z = (model->dataloglk + logprior) - rexp(1);

        u = unif_rand();
        L = x - w*u;
        R = L + w;

        model_setbeta(j, L, model);
        while (L > 1e-8 && ((model_dataloglk(model) + logprior) > z)) {
            L -= w;
            model_setbeta(j, L, model);
        }

        model_setbeta(j, R, model);
        while (R < 10 && ((model_dataloglk(model) + logprior) > z)) {
            R += w;
            model_setbeta(j, R, model);
        }

        L = (L < 1e-8) ? 1e-8 : L;
        R = (R > 10) ? 10 : R;

        Lbar = L;
        Rbar = R;
        do {
            u = unif_rand();
            x1 = Lbar + u * (Rbar - Lbar);

            if (x1 < x)
                Lbar = x1;
            else
                Rbar = x1;

            model_setbeta(j, x1, model);
            model->dataloglk = model_dataloglk(model);

        } while ((model->dataloglk + logprior) < z);
    }
}



static void mcmc_sample(struct mcmc *mcmc)
{
    Rconnection con;
    con = R_GetConnection(mcmc->outputConn);
    con->write(&mcmc->model->dataloglk, sizeof(double), 1, con);
    con->write(&mcmc->model->treeloglk, sizeof(double), 1, con);
    con->write(&mcmc->model->alpha, sizeof(double), 1, con);
    con->write(mcmc->model->beta, sizeof(double), mcmc->model->nrcat, con);
    con->write(mcmc->model->stateid_r, sizeof(int), mcmc->model->phy->ntip, con);
}


static void mcmc_step(struct mcmc *mcmc)
{
    double u = unif_rand() * mcmc->update_weight;

    u -= mcmc->update_node;

    if (u <= 0)
        model_update_nodes(mcmc->model);
    else {
        u -= mcmc->update_alpha;
        if (u <= 0)
            model_update_alpha(mcmc->tune_alpha, mcmc->model);
        else
            model_update_beta(mcmc->tune_beta, mcmc->model);
    }
}


static void mcmc_run(struct mcmc *mcmc)
{
    int i;

    mcmc_sample(mcmc);

    for (i = 0; i < mcmc->niter; ++i) {
        mcmc_step(mcmc);

        if (((i+1) % mcmc->sample_freq) == 0)
            mcmc_sample(mcmc);

        if (((i+1) % 1048576) == 0)
            R_CheckUserInterrupt();
    }
}


void mkdmm_model_free(SEXP model)
{
    model_free((struct model *)R_ExternalPtrAddr(model));
}


SEXP mkdmm_model_init(
    SEXP rtree,
    SEXP data,
    SEXP nrcat,
    SEXP nrcls,
    SEXP stateid,
    SEXP alpha,
    SEXP beta)
{
    SEXP exptr;
    struct model *model;

    model = model_init(
        INTEGER(nrcat)[0],
        INTEGER(nrcls)[0],
        INTEGER(getAttrib(data, R_DimSymbol))[0],
        INTEGER(data),
        INTEGER(stateid),
        REAL(alpha)[0],
        REAL(beta),
        (struct phy *)R_ExternalPtrAddr(rtree));

    exptr = PROTECT(R_MakeExternalPtr(model, R_NilValue, R_NilValue));
    R_RegisterCFinalizer(exptr, &mkdmm_model_free);

    UNPROTECT(1);
    return exptr;
}


SEXP mkdmm_mcmc_run(
    SEXP model,
    SEXP niter,
    SEXP thin,
    SEXP tunealpha,
    SEXP tunebeta,
    SEXP updatenode,
    SEXP updatealpha,
    SEXP updatebeta,
    SEXP outputConn
){
    struct mcmc mcmc;

    if (R_CONNECTIONS_VERSION != 1)
        error("mkdmm package code was written for R_CONNECTIONS_VERSION=1 "
            "but this version of R is running R_CONNECTIONS_VERSION=%d", R_CONNECTIONS_VERSION);

    mcmc.niter          = INTEGER(niter)[0];
    mcmc.sample_freq    = INTEGER(thin)[0];
    mcmc.update_node    = REAL(updatenode)[0];
    mcmc.update_alpha   = REAL(updatealpha)[0];
    mcmc.update_beta    = REAL(updatebeta)[0];
    mcmc.tune_alpha     = REAL(tunealpha)[0];
    mcmc.tune_beta      = REAL(tunebeta)[0];
    mcmc.update_weight  = mcmc.update_node + mcmc.update_alpha + mcmc.update_beta;
    mcmc.model          = (struct model *)R_ExternalPtrAddr(model);
    mcmc.outputConn     = outputConn;

    GetRNGstate();

    mcmc_run(&mcmc);

    PutRNGstate();

    return R_NilValue;
}


/* Marginal ancestral state reconstruction */

static double asr_compute(
    struct state *j,
    struct node *node,
    struct model *model)
{
    struct statelist *sl = &model->sl;
    struct state *i;
    struct state *non = sl->tail;
    int r = sl->len - 1;  // number of analog states
    int r_max = sl->len_max;

    double lk = 0;

    struct node *anc;
    struct node *sib;

    anc = node->anc;

    if (r == r_max) {
        // there are zero non-analog states
        r = r_max - 1;
    }

    if (anc) {
        sib = anc->lfdesc;

        if (sib == node)
            sib = sib->next;

        if (j != non) {
            lk = (r_max - r) * UCLK(non, anc) * SCLK(non, sib) * model->pij * DCLK(j, node);
        } else {
            lk = UCLK(non, anc) * SCLK(non, sib) * model->pii * DCLK(j, node)
                + (r_max - r - 1) * UCLK(non, anc) * SCLK(non, sib) * model->pij * DCLK(j, node);
        }

        for (i = sl->head; i != non; i = i->next) {
            if (i != j)
                lk += UCLK(i, anc) * SCLK(i, sib) * model->pij * DCLK(j, node);
            else
                lk += UCLK(i, anc) * SCLK(i, sib) * model->pii * DCLK(j, node);
        }

    } else {

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


SEXP mkdmm_marginal_asr(SEXP model)
{
    int i;
    int r;
    int ntip;
    struct statelist *sl;
    struct state *j;
    double *asr;

    struct model *m;
    struct node *node;

    SEXP ASR;

    m = (struct model *)R_ExternalPtrAddr(model);
    ntip = m->phy->ntip;
    sl = &m->sl;
    r = sl->len - 1;

    if (r < sl->len_max)
        ++r;

    ASR = PROTECT(allocMatrix(REALSXP, r, ntip - 1));
    asr = REAL(ASR);

    phy_traverse_prepare(m->phy, m->phy->root, INTERNAL_NODES_ONLY, POSTORDER);

    while ((node = phy_traverse_step(m->phy)) != 0) {

        for (j = sl->head, i = 0; j != NULL; j = j->next, ++i)
            asr[i + r * (node->index - ntip)] = asr_compute(j, node, m);

        asr_normalize(r, asr + r * (node->index - ntip));
    }

    UNPROTECT(1);
    return ASR;
}


SEXP mkdmm_posterior_multinomial(SEXP model)
{
    int j;
    int k;
    int r;
    struct statelist *sl;
    struct state *i;

    struct model *m;
    struct node *node;

    SEXP result;
    double *d;

    m = (struct model *)R_ExternalPtrAddr(model);
    sl = &m->sl;

    r = sl->len - 1;

    if (r < sl->len_max)
        ++r;

    result = PROTECT(allocMatrix(REALSXP, r, m->nrcat));
    d = REAL(result);

    for (i = sl->head, k = 0; i != NULL; i = i->next, ++k) {
        for (j = 0; j < m->nrcat; ++j)
            d[k + j*r] = i->count[j] + m->beta[j];
    }

    UNPROTECT(1);
    return result;
}

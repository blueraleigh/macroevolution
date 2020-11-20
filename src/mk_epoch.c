#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/BLAS.h>
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
#define DCLK(state, node) model->dclk[(state) + (node)->index * model->nstate]
#define SCLK(state, node) model->sclk[(state) + (node)->index * model->nstate]
#define UCLK(state, node) model->uclk[(state) + (node)->index * model->nstate]

/* Macro to fetch the rate matrix for a given epoch */
#define RATE(epoch) \
    model->rate->rate + (epoch) * model->nstate * model->nstate

/* Macro for accessing transition probability */
#define PROB(i, j) model->pij[(i) + (j) * model->nstate]


/* Structure to facilitate iteration over epoch durations on each branch */
struct epoch {
    /* Number of epochs on each branch */
    int *m;

    /* Offset of each node in the index and len vectors (below) */
    int *ofs;

    /* Indices of the epochs on each branch, ordered from most recent to most
    ** ancient */
    int *index;

    /* Length of time lineage spends in each epoch. */
    double *len;

    /* State information for the iterator. Use the following idiom to
    ** to iterate over epochs on a branch
    **
    ** int i;
    ** double t;
    ** epoch_prepare(e, node);
    ** while (epoch_step(e)) {
    **      i = e->i;     // epoch index
    **      t = e->t;     // length of time spent in that epoch
    ** }
    */
    int cursor;

    int stop;

    /* Direction of iteration. 1 means present to past, -1 means past to present */
    int dir;

    int i;

    double t;
};


static void epoch_prepare(struct epoch *e, struct node *node, int dir)
{
    if (dir > 0) {
        e->dir = 1;
        e->cursor = e->ofs[node->index];
        e->stop = e->cursor + e->m[node->index];
    } else {
        e->dir = -1;
        e->cursor = e->ofs[node->index] + e->m[node->index] - 1;
        e->stop = e->ofs[node->index];
    }
}


static int epoch_step(struct epoch *e)
{
    if (e->dir > 0) {
        if (e->cursor < e->stop) {
            e->i = e->index[e->cursor];
            e->t = e->len[e->cursor++];
            return 1;
        }
        return 0;
    } else {
        if (e->cursor >= e->stop) {
            e->i = e->index[e->cursor];
            e->t = e->len[e->cursor--];
            return 1;
        }
        return 0;
    }
    return 0;
}


static void epoch_free(struct epoch *e)
{
    if (e) {
        free(e->m);
        free(e->ofs);
        free(e->index);
        free(e->len);
        // e itself is never malloc'd because it is
        // a struct member of model (below)
    }
}


/* Epoch rate matrix structure */
struct rate {
    /* number of rate matrices (one for each epoch) */
    int nmat;

    /* number of non-zero transition rate params */
    int npar;

    /* number of states */
    int nstate;

    /* layout[i + j*nstate + k*nstate*nstate] = parameter index for the
    ** i -> j transition rate in rate matrix k. These indices are 1-based,
    ** so a value of 0 means that transition rate is constrained to be 0.
    ** Otherwise, this value minus 1 is the index of a rate parameter passed
    ** to the log likelihood function */
    int *layout;

    /* rate[i + j*nstate + k*nstate*nstate] = i -> j transition rate
    ** in rate matrix k */
    double *rate;
};


static struct rate *rate_init(int nmat, int npar, int nstate,
    int *layout)
{
#define PARAM(i, j, k) \
    rate->layout[(i) + (j)*rate->nstate + (k)*rate->nstate*rate->nstate]

    int i;
    int j;
    int k;
    int n;
    int nelem = nmat*nstate*nstate;

    struct rate *rate;

    rate = malloc(sizeof(struct rate));

    rate->nmat = nmat;
    rate->npar = npar;
    rate->nstate = nstate;
    rate->layout = calloc(nelem, sizeof(int));
    rate->rate = calloc(nelem, sizeof(double));

    memcpy(rate->layout, layout, nelem*sizeof(int));

    return rate;
}


static void rate_set(double *par, struct rate *rate)
{
#define RATE_ELEM(i, j, k) \
    rate->rate[(i) + (j)*rate->nstate + (k)*rate->nstate*rate->nstate]

    int i;
    int j;
    int k;
    int n;

    memset(rate->rate, 0, rate->nmat*rate->nstate*rate->nstate*sizeof(double));

    for (k = 0; k < rate->nmat; ++k) {
        for (i = 0; i < (rate->nstate-1); ++i) {
            for (j = (i+1); j < rate->nstate; ++j) {
                if (PARAM(i, j, k)) {
                    n = PARAM(i, j, k) - 1;
                    RATE_ELEM(i, j, k) = par[n];
                    RATE_ELEM(i, i, k) -= par[n];
                }
                if (PARAM(j, i, k)) {
                    n = PARAM(j, i, k) - 1;
                    RATE_ELEM(j, i, k) = par[n];
                    RATE_ELEM(j, j, k) -= par[n];
                }
            }
        }
    }
}


void mkepoch_rate_free(SEXP rate)
{
    struct rate *r = (struct rate *)R_ExternalPtrAddr(rate);
    if (r) {
        free(r->layout);
        free(r->rate);
        free(r);
    }
    R_ClearExternalPtr(rate);
}


SEXP mkepoch_rate_init(SEXP nmat, SEXP npar, SEXP nstate, SEXP layout)
{
    SEXP exptr;
    SEXP dims;
    SEXP lattr;

    int nelem = INTEGER(nmat)[0] * INTEGER(nstate)[0] * INTEGER(nstate)[0];
    struct rate *rate;

    rate = rate_init(
        INTEGER(nmat)[0],
        INTEGER(npar)[0],
        INTEGER(nstate)[0],
        INTEGER(layout)
    );

    exptr = PROTECT(R_MakeExternalPtr(rate, R_NilValue, R_NilValue));
    R_RegisterCFinalizer(exptr, &mkepoch_rate_free);

    dims = PROTECT(allocVector(INTSXP, 3));
    INTEGER(dims)[0] = INTEGER(nstate)[0];
    INTEGER(dims)[1] = INTEGER(nstate)[0];
    INTEGER(dims)[2] = INTEGER(nmat)[0];

    lattr = PROTECT(allocArray(INTSXP, dims));
    memcpy(INTEGER(lattr), rate->layout, nelem * sizeof(int));

    setAttrib(exptr, install("layout"), lattr);

    UNPROTECT(3);
    return exptr;
}


struct model {
    int nstate;

    /* Number of epochs */
    int nepochs;

    /* Epoch boundaries. Ordered from present to past. i.e., the oldest epoch
    ** comes last. One rate matrix for each epoch. */
    double *epochs;

    /* downpass conditional (scaled) likelihoods */
    double *dclk;

    /* downpass conditional (scaled) stem likelihoods */
    double *sclk;

    /* uppass conditional (scaled) likelihoods */
    double *uclk;

    /* downpass scaling factor exponent for each node */
    int *lzd;

    /* uppass scaling factor exponent for each node */
    int *lzu;

    /* workspace for storing transition probs from matrix exponential */
    double *pij;

    /* workspace for storing stem likelihoods for multi-epoch branches */
    double *init;

    void *mem;

    double loglk;

    struct phy *phy;

    struct epoch e;

    struct rate *rate;
};


BLAS_extern void
F77_NAME(dgpadm)(int*, int*, double*, double*, int*,
    double*, int*, int*, int*, int*, int*);

/* Approximate the matrix exponential exp(Rt) of matrix R (t is a scalar)
** using the Padé approximation and scaling-and-squaring method implemented
** in the Expokit fortran library function dgpadm. */
static void matrix_exponential(double *R, double t, int nrow, double *matexp)
{
    int ideg = 6;
    int m = nrow;
    int lwsp = 4*nrow*nrow + ideg + 1;
    int ipiv[nrow];
    int iexph;
    int ns;
    int iflag;

    double wsp[lwsp];
    double alpha = t;

    F77_CALL(dgpadm)(&ideg, &m, &alpha, R, &m, wsp, &lwsp, ipiv, &iexph, &ns,
        &iflag);

    memcpy(matexp, wsp+iexph-1, nrow * nrow * sizeof(double));
}


static void branch_downpass(struct node *node, struct model *model)
{
    int i;
    int j;
    double *rate;

    // do the first epoch
    epoch_prepare(&model->e, node, 1);
    epoch_step(&model->e);

    rate = RATE(model->e.i);
    matrix_exponential(rate, model->e.t, model->nstate, model->pij);

    for (i = 0; i < model->nstate; ++i)
        SCLK(i, node) = PROB(i, i) * DCLK(i, node);

    for (i = 0; i < (model->nstate-1); ++i) {
        for (j = (i+1); j < model->nstate; ++j) {
            SCLK(i, node) += PROB(i, j) * DCLK(j, node);
            SCLK(j, node) += PROB(j, i) * DCLK(i, node);
        }
    }

    // step through any remaining epochs, using the stem likelihoods from the
    // previous epoch as new "downpass" likelihoods
    while (epoch_step(&model->e)) {
        memcpy(model->init, &(SCLK(0, node)), model->nstate * sizeof(double));

        rate = RATE(model->e.i);
        matrix_exponential(rate, model->e.t, model->nstate, model->pij);

        for (i = 0; i < model->nstate; ++i)
            SCLK(i, node) = PROB(i, i) * model->init[i];

        for (i = 0; i < (model->nstate-1); ++i) {
            for (j = (i+1); j < model->nstate; ++j) {
                SCLK(i, node) += PROB(i, j) * model->init[j];
                SCLK(j, node) += PROB(j, i) * model->init[i];
            }
        }
    }

}


static void node_downpass(struct node *node, struct model *model)
{
    int i;
    int scale;
    struct node *lfdesc;
    struct node *rtdesc;

    lfdesc = node->lfdesc;
    rtdesc = lfdesc->next;

    model->lzd[node->index] = model->lzd[lfdesc->index]
        + model->lzd[rtdesc->index];

    for (i = 0; i < model->nstate; ++i)
        DCLK(i, node) = SCLK(i, lfdesc) * SCLK(i, rtdesc);

    scale = 1;
    for (i = 0; scale && (i < model->nstate); ++i)
        scale = (DCLK(i, node) < minlikelihood) && (DCLK(i, node) > minusminlikelihood);

    if (scale) {
        for (i = 0; i < model->nstate; ++i)
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
    int i;
    int j;
    double *rate;

    struct node *anc;
    struct node *sib;

    anc = node->anc;

    sib = anc->lfdesc;

    if (sib == node)
        sib = sib->next;

    epoch_prepare(&model->e, node, -1);
    epoch_step(&model->e);

    rate = RATE(model->e.i);
    matrix_exponential(rate, model->e.t, model->nstate, model->pij);

    for (j = 0; j < model->nstate; ++j)
        UCLK(j, node) = UCLK(j, anc) * SCLK(j, sib) * PROB(j, j);

    for (j = 0; j < (model->nstate - 1); ++j) {
        for (i = (j+1); i < model->nstate; ++i) {
            UCLK(j, node) += UCLK(i, anc) * SCLK(i, sib) * PROB(i, j);
            UCLK(i, node) += UCLK(j, anc) * SCLK(j, sib) * PROB(j, i);
        }
    }

    // step through any remaining epochs, using the uppass likelihoods from the
    // previous epoch as new ancestral uppass likelihoods and dropping term from
    // stem likelihoods of sib as they are already accounted for
    while (epoch_step(&model->e)) {
        memcpy(model->init, &(UCLK(0, node)), model->nstate * sizeof(double));

        rate = RATE(model->e.i);
        matrix_exponential(rate, model->e.t, model->nstate, model->pij);

        for (j = 0; j < model->nstate; ++j)
            UCLK(j, node) = model->init[j] * PROB(j, j);

        for (j = 0; j < (model->nstate - 1); ++j) {
            for (i = (j+1); i < model->nstate; ++i) {
                UCLK(j, node) += model->init[i] * PROB(i, j);
                UCLK(i, node) += model->init[j] * PROB(j, i);
            }
        }
    }

    model->lzu[node->index] = model->lzu[anc->index] + model->lzd[sib->index];

    scale = 1;
    for (j = 0; scale && (j < model->nstate); ++j)
        scale = (UCLK(j, node) < minlikelihood) && (UCLK(j, node) > minusminlikelihood);

    if (scale) {
        for (j = 0; j < model->nstate; ++j)
            UCLK(j, node) *= twotothe256;
        model->lzu[node->index] += 1;
    }
}


// Note that root_loglk needs to have been called first to initialize
// root node uppass
static void model_uppass(struct model *model)
{
    struct node *node;

    phy_traverse_prepare(model->phy, model->phy->root, ALL_NODES, PREORDER);

    // step past root
    phy_traverse_step(model->phy);

    while ((node = phy_traverse_step(model->phy)) != 0)
        node_uppass(node, model);
}


static double root_loglk(struct model *model)
{
    int i;
    double norm = 0;
    double lk = 0;

    for (i = 0; i < model->nstate; ++i)
        norm += DCLK(i, model->phy->root);

    for (i = 0; i < model->nstate; ++i) {
        UCLK(i, model->phy->root) = DCLK(i, model->phy->root) / norm;
        lk += DCLK(i, model->phy->root) * UCLK(i, model->phy->root);
    }

    return log(lk) + model->lzd[model->phy->root->index] * log(minlikelihood);
}


static double model_treeloglk(double *par, struct model *model)
{
    rate_set(par, model->rate);
    model_downpass(model);
    return root_loglk(model);
}


static void model_free(struct model *model)
{
    if (model) {
        epoch_free(&model->e);
        free(model->mem);
        free(model);
    }
}


static void model_init_epochs(
    struct model *model,
    double *cage,       // crown age (younger)
    double *sage        // stem age (older)
){
    int i;
    int A;
    int B;
    int C;
    int D;
    int nepochs;
    int n = 0;

    double start;
    double stop;

    double *epochs;
    struct phy *phy;
    struct node *node;

    nepochs = model->nepochs;
    epochs = model->epochs;
    phy = model->phy;

    phy_traverse_prepare(phy, phy->root, ALL_NODES, POSTORDER);
    while ((node = phy_traverse_step(phy)) != 0) {

        model->e.ofs[node->index] = n;

        // count the number of epochs on this branch
        for (i = 0; i < nepochs; ++i) {
            // spans epoch i
            A = (sage[node->index] >= epochs[i+1]) &&
                    (cage[node->index] <= epochs[i]);
            // straddles epoch i and i+1
            B = (sage[node->index] >= epochs[i+1]) &&
                    (cage[node->index] <= epochs[i+1] && cage[node->index] >= epochs[i]);
            // straddles epoch i-1, i
            C = (sage[node->index] <= epochs[i+1] && sage[node->index] >= epochs[i]) &&
                    (cage[node->index] <= epochs[i]);
            // contained within epoch i
            D = (sage[node->index] <= epochs[i+1]) &&
                    (cage[node->index] >= epochs[i]);
            if (A || B || C || D) {
                model->e.m[node->index] += 1;
                ++n;
            }
        }
    }

    model->e.index = calloc(n, sizeof(int));
    model->e.len = calloc(n, sizeof(double));

    phy_traverse_prepare(phy, phy->root, ALL_NODES, POSTORDER);
    while ((node = phy_traverse_step(phy)) != 0) {
        n = model->e.ofs[node->index];
        start = cage[node->index];
        for (i = 0; i < nepochs; ++i) {
            A = (sage[node->index] >= epochs[i+1]) &&
                    (cage[node->index] <= epochs[i]);
            B = (sage[node->index] >= epochs[i+1]) &&
                    (cage[node->index] <= epochs[i+1] && cage[node->index] >= epochs[i]);
            C = (sage[node->index] <= epochs[i+1] && sage[node->index] >= epochs[i]) &&
                    (cage[node->index] <= epochs[i]);
            D = (sage[node->index] <= epochs[i+1]) &&
                    (cage[node->index] >= epochs[i]);
            if (A || B || C || D) {
                if (C || D)
                    stop = sage[node->index];
                else
                    stop = epochs[i+1];
                model->e.index[n] = i;
                model->e.len[n++] = stop - start;
                start = stop;
            }
        }
    }
}


static struct model *model_init(int nstate, int nepochs, double *epochs,
    double *dclk, struct phy *phy, struct rate *rate, double *cage, double *sage)
{
    size_t nbytes;
    struct model *model;

    model = malloc(sizeof(struct model));

    if (!model)
        error("unable to allocate memory for Mk model");

    nbytes = 2 * phy->nnode * sizeof(int) +
        (3*nstate*phy->nnode + nstate*nstate + nstate + nepochs+1) * sizeof(double);

    model->mem = malloc(nbytes);
    memset(model->mem, 0, nbytes);

    if (!model->mem) {
        error("unable to allocate memory for Mk model");
        free(model);
    }

    model->nstate = nstate;
    model->nepochs = nepochs;
    model->lzd = (int *)(model->mem);
    model->lzu = model->lzd + phy->nnode;
    model->dclk = (double *)(model->lzu + phy->nnode);
    model->sclk = model->dclk + nstate*phy->nnode;
    model->uclk = model->sclk + nstate*phy->nnode;
    model->pij = model->uclk + nstate*phy->nnode;
    model->init = model->pij + nstate*nstate;
    model->epochs = model->init + nstate;
    model->loglk = 0;
    model->phy = phy;
    model->rate = rate;
    model->e.m = calloc(phy->nnode, sizeof(int));
    model->e.ofs = calloc(phy->nnode, sizeof(int));
    model->e.index = NULL;
    model->e.len = NULL;

    memcpy(model->epochs, epochs, (nepochs+1) * sizeof(double));
    memcpy(model->dclk, dclk, nstate * phy->ntip * sizeof(double));

    model_init_epochs(model, cage, sage);

    return model;
}


void mkepoch_model_free(SEXP model)
{
    model_free((struct model *)R_ExternalPtrAddr(model));
    R_ClearExternalPtr(model);
}


SEXP mkepoch_model_init(SEXP epochs, SEXP clk, SEXP tree, SEXP rate,
    SEXP crownage, SEXP stemage)
{
    SEXP exptr;
    struct model *model;
    struct rate *r;

    r = (struct rate *)R_ExternalPtrAddr(rate);

    model = model_init(
        r->nstate,
        r->nmat,
        REAL(epochs),
        REAL(clk),
        (struct phy *)R_ExternalPtrAddr(tree),
        r,
        REAL(crownage),
        REAL(stemage)
    );

    exptr = PROTECT(R_MakeExternalPtr(model, R_NilValue, R_NilValue));
    R_RegisterCFinalizer(exptr, &mkepoch_model_free);

    UNPROTECT(1);
    return exptr;
}


SEXP mkepoch_loglk(SEXP par, SEXP model)
{
    struct model *m = (struct model *)R_ExternalPtrAddr(model);

    if (LENGTH(par) != m->rate->npar)
        error("Expected %d parameters, received %d", m->rate->npar, LENGTH(par));

    return ScalarReal(model_treeloglk(REAL(par), m));
}



/* Marginal ancestral state reconstruction */

static double asr_compute(
    int j,
    struct node *node,
    struct model *model)
{
    int i;
    double lk = 0;
    double *rate;

    struct node *anc;
    struct node *sib;

    double *init = &(DCLK(0, node));

    anc = node->anc;

    if (anc) {
        sib = anc->lfdesc;

        if (sib == node)
            sib = sib->next;

        if (model->e.m[node->index] > 1) {
            branch_downpass(node, model);
            // model->init workspace now contains the stem likelihoods at the
            // epoch boundary just tipward of the most rootward epoch on this
            // branch. we use these in lieu of DCLK at tipward head of branch
            init = model->init;
        }

        epoch_prepare(&model->e, node, -1);
        epoch_step(&model->e);

        rate = RATE(model->e.i);
        matrix_exponential(rate, model->e.t, model->nstate, model->pij);

        for (i = 0; i < model->nstate; ++i)
            lk += UCLK(i, anc) * SCLK(i, sib) * PROB(i, j) * init[j];

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


SEXP mkepoch_marginal_asr(SEXP par, SEXP model)
{
    int i;
    int ntip;
    double *asr;

    struct model *m;
    struct node *node;

    SEXP ASR;

    m = (struct model *)R_ExternalPtrAddr(model);
    ntip = m->phy->ntip;

    if (LENGTH(par) != m->rate->npar)
        error("Expected %d parameters, received %d", m->rate->npar, LENGTH(par));

    rate_set(REAL(par), m->rate);
    model_downpass(m);
    root_loglk(m);
    model_uppass(m);

    ASR = PROTECT(allocMatrix(REALSXP, m->nstate, ntip - 1));
    asr = REAL(ASR);

    phy_traverse_prepare(m->phy, m->phy->root, INTERNAL_NODES_ONLY, POSTORDER);

    while ((node = phy_traverse_step(m->phy)) != 0) {

        for (i = 0; i < m->nstate; ++i)
            asr[i + m->nstate * (node->index - ntip)] = asr_compute(i, node, m);

        asr_normalize(m->nstate, asr + m->nstate * (node->index - ntip));
    }

    UNPROTECT(1);
    return ASR;
}


SEXP mkepoch_simulate(SEXP par, SEXP m, SEXP nsims, SEXP rootp)
{
    int i;
    int j;
    int nsim = INTEGER(nsims)[0];
    int *nodestate;
    double u;
    double t;
    double lam;
    double norm;

    int pstate;

    struct model *model;
    struct node *d;
    struct node *p;
    struct rate *rate;

    SEXP NodeState;

    model = (struct model *)R_ExternalPtrAddr(m);

    if (LENGTH(par) != model->rate->npar)
        error("Expected %d parameters, received %d", model->rate->npar,
            LENGTH(par));

    rate = model->rate;
    rate_set(REAL(par), rate);

    NodeState = PROTECT(allocMatrix(INTSXP, model->phy->nnode, nsim));
    nodestate = INTEGER(NodeState);

    phy_traverse_prepare(model->phy, model->phy->root, ALL_NODES, PREORDER);

    while ((d = phy_traverse_step(model->phy)) != 0)
    {
        for (i = 0; i < nsim; ++i)
        {
            if (d->anc)
            {
                p = d->anc;
                pstate = nodestate[p->index + i * model->phy->nnode] - 1;
                epoch_prepare(&model->e, d, -1);
                epoch_step(&model->e);
                do {
                    t = model->e.t;
                    do {
                        lam = -1*RATE_ELEM(pstate, pstate, model->e.i);
                        t -= rexp(1 / lam);
                        if (t > 0) {
                            u = unif_rand() * lam;
                            for (j = (pstate+1); j < model->nstate; ++j) {
                                u -= RATE_ELEM(pstate, j, model->e.i);
                                if (u < 0) {
                                    pstate = j;
                                    break;
                                }
                            }
                            if (u > 0) {
                                for (j = (pstate-1); j >= 0; --j) {
                                    u -= RATE_ELEM(pstate, j, model->e.i);
                                    if (u < 0) {
                                        pstate = j;
                                        break;
                                    }
                                }
                            }
                        }
                    } while (t > 0);
                } while (epoch_step(&model->e));
            }
            else
            {
                if (rootp != R_NilValue)
                {
                    if (TYPEOF(rootp) != REALSXP) {
                        UNPROTECT(1);
                        error("root.p is not a numeric vector");
                    }
                    if (LENGTH(rootp) != model->nstate) {
                        UNPROTECT(1);
                        error(
                            "length of root.p does not equal number of states");
                    }
                    norm = 0;
                    for (j = 0; j < model->nstate; ++j)
                        norm += REAL(rootp)[j];
                    u = unif_rand() * norm;
                    for (j = 0; j < model->nstate; ++j) {
                        u -= REAL(rootp)[j];
                        if (u < 0) {
                            pstate = j;
                            break;
                        }
                    }
                }
                else
                {
                    t = 10000;
                    pstate = 0;
                    do {
                        lam = -1 * RATE_ELEM(pstate, pstate, model->rate->nmat-1);
                        t -= rexp(1 / lam);
                        if (t > 0) {
                            u = unif_rand() * lam;
                            for (j = (pstate+1); j < model->nstate; ++j) {
                                u -= RATE_ELEM(pstate, j, model->rate->nmat-1);
                                if (u < 0) {
                                    pstate = j;
                                    break;
                                }
                            }
                            if (u > 0) {
                                for (j = (pstate-1); j >= 0; --j) {
                                    u -= RATE_ELEM(pstate, j, model->rate->nmat-1);
                                    if (u < 0) {
                                        pstate = j;
                                        break;
                                    }
                                }
                            }
                        }
                    } while (t > 0);
                }
            }
            nodestate[d->index + i * model->phy->nnode] = pstate + 1;
        }
    }
    UNPROTECT(1);
    return NodeState;
}

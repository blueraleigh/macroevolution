#include "rcm.h"


/* hyperparameter on Dirichlet prior, shared by all states */
struct hyperparam {
    double betasum;
    double lgamma_betasum;
    double *beta;
    double *lgamma_beta;
};


/* Sparse representation of count data */
struct rcm_data {
    /* Number of species in the dataset */
    int n;

    /* Number of categories in the multinomial */
    int p;

    /* Records position in the data array */
    int cursor;

    /* Records the end of data position */
    int stop;

    /* Records the current count */
    int cnt;

    /* Records the current category */
    int cat;

    /* Records current total number of observations */
    int sz;

    /* Number of resource categories used
    ** by each species */
    int *m;

    /* Total number of observations for each species */
    int *tsz;

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
static int data_prepare(int id, struct rcm_data *data)
{
    if (id < data->n)
    {
        data->sz = data->tsz[id];
        data->cursor = data->ofs[id];
        data->stop = data->cursor + 2 * data->m[id];
        return 1;
    }
    return 0;
}


/* Returns >0 if there is data, 0 otherwise */
static int data_step(struct rcm_data *data)
{
    if (data->cursor < data->stop)
    {
        data->cnt = data->data[data->cursor++];
        data->cat = data->data[data->cursor++];
        return 1;
    }
    return 0;
}


static struct rcm_data *data_alloc(int n, int *data, struct rcm *model)
{
    int i;
    int j;
    int m;
    int sz;
    int ofs;

    struct rcm_data *d = malloc(sizeof(struct rcm_data));

    d->n       = model->phy->ntip;
    d->p       = model->p;
    d->m       = calloc(model->phy->ntip, sizeof(int));
    d->tsz     = calloc(model->phy->ntip, sizeof(int));
    d->ofs     = calloc(model->phy->ntip, sizeof(int));
    d->data    = calloc(2 * n, sizeof(int));

    if (n)
    {
        // n will be 0 if matrix(0L, 0L, 3L) passed from R, in which case
        // using the data arg here is undefined behavior
        ofs = 0;
        m = 0;
        sz = 0;
        j = data[0];
        d->ofs[j] = 0;
        for (i = 0; i < n; ++i)
        {
            if (data[i] != j)
            {
                d->m[j] = m;
                d->tsz[j] = sz;
                j = data[i];
                d->ofs[j] = ofs;
                m = 0;
                sz = 0;
            }
            d->data[ofs++] = data[i + 2 * n];
            d->data[ofs++] = data[i + 1 * n];
            sz += data[i + 2 * n];
            ++m;
        }
        // don't forget the last species
        d->m[j] = m;
        d->tsz[j] = sz;
    }

    return d;
}


static void data_free(struct rcm_data *data)
{
    free(data->m);
    free(data->tsz);
    free(data->ofs);
    free(data->data);
    free(data);
}


struct rcm_stat {
    /* number of multinomial observations */
    int n;

    /* dimension of each (complete) observation */
    int p;

    /* aggregated counts */
    int *count;

    struct hyperparam *h;
};


static void stat_free(struct rcm_stat *stat)
{
    free(stat->count);
    free(stat);
}


static struct rcm_stat *stat_alloc(int p, void *hyperparam)
{
    struct rcm_stat *stat = malloc(sizeof(struct rcm_stat));

    if (!stat)
        error("failed to allocate memory for running stat struct.");
    stat->n = 0;
    stat->p = p;
    stat->count = calloc(p, sizeof(int));
    if (!stat->count)
    {
        free(stat);
        error("failed to allocate memory for running stat struct.");
    }
    stat->h = (struct hyperparam *)(hyperparam);

    return stat;
}


/* Returns predictive log likelihood that observation belongs to state */
static double stat_push(int index, struct rcm_data *data, struct rcm_stat *stat)
{
    int j;
    int sz;
    int x;
    double plnl = 0;

    struct hyperparam *h = (struct hyperparam *)(stat->h);

    sz = data->tsz[index];

    if (sz)
    {
        plnl += lgammafn(h->betasum + stat->n) -
            lgammafn(h->betasum + stat->n + sz);
    }

    data_prepare(index, data);
    while (data_step(data))
    {
        j = data->cat;
        x = data->cnt;
        plnl += lgammafn(stat->count[j] + x + h->beta[j])
            - lgammafn(stat->count[j] + h->beta[j]);
        stat->count[j] += x;
    }

    stat->n += sz;

    return plnl;
}


static void stat_push0(int index, struct rcm_data *data, struct rcm_stat *stat)
{
    int j;
    int sz;
    int x;

    sz = data->tsz[index];

    data_prepare(index, data);
    while (data_step(data))
    {
        j = data->cat;
        x = data->cnt;
        stat->count[j] += x;
    }

    stat->n += sz;
}


static void stat_pop(int index, struct rcm_data *data, struct rcm_stat *stat)
{
    int j;
    int sz;
    int x;

    sz = data->tsz[index];
    stat->n -= sz;

    data_prepare(index, data);
    while (data_step(data))
    {
        j = data->cat;
        x = data->cnt;
        stat->count[j] -= x;
    }
}


static double stat_loglk(struct rcm_stat *stat, struct rcm_data *data)
{
    int j;
    double loglk = 0;

    struct hyperparam *h = (struct hyperparam *)(stat->h);

    for (j = 0; j < stat->p; ++j)
    {
        loglk += lgammafn(stat->count[j] + h->beta[j]) -
            h->lgamma_beta[j];
    }

    loglk += h->lgamma_betasum - lgammafn(stat->n + h->betasum);

    return loglk;
}


void rcm_dmm_free(SEXP model)
{
    struct rcm *m = (struct rcm *)R_ExternalPtrAddr(model);
    struct hyperparam *h = (struct hyperparam *)(m->sl.hyperparam);
    rcm_clear(m);
    data_free(m->data);
    free(h->beta);
    free(h->lgamma_beta);
    free(h);
    free(m);
    R_ClearExternalPtr(model);
}


SEXP rcm_dmm_model_init(
    SEXP rtree,
    SEXP data,
    SEXP p,
    SEXP r,
    SEXP rate,
    SEXP beta,
    SEXP stateid)
{
    int j;
    SEXP exptr;
    struct rcm *model;
    struct hyperparam *h;

    model = rcm_init_start(
        INTEGER(p)[0],
        INTEGER(r)[0],
        REAL(rate)[0],
        (struct phy *)R_ExternalPtrAddr(rtree));

    if (!model)
        error("model allocation failed.");

    model->push = &stat_push;
    model->push0 = &stat_push0;
    model->pop = &stat_pop;
    model->data = data_alloc(INTEGER(getAttrib(data, R_DimSymbol))[0],
        INTEGER(data), model);

    h = malloc(sizeof(struct hyperparam));

    if (!h)
        error("memory allocation failure");

    h->betasum = 0;
    h->lgamma_betasum = 0;
    h->beta = calloc(INTEGER(p)[0], sizeof(double));
    h->lgamma_beta = calloc(INTEGER(p)[0], sizeof(double));

    if (!h->beta || !h->lgamma_beta)
        error("memory allocation failure");

    for (j = 0; j < INTEGER(p)[0]; ++j)
    {
        h->betasum += REAL(beta)[j];
        h->beta[j] = REAL(beta)[j];
        h->lgamma_beta[j] = lgammafn(REAL(beta)[j]);
    }
    h->lgamma_betasum = lgammafn(h->betasum);

    model->sl.hyperparam = (void *)(h);

    model->sl.stat_loglk = &stat_loglk;
    model->sl.stat_alloc = &stat_alloc;
    model->sl.stat_free = &stat_free;

    rcm_init_states(INTEGER(stateid), model);

    model->dataloglk = rcm_dataloglk(model);
    model->treeloglk = rcm_treeloglk(model);

    exptr = PROTECT(R_MakeExternalPtr(model, R_NilValue, R_NilValue));
    R_RegisterCFinalizer(exptr, &rcm_dmm_free);

    UNPROTECT(1);
    return exptr;
}


SEXP rcm_dmm_posterior_multinomial(SEXP model)
{
    int j;
    int k;
    int r;
    struct rcm_statelist *sl;
    struct rcm_state *i;

    struct rcm *m;
    struct node *node;
    struct hyperparam *h;

    SEXP result;
    double *d;

    m = (struct rcm *)R_ExternalPtrAddr(model);
    sl = &m->sl;

    r = sl->len - 1;

    if (r < sl->r)
        ++r;

    h = (struct hyperparam *)(sl->hyperparam);

    result = PROTECT(allocMatrix(REALSXP, r, m->p));
    d = REAL(result);

    for (i = sl->head, k = 0; i != NULL; i = i->next, ++k)
    {
        for (j = 0; j < m->p; ++j)
            d[k + j*r] = i->stat->count[j] + h->beta[j];
    }

    UNPROTECT(1);
    return result;
}


/* Additively decompose trait diversity over nodes sensu
** Pavoine et al. 2010. Ecological Monographs 80(3): 485-207
** where
**
**    tipstate_arr is an integer vector of state id's for the terminals
**
**    prob_arr is a matrix of multinomials corresponding to each state
**
**    perm_arr is a matrix each row containing a permutation order for the
**      rows in prob_arr
**
**    dist_arr is a distance matrix between the rows of prob_arr
**
**    rtree is a pointer to the phylogeny
*/
SEXP rcm_dmm_decomp(SEXP tipstate_arr, SEXP prob_arr, SEXP perm_arr,
    SEXP dist_arr, SEXP rtree)
{
    int i;
    int j;
    int k;
    int n;
    int p;
    int l;
    int z;
    double f0;
    double f1;
    double f2;
    double d;
    double *dij = REAL(dist_arr);
    double *prob = REAL(prob_arr);
    int *tipstate = INTEGER(tipstate_arr);
    int *perm = INTEGER(perm_arr);

    struct node *node;
    struct phy *phy = (struct phy *)R_ExternalPtrAddr(rtree);

    int *tip1 = calloc(phy->nnode, sizeof(int));
    int *tip2 = calloc(phy->nnode, sizeof(int));
    int *ndesc = calloc(phy->nnode, sizeof(int));
    double *h = calloc(phy->nnode, sizeof(double));

    SEXP W;
    double *w;

    // number of states
    p = INTEGER(getAttrib(prob_arr, R_DimSymbol))[0];

    // number of resource categories
    k = INTEGER(getAttrib(prob_arr, R_DimSymbol))[1];

    //  perm_arr is an n by p matrix
    n = INTEGER(getAttrib(perm_arr, R_DimSymbol))[0];

    W = PROTECT(allocMatrix(REALSXP, phy->nnode - phy->ntip, n));
    w = REAL(W);

    for (i = 0; i < phy->ntip; ++i)
        ndesc[i] = 1;

    for (i = phy->root->index; i < phy->nnode; ++i)
    {
        node = phy_getnode_with_index(phy, i);
        tip2[i] = node->lastvisit->index;
        phy_traverse_prepare(phy, node, ALL_NODES, PREORDER);
        while ((node = phy_traverse_step(phy)) != 0)
        {
            if (node->ndesc == 0)
            {
                tip1[i] = node->index;
                break;
            }
        }
    }


    for (l = 0; l < n; ++l)
    {
        memset(h, 0, phy->nnode * sizeof(double));
        phy_traverse_prepare(phy, phy->root, INTERNAL_NODES_ONLY, POSTORDER);

        while ((node = phy_traverse_step(phy)) != 0)
        {
            ndesc[node->index] = ndesc[node->lfdesc->index] +
                ndesc[node->lfdesc->next->index];

            // compute expected distance btwn two terminals drawn
            // at random from among the set descended from node
            for (i = tip1[node->index]; i < tip2[node->index]; ++i)
            {
                for (j = (i+1); j <= tip2[node->index]; ++j)
                {
                    d = dij[perm[l + tipstate[i] * n]
                                + perm[l + tipstate[j] * n] * p];
                    h[node->index] += (d * d) / (double)(ndesc[node->index] * ndesc[node->index]);
                }
            }

            f0 = ndesc[node->index] / (double)(phy->ntip);
            f1 = ndesc[node->lfdesc->index] / (double)(ndesc[node->index]);
            f2 = ndesc[node->lfdesc->next->index] / (double)(ndesc[node->index]);

            w[(node->index - phy->ntip) + l * (phy->nnode - phy->ntip)] =
                f0 * (h[node->index] - (
                    f1 * h[node->lfdesc->index]
                        + f2 * h[node->lfdesc->next->index]));
        }
    }

    free(tip1);
    free(tip2);
    free(ndesc);
    free(h);

    UNPROTECT(1);
    return W;
}


/* Computes the optimal transport matrix between two discrete
** probability distributions using the Sinkhorn-Knopp algorithm
** where
**
**    N is the number of categories in each distribution
**    r and c are the two distributions
**    P is the optimal transport matrix
**    M is the cost matrix
**
** Under a 0,1 cost scheme the optimal cost (minimum cost) to
** transfrom r into c is equal 1/2 the L1 norm, and this function
** computes the transport matrix P that achieves this. P has the
** property that row sums equal r and column sums equal c, and
** each i,j element can be thought of as how much of r[i] is
** transferred to c[j].
*/
void wasserstein(int *N, double *r, double *c, double *P, double *M)
{
    int i;
    int j;
    int n = *N;
    int nelem = n*n;

    // regularization
    double lambda = 10;

    double norm = 0;
    double diff;
    double maxdiff;

    // rowsums
    double u[n];

    // colsums
    double v[n];

    double wkspc[n];

    for (i = 0; i < nelem; ++i)
    {
        P[i] = exp(-lambda * M[i]);
        norm += P[i];
    }

    for (i = 0; i < nelem; ++i)
        P[i] /= norm;

    // P now sums to 1

    #define do_rowsums(u)               \
    do {                                \
        for (i = 0; i < n; ++i)         \
        {                               \
            (u)[i] = 0;                   \
            for (j = 0; j < n; ++j)     \
                (u)[i] += P[i + j*n];     \
        }                               \
    } while (0)

    #define do_colsums                  \
    do {                                \
        for (j = 0; j < n; ++j)         \
        {                               \
            v[j] = 0;                   \
            for (i = 0; i < n; ++i)     \
                v[j] += P[i + j*n];     \
        }                               \
    } while (0)

    #define scale_rowsums                       \
    do {                                        \
        for (j = 0; j < n; ++j)                 \
        {                                       \
            for (i = 0; i < n; ++i)             \
            {                                   \
                if (u[i] > 0)                   \
                    P[i + j*n] *= r[i]/u[i];    \
            }                                   \
        }                                       \
    } while (0)

    #define scale_colsums                       \
    do {                                        \
        for (i = 0; i < n; ++i)                 \
        {                                       \
            for (j = 0; j < n; ++j)             \
            {                                   \
                if (v[j] > 0)                   \
                    P[i + j*n] *= c[j]/v[j];    \
            }                                   \
        }                                       \
    } while (0)

    #define do_maxdiff                      \
    do {                                    \
        for (i = 0; i < n; ++i)             \
        {                                   \
            diff = fabs(wkspc[i] - u[i]);   \
            if (diff > maxdiff)             \
                maxdiff = diff;             \
        }                                   \
    } while(0)

    do {
        maxdiff = 0;
        do_rowsums(u);
        scale_rowsums;
        do_colsums;
        scale_colsums;
        do_rowsums(wkspc);
        do_maxdiff;
    } while (maxdiff > 1e-8);
}


/* Computes the average optimal transport matrix on each branch
** of the phylogeny where
**
**   pij is the result from rcm_pij giving the probability of
**   each ancestor-descendant state condition for each branch
**
**   prob is the matrix of multinomials for each state
**
**   cost is the cost matrix for tranforming a unit of category i
**   into a unit of category j
**
**   rtree is the phylogeny
**
** Note that each state is expected to stored as a column
** in prob whereas they are stored as rows in the R object
** that summarizes a posterior sample, meaning that matrix
** will need to be transposed prior to passing to this
** function. */
SEXP rcm_dmm_doflux(SEXP pij, SEXP prob, SEXP cost, SEXP rtree)
{
    int i;
    int j;
    int k;
    int ndim;
    int nstate;
    double weight;
    double *M;
    double *Pr;
    double *Pij;
    double *flux;
    struct phy *phy;
    struct node *d;

    SEXP Flux;

    phy = (struct phy *)R_ExternalPtrAddr(rtree);

    // number of multinomials
    nstate = INTEGER(getAttrib(prob, R_DimSymbol))[1];

    // number of categories in each multinomial
    ndim = INTEGER(getAttrib(prob, R_DimSymbol))[0];

    Pr = REAL(prob);

    M = REAL(cost);

    phy_traverse_prepare(phy, phy->root, ALL_NODES, PREORDER);
    phy_traverse_step(phy);

    Flux = PROTECT(allocVector(VECSXP, phy->nnode));

    SET_VECTOR_ELT(Flux, phy->root->index, allocMatrix(REALSXP, ndim, ndim));
    memset(REAL(VECTOR_ELT(Flux, phy->root->index)),
        0, ndim*ndim*sizeof(double));

    double U[ndim * ndim];

    while ((d = phy_traverse_step(phy)) != 0)
    {
        SET_VECTOR_ELT(Flux, d->index, allocMatrix(REALSXP, ndim, ndim));

        Pij = REAL(VECTOR_ELT(pij, d->index));
        flux = REAL(VECTOR_ELT(Flux, d->index));

        memset(flux, 0, ndim*ndim*sizeof(double));

        weight = 0;

        for (i = 0; i < (nstate-1); ++i)
        {
            for (j = (i+1); j < nstate; ++j)
            {
                if (Pij[i + j * nstate] > 1e-3)
                {
                    weight += Pij[i + j * nstate];
                    wasserstein(&ndim, Pr + i*ndim, Pr + j*ndim, U, M);
                    for (k = 0; k < ndim*ndim; ++k)
                        flux[k] += Pij[i + j * nstate] * U[k];
                }

                if (Pij[j + i * nstate] > 1e-3)
                {
                    weight += Pij[j + i * nstate];
                    wasserstein(&ndim, Pr + j*ndim, Pr + i*ndim, U, M);
                    for (k = 0; k < ndim*ndim; ++k)
                        flux[k] += Pij[j + i * nstate] * U[k];
                }
            }
        }

        if (weight > 0)
        {
            for (k = 0; k < ndim*ndim; ++k)
                flux[k] /= weight;
        }
    }

    UNPROTECT(1);
    return Flux;
}

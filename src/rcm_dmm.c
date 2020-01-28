#include "rcm.h"

/* hyperparameter on Dirichlet prior, shared by all states */
static double BETASUM = 0;
static double LGAMMA_BETASUM = 0;
static double *BETA = 0;
static double *LGAMMA_BETA = 0;

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
};


static void stat_free(struct rcm_stat *stat)
{
    free(stat->count);
    free(stat);
}


static struct rcm_stat *stat_alloc(int p)
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

    return stat;
}


/* Returns predictive log likelihood that observation belongs to state */
static double stat_push(int index, struct rcm_data *data, struct rcm_stat *stat)
{
    int j;
    int sz;
    int x;
    double plnl = 0;

    sz = data->tsz[index];

    if (sz)
    {
        plnl += lgammafn(BETASUM + stat->n) -
            lgammafn(BETASUM + stat->n + sz);
    }

    data_prepare(index, data);
    while (data_step(data))
    {
        j = data->cat;
        x = data->cnt;
        plnl += lgammafn(stat->count[j] + x + BETA[j])
            - lgammafn(stat->count[j] + BETA[j]);
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

    for (j = 0; j < stat->p; ++j)
    {
        loglk += lgammafn(stat->count[j] + BETA[j]) -
            LGAMMA_BETA[j];
    }

    loglk += LGAMMA_BETASUM - lgammafn(stat->n + BETASUM);

    return loglk;
}


void rcm_dmm_free(SEXP model)
{
    struct rcm *m = (struct rcm *)R_ExternalPtrAddr(model);
    rcm_clear(m);
    data_free(m->data);
    free(m);
    R_ClearExternalPtr(model);
    free(BETA);
    free(LGAMMA_BETA);
    BETA = 0;
    LGAMMA_BETA = 0;
    BETASUM = 0;
    LGAMMA_BETASUM = 0;
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

    model = rcm_init_start(
        INTEGER(p)[0],
        INTEGER(r)[0],
        REAL(rate)[0],
        (struct phy *)R_ExternalPtrAddr(rtree));

    model->push = &stat_push;
    model->push0 = &stat_push0;
    model->pop = &stat_pop;
    model->data = data_alloc(INTEGER(getAttrib(data, R_DimSymbol))[0],
        INTEGER(data), model);

    BETASUM = 0;
    LGAMMA_BETASUM = 0;
    BETA = calloc(INTEGER(p)[0], sizeof(double));
    LGAMMA_BETA = calloc(INTEGER(p)[0], sizeof(double));
    for (j = 0; j < INTEGER(p)[0]; ++j)
    {
        BETASUM += REAL(beta)[j];
        BETA[j] = REAL(beta)[j];
        LGAMMA_BETA[j] = lgammafn(REAL(beta)[j]);
    }
    LGAMMA_BETASUM = lgammafn(BETASUM);

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

    SEXP result;
    double *d;

    m = (struct rcm *)R_ExternalPtrAddr(model);
    sl = &m->sl;

    r = sl->len - 1;

    if (r < sl->r)
        ++r;

    result = PROTECT(allocMatrix(REALSXP, r, m->p));
    d = REAL(result);

    for (i = sl->head, k = 0; i != NULL; i = i->next, ++k)
    {
        for (j = 0; j < m->p; ++j)
            d[k + j*r] = i->stat->count[j] + BETA[j];
    }

    UNPROTECT(1);
    return result;
}


// additively decompose trait diversity over nodes
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

            // compute quadratic entropy
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




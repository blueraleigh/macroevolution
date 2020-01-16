#include "rcm.h"

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


/*
** stat_push and stat_pop are used to efficiently update the mean and
** error sum of squares of a set of multivariate count vectors. We define
** the mean proportion for the j-th component to be p_j =
**
**        ____               /  ____
**        \    (n_ij * n_i) /   \     (n_i * n_i)
**        /___             /    /___
**             i                     i
**
** Note that n_i is the sample size for the i-th observation, which is
** just the sum of n_ij over j. See
** Journal of the Royal Statistical Society. Series C (Applied Statistics), Vol. 38, No. 1 (1989), pp. 71-80
** for a justification for using the weighted proportion rather than the naive
** unweighted MLE for the proportion.
**
** The error sum of squares for the j-th component is
**
**   _____
**   \                        2
**    \     (n_ij - n_i * p_j)
**    /
**   /____
**         i
**
** The following R code can be used to verify that the push and pop
** algorithms are correct

push = function(x, stat) {
    sz = sum(x)
    stat$n = stat$n + sz*sz
    for (i in 1:3)
    {
        oldmu = stat$mu[i]
        stat$mu[i] = stat$mu[i] + (sz*x[i] - oldmu*sz*sz) / stat$n
        stat$ss[i] = stat$ss[i] + (x[i] - sz*oldmu) * (x[i] - sz*stat$mu[i])
    }
}

pop = function(x, stat) {
    sz = sum(x)
    stat$n = stat$n - sz*sz
    for (i in 1:3)
    {
        oldmu = stat$mu[i]
        stat$mu[i] = stat$mu[i] + (oldmu*sz*sz - sz*x[i]) / stat$n
        stat$ss[i] = stat$ss[i] + (sz*oldmu - x[i]) * (x[i] - sz*stat$mu[i])
    }
}

do.test = function() {
    stat = new.env()
    stat$n = 0
    stat$mu = c(0, 0, 0)
    stat$ss = c(0, 0, 0)

    obs = matrix(c(0, 2, 3, 1, 2, 0, 5, 0, 1), 3, 3, byrow=TRUE)

    push(obs[1, ], stat)
    push(obs[2, ], stat)

    mu.expected = colSums(sweep(obs[1:2, ], 1, rowSums(obs[1:2, ]), "*")) /
        sum(rowSums(obs[1:2, ])^2)

    ss.expected = numeric(3)
    ss.expected[1] = (0 - 5*stat$mu[1])^2 + (1 - 3*stat$mu[1])^2
    ss.expected[2] = (2 - 5*stat$mu[2])^2 + (2 - 3*stat$mu[2])^2
    ss.expected[3] = (3 - 5*stat$mu[3])^2 + (0 - 3*stat$mu[3])^2

    A = all.equal(stat$mu, mu.expected)
    B = all.equal(stat$ss, ss.expected)

    push(c(5, 0, 1), stat)
    pop(c(5, 0, 1), stat)

    C = all.equal(stat$mu, mu.expected)
    D = all.equal(stat$ss, ss.expected)

    return(all(A, B, C, D))
}
**
*/
/*
static double stat_push(int index, struct rcm_data *data, struct rcm_stat *stat)
{
    int j;
    int sz;
    double x;
    double oldmu;
    double plnl = 0;

    if (data->tsz[index])
    {
        plnl = stat->tss;
        sz = data->tsz[index];
        stat->n += sz * sz;
        for (j = 0; j < data->p; ++j)
        {
            x = (double)(data->data[index + j * data->n]);
            oldmu = stat->mu[j];
            stat->tss -= stat->ss[j];
            stat->mu[j] += (x*sz - oldmu*sz*sz) / stat->n;
            stat->ss[j] += (x - oldmu*sz) * (x - stat->mu[j]*sz);
            stat->tss += stat->ss[j];
        }
        plnl -= stat->tss;
    }

    return plnl;
}


static void stat_push0(int index, struct rcm_data *data, struct rcm_stat *stat)
{
    stat_push(index, data, stat);
}


static void stat_pop(int index, struct rcm_data *data, struct rcm_stat *stat)
{
    int j;
    int sz;
    double x;
    double oldmu;

    if (data->tsz[index])
    {
        sz = data->tsz[index];
        stat->n -= sz * sz;
        if (stat->n == 0)
        {
            memset(stat->mu, 0, stat->p * sizeof(double));
            memset(stat->ss, 0, stat->p * sizeof(double));
            stat->tss = 0;
        }
        else
        {
            for (j = 0; j < data->p; ++j)
            {
                x = (double)(data->data[index + j * data->n]);
                oldmu = stat->mu[j];
                stat->tss -= stat->ss[j];
                stat->mu[j] += (oldmu*sz*sz - x*sz) / stat->n;
                stat->ss[j] += (oldmu*sz - x) * (x - stat->mu[j]*sz);
                stat->tss += stat->ss[j];
            }
        }
    }
}
*/


struct rcm_stat {
    /* number of multinomial observations */
    int n;

    /* dimension of each (complete) observation */
    int p;

    /* aggregated counts */
    int *count;

    /* sum of all elements in beta vector */
    double betasum;

    /* log gamma(betasum) */
    double lgamma_betasum;

    /* Dirichlet hyperparameter for each state */
    double *beta;

    /* log gamma(beta[j]) for each j */
    double *lgamma_beta;
};


static void stat_free(struct rcm_stat *stat)
{
    free(stat->count);
    free(stat->beta);
    free(stat->lgamma_beta);
    free(stat);
}


static struct rcm_stat *stat_alloc(int p)
{
    struct rcm_stat *stat = malloc(sizeof(struct rcm_stat));
    if (!stat)
        error("failed to allocate memory for running stat struct.");
    stat->n = 0;
    stat->p = p;
    stat->betasum = 0;
    stat->lgamma_betasum = 0;
    stat->count = calloc(p, sizeof(int));
    stat->beta = calloc(p, sizeof(double));
    stat->lgamma_beta = calloc(p, sizeof(double));
    if (!stat->count || !stat->beta || !stat->lgamma_beta)
    {
        free(stat);
        error("failed to allocate memory for running stat struct.");
    }

    /* for now beta is fixed to vector of 1's */
    int i;
    for (i = 0; i < stat->p; ++i)
    {
        stat->betasum += 1;
        stat->beta[i] = 1;
        stat->lgamma_beta[i] = lgammafn(1);
    }
    stat->lgamma_betasum = lgammafn(stat->betasum);

    return stat;
}


/* Returns predictive log likelihood that observation belongs to state */
static double stat_push(int index, struct rcm_data *data, struct rcm_stat *stat)
{
    int j;
    int sz;
    int x;
    double *beta = stat->beta;
    double plnl = 0;

    sz = data->tsz[index];

    if (sz)
    {
        plnl += lgammafn(stat->betasum + stat->n) -
            lgammafn(stat->betasum + stat->n + sz);
    }

    data_prepare(index, data);
    while (data_step(data))
    {
        j = data->cat;
        x = data->cnt;
        plnl += lgammafn(stat->count[j] + x + beta[j])
            - lgammafn(stat->count[j] + beta[j]);
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
        loglk += lgammafn(stat->count[j] + stat->beta[j]) -
            stat->lgamma_beta[j];
    }

    loglk += stat->lgamma_betasum - lgammafn(stat->n + stat->betasum);

    return loglk;
}


void rcm_dmm_free(SEXP model)
{
    struct rcm *m = (struct rcm *)R_ExternalPtrAddr(model);
    rcm_clear(m);
    data_free(m->data);
    free(m);
    R_ClearExternalPtr(model);
}


SEXP rcm_dmm_model_init(
    SEXP rtree,
    SEXP data,
    SEXP p,
    SEXP r,
    SEXP rate,
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
            d[k + j*r] = i->stat->count[j] + i->stat->beta[j];
    }

    UNPROTECT(1);
    return result;
}

#include "rcm.h"


struct rcm_data {
    /* Number of species in the dataset */
    int n;

    /* Number of categories in the multinomial */
    int p;

    /* Total number of observations for each species */
    int *tsz;

    /* An n by p matrix of count data for each species */
    int *data;
};


static struct rcm_data *data_alloc(int n, int *data, struct rcm *model)
{
    int i;
    int j;

    struct rcm_data *d = malloc(sizeof(struct rcm_data));

    d->n       = model->phy->ntip;
    d->p       = model->p;
    d->tsz     = calloc(model->phy->ntip, sizeof(int));
    d->data    = calloc(d->n * d->p, sizeof(int));

    memcpy(d->data, data, d->n * d->p * sizeof(int));

    for (i = 0; i < d->n; ++i)
    {
        for (j = 0; j < d->p; ++j)
            d->tsz[i] += d->data[i + j * d->n];
    }

    return d;
}


static void data_free(struct rcm_data *data)
{
    free(data->tsz);
    free(data->data);
    free(data);
}


struct rcm_stat {
    /* number of observations squared */
    int n;

    /* number of categories for each observation */
    int p;

    /* mean proportion vector for cluster */
    double* mu;

    /* sum of squares within cluster (by each category) */
    double *ss;

    /* total sum of squares within cluster */
    double tss;
};


static void stat_free(struct rcm_stat *stat)
{
    free(stat->mu);
    free(stat->ss);
    free(stat);
}


static struct rcm_stat *stat_alloc(int p, void *hyperparam)
{
    struct rcm_stat *stat = malloc(sizeof(struct rcm_stat));
    if (!stat)
        error("failed to allocate memory for running stat struct.");
    stat->n = 0;
    stat->p = p;
    stat->tss = 0;
    stat->mu = calloc(p, sizeof(double));
    stat->ss = calloc(p, sizeof(double));
    if (!stat->mu || !stat->ss)
    {
        free(stat);
        error("failed to allocate memory for running stat struct.");
    }

    return stat;
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
static double stat_push(int index, struct rcm_data *data, struct rcm_stat *stat)
{
    int j;
    int sz;
    double x;
    double oldmu;

    /* total sum of squares without contribution from new species */
    double tss = stat->tss;

    if (data->tsz[index])
    {
        sz = data->tsz[index];
        stat->n += sz * sz;
        for (j = 0; j < data->p; ++j)
        {
            x = (double)(data->data[index + j * data->n]);
            oldmu = stat->mu[j];
            stat->tss -= stat->ss[j];
            stat->mu[j] += (x*sz - oldmu*sz*sz) / (double)(stat->n);
            stat->ss[j] += (x - oldmu*sz) * (x - stat->mu[j]*sz);
            stat->tss += stat->ss[j];
        }
    }

    /* predictive log likelihood is the negative of the increase in sum of
    ** squares resulting from addition of new species to cluster. so,
    ** observations are preferentially added to clusters where they are
    ** closest to the mean of the cluster. */
    return -1 * stat->tss + tss;
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
                stat->mu[j] += (oldmu*sz*sz - x*sz) / (double)(stat->n);
                stat->ss[j] += (oldmu*sz - x) * (x - stat->mu[j]*sz);
                stat->tss += stat->ss[j];
            }
        }
    }
}


/*
** We directly assign the marginal likelhood of count data within each
** cluster to be exp(-psi), where psi is defined to be the error sum of
** squares.
*/
static double stat_loglk(struct rcm_stat *stat, struct rcm_data *data)
{
    return -1 * stat->tss;
}


void rcm_kmeans_counts_free(SEXP model)
{
    struct rcm *m = (struct rcm *)R_ExternalPtrAddr(model);
    rcm_clear(m);
    data_free(m->data);
    free(m);
    R_ClearExternalPtr(model);
}


SEXP rcm_kmeans_counts_model_init(
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
    R_RegisterCFinalizer(exptr, &rcm_kmeans_counts_free);

    UNPROTECT(1);
    return exptr;
}

#include "rcm.h"

/* Sparse representation of multivariate continuous data */
struct rcm_data {
    /* Number of species in the dataset */
    int n;

    /* Records position in the p array */
    int cursor1;

    /* Records position in the idx, data arrays */
    int cursor2;

    /* Records the end of data position */
    int stop;

    /* Number of individuals measured for each
    ** species */
    int *m;

    /* Number of non-missing measurements for each
    ** individual of each species. */
    int *p;

    /* Beginning of data for each species
    ** in the p, idx, and data arrays */
    int *ofs;

    /* Measurement indices for each datum in data. i.e.,
    ** if we have (NA,1.34,0.54,NA,3.42) for an idividual
    ** we store (1,2,4) in the idx array and (1.34,0.54,3.42)
    ** in the data array. */
    int *idx;

    /* An array of measurement data for each species.
    ** Observations for species are contiguous. */
    double *data;

    /* record positional info for current measurement vector */
    struct {
        int len;
        int *j;
        double *x;
    } datum;
};


static void data_free(struct rcm_data *data)
{
    free(data->m);
    free(data->p);
    free(data->ofs);
    free(data->idx);
    free(data->data);
    free(data);
}


/*
** Data is transformed on the R side to a four column matrix with column
** structure,
**
** tip index, individual index, measurement index, measurement value
*/
static struct rcm_data *data_alloc(SEXP data, struct rcm *model)
{
    int i;
    int m;
    int sz;
    int sid;
    int iid;
    int ofs;
    int nnz = INTEGER(getAttrib(data, R_DimSymbol))[0];

    double *measure = REAL(data);

    struct rcm_data *d = malloc(sizeof(struct rcm_data));
    d->n       = model->phy->ntip;
    d->m       = calloc(model->phy->ntip, sizeof(int));
    d->ofs     = calloc(model->phy->ntip, sizeof(int));
    d->p       = calloc(nnz, sizeof(int));
    d->idx     = calloc(nnz, sizeof(int));
    d->data    = calloc(nnz, sizeof(double));

    if (nnz)
    {
        m = 0;
        ofs = 0;
        sid = (int)(measure[0]);
        iid = (int)(measure[0 + 1*nnz]);
        d->ofs[sid] = 0;
        for (i = 0; i < nnz; ++i)
        {
            if ((int)(measure[i + 1*nnz]) != iid || (int)(measure[i]) != sid)
            {
                d->p[d->ofs[sid] + m] = sz;
                iid = (int)(measure[i + 1*nnz]);
                sz = 0;
                ++m;
            }
            if ((int)(measure[i]) != sid)
            {
                d->m[sid] = m;
                sid = (int)(measure[i]);
                d->ofs[sid] = ofs;
                m = 0;
                sz = 0;
            }
            d->idx[ofs] = (int)(measure[i + 2*nnz]);
            d->data[ofs++] = measure[i + 3*nnz];
            ++sz;
        }
        d->p[d->ofs[sid] + m] = sz;
        d->m[sid] = m;
    }

    return d;
}


static int data_prepare(int id, struct rcm_data *data)
{
    if (id < data->n)
    {
        data->cursor1 = data->cursor2 = data->ofs[id];
        data->stop = data->cursor1 + data->m[id];
        return 1;
    }
    return 0;
}


/* Loop over individual data vectors
**
** Returns >0 if there is data, 0 otherwise.
** So, to loop over the data for a particular
** species use the following idiom,
**
** data_prepare(tip->index, data);
** while (data_step(data));
**     // do something
** }
**
*/
static int data_step(struct rcm_data *data)
{
    if (data->cursor1 < data->stop)
    {
        data->datum.len = data->p[data->cursor1];
        data->datum.j = &data->idx[data->cursor2];
        data->datum.x = &data->data[data->cursor2];
        data->cursor2 += data->p[data->cursor1++];
        return 1;
    }
    return 0;
}


struct rcm_stat {
    /* number of observations */
    int n;

    /* dimension of each (complete) observation */
    int p;

    /* average of each dimension */
    double* mu;

    /* sum of squares of each dimension */
    double *ss;

    /* total sum of squares */
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


static double stat_push(int index, struct rcm_data *data, struct rcm_stat *stat)
{
    int i;
    int j;
    double x;
    double oldmu;

    /* total sum of squares without contribution from new species */
    double tss = stat->tss;

    data_prepare(index, data);

    while (data_step(data))
    {
        stat->n += 1;
        if (stat->n == 1)
        {
            for (i = 0; i < data->datum.len; ++i)
            {
                j = data->datum.j[i];
                x = data->datum.x[i];
                stat->mu[j] = x;
                stat->ss[j] = 0;
            }
            stat->tss = 0;
        }
        else
        {
            for (i = 0; i < data->datum.len; ++i)
            {
                j = data->datum.j[i];
                x = data->datum.x[i];
                oldmu = stat->mu[j];
                stat->tss -= stat->ss[j];
                stat->mu[j] += (x - oldmu) / stat->n;
                stat->ss[j] += (x - oldmu) * (x - stat->mu[j]);
                stat->tss += stat->ss[j];
            }
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
    int i;
    int j;
    double x;
    double oldmu;

    data_prepare(index, data);

    while (data_step(data))
    {
        stat->n -= 1;
        if (stat->n == 0)
        {
            memset(stat->mu, 0, stat->p * sizeof(double));
            memset(stat->ss, 0, stat->p * sizeof(double));
            stat->tss = 0;
        }
        else
        {
            for (i = 0; i < data->datum.len; ++i)
            {
                j = data->datum.j[i];
                x = data->datum.x[i];
                oldmu = stat->mu[j];
                stat->tss -= stat->ss[j];
                stat->mu[j] += (oldmu - x) / stat->n;
                stat->ss[j] += (oldmu - x) * (x - stat->mu[j]);
                stat->tss += stat->ss[j];
            }
        }
    }
}


static double stat_loglk(struct rcm_stat *stat, struct rcm_data *data)
{
    return -1 * stat->tss;
}


void rcm_kmeans_cont_free(SEXP model)
{
    struct rcm *m = (struct rcm *)R_ExternalPtrAddr(model);
    rcm_clear(m);
    data_free(m->data);
    free(m);
    R_ClearExternalPtr(model);
}


SEXP rcm_kmeans_cont_model_init(
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
    model->data = data_alloc(data, model);

    model->sl.stat_loglk = &stat_loglk;
    model->sl.stat_alloc = &stat_alloc;
    model->sl.stat_free = &stat_free;

    rcm_init_states(INTEGER(stateid), model);

    model->dataloglk = rcm_dataloglk(model);
    model->treeloglk = rcm_treeloglk(model);

    exptr = PROTECT(R_MakeExternalPtr(model, R_NilValue, R_NilValue));
    R_RegisterCFinalizer(exptr, &rcm_kmeans_cont_free);

    UNPROTECT(1);
    return exptr;
}

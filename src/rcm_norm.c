#include "rcm.h"
#include "R_ext/Lapack.h"


// modified from function det_ge_real, lines 1260-1316 in
// src/modules/lapack/Lapack.c in R source tree
static double logdet(int n, double *Ain) {
    if (n == 1)
        return log(Ain[0]);
    int i;
    int info;
    int jpvt[n];
    double dii;
    double modulus = 0.0;
    double A[n*n];
    memset(jpvt, 0, n * sizeof(int));
    memcpy(A, Ain, n*n*sizeof(double));
    F77_CALL(dgetrf)(&n, &n, A, &n, jpvt, &info);
    if (info != 0)
        error("Lapack dgetrf() did not return success");
    for (i = 0; i < n; i++) {
        dii = A[i + i * n]; /* i-th diagonal element */
        modulus += log(dii < 0 ? -dii : dii);
    }
    return modulus;
}


static double multigammaln(double x, int d) {
    int i;
    double val = 0.25 * d * (d-1) * log(M_PI);
    for (i = 1; i <= d; ++i) {
        val += lgammafn( x + (1-i)/2 );
    }
    return val;
}


static double calc_log_z(
    double p, double nu, double kappa, double *mu, double *lambda)
{
    return M_LN2*(nu*p/2)
        + (p/2)*log(M_2PI/kappa)
        + multigammaln(nu/2, p)
        - (nu/2)*logdet(p, lambda);
}

/* hyperparameter on Normal-inverse Wishart prior, shared by all states */
struct hyperparam {
    /* prior degrees of freedom */
    double nu0;
    /* prior number of observations */
    double kappa0;
    /* prior mean vector */
    double *mu0;
    /* prior sum of squares matrix */
    double *lambda0;
    /* log normalization constant */
    double logz0;
};


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
        int len;     // number of measurements
        int *j;      // measurement indices
        double *x;   // measurement values
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
        m = 0;  /* number of individuals measured for a species */
        sz = 0; /* number of measurements for an individual */
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
            d->idx[ofs] = (int)(measure[i + 2*nnz]); /* character index */
            d->data[ofs++] = measure[i + 3*nnz];     /* character value */
            ++sz;
        }
        // the two if branches inside the above loop will never be
        // entered when we reach the last row of data, so the
        // following lines are necessary to pick up the last species
        // or individual
        d->p[d->ofs[sid] + m] = sz;
        d->m[sid] = ++m;
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
**     // loop over measurement vector for current datum
**     for (int i = 0; i < data->datum.len; ++i) {
**          // measurement index
**          int j = data->datum.j[i];
**          // measurement value
**          double x = data->datum.x[i];
**     }
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

    /* trait means for each dimension */
    double *xbar;

    /* sum of squares matrix */
    double *S;

    /* Normal-inverse Wishart hyperparameters */
    struct hyperparam *h;
};


static void stat_free(struct rcm_stat *stat)
{
    free(stat->xbar);
    free(stat->S);
    free(stat);
}


static struct rcm_stat *stat_alloc(int p, void *hyperparam)
{
    struct rcm_stat *stat = malloc(sizeof(struct rcm_stat));
    if (!stat)
        error("failed to allocate memory for running stat struct.");
    stat->n = 0;
    stat->p = p;
    stat->xbar = calloc(p, sizeof(double));
    stat->S = calloc(p*p, sizeof(double));
    if (!stat->xbar || !stat->S)
    {
        free(stat);
        error("failed to allocate memory for running stat struct.");
    }
    stat->h = (struct hyperparam *)(hyperparam);

    return stat;
}


/* marginal log likelihood of data in a state */
static double stat_loglk(struct rcm_stat *stat, struct rcm_data *data)
{
    int i;
    int j;
    int p = stat->p;
    double nu;
    double kappa;
    double logz;
    double loglk = 0;
    double mu[p];
    double dt[p];
    double lambda[p * p];
    double back[p * p];

    struct hyperparam *h = (struct hyperparam *)(stat->h);

    nu = h->nu0 + stat->n;
    kappa = h->kappa0 + stat->n;

    for (j = 0; j < p; ++j) {
        mu[j] = (h->kappa0 * h->mu0[j] + stat->n * stat->xbar[j]) / kappa;
        dt[j] = stat->xbar[j] - h->mu0[j];
    }

    for (i = 0; i < p; ++i) {
        for (j = 0; j <= i; ++j) {
            back[i + j * p] = dt[i] * dt[j];
            back[j + i * p] = back[i + j * p];
        }
    }

    for (j = 0; j < p*p; ++j) {
        lambda[j] = h->lambda0[j]
            + stat->S[j] + (h->kappa0*stat->n/kappa)*back[j];
    }
    logz = calc_log_z((double)p, nu, kappa, mu, lambda);

    loglk = logz - h->logz0 - M_LN_2PI * (stat->n * p / 2);
    return loglk;
}


static void stat_push0(int index, struct rcm_data *data, struct rcm_stat *stat)
{
    int i;
    int j;
    int k;
    int p;
    double x;
    double y;
    double oldmu[stat->p];
    memcpy(oldmu, stat->xbar, stat->p * sizeof(double));

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
                stat->xbar[j] = x;
            }
        }
        else
        {
            for (i = 0; i < data->datum.len; ++i)
            {
                j = data->datum.j[i];
                x = data->datum.x[i];
                stat->xbar[j] += (x - oldmu[j]) / stat->n;
            }
            for (i = 0; i < data->datum.len; ++i) {
                x = data->datum.x[i];
                p = data->datum.j[i];
                for (k = 0; k <= i; ++k) {
                    y = data->datum.x[k];
                    j = data->datum.j[k];
                    stat->S[p+j*stat->p] += (x - oldmu[p]) * (y - stat->xbar[j]);
                    stat->S[j+p*stat->p] = stat->S[p+j*stat->p];
                }
            }
        }
    }
}


static double stat_push(int index, struct rcm_data *data, struct rcm_stat *stat)
{
    double l0;
    double l1;
    l0 = stat_loglk(stat, data);
    stat_push0(index, data, stat);
    l1 = stat_loglk(stat, data);
    return l1 - l0;
}


static void stat_pop(int index, struct rcm_data *data, struct rcm_stat *stat)
{
    int i;
    int j;
    int k;
    int p;
    double x;
    double y;
    double oldmu[stat->p];
    memcpy(oldmu, stat->xbar, stat->p * sizeof(double));

    data_prepare(index, data);

    while (data_step(data))
    {
        stat->n -= 1;
        if (stat->n == 0)
        {
            memset(stat->xbar, 0, stat->p * sizeof(double));
            memset(stat->S, 0, stat->p * stat->p * sizeof(double));
        }
        else
        {
            for (i = 0; i < data->datum.len; ++i)
            {
                j = data->datum.j[i];
                x = data->datum.x[i];
                stat->xbar[j] += (oldmu[j] - x) / stat->n;
            }
            for (i = 0; i < data->datum.len; ++i) {
                x = data->datum.x[i];
                p = data->datum.j[i];
                for (k = 0; k <= i; ++k) {
                    y = data->datum.x[k];
                    j = data->datum.j[k];
                    stat->S[p+j*stat->p] += (oldmu[p] - x) * (y - stat->xbar[j]);
                    stat->S[j+p*stat->p] = stat->S[p+j*stat->p];
                }
            }
        }
    }
}


void rcm_norm_free(SEXP model)
{
    struct rcm *m = (struct rcm *)R_ExternalPtrAddr(model);
    struct hyperparam *h = (struct hyperparam *)(m->sl.hyperparam);
    rcm_clear(m);
    data_free(m->data);
    free(h->mu0);
    free(h->lambda0);
    free(h);
    free(m);
    R_ClearExternalPtr(model);
}


SEXP rcm_norm_model_init(
    SEXP rtree,
    SEXP data,
    SEXP p,
    SEXP r,
    SEXP integrate_brlen,
    SEXP rate,
    SEXP nu,
    SEXP kappa,
    SEXP mu,
    SEXP lambda,
    SEXP stateid)
{
    int j;
    SEXP exptr;
    struct rcm *model;
    struct hyperparam *h;

    model = rcm_init_start(
        INTEGER(p)[0],
        INTEGER(r)[0],
        INTEGER(integrate_brlen)[0],
        REAL(rate)[0],
        (struct phy *)R_ExternalPtrAddr(rtree));

    if (!model)
        error("model allocation failed.");

    model->push = &stat_push;
    model->push0 = &stat_push0;
    model->pop = &stat_pop;
    model->data = data_alloc(data, model);

    h = malloc(sizeof(struct hyperparam));

    if (!h)
        error("memory allocation failure");

    h->mu0 = calloc(INTEGER(p)[0], sizeof(double));
    h->lambda0 = calloc(INTEGER(p)[0]*INTEGER(p)[0], sizeof(double));

    if (!h->mu0 || !h->lambda0)
        error("memory allocation failure");

    h->nu0 = REAL(nu)[0];
    h->kappa0 = REAL(kappa)[0];
    memcpy(h->mu0, REAL(mu), INTEGER(p)[0]*sizeof(double));
    memcpy(h->lambda0, REAL(lambda), INTEGER(p)[0]*INTEGER(p)[0]*sizeof(double));
    h->logz0 = calc_log_z((double)(INTEGER(p)[0]), h->nu0, h->kappa0, h->mu0, h->lambda0);

    model->sl.hyperparam = (void *)(h);

    model->sl.stat_loglk = &stat_loglk;
    model->sl.stat_alloc = &stat_alloc;
    model->sl.stat_free = &stat_free;

    rcm_init_states(INTEGER(stateid), model);

    model->dataloglk = rcm_dataloglk(model);
    model->treeloglk = (INTEGER(r)[0] > 1) ? rcm_treeloglk(model) : 0;

    exptr = PROTECT(R_MakeExternalPtr(model, R_NilValue, R_NilValue));
    R_RegisterCFinalizer(exptr, &rcm_norm_free);

    UNPROTECT(1);
    return exptr;
}


SEXP rcm_norm_posterior_normal(SEXP model)
{
    int j;
    int k;
    int l;
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

    result = PROTECT(allocMatrix(REALSXP, r, 2+m->p+m->p*m->p));
    d = REAL(result);

    double dt[m->p];
    double back[m->p * m->p];

    for (i = sl->head, k = 0; i != NULL; i = i->next, ++k)
    {
        d[k + 0*r] = i->stat->n + h->kappa0;
        d[k + 1*r] = i->stat->n + h->nu0;
        for (j = 0; j < m->p; ++j) {
            d[k + (j+2)*r] = (h->kappa0 * h->mu0[j] +
                i->stat->n * i->stat->xbar[j]) / d[k + 0*r];

            dt[j] = i->stat->xbar[j] - h->mu0[j];
        }

        for (j = 0; j < m->p; ++j) {
            for (l = 0; l <= j; ++l) {
                back[j + l * m->p] = dt[j] * dt[l];
                back[l + j * m->p] = back[j + l * m->p];
            }
        }

        for (j = 0; j < m->p*m->p; ++j) {
            d[k + (j+2+m->p)*r] = h->lambda0[j] +
                i->stat->S[j] + (h->kappa0*i->stat->n/d[k + 0*r])*back[j];
        }
    }

    UNPROTECT(1);
    return result;
}

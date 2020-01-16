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
#define DCLK(state, node) model->dclk[(state) + (node)->index * 2]
#define SCLK(state, node) model->sclk[(state) + (node)->index * 2]
#define UCLK(state, node) model->uclk[(state) + (node)->index * 2]

struct model {

    double eps;

    double tau;

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

    void *mem;

    double loglk;

    struct phy *phy;
};


static void branch_downpass(struct node *node, struct model *model)
{
    int i;
    int j;
    double p01;
    double p10;
    double p00;
    double p11;

    double eps = model->eps;
    double tau = model->tau;
    double D = exp(-tau * node->brlen);

    p01 = (1 - D) * (eps) / (1+eps);
    p10 = (1 - D) / (1 + eps);
    p00 = (1 / (1 + eps)) + D * (eps / (1 + eps));
    p11 = (eps / (1 + eps)) + D / (1 + eps);

    SCLK(0, node) =  p00 * DCLK(0, node);
    SCLK(1, node) =  p11 * DCLK(1, node);
    SCLK(0, node) += p01 * DCLK(1, node);
    SCLK(1, node) += p10 * DCLK(0, node);
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

    for (i = 0; i < 2; ++i)
        DCLK(i, node) = SCLK(i, lfdesc) * SCLK(i, rtdesc);

    scale = 1;
    for (i = 0; scale && (i < 2); ++i)
        scale = (DCLK(i, node) < minlikelihood) && (DCLK(i, node) > minusminlikelihood);

    if (scale) {
        for (i = 0; i < 2; ++i)
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
    double p01;
    double p10;
    double p00;
    double p11;

    double eps = model->eps;
    double tau = model->tau;
    double D = exp(-tau * node->brlen);

    p01 = (1 - D) * (eps) / (1+eps);
    p10 = (1 - D) / (1 + eps);
    p00 = (1 / (1 + eps)) + D * (eps / (1 + eps));
    p11 = (eps / (1 + eps)) + D / (1 + eps);

    struct node *anc;
    struct node *sib;

    anc = node->anc;

    sib = anc->lfdesc;

    if (sib == node)
        sib = sib->next;

    UCLK(0, node)  = UCLK(0, anc) * SCLK(0, sib) * p00;
    UCLK(1, node)  = UCLK(1, anc) * SCLK(1, sib) * p11;
    UCLK(0, node) += UCLK(1, anc) * SCLK(1, sib) * p10;
    UCLK(1, node) += UCLK(0, anc) * SCLK(0, sib) * p01;

    model->lzu[node->index] = model->lzu[anc->index] + model->lzd[sib->index];

    scale = 1;
    for (j = 0; scale && (j < 2); ++j)
        scale = (UCLK(j, node) < minlikelihood) && (UCLK(j, node) > minusminlikelihood);

    if (scale) {
        for (j = 0; j < 2; ++j)
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
    //double norm = 0;
    double lk = 0;

    //for (i = 0; i < 2; ++i)
    //    norm += DCLK(i, model->phy->root);

    //for (i = 0; i < 2; ++i) {
    //    UCLK(i, model->phy->root) = DCLK(i, model->phy->root) / norm;
    //    lk += DCLK(i, model->phy->root) * UCLK(i, model->phy->root);
    //}

    UCLK(0, model->phy->root) = 1 / (1 + model->eps);
    lk += DCLK(0, model->phy->root) * UCLK(0, model->phy->root);
    UCLK(1, model->phy->root) = model->eps / (1 + model->eps);
    lk += DCLK(1, model->phy->root) * UCLK(1, model->phy->root);

    return log(lk) + model->lzd[model->phy->root->index] * log(minlikelihood);
}


static double model_treeloglk(double *par, struct model *model)
{
    model->eps = par[0];
    model->tau = par[1];
    model_downpass(model);
    return root_loglk(model);
}


static void model_free(struct model *model)
{
    if (model) {
        free(model->mem);
        free(model);
    }
}


static struct model *model_init(double *dclk, struct phy *phy)
{
    size_t nbytes;
    struct model *model;

    model = malloc(sizeof(struct model));

    if (!model)
        error("unable to allocate memory for Mk model");

    nbytes = 2 * phy->nnode * sizeof(int) + 3 * 2 * phy->nnode * sizeof(double);

    model->mem = malloc(nbytes);
    memset(model->mem, 0, nbytes);

    if (!model->mem) {
        error("unable to allocate memory for Mk model");
        free(model);
    }

    model->eps = 1;
    model->tau = 0;
    model->lzd = (int *)(model->mem);
    model->lzu = model->lzd + phy->nnode;
    model->dclk = (double *)(model->lzu + phy->nnode);
    model->sclk = model->dclk + 2 * phy->nnode;
    model->uclk = model->sclk + 2 * phy->nnode;
    model->loglk = 0;
    model->phy = phy;

    memcpy(model->dclk, dclk, 2 * phy->ntip * sizeof(double));

    return model;
}


void mk2_model_free(SEXP model)
{
    model_free((struct model *)R_ExternalPtrAddr(model));
    R_ClearExternalPtr(model);
}


SEXP mk2_model_init(SEXP clk, SEXP tree)
{
    SEXP exptr;
    struct model *model;

    model = model_init(REAL(clk), (struct phy *)R_ExternalPtrAddr(tree));

    exptr = PROTECT(R_MakeExternalPtr(model, R_NilValue, R_NilValue));
    R_RegisterCFinalizer(exptr, &mk2_model_free);

    UNPROTECT(1);
    return exptr;
}


SEXP mk2_loglk(SEXP par, SEXP model)
{
    struct model *m = (struct model *)R_ExternalPtrAddr(model);

    if (LENGTH(par) != 2)
        error("Expected 2 parameters, received %d", LENGTH(par));

    return ScalarReal(model_treeloglk(REAL(par), m));
}


static double grad1(double tau, struct node *node, struct model *model)
{
    int i;
    int j;
    double d;
    double g;
    double lk;
    struct node *anc;
    struct node *sib;

    double p01;
    double p10;
    double p00;
    double p11;

    double D = exp(-tau * node->brlen);

    p01 = (1 - D) / 2;
    p10 = (1 - D) / 2;
    p00 = 0.5 + D * 0.5;
    p11 = 0.5 + D * 0.5;

    anc = node->anc;
    sib = anc->lfdesc;

    if (sib == node)
        sib = sib->next;

    d = (node->brlen * D) / 2;
    g = UCLK(0, anc) * SCLK(0, sib) * -d * DCLK(0, node);
    g += UCLK(0, anc) * SCLK(0, sib) * d * DCLK(1, node);
    g += UCLK(1, anc) * SCLK(1, sib) * -d * DCLK(1, node);
    g += UCLK(1, anc) * SCLK(1, sib) * d * DCLK(0, node);

    lk = UCLK(0, anc) * SCLK(0, sib) * p00 * DCLK(0, node);
    lk += UCLK(0, anc) * SCLK(0, sib) * p01 * DCLK(1, node);
    lk += UCLK(1, anc) * SCLK(1, sib) * p11 * DCLK(1, node);
    lk += UCLK(1, anc) * SCLK(1, sib) * p10 * DCLK(0, node);

    return g / lk;
}


static void grad2(double eps, double tau, double *g,
    struct node *node, struct model *model)
{
    int i;
    int j;
    double d;
    double lk;
    double D;
    struct node *anc;
    struct node *sib;

    double p01;
    double p10;
    double p00;
    double p11;

    if (!node->anc)
    {
        // handle root separately
        d = 1 / ((1+eps)*(1+eps));
        g[0] = -d * DCLK(0, node) + d * DCLK(1, node);
        g[0] /= (1/(1+eps))*DCLK(0, node) + (eps/(1+eps))*DCLK(1, node);
        g[1] = 0;
        return;
    }

    D = exp(-tau * node->brlen);

    p01 = (1 - D) * (eps) / (1+eps);
    p10 = (1 - D) / (1 + eps);
    p00 = (1 / (1 + eps)) + D * (eps / (1 + eps));
    p11 = (eps / (1 + eps)) + D / (1 + eps);

    memset(g, 0, 2 * sizeof(double));

    anc = node->anc;
    sib = anc->lfdesc;

    if (sib == node)
        sib = sib->next;

    // g[0] deriv wrt eps
    // g[1] deriv wrt tau

    d = (1 - D) / ((1 + eps) * (1 + eps));
    g[0] += UCLK(0, anc) * SCLK(0, sib) * -d * DCLK(0, node);
    g[0] += UCLK(1, anc) * SCLK(1, sib) * d * DCLK(1, node);
    g[0] += UCLK(0, anc) * SCLK(0, sib) * d * DCLK(1, node);
    g[0] += UCLK(1, anc) * SCLK(1, sib) * -d * DCLK(0, node);

    d = (node->brlen * eps * D) / (1 + eps);
    g[1] += UCLK(0, anc) * SCLK(0, sib) * -d * DCLK(0, node);
    g[1] += UCLK(0, anc) * SCLK(0, sib) * d * DCLK(1, node);
    d = (node->brlen * D) / (1 + eps);
    g[1] += UCLK(1, anc) * SCLK(1, sib) * -d * DCLK(1, node);
    g[1] += UCLK(1, anc) * SCLK(1, sib) * d * DCLK(0, node);

    lk = UCLK(0, anc) * SCLK(0, sib) * p00 * DCLK(0, node);
    lk += UCLK(1, anc) * SCLK(1, sib) * p11 * DCLK(1, node);
    lk += UCLK(0, anc) * SCLK(0, sib) * p01 * DCLK(1, node);
    lk += UCLK(1, anc) * SCLK(1, sib) * p10 * DCLK(0, node);

    g[0] /= lk;
    g[1] /= lk;
}


static double hessian1(double tau, struct node *node, struct model *model)
{
    int i;
    int j;
    double d;
    double g;
    double A;
    struct node *anc;
    struct node *sib;

    anc = node->anc;
    sib = anc->lfdesc;

    if (sib == node)
        sib = sib->next;

    A = exp(-tau * node->brlen);

    d = (node->brlen * node->brlen * A) / 2;
    g = UCLK(0, anc) * SCLK(0, sib) * d * DCLK(0, node);
    g += UCLK(0, anc) * SCLK(0, sib) * -d * DCLK(1, node);
    g += UCLK(1, anc) * SCLK(1, sib) * d * DCLK(1, node);
    g += UCLK(1, anc) * SCLK(1, sib) * -d * DCLK(0, node);

    return g;
}


static void hessian2(double eps, double tau, double *g,
    struct node *node, struct model *model)
{
    int i;
    int j;
    double d;
    double A;
    struct node *anc;
    struct node *sib;

    memset(g, 0, 4 * sizeof(double));

    anc = node->anc;
    sib = anc->lfdesc;

    if (sib == node)
        sib = sib->next;

    // g[0] second deriv wrt eps
    // g[1] and g[2] deriv wrt to eps and tau
    // g[3] second deriv wrt tau

    A = exp(-tau * node->brlen);

    d = 2 * (1 - A) / ((1 + eps) * (1 + eps) * (1 + eps));
    g[0] += UCLK(0, anc) * SCLK(0, sib) * d * DCLK(0, node);
    g[0] += UCLK(1, anc) * SCLK(1, sib) * -d * DCLK(1, node);
    g[0] += UCLK(0, anc) * SCLK(0, sib) * -d * DCLK(1, node);
    g[0] += UCLK(1, anc) * SCLK(1, sib) * d * DCLK(0, node);

    d = (node->brlen * A) / ((1 + eps) * (1 + eps));
    g[1] += UCLK(0, anc) * SCLK(0, sib) * -d * DCLK(0, node);
    g[1] += UCLK(1, anc) * SCLK(1, sib) * d * DCLK(1, node);
    g[1] += UCLK(0, anc) * SCLK(0, sib) * d * DCLK(1, node);
    g[1] += UCLK(1, anc) * SCLK(1, sib) * -d * DCLK(0, node);
    g[2] = g[1];

    d = (node->brlen * node->brlen * eps * A) / (1 + eps);
    g[3] += UCLK(0, anc) * SCLK(0, sib) * d * DCLK(0, node);
    g[3] += UCLK(0, anc) * SCLK(0, sib) * -d * DCLK(1, node);
    d = (node->brlen * node->brlen * A) / (1 + eps);
    g[3] += UCLK(1, anc) * SCLK(1, sib) * d * DCLK(1, node);
    g[3] += UCLK(1, anc) * SCLK(1, sib) * -d * DCLK(0, node);
}


SEXP mk2_grad1(SEXP par, SEXP m)
{
    SEXP grad;
    double *g;
    struct node *node;
    struct model *model = (struct model *)R_ExternalPtrAddr(m);
    model_treeloglk(REAL(par), model);
    model_uppass(model);

    grad = PROTECT(allocVector(REALSXP, model->phy->nnode));
    g = REAL(grad);

    memset(g, 0, model->phy->nnode * sizeof(double));

    phy_traverse_prepare(model->phy, model->phy->root, ALL_NODES, PREORDER);
    phy_traverse_step(model->phy);

    while ((node = phy_traverse_step(model->phy)) != 0)
    {
        g[node->index] = grad1(model->tau, node, model);
    }


    UNPROTECT(1);
    return grad;
}


SEXP mk2_grad2(SEXP par, SEXP m)
{
    SEXP grad;
    double *g;
    struct node *node;
    struct model *model = (struct model *)R_ExternalPtrAddr(m);
    model_treeloglk(REAL(par), model);
    model_uppass(model);

    grad = PROTECT(allocMatrix(REALSXP, 2, model->phy->nnode));
    g = REAL(grad);

    memset(g, 0, model->phy->nnode * sizeof(double));

    phy_traverse_prepare(model->phy, model->phy->root, ALL_NODES, PREORDER);

    while ((node = phy_traverse_step(model->phy)) != 0)
    {
        grad2(model->eps, model->tau, g + 2 * node->index, node, model);
    }


    UNPROTECT(1);
    return grad;
}


SEXP mk2_uclk(SEXP par, SEXP m)
{
    int i;
    SEXP p;
    struct node *node;
    struct model *model = (struct model *)R_ExternalPtrAddr(m);
    model_treeloglk(REAL(par), model);
    model_uppass(model);

    p = PROTECT(allocMatrix(REALSXP, 2, model->phy->ntip));

    for (i = 0; i < model->phy->ntip; ++i)
    {
        node = phy_getnode_with_index(model->phy, i);
        REAL(p)[0 + i * 2] = log(UCLK(0, node)) + model->lzu[i] * log(minlikelihood);
        REAL(p)[1 + i * 2] = log(UCLK(1, node)) + model->lzu[i] * log(minlikelihood);
    }

    UNPROTECT(1);
    return p;
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

    double p01;
    double p10;
    double p00;
    double p11;

    double eps = model->eps;
    double tau = model->tau;
    double D = exp(-tau * node->brlen);

    p01 = (1 - D) * (eps) / (1+eps);
    p10 = (1 - D) / (1 + eps);
    p00 = (1 / (1 + eps)) + D * (eps / (1 + eps));
    p11 = (eps / (1 + eps)) + D / (1 + eps);

    anc = node->anc;

    if (anc) {
        sib = anc->lfdesc;

        if (sib == node)
            sib = sib->next;

        if (j == 0)
        {
            lk = UCLK(0, anc) * SCLK(0, sib) * p00 * init[j];
            lk += UCLK(1, anc) * SCLK(1, sib) * p10 * init[j];
        }
        else
        {
            lk = UCLK(0, anc) * SCLK(0, sib) * p01 * init[j];
            lk += UCLK(1, anc) * SCLK(1, sib) * p11 * init[j];
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


SEXP mk2_marginal_asr(SEXP par, SEXP model)
{
    int i;
    int ntip;
    double *asr;

    struct model *m;
    struct node *node;

    SEXP ASR;

    m = (struct model *)R_ExternalPtrAddr(model);
    ntip = m->phy->ntip;

    if (LENGTH(par) != 2)
        error("Expected 2 parameters, received %d", LENGTH(par));

    model_treeloglk(REAL(par), m);
    model_uppass(m);

    ASR = PROTECT(allocMatrix(REALSXP, 2, ntip - 1));
    asr = REAL(ASR);

    phy_traverse_prepare(m->phy, m->phy->root, INTERNAL_NODES_ONLY, POSTORDER);

    while ((node = phy_traverse_step(m->phy)) != 0) {

        for (i = 0; i < 2; ++i)
            asr[i + 2 * (node->index - ntip)] = asr_compute(i, node, m);

        asr_normalize(2, asr + 2 * (node->index - ntip));
    }

    UNPROTECT(1);
    return ASR;
}



SEXP mk2_perr(SEXP par, SEXP m)
{
    int i;
    double *e;
    double mx;
    double loglk;
    double err = 0;
    SEXP p;
    struct node *node;
    struct model *model = (struct model *)R_ExternalPtrAddr(m);
    loglk = model_treeloglk(REAL(par), model);
    model_uppass(model);

    p = PROTECT(allocMatrix(REALSXP, 2, model->phy->ntip));
    e = REAL(p);

    for (i = 0; i < model->phy->ntip; ++i)
    {
        node = phy_getnode_with_index(model->phy, i);
        e[0 + i * 2] = log(UCLK(0, node)) + model->lzu[i] * log(minlikelihood);
        e[1 + i * 2] = log(UCLK(1, node)) + model->lzu[i] * log(minlikelihood);

        mx = e[0 + i * 2] > e[1 + i * 2] ? e[0 + i * 2] : e[1 + i * 2];

        e[0 + i * 2] -= mx;
        e[1 + i * 2] -= mx;

        mx = e[0 + i * 2] = exp(e[0 + i * 2]);
        mx += e[1 + i * 2] = exp(e[1 + i * 2]);

        e[0 + i * 2] /= mx;
        e[1 + i * 2] /= mx;

        if (((DCLK(0, node) - 1) < 1e-9) && ((DCLK(0, node) - 1) > -1e-9))
            err += e[1 + i * 2];
        else
            err += e[0 + i * 2];
    }

    UNPROTECT(1);

    return ScalarReal(loglk - err);
}


SEXP mk2_mcmc_slice(SEXP n, SEXP p, SEXP win, SEXP m)
{
    int i;
    int j;
    int niter;
    double L;
    double R;
    double Lbar;
    double Rbar;
    double x;
    double x1;
    double u;
    double w;
    double z;
    double loglk;
    double *par;
    double *res;
    struct model *model = (struct model *)R_ExternalPtrAddr(m);
    SEXP result;

    niter = INTEGER(n)[0];
    w = REAL(win)[0];
    par = REAL(p);

    result = PROTECT(allocMatrix(REALSXP, niter, 3));
    res = REAL(result);

    loglk = model_treeloglk(par, model);

    GetRNGstate();

    for (i = 0; i < niter; ++i)
    {
        for (j = 0; j < 2; ++j)
        {
            x = par[j];

            // step 1
            z = loglk - rexp(1);

            // step 2
            u = unif_rand();
            L = x - w * u;
            R = L + w;

            par[j] = L;
            while (L > 0 && model_treeloglk(par, model) > z) {
                L -= w;
                par[j] = L;
            }

            par[j] = R;
            while (R < 10 && model_treeloglk(par, model) > z) {
                R += w;
                par[j] = R;
            }

            L = (L < 0) ? 1e-8 : L;
            R = (R > 10) ? 10 : R;

            // step 3
            Lbar = L;
            Rbar = R;

            do {
                u = unif_rand();
                x1 = Lbar + u * (Rbar - Lbar);

                if (x1 < x)
                    Lbar = x1;
                else
                    Rbar = x1;

                par[j] = x1;
                loglk = model_treeloglk(par, model);
            } while (loglk < z);
        }
        res[i + 0 * niter] = loglk;
        res[i + 1 * niter] = par[0];
        res[i + 2 * niter] = par[1];
    }

    PutRNGstate();

    UNPROTECT(1);
    return result;
}


static double pij(int i, int j, double *par, double brlen)
{
    int xtype = i + j;
    double eps = par[0];
    double tau = par[1];
    double D = exp(-tau * brlen);

    switch (xtype)
    {
        case 0:
            return 1 / (1 + eps) + (eps / (1 + eps)) * D;
        case 2:
            return (eps / (1 + eps)) + (1 / (1 + eps)) * D;
        case 1:
            return (i == 1) ? (1 - D) / (1 + eps) : eps * (1 - D) / (1 + eps);
    }

    return 0;
}


SEXP mk2_loglk_conditional(SEXP p, SEXP state, SEXP rtree)
{
    int pstate;
    int cstate;
    double loglk = 0;
    int *nodestate = INTEGER(state);
    double *par = REAL(p);
    struct phy *phy = (struct phy *)R_ExternalPtrAddr(rtree);

    struct node *d;
    struct node *node;

    phy_traverse_prepare(phy, phy->root, INTERNAL_NODES_ONLY,
        POSTORDER);

    while ((node = phy_traverse_step(phy)) != 0)
    {
        pstate = nodestate[node->index];
        for (d = node->lfdesc; d != 0; d = d->next)
        {
            cstate = nodestate[d->index];
            loglk += log(pij(pstate, cstate, par, d->brlen));
        }
    }

    if (nodestate[phy->root->index])
        loglk += log(par[0] / (1+par[0]));
    else
        loglk += log(1 / (1 + par[0]));

    return ScalarReal(loglk);
}


static void pij_grad(int i, int j, double *par, double *grad, double brlen)
{
    int xtype = i + j;
    double eps = par[0];
    double tau = par[1];
    double D = exp(-tau * brlen);

    switch (xtype)
    {
        case 0:
            grad[0] -= (1 - D) / ((1 + eps) * (1 + eps * D));
            grad[1] -= (brlen * eps * D) / (1 + eps * D);
            break;
        case 2:
            grad[0] += (1 - D) / ((1 + eps) * (eps + D));
            grad[1] -= (brlen * D) / (eps + D);
            break;
        case 1:
            if (i == 1)
                grad[0] -= 1 / (1 + eps);
            else
                grad[0] += 1 / ((1+eps)*eps);
            grad[1] += (brlen * D) / (1 - D);
        break;
    }
}


SEXP mk2_gradient_conditional(SEXP p, SEXP state, SEXP rtree)
{
    int pstate;
    int cstate;
    int *nodestate = INTEGER(state);
    double *par = REAL(p);
    struct phy *phy = (struct phy *)R_ExternalPtrAddr(rtree);

    struct node *d;
    struct node *node;

    SEXP grad = PROTECT(allocVector(REALSXP, 2));
    double *g = REAL(grad);
    g[0] = g[1] = 0;

    phy_traverse_prepare(phy, phy->root, INTERNAL_NODES_ONLY,
        POSTORDER);

    while ((node = phy_traverse_step(phy)) != 0)
    {
        pstate = nodestate[node->index];
        for (d = node->lfdesc; d != 0; d = d->next)
        {
            cstate = nodestate[d->index];
            pij_grad(pstate, cstate, par, g, node->brlen);
        }
    }

    if (nodestate[phy->root->index])
        g[0] += (1+par[0]) / par[0];
    else
        g[0] -= 1 / (1+par[0]);

    UNPROTECT(1);
    return grad;
}

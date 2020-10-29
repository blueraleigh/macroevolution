#include <R.h>
#include <Rinternals.h>
#include "phy.h"

/* Calculate node and stem state costs. */
static void sankoff_downpass(struct phy *phy, int r, double *g, double *h,
    double *cost)
{
    int i;
    int j;
    int k;
    double t;
    double min_t;
    struct node *node;
    struct node *desc;

#define C(i,j) cost[(i) + r * (j)]
#define G(i,j) g[(i) + phy->nnode * (j)]
#define H(i,j) h[(i) + phy->nnode * (j)]

    for (k = 0; k < phy->ntip; ++k) {
        node = phy_getnode_with_index(phy, k);
        for (i = 0; i < r; ++i) {
            min_t = R_PosInf;
            for (j = 0; j < r; ++j) {
                t = C(i, j) + G(node->index, j);
                if (t < min_t)
                    min_t = t;
            }
            H(node->index, i) = min_t;
        }
    }

    phy_traverse_prepare(phy, phy->root, INTERNAL_NODES_ONLY, POSTORDER);

    while ((node = phy_traverse_step(phy)) != 0)
    {
        for (j = 0; j < r; ++j) {
            desc = node->lfdesc;
            G(node->index, j) = H(desc->index, j);
            for (desc = desc->next; desc; desc = desc->next)
                G(node->index, j) += H(desc->index, j);
        }

        for (i = 0; i < r; ++i) {
            min_t = R_PosInf;
            for (j = 0; j < r; ++j) {
                t = C(i, j) + G(node->index, j);
                if (t < min_t)
                    min_t = t;
            }
            H(node->index, i) = min_t;
        }
    }
}


/* Calculate final state costs */
static void sankoff_uppass(struct phy *phy, int r, double *g, double *h,
    double *f, double *cost)
{
    int i;
    int j;
    int focal;
    int parent;
    double t;
    double min_t;
    struct node *node;

#define F(i,j) f[(i) + phy->nnode * (j)]

    for (j = 0; j < r; ++j)
        F(phy->root->index, j) = G(phy->root->index, j);

    phy_traverse_prepare(phy, phy->root, ALL_NODES, PREORDER);
    phy_traverse_step(phy);
    while ((node = phy_traverse_step(phy)) != 0)
    {
        focal = node->index;
        parent = node->anc->index;

        for (i = 0; i < r; ++i) {
            min_t = R_PosInf;
            for (j = 0; j < r; ++j) {
                t = (F(parent, j) - H(focal, j)) + C(j, i) + G(focal, i);
                if (t < min_t)
                    min_t = t;
            }
            F(focal, i) = min_t;
        }
    }
}


/* Count the number of MPR histories */
static double sankoff_count(struct phy *phy, int r, double *g, double *h,
    double *f, double *cost)
{
    int i;
    int j;
    int k;
    double np;
    double L;
    double nbr = 0;
    double nh[phy->nnode * r];
    struct node *node;
    struct node *child;

    memset(nh, 0, (phy->nnode * r) * sizeof(double));

#define NH(i, j) nh[(i) + (j) * phy->nnode]

    for (i = 0; i < phy->ntip; ++i) {
        node = phy_getnode_with_index(phy, i);
        for (j = 0; j < r; ++j) {
            if (G(i, j) < R_PosInf)
                NH(i, j) = 1;
            else
                NH(i, j) = 0;
        }
    }

    L = R_PosInf;
    for (j = 0; j < r; ++j) {
        if (F(phy->root->index, j) < L)
            L = F(phy->root->index, j);
    }

    phy_traverse_prepare(phy, phy->root, INTERNAL_NODES_ONLY, POSTORDER);

    while ((node = phy_traverse_step(phy)) != 0) {
        for (j = 0; j < r; ++j) {
            if (F(node->index, j) > L)
                continue;
            NH(node->index, j) = 1;
        }
        for (child = node->lfdesc; child; child = child->next) {
            for (j = 0; j < r; ++j) {
                if (F(node->index, j) > L)
                    continue;
                np = 0;
                for (k = 0; k < r; ++k) {
                    if ((F(node->index, j) - H(child->index, j)
                            + C(j, k) + G(child->index, k)) > L)
                        continue;
                    np += NH(child->index, k);
                }
                NH(node->index, j) *= np;
            }
        }
    }

    for (j = 0; j < r; ++j)
        nbr += NH(phy->root->index, j);

    return nbr;
}


/* Calculate the set of MPR state-to-state transitions for each branch. */
static void sankoff_transitions(struct phy *phy, int r, double *g, double *h,
    double *f, double *cost)
{
    int i;
    int j;
    int k;
    double L;
    struct node *node;
    struct node *child;

    L = R_PosInf;
    for (j = 0; j < r; ++j) {
        if (F(phy->root->index, j) < L)
            L = F(phy->root->index, j);
    }

    phy_traverse_prepare(phy, phy->root, INTERNAL_NODES_ONLY, POSTORDER);

    while ((node = phy_traverse_step(phy)) != 0) {
        for (child = node->lfdesc; child; child = child->next) {
            for (j = 0; j < r; ++j) {
                if (F(node->index, j) > L)
                    continue;
                for (k = 0; k < r; ++k) {
                    if ((F(node->index, j) - H(child->index, j)
                            + C(j, k) + G(child->index, k)) > L)
                        continue;
                    // j -> k is a transition for this branch in at least
                    // one MPR history
                }
            }
        }
    }
}


static int sankoff_choose(int r, int stride, double L, double *f, int only_mpr)
{
    int i;
    int j = 0;
    double g;
    double w;
    double max_w = R_NegInf;
    for (i = 0; i < r; ++i, f += stride) {
        if (*f > L && only_mpr)
            continue;
        w = L - *f;
        g = -log(-log(unif_rand()));
        if ((g + w) > max_w) {
            max_w = g + w;
            j = i;
        }
    }
    return j;
}


/* Sample histories of character evolution */
static void sankoff_sample(struct phy *phy, int r, double *g, double *h,
    double *f, double *cost, int nsample, int *nodestate, int only_mpr)
{
    int i;
    int j;
    int k;
    double L;
    double w[r];
    struct node *node;

    L = R_PosInf;
    for (j = 0; j < r; ++j) {
        if (F(phy->root->index, j) < L)
            L = F(phy->root->index, j);
    }

#define NS(i,j) nodestate[(i)+(j)*phy->nnode]

    for (k = 0; k < nsample; ++k)
        NS(phy->root->index, k) = sankoff_choose(r, phy->nnode, L, f, only_mpr);

    phy_traverse_prepare(phy, phy->root, ALL_NODES, PREORDER);
    phy_traverse_step(phy);

    while ((node = phy_traverse_step(phy)) != 0) {
        for (k = 0; k < nsample; ++k) {
            i = NS(node->anc->index, k);
            for (j = 0; j < r; ++j)
                w[j] = F(node->anc->index, i) - H(node->index, i)
                            + C(i, j) + G(node->index, j);
            NS(node->index, k) = sankoff_choose(
                r, 1, F(node->anc->index, i), w, only_mpr);
        }
    }
}


SEXP do_sankoff_downpass(SEXP rtree, SEXP r, SEXP g, SEXP h, SEXP cost)
{
    sankoff_downpass((struct phy *)R_ExternalPtrAddr(rtree), INTEGER(r)[0],
        REAL(g), REAL(h), REAL(cost));
    return R_NilValue;
}


SEXP do_sankoff_uppass(SEXP rtree, SEXP r, SEXP g, SEXP h, SEXP f, SEXP cost)
{
    sankoff_uppass((struct phy *)R_ExternalPtrAddr(rtree), INTEGER(r)[0],
        REAL(g), REAL(h), REAL(f), REAL(cost));
    return R_NilValue;
}


SEXP do_sankoff_count(SEXP rtree, SEXP r, SEXP g, SEXP h, SEXP f, SEXP cost)
{
    return ScalarReal(sankoff_count((struct phy *)R_ExternalPtrAddr(rtree),
        INTEGER(r)[0], REAL(g), REAL(h), REAL(f), REAL(cost)));
}


SEXP do_sankoff_sample(SEXP rtree, SEXP r, SEXP g, SEXP h, SEXP f, SEXP cost,
    SEXP nsample, SEXP nodestate, SEXP only_mpr)
{
    GetRNGstate();
    sankoff_sample((struct phy *)R_ExternalPtrAddr(rtree), INTEGER(r)[0],
        REAL(g), REAL(h), REAL(f), REAL(cost), INTEGER(nsample)[0],
        INTEGER(nodestate), INTEGER(only_mpr)[0]);
    PutRNGstate();
    return R_NilValue;
}

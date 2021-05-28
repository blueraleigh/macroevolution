#include <R.h>
#include <Rinternals.h>
#include "phy.h"

/*
** Fitch parsimony routines. We represent state sets
** as bit masks using 32 bit signed integers. For example
** the state set {1, 3, 7} is represented as
**
** 00000000000000000000000001000101
**
** corresponding to the integer value 138.
**
** Consequently, these routines can only handle a maximum
** of 31 states.
**
** States are numbered starting from 1 on the R side but
** starting from 0 on the C side.
*/


static int stateset_contains(int state, int stateset) {
    return (stateset & (1<<state)) != 0 ? 1 : 0;
}


static int stateset_add(int state, int stateset) {
    return stateset | 1<<state;
}


static int stateset_remove(int state, int stateset) {
    return stateset & ~(1<<state);
}


/*
** Perform the Fitch parsimony downpass algorithm
**
** @param ddata binary vector to hold the downpass state sets
**  at each internal node. The i-th position will correspond
**  to the downpass state set for the node with i-th index.
**  i.e. This is a r x phy->nnode matrix stored in column major
**  as a flat array. This same convention applies to all other
**  variables in this file with names like udata, dcost, ucost.
** @param pscores vector of parsimony scores. Each value is
**  the number of parsimony changes occurring in the subclade
**  rooted at the node with the corresponding index number.
** @param r the number of character states
*/
static void fitch_downpass(struct phy *phy, int *g, int *pscores, int r)
{
    int parent;
    int lfchild;
    int rtchild;
    int pscore;
    struct node *node;

    memset(pscores, 0, phy->nnode * sizeof(int));

    phy_traverse_prepare(phy, phy->root, INTERNAL_NODES_ONLY, POSTORDER);

    while ((node = phy_traverse_step(phy)) != 0)
    {
        pscore = 0;
        parent = node->index;
        lfchild = node->lfdesc->index;
        rtchild = node->lfdesc->next->index;

        g[parent] = g[lfchild] & g[rtchild];

        if (!g[parent]) {
            g[parent] = g[lfchild] | g[rtchild];
            ++pscore;
        }

        pscores[parent] += pscores[lfchild];
        pscores[parent] += pscores[rtchild];
        pscores[parent] += pscore;
    }
}


/*
** Perform the Fitch uppass algorithm
**
** @param udata binary vector to hold the uppass state sets
**  at each internal node. The i-th position will correspond
**  to the uppass state set for the node with i-th index
** @param ddata the downpass state sets resulting from the
**  Fitch downpass algorithm
*/
static void fitch_uppass(struct phy *phy, int *f, int *g, int r)
{
    int focal;
    int parent;
    int lfchild;
    int rtchild;
    struct node *node;

    f[phy->root->index] = g[phy->root->index];

    phy_traverse_prepare(phy, phy->root, INTERNAL_NODES_ONLY, PREORDER);
    phy_traverse_step(phy);

    while ((node = phy_traverse_step(phy)) != 0)
    {
        focal = node->index;
        parent = node->anc->index;
        lfchild = node->lfdesc->index;
        rtchild = node->lfdesc->next->index;

        if (f[parent] == (g[focal] & f[parent])) {
            f[focal] = f[parent];
        } else if ((g[lfchild] & g[rtchild]) == 0) {
            f[focal] = g[focal] | f[parent];
        } else {
            f[focal] = g[lfchild] | g[rtchild];
            f[focal] &= f[parent];
            f[focal] |= g[focal];
        }
    }
}


SEXP do_fitch_mpr2(SEXP rtree, SEXP r, SEXP up, SEXP down, SEXP pscore)
{
    struct phy *phy = (struct phy *)R_ExternalPtrAddr(rtree);
    int *ddata = INTEGER(down);
    int *udata = INTEGER(up);
    fitch_downpass(phy, ddata, INTEGER(pscore), *INTEGER(r));
    fitch_uppass(phy, udata, ddata, *INTEGER(r));

    SEXP ans = PROTECT(allocMatrix(INTSXP, phy->nnode, *INTEGER(r)));

    for (int i = 0; i < phy->nnode; ++i) {
        for (int j = 0; j < *INTEGER(r); ++j) {
            if (stateset_contains(j, udata[i]))
                INTEGER(ans)[i + j*phy->nnode] = 1;
            else
                INTEGER(ans)[i + j*phy->nnode] = 0;
        }
    }

    UNPROTECT(1);
    return ans;
}



/*
** Uniformly sample a 0-valued index from an r-length
** binary vector (reservoir sampling)
*/
static int choose_state(int r, int stateset)
{
    int i;
    int j;
    double k = 1;
    for (i = 0; i < r; ++i) {
        if (stateset_contains(i, stateset)) {
            if (unif_rand() < 1/k)
                j = i;
            k += 1;
        }
    }
    return j;
}


static void fitch_history(struct phy *phy, int *nodestate, int *f,
    int *g, int r, int nsample)
{
    int k;
    int n = phy->nnode - phy->ntip;
    int focal;
    int parent;
    int cstate;
    int pstate;

    struct node *node;

    for (k = 0; k < nsample; ++k)
    {
        nodestate[phy->root->index + k * phy->nnode] =
            choose_state(r, f[phy->root->index]) + 1;

        phy_traverse_prepare(phy, phy->root, INTERNAL_NODES_ONLY, PREORDER);
        phy_traverse_step(phy);

        while ((node = phy_traverse_step(phy)) != NULL)
        {
            focal = node->index;
            parent = node->anc->index;
            pstate = nodestate[parent + k * phy->nnode] - 1;

            if (stateset_contains(pstate, g[focal]))
            {
                cstate = pstate;
            }
            else
            {
                if (stateset_contains(pstate, f[focal]))
                {
                    // temporarily modify the downpass state set
                    // to include the parent state, which is a
                    // valid MPR state-to-state transition in this case.
                    g[focal] = stateset_add(pstate, g[focal]);
                    cstate = choose_state(r, g[focal]);
                    g[focal] = stateset_remove(pstate, g[focal]);
                }
                else
                {
                    cstate = choose_state(r, g[focal]);
                }
            }
            nodestate[focal + k*phy->nnode] = cstate + 1;
        }
    }
}


/*
** Count the number of MPR reconstructions
**
** Algorithm modified from Mesquite project
** https://github.com/MesquiteProject/MesquiteCore/blob/master/Source/mesquite/parsimony/lib/MPRProcessor.java
*/
static double fitch_count(struct phy *phy, int *f, int *g, int r)
{
    int i;
    int j;
    int k;
    int parent;
    int child;
    double np;
    double nbr = 0;
    double nh[phy->nnode * r];
    struct node *d;
    struct node *node;

    memset(nh, 0, (phy->nnode * r) * sizeof(double));

    for (i = 0; i < phy->ntip; ++i) {
        for (j = 0; j < r; ++j) {
            if (stateset_contains(j, g[i]))
                nh[i + j * phy->nnode] = 1;
        }
    }

    phy_traverse_prepare(phy, phy->root, INTERNAL_NODES_ONLY, POSTORDER);

    while ((node = phy_traverse_step(phy)) != 0) {
        parent = node->index;
        for (j = 0; j < r; ++j) {
            if (stateset_contains(j, f[parent]))
                nh[parent + j * phy->nnode] = 1;
        }
        for (d = node->lfdesc; d; d = d->next) {
            child = d->index;
            for (j = 0; j < r; ++j) {
                if (stateset_contains(j, f[parent])) {
                    np = 0;
                    if (stateset_contains(j, g[child])) {
                        np += nh[child + j * phy->nnode];
                    } else {
                        for (k = 0; k < r; ++k) {
                            if (stateset_contains(k, g[child]) ||
                                stateset_contains(k, f[child]))
                            {
                                np += nh[child + k * phy->nnode];
                            }
                        }
                    }
                    nh[parent + j * phy->nnode] *= np;
                }
            }
        }
    }

    for (j = 0; j < r; ++j)
        nbr += nh[phy->root->index + j * phy->nnode];

    return nbr;
}


SEXP do_fitch_count3(SEXP rtree, SEXP up, SEXP down, SEXP r)
{
    struct phy *phy = (struct phy *)R_ExternalPtrAddr(rtree);
    int *ddata = INTEGER(down);
    int *udata = INTEGER(up);
    return ScalarReal(fitch_count(phy, udata, ddata, *INTEGER(r)));
}

#include <R.h>
#include <Rinternals.h>
#include "phy.h"


/* (inefficient) Parsimony routines for a single r-state characer */


/*
** Fitch parsimony
*/

/*
** Perform a union and intersection operation on
** two r-length binary vectors
**
** @param r number of character states
** @param a left binary vector
** @param b right binary vector
** @param x vector to hold intersection of a and b
** @param u vector to hold union of a and b
** @return 1 if we can take an intersection, 0 otherwise
*/
static int xoru(int r, int *a, int *b, int *x, int *u)
{
    int i;
    int xtype;
    int n = 0;
    for (i = 0; i < r; ++i)
    {
        xtype = a[i] + b[i];
        switch (xtype)
        {
            case 0:
                x[i] = u[i] = 0;
                break;
            case 1:
                x[i] = 0;  // intersection
                u[i] = 1;  // union
                break;
            case 2:
                x[i] = u[i] = 1;
                n = 1;
                break;
        }
    }
    return n;
}


/*
** Determine if two r-length binary vectors are identical,
** disjoint or somewhere in between
**
** @return 1 if identical, 0 if disjoint, -1 if in between
*/
static int iddj(int r, int *a, int *b)
{
    int i;
    int size_a = 0, size_b = 0, id = 0;
    for (i = 0; i < r; ++i)
    {
        if (a[i] && b[i])
            id += 1;
        if (a[i])
            size_a += 1;
        if (b[i])
            size_b += 1;
    }
    if (id == 0)
        return 0;
    if ((size_a == size_b) && (id == size_a))
        return 1;
    return -1;
}


/*
** Uniformly sample a 1-valued index from an r-length
** binary vector (reservoir sampling)
*/
static int choose_state(int r, int *a)
{
    int i;
    int j;
    double k = 1;
    for (i = 0; i < r; ++i)
    {
        if (a[i])
        {
            if (unif_rand() < 1/k)
                j = i;
            k += 1;
        }
    }
    return j;
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
static void fitch_downpass(struct phy *phy, int *ddata, int *pscores, int r)
{
    int parent;
    int lfchild;
    int rtchild;
    int pscore;
    int x[r];
    int u[r];
    struct node *node;

    memset(pscores, 0, phy->nnode * sizeof(int));

    phy_traverse_prepare(phy, phy->root, INTERNAL_NODES_ONLY, POSTORDER);

    while ((node = phy_traverse_step(phy)) != 0)
    {
        pscore = 0;
        parent = node->index;
        lfchild = node->lfdesc->index;
        rtchild = node->lfdesc->next->index;

        if (xoru(r, ddata + lfchild*r, ddata + rtchild*r, x, u))
        {
            memcpy(ddata + parent*r, x, r * sizeof(int));
        }
        else
        {
            pscore += 1;
            memcpy(ddata + parent*r, u, r * sizeof(int));
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
static void fitch_uppass(struct phy *phy, int *udata, int *ddata, int r)
{
    int focal;
    int parent;
    int lfchild;
    int rtchild;
    int diff;
    int x[r];
    int u[r];
    struct node *node;

    memcpy(udata, ddata, (phy->nnode * r) * sizeof(int));

    phy_traverse_prepare(phy, phy->root, INTERNAL_NODES_ONLY, PREORDER);
    phy_traverse_step(phy);

    while ((node = phy_traverse_step(phy)) != 0)
    {
        focal = node->index;
        parent = node->anc->index;
        lfchild = node->lfdesc->index;
        rtchild = node->lfdesc->next->index;

        xoru(r, ddata + focal*r, udata + parent*r, x, u);

        if (iddj(r, udata + parent*r, x) == 1)
        {
            memcpy(udata + focal*r, x, r * sizeof(int));
        }
        else if (iddj(r, ddata + lfchild*r, ddata + rtchild*r) == 0)
        {
            memcpy(udata + focal*r, u, r * sizeof(int));
        }
        else
        {
            xoru(r, ddata + lfchild*r, ddata + rtchild*r, x, u);
            xoru(r, udata + parent*r, u, x, u);
            xoru(r, ddata + focal*r, x, x, u);
            memcpy(udata + focal*r, u, r * sizeof(int));
        }
    }
}


SEXP do_fitch_pscore(SEXP rtree, SEXP ddata, SEXP pscore)
{
    int r = INTEGER(getAttrib(ddata, R_DimSymbol))[0];
    struct phy *phy = (struct phy *)R_ExternalPtrAddr(rtree);
    fitch_downpass(phy, INTEGER(ddata), INTEGER(pscore), r);
    return R_NilValue;
}


SEXP do_fitch_mpr(SEXP rtree, SEXP up, SEXP down, SEXP pscore)
{
    int r = INTEGER(getAttrib(down, R_DimSymbol))[0];
    struct phy *phy = (struct phy *)R_ExternalPtrAddr(rtree);
    int *ddata = INTEGER(down);
    int *udata = INTEGER(up);
    fitch_downpass(phy, ddata, INTEGER(pscore), r);
    fitch_uppass(phy, udata, ddata, r);
    return R_NilValue;
}


/*
** Count the number of MPR reconstructions
**
** Algorithm modified from Mesquite project
** https://github.com/MesquiteProject/MesquiteCore/blob/master/Source/mesquite/parsimony/lib/MPRProcessor.java
*/
static double fitch_count(struct phy *phy, int *udata, int *ddata, int r)
{
    int i;
    int j;
    int k;
    int parent;
    int child;
    double np;
    double nbr = 0;
    double nh[phy->nnode * r];
    struct node *node;

    memset(nh, 0, (phy->nnode * r) * sizeof(double));

    for (i = 0; i < phy->ntip; ++i)
    {
        for (j = 0; j < r; ++j)
        {
            if (ddata[j + i * r])
                nh[i + j * phy->nnode] = 1;
        }
    }

    phy_traverse_prepare(phy, phy->root, INTERNAL_NODES_ONLY, POSTORDER);

    while ((node = phy_traverse_step(phy)) != 0)
    {
        parent = node->index;
        for (j = 0; j < r; ++j)
        {
            if (udata[j + parent * r])
                nh[parent + j * phy->nnode] = 1;
        }
        child = node->lfdesc->index;
        for (j = 0; j < r; ++j)
        {
            if (udata[j + parent * r])
            {
                np = 0;
                if (ddata[j + child * r])
                {
                    np += nh[child + j * phy->nnode];
                }
                else
                {
                    for (k = 0; k < r; ++k)
                    {
                        if (ddata[k + child * r] || udata[k + child * r])
                            np += nh[child + k * phy->nnode];
                    }
                }
                nh[parent + j * phy->nnode] *= np;
            }
        }
        child = node->lfdesc->next->index;
        for (j = 0; j < r; ++j)
        {
            if (udata[j + parent * r])
            {
                np = 0;
                if (ddata[j + child * r])
                {
                    np += nh[child + j * phy->nnode];
                }
                else
                {
                    for (k = 0; k < r; ++k)
                    {
                        if (ddata[k + child * r] || udata[k + child * r])
                            np += nh[child + k * phy->nnode];
                    }
                }
                nh[parent + j * phy->nnode] *= np;
            }
        }
    }

    for (j = 0; j < r; ++j)
        nbr += nh[phy->root->index + j * phy->nnode];

    return nbr;
}


/*
** Sample histories of character evolution from a maximum parsimony reconstruction
** of ancestral states.
**
** @param nodestate a matrix to store the sampled states at each node in each sample
** @param nsample the number of samples to take
*/
static void fitch_history(struct phy *phy, int *nodestate, int *udata,
    int *ddata, int r, int nsample)
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
            choose_state(r, udata + r * phy->root->index);

        phy_traverse_prepare(phy, phy->root, INTERNAL_NODES_ONLY, PREORDER);
        phy_traverse_step(phy);

        while ((node = phy_traverse_step(phy)) != NULL)
        {
            focal = node->index;
            parent = node->anc->index;
            pstate = nodestate[parent + k * phy->nnode];

            if (ddata[pstate + focal*r])
            {
                cstate = pstate;
            }
            else
            {
                if (udata[pstate + focal*r])
                {
                    // temporarily modify the downpass state set
                    // to include the parent state, which is a
                    // valid MPR state-to-state transition in this case.
                    ddata[pstate + r*focal] = 1;
                    cstate = choose_state(r, ddata + r*focal);
                    ddata[pstate + r*focal] = 0;
                }
                else
                {
                    cstate = choose_state(r, ddata + r*focal);
                }
            }
            nodestate[focal + k*phy->nnode] = cstate;
        }
    }
}


SEXP do_fitch_history(SEXP rtree, SEXP nodestate, SEXP up, SEXP down,
    SEXP nsample)
{
    GetRNGstate();
    int r = INTEGER(getAttrib(down, R_DimSymbol))[0];
    struct phy *phy = (struct phy *)R_ExternalPtrAddr(rtree);
    int *ddata = INTEGER(down);
    int *udata = INTEGER(up);
    int *state = INTEGER(nodestate);
    fitch_history(phy, state, udata, ddata, r, INTEGER(nsample)[0]);
    PutRNGstate();
    return R_NilValue;
}


SEXP do_fitch_count(SEXP rtree, SEXP up, SEXP down)
{
    int r = INTEGER(getAttrib(down, R_DimSymbol))[0];
    struct phy *phy = (struct phy *)R_ExternalPtrAddr(rtree);
    int *ddata = INTEGER(down);
    int *udata = INTEGER(up);
    return ScalarReal(fitch_count(phy, udata, ddata, r));
}


static void fitch_count2(int r, int *count, int *uppass, int *downpass,
    struct phy *phy)
{
    int i;
    int j;
    struct node *node;

    phy_traverse_prepare(phy, phy->root, ALL_NODES, PREORDER);
    phy_traverse_step(phy);

    while ((node = phy_traverse_step(phy)) != 0)
    {
        for (i = 0; i < r; ++i)
        {
            if (uppass[i + r * node->anc->index])
            {
                // state is in uppass set of parent
                if (downpass[i + r * node->index])
                {
                    // state is in downpass set
                    count[i + i * r] += 1;
                }
                else
                {
                    for (j = 0; j < r; ++j)
                    {
                        if (downpass[j + r * node->index])
                            count[i + j * r] += 1;
                    }
                    if (uppass[i + r * node->index])
                        count[i + i * r] += 1;
                }
            }
        }
    }
}


SEXP do_fitch_count2(SEXP down, SEXP up, SEXP pscore, SEXP rtree)
{
    int r;
    struct phy *phy;

    r = INTEGER(getAttrib(down, R_DimSymbol))[0];
    phy = (struct phy *)R_ExternalPtrAddr(rtree);

    SEXP count = PROTECT(allocMatrix(INTSXP, r, r));
    memset(INTEGER(count), 0, r * r * sizeof(int));

    fitch_downpass(phy, INTEGER(down), INTEGER(pscore), r);
    fitch_uppass(phy, INTEGER(up), INTEGER(down), r);
    fitch_count2(r, INTEGER(count), INTEGER(up), INTEGER(down), phy);

    UNPROTECT(1);
    return count;
}


SEXP do_fitch_phat(SEXP nstate, SEXP h, SEXP rtree)
{
    int i;
    int r;
    int k;
    int l;
    int n;
    int nnode;
    int *nodestate;
    double *phat;
    struct phy *phy;
    struct node *node;
    struct node *anc;

    SEXP P;

    phy = (struct phy *)R_ExternalPtrAddr(rtree);
    nnode = phy->nnode;
    r = INTEGER(nstate)[0];
    n = INTEGER(getAttrib(h, R_DimSymbol))[1];
    nodestate = INTEGER(h);

    double rowsum[r];

    P = PROTECT(allocMatrix(REALSXP, r, r));
    phat = REAL(P);
    memset(phat, 0, r * r * sizeof(double));
    memset(rowsum, 0, r * sizeof(double));

    phy_traverse_prepare(phy, phy->root, ALL_NODES, PREORDER);
    phy_traverse_step(phy);

    while ((node = phy_traverse_step(phy)) != 0)
    {
        anc = node->anc;

        for (i = 0; i < n; ++i)
        {
            k = nodestate[anc->index + i * nnode];
            l = nodestate[node->index + i * nnode];
            phat[k + l * r] += 1 / (double)n;
            rowsum[k] += 1 / (double)n;
        }
    }
/*
    for (i = 0; i < r; ++i)
    {
        for (k = 0; k < r; ++k)
            phat[i + k * r] /= rowsum[i];
    }
*/
    UNPROTECT(1);
    return P;
}



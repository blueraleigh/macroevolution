#include <R.h>
#include <Rinternals.h>
#include "phy.h"


SEXP phyr_plot_ctree(SEXP rtree, SEXP ages, SEXP direction)
{
    struct phy *phy = (struct phy *)R_ExternalPtrAddr(rtree);
    int i;
    int nnode = phy->nnode;
    int ntip = phy->ntip;
    int d = INTEGER(direction)[0];
    SEXP coord = PROTECT(allocMatrix(REALSXP, nnode, 4));
    SEXP bar = PROTECT(allocMatrix(REALSXP, (nnode-ntip),  4));
    double *segs = REAL(coord);
    double *bars = REAL(bar);
    double *age = REAL(ages);
    double a, b;
    struct node *desc;
    struct node *node;
    double ypos = ntip;
    double maxage = 0;
    for (i = 0; i < phy->nnode; ++i) {
        if (age[i] > maxage)
            maxage = age[i];
    }

    phy_traverse_prepare(phy, phy->root, ALL_NODES, POSTORDER);

    while ((node = phy_traverse_step(phy)) != 0) {
        i = node->index;
        switch (d) {
            case 0:
                segs[i + 0 * nnode] = age[i];                   // x0
                segs[i + 1 * nnode] = age[i] - node->brlen;     // x1
                if (!node->ndesc) {
                    segs[i + 2 * nnode] = ypos;                 // y0
                    segs[i + 3 * nnode] = ypos--;               // y1
                } else {
                    a = segs[node->lfdesc->index + 2 * nnode];
                    for (desc = node->lfdesc; desc->next != 0; desc = desc->next) {}
                    b = segs[desc->index + 2 * nnode];
                    segs[i + 2 * nnode] = (a + b) / 2;
                    segs[i + 3 * nnode] = (a + b) / 2;
                    bars[(i-ntip) + 0 * (nnode-ntip)] = age[i];
                    bars[(i-ntip) + 1 * (nnode-ntip)] = age[i];
                    bars[(i-ntip) + 2 * (nnode-ntip)] = a;
                    bars[(i-ntip) + 3 * (nnode-ntip)] = b;
                }
                break;
            case 1:
                segs[i + 0 * nnode] = maxage - age[i];
                segs[i + 1 * nnode] = maxage - age[i] + node->brlen;
                if (!node->ndesc) {
                    segs[i + 2 * nnode] = ypos;
                    segs[i + 3 * nnode] = ypos--;
                } else {
                    a = segs[node->lfdesc->index + 2 * nnode];
                    for (desc = node->lfdesc; desc->next != 0; desc = desc->next) {}
                    b = segs[desc->index + 2 * nnode];
                    segs[i + 2 * nnode] = (a + b) / 2;
                    segs[i + 3 * nnode] = (a + b) / 2;
                    bars[(i-ntip) + 0 * (nnode-ntip)] = maxage - age[i];
                    bars[(i-ntip) + 1 * (nnode-ntip)] = maxage - age[i];
                    bars[(i-ntip) + 2 * (nnode-ntip)] = a;
                    bars[(i-ntip) + 3 * (nnode-ntip)] = b;
                }
                break;
            case 2:
                segs[i + 2 * nnode] = age[i];
                segs[i + 3 * nnode] = age[i] - node->brlen;
                if (!node->ndesc) {
                    segs[i + 0 * nnode] = ypos;
                    segs[i + 1 * nnode] = ypos--;
                } else {
                    a = segs[node->lfdesc->index + 0 * nnode];
                    for (desc = node->lfdesc; desc->next != 0; desc = desc->next) {}
                    b = segs[desc->index + 0 * nnode];
                    segs[i + 0 * nnode] = (a + b) / 2;
                    segs[i + 1 * nnode] = (a + b) / 2;
                    bars[(i-ntip) + 2 * (nnode-ntip)] = age[i];
                    bars[(i-ntip) + 3 * (nnode-ntip)] = age[i];
                    bars[(i-ntip) + 0 * (nnode-ntip)] = a;
                    bars[(i-ntip) + 1 * (nnode-ntip)] = b;
                }
                break;
            case 3:
                segs[i + 2 * nnode] = maxage - age[i];
                segs[i + 3 * nnode] = maxage - age[i] + node->brlen;
                if (!node->ndesc) {
                    segs[i + 0 * nnode] = ypos;
                    segs[i + 1 * nnode] = ypos--;
                } else {
                    a = segs[node->lfdesc->index + 0 * nnode];
                    for (desc = node->lfdesc; desc->next != 0; desc = desc->next) {}
                    b = segs[desc->index + 0 * nnode];
                    segs[i + 0 * nnode] = (a + b) / 2;
                    segs[i + 1 * nnode] = (a + b) / 2;
                    bars[(i-ntip) + 2 * (nnode-ntip)] = maxage - age[i];
                    bars[(i-ntip) + 3 * (nnode-ntip)] = maxage - age[i];
                    bars[(i-ntip) + 0 * (nnode-ntip)] = a;
                    bars[(i-ntip) + 1 * (nnode-ntip)] = b;
                }
                break;
        }
    }
    SEXP ret = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ret, 0, coord);
    SET_VECTOR_ELT(ret, 1, bar);
    UNPROTECT(3);
    return ret;
}


SEXP phyr_plot_ptree(SEXP rtree, SEXP step)
{
    struct phy *phy = (struct phy *)R_ExternalPtrAddr(rtree);
    int i;
    int z = 0;
    int nnode = phy->nnode;
    SEXP theta = PROTECT(allocMatrix(REALSXP, nnode, 3));
    double *th = REAL(theta);
    double a, b, vstep = REAL(step)[0];
    struct node *node;
    struct node *desc;

    phy_traverse_prepare(phy, phy->root, ALL_NODES, POSTORDER);

    while ((node = phy_traverse_step(phy)) != 0) {
        i = node->index;
        if (!node->ndesc) {
            th[i + 0 * nnode] = vstep * z++;
            th[i + 1 * nnode] = 0;
            th[i + 2 * nnode] = 0;
        } else {
            a = th[node->lfdesc->index];
            for (desc = node->lfdesc; desc->next != 0; desc = desc->next) {}
            b = th[desc->index];
            th[i + 0 * nnode] = (a + b) / 2;
            th[i + 1 * nnode] = a;
            th[i + 2 * nnode] = b;
        }
    }
    UNPROTECT(1);
    return theta;
}


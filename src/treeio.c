#include <R.h>
#include <Rinternals.h>
#include "phy.h"


void phyr_tree_free(SEXP rtree)
{
    if (TYPEOF(rtree) == NILSXP)
        return;
    struct phy *phy = (struct phy*)R_ExternalPtrAddr(rtree);
    phy_free(phy);
    R_ClearExternalPtr(rtree);
}


SEXP phyr_read_newick(SEXP newick)
{
    SEXP rtree;
    struct phy *phy = phy_read(CHAR(STRING_ELT(newick, 0)));
    if (phy) {
        rtree = PROTECT(R_MakeExternalPtr(phy, R_NilValue, R_NilValue));
        R_RegisterCFinalizer(rtree, &phyr_tree_free);
        setAttrib(rtree, install("root"), ScalarInteger(phy->root->index+1));
        setAttrib(rtree, install("Ntip"), ScalarInteger(phy->ntip));
        setAttrib(rtree, install("Nnode"), ScalarInteger(phy->nnode));
        UNPROTECT(1);
        return rtree;
    }
    error(phy_errmsg());
}


SEXP phyr_write_newick(SEXP rtree)
{
    struct phy *phy = (struct phy*)R_ExternalPtrAddr(rtree);
    char *newick = phy_write(phy);
    SEXP ret = PROTECT(allocVector(STRSXP, 1));
    SET_STRING_ELT(ret, 0, mkChar(newick));
    free(newick);
    UNPROTECT(1);
    return ret;
}


SEXP phyr_tiplabels(SEXP rtree)
{
    int i;
    SEXP tiplabel;
    struct phy *phy = (struct phy*)R_ExternalPtrAddr(rtree);

    tiplabel = PROTECT(allocVector(STRSXP, phy->ntip));

    for (i = 0; i < phy->ntip; ++i)
        SET_STRING_ELT(tiplabel, i, mkChar(phy->nodes[phy->vseq[i]]->lab));

    UNPROTECT(1);
    return tiplabel;
}


SEXP phyr_node_notes(SEXP rtree)
{
    int i;
    SEXP notes;
    struct phy *phy = (struct phy*)R_ExternalPtrAddr(rtree);

    notes = PROTECT(allocVector(STRSXP, phy->nnode));

    for (i = 0; i < phy->nnode; ++i) {
        if (phy->nodes[phy->vseq[i]]->note)
            SET_STRING_ELT(notes, i, mkChar(phy->nodes[phy->vseq[i]]->note));
        else
            SET_STRING_ELT(notes, i, mkChar(""));
    }

    UNPROTECT(1);
    return notes;
}


SEXP phyr_node_brlens(SEXP rtree)
{
    int i;
    SEXP brlen;
    struct phy *phy = (struct phy*)R_ExternalPtrAddr(rtree);

    brlen = PROTECT(allocVector(REALSXP, phy->nnode));

    for (i = 0; i < phy->nnode; ++i)
        REAL(brlen)[i] = phy->nodes[phy->vseq[i]]->brlen;

    UNPROTECT(1);
    return brlen;
}


SEXP phyr_node_ages(SEXP rtree)
{
    int i;
    double *node_age;
    SEXP age;
    struct phy *phy = (struct phy*)R_ExternalPtrAddr(rtree);
    struct node *node;

    age = PROTECT(allocVector(REALSXP, phy->nnode));
    node_age = REAL(age);

    for (i = 0; i < phy->nnode; ++i) {
        node = phy->nodes[phy->vseq[i]];
        node_age[i] = 0.0;
        while (node != NULL) {
            node_age[i] += node->brlen;
            node = node->anc;
        }
    }

    UNPROTECT(1);
    return age;
}


SEXP phyr_ancestors(SEXP rtree, SEXP node)
{
    int i;
    int sz;
    int end = 0;
    int cnt = 0;
    int nanc = 0;
    struct phy *phy = (struct phy*)R_ExternalPtrAddr(rtree);
    struct node *p = phy->nodes[phy->vseq[INTEGER(node)[0]-1]];

    SEXP buf;
    SEXP root = PROTECT(list1(buf = allocVector(VECSXP, 1000)));
    SEXP tail = root;

    while (p != NULL) {
        nanc++;
        SET_VECTOR_ELT(buf, cnt++, ScalarInteger(p->index+1));
        if (cnt == 1000) {
            tail = SETCDR(tail, list1(buf = allocVector(VECSXP, 1000)));
            cnt = 0;
        }
        p = p->anc;
    }

    SEXP ret = PROTECT(allocVector(INTSXP, nanc));

    while (root != R_NilValue) {
        sz = CDR(root) == R_NilValue ? cnt : 1000;
        for (i = 0; i < sz; ++i)
            INTEGER(ret)[end++] = INTEGER(VECTOR_ELT(CAR(root), i))[0];
        root = CDR(root);
    }

    UNPROTECT(2);
    return ret;
}


SEXP phyr_children(SEXP rtree, SEXP node)
{
    int i = 0;
    struct phy *phy = (struct phy*)R_ExternalPtrAddr(rtree);
    struct node *p = phy->nodes[phy->vseq[INTEGER(node)[0]-1]];
    SEXP ret = PROTECT(allocVector(INTSXP, p->ndesc));
    for (p = p->lfdesc; p != 0; p = p->next)
        INTEGER(ret)[i++] = p->index + 1;
    UNPROTECT(1);
    return ret;
}


SEXP phyr_descendants(SEXP rtree, SEXP node, SEXP visit, SEXP order)
{
    int i;
    int sz;
    int end = 0;
    int cnt = 0;
    int ndesc = 0;
    struct node *d;
    struct phy *phy = (struct phy*)R_ExternalPtrAddr(rtree);

    SEXP buf;
    SEXP root = PROTECT(list1(buf = allocVector(VECSXP, 1000)));
    SEXP tail = root;

    phy_traverse_prepare(phy, phy->nodes[phy->vseq[INTEGER(node)[0]-1]],
        INTEGER(visit)[0], INTEGER(order)[0]);

    while ((d = phy_traverse_step(phy)) != 0) {
        ndesc++;
        SET_VECTOR_ELT(buf, cnt++, ScalarInteger(d->index+1));
        if (cnt == 1000) {
            tail = SETCDR(tail, list1(buf = allocVector(VECSXP, 1000)));
            cnt = 0;
        }
    }

    SEXP descendants = PROTECT(allocVector(INTSXP, ndesc));

    while (root != R_NilValue) {
        sz = CDR(root) == R_NilValue ? cnt : 1000;
        for (i = 0; i < sz; ++i)
            INTEGER(descendants)[end++] = INTEGER(VECTOR_ELT(CAR(root), i))[0];
        root = CDR(root);
    }

    UNPROTECT(2);
    return descendants;
}


SEXP phyr_extract_clade(SEXP rtree, SEXP node)
{
    SEXP rclade;
    struct phy *phy = (struct phy *)R_ExternalPtrAddr(rtree);
    struct phy *clade = phy_extract_clade(
        phy->nodes[phy->vseq[INTEGER(node)[0]-1]]);
    if (clade) {
        rclade = PROTECT(R_MakeExternalPtr(clade, R_NilValue, R_NilValue));
        R_RegisterCFinalizer(rclade, &phyr_tree_free);
        setAttrib(rclade, install("root"), ScalarInteger(clade->root->index+1));
        setAttrib(rclade, install("Ntip"), ScalarInteger(clade->ntip));
        setAttrib(rclade, install("Nnode"), ScalarInteger(clade->nnode));
        UNPROTECT(1);
        return rclade;
    }
    error(phy_errmsg());
}


SEXP phyr_extract_subtree(SEXP rtree, SEXP ntip, SEXP tips)
{
    SEXP rsubtree;
    int i;
    int Ntip = INTEGER(ntip)[0];
    struct phy *phy = (struct phy *)R_ExternalPtrAddr(rtree);
    struct node *nodes[Ntip];

    for (i = 0; i < Ntip; ++i)
        nodes[i] = phy->nodes[phy->vseq[INTEGER(tips)[i]-1]];

    struct phy *subtree = phy_extract_subtree(Ntip, nodes, phy);
    if (subtree) {
        rsubtree = PROTECT(R_MakeExternalPtr(subtree, R_NilValue, R_NilValue));
        R_RegisterCFinalizer(rsubtree, &phyr_tree_free);
        setAttrib(rsubtree, install("root"), ScalarInteger(subtree->root->index+1));
        setAttrib(rsubtree, install("Ntip"), ScalarInteger(subtree->ntip));
        setAttrib(rsubtree, install("Nnode"), ScalarInteger(subtree->nnode));
        UNPROTECT(1);
        return rsubtree;
    }
    error(phy_errmsg());
}


SEXP phyr_ladderize(SEXP rtree, SEXP ndesc)
{
    int i = 0, j = 0, k = 0;
    int *n;
    struct node *p;
    struct node *q;
    struct node *lfdesc;
    struct node *rtdesc;
    struct phy *phy = (struct phy *)R_ExternalPtrAddr(rtree);

    SEXP perm = PROTECT(allocVector(INTSXP, phy->nnode));

    n = INTEGER(ndesc);

    phy_traverse_prepare(phy, phy->root, INTERNAL_NODES_ONLY, PREORDER);

    while ((p = phy_traverse_step(phy)) != 0) {
        lfdesc = p->lfdesc;
        rtdesc = lfdesc->next;

        if (n[lfdesc->index] < n[rtdesc->index]) {
            rtdesc->next = lfdesc;
            rtdesc->prev = NULL;
            lfdesc->next = NULL;
            lfdesc->prev = rtdesc;
            p->lfdesc = rtdesc;
        }
    }

    p = phy->root;
    while (p) {
        phy->nodes[i] = p;
        if (p->ndesc)
        {
            INTEGER(perm)[phy->ntip + j] = p->index + 1;
            p->index = phy->ntip + j;
            phy->inodes[j++] = p;
        }
        else
        {
            INTEGER(perm)[k] = p->index + 1;
            p->index = k++;
        }
        phy->vseq[p->index] = i++;

        if (p->lfdesc)
            p = p->lfdesc;
        else if (p->next)
            p = p->next;
        else {
            q = p;
            while (p->anc != 0 && p->next == 0) {
                p = p->anc;
                p->lastvisit = q;
            }
            p = p->next;
        }
    }

    UNPROTECT(1);
    return perm;
}


SEXP phyr_rotate(SEXP rtree, SEXP index)
{
    int i;
    int j;
    int k;
    int n;
    struct node *p;
    struct node *q;
    struct node *node;
    struct phy *phy;

    n = LENGTH(index);

    phy = (struct phy *)R_ExternalPtrAddr(rtree);

    for (i = 0; i < n; ++i) {
        node = phy_getnode_with_index(phy, INTEGER(index)[i]);

        p = node->lfdesc;
        q = p->next;

        p->next = NULL;
        p->prev = q;
        q->next = p;
        q->prev = NULL;

        node->lfdesc = q;
    }

    i = 0;
    j = 0;
    k = 0;

    p = phy->root;
    while (p) {
        phy->nodes[i] = p;
        if (p->ndesc) {
            p->index = phy->ntip + j;
            phy->inodes[j++] = p;
        } else
            p->index = k++;
        phy->vseq[p->index] = i++;

        if (p->lfdesc)
            p = p->lfdesc;
        else if (p->next)
            p = p->next;
        else {
            q = p;
            while (p->anc != 0 && p->next == 0) {
                p = p->anc;
                p->lastvisit = q;
            }
            p = p->next;
        }
    }

    return R_NilValue;
}

#ifndef PHY_H
#define PHY_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#ifdef __cplusplus
extern "C" {
#endif

#define LIBPHY_VERSION "1.0.0"

/* Macros used for tree traversals. */
#define PREORDER 0
#define POSTORDER 1
#define ALL_NODES 0
#define INTERNAL_NODES_ONLY 1

struct node {
    /* Terminal nodes are numbered 0 to N-1, where
    ** N is the number of terminal nodes in the
    ** phylogeny. Internal nodes are numbered N
    ** up to one less than the number of nodes in
    ** the phylogeny. */
    int index;

    /* Number of immediate descendants */
    int ndesc;

    /* Name of node */
    char *lab;

    /* First immediate descendant. The other descendants
    ** are stored as a doubly-linked list with this descendant
    ** as the head. The next descendant is found as lfdesc->next,
    ** the next as lfdesc->next->next, and so on.
    */
    struct node *lfdesc;

    /* Next sibling node. Sibling nodes share the same ancestor. */
    struct node *next;

    /* Previous sibling node */
    struct node *prev;

    /* Ancestral node */
    struct node *anc;

    /* This is the last node visited in a preorder traversal of
    ** the subtree rooted at this node (null if this node
    ** is a terminal node). It is always a terminal node if not null. */
    struct node *lastvisit;

    /* The length of the branch that leads to this node */
    double brlen;

    /* Arbitrary data for client programs using this library */
    void *data;

    /* Function pointer to free data held by node */
    void (*data_free)(void *);
};


struct phy {
    /* Number of terminal nodes in phylogeny */
    int ntip;

    /* Number of nodes in phylogeny */
    int nnode;

    /* Root node of phylogeny */
    struct node *root;

    /* Array of all nodes arranged in a preorder traversal */
    struct node **nodes;

    /* Array of internal nodes arranged in a preorder traversal.
    ** Note that as a consequence of how nodes are numbered, the
    ** index of an internal node minus the number of tips corresponds
    ** to its position in this array. */
    struct node **inodes;

    /* Preorder visitation sequence of node indices.
    ** For example, vseq[i] will return the position of
    ** the node having index i in the nodes array. */
    int *vseq;

    /* State information used during tree traversals.
    ** The canonical way to perform a tree traversal from
    ** a particular node is like so,
    **
    ** phy_traverse_prepare(phy, node, ALL_NODES, PREORDER);
    ** while ((node = phy_traverse_step(phy)) != 0) {
    **      // perform some operations on node
    ** }
    **
    ** Which will visit all nodes in the tree in preorder
    ** traversal sequence. Other permissible values are,
    **
    **      phy_traverse_prepare(phy, node, ALL_NODES, POSTORDER);
    **      phy_traverse_prepare(phy, node, INTERNAL_NODES_ONLY, PREORDER);
    **      phy_traverse_prepare(phy, node, INTERNAL_NODES_ONLY, POSTORDER);
    **
    ** Attempting to perform a traversal by using the
    ** loop construct without calling phy_traverse_prepare
    ** will break from the loop immediately.
    */
    struct {
        int visit;
        int order;
        int cursor;
        struct node *begin;
        struct node *end;
        struct node *next;
    } t;
};

struct phy *phy_read(const char *newick);
char *phy_write(struct phy *phy);
struct phy *phy_fread(FILE *in);
void phy_fwrite(struct phy *phy, FILE *out);
void phy_free(struct phy *phy);
void phy_node_add_data(struct node *node, void *data, void (*data_free)(void *));
void phy_traverse_prepare(struct phy *phy, struct node *node, int visit, int order);
struct node *phy_traverse_step(struct phy *phy);
int phy_isbinary(struct phy *phy);
struct node *phy_getnode_with_index(struct phy *phy, int index);
struct node *phy_getnode_with_label(struct phy *phy, const char *label);
struct phy *phy_extract_clade(struct node *node);
struct phy *phy_extract_subtree(int ntip, struct node **tips, struct phy *phy);
const char *phy_errmsg();

#ifdef __cplusplus
}
#endif

#endif /* PHY_H */

#include "phy.h"

#define PHY_ERR1 "cannot allocate memory"
#define PHY_ERR2 "encountered unexpected character in Newick string node label/branch length"
#define PHY_ERR3 "detected unifurcation in Newick string"
#define PHY_ERR4 "malformed Newick string"


static int phy_errno = 0;


static struct node *node_new()
{
    struct node *node = malloc(sizeof(struct node));
    if (!node) {
        phy_errno = 1;
        return NULL;
    }
    node->index = -1;
    node->ndesc = 0;
    node->lab = 0;
    node->note = 0;
    node->lfdesc = 0;
    node->next = 0;
    node->prev = 0;
    node->anc = 0;
    node->lastvisit = 0;
    node->brlen = 0;
    node->data = 0;
    node->data_free = 0;
    return node;
}


static void node_free(struct node *node)
{
    if (node) {
        if (node->data && node->data_free)
            node->data_free(node->data);
        free(node->lab);
        free(node->note);
        free(node);
    }
}


/*
** This function is called on the root node whenever
** an error is encountered during the process of
** building the phylogeny.
*/
static void cleanup(struct node *root) {
    struct node *q, *r, *p = root;
    while (p) {
        if (p->lfdesc)
            p = p->lfdesc;
        else if (p->next)
            p = p->next;
        else {
            // on entry p is a terminal node that marks a clade
            // boundary. this will be the last node visited in a
            // preorder traversal of the subtree rooted at each node
            // on the path back from this node up to and including the
            // first encountered node with a ->next sibling
            q = p->prev;
            while (q != 0) {
                r = q;
                q = q->prev;
                free(r->lab);
                free(r->note);
                free(r);
                r = 0;
            }
            while (p->anc != 0 && p->next == 0) {
                q = p;
                p = p->anc;
                free(q->lab);
                free(q->note);
                free(q);
                q = 0;
            }
            q = p;
            p = p->next;
            free(q->lab);
            free(q->note);
            free(q);
            q = 0;
        }
    }
}


static void node_add_child(struct node *parent, struct node *child)
{
    struct node *r;
    switch (parent->ndesc) {
        case 0:
            parent->lfdesc = child;
            break;
        case 1:
            parent->lfdesc->next = child;
            child->prev = parent->lfdesc;
            break;
        case 2:
            parent->lfdesc->next->next = child;
            child->prev = parent->lfdesc->next;
            break;
        case 3:
            parent->lfdesc->next->next->next = child;
            child->prev = parent->lfdesc->next->next;
            break;
        case 4:
            parent->lfdesc->next->next->next->next = child;
            child->prev = parent->lfdesc->next->next->next;
            break;
        case 5:
            parent->lfdesc->next->next->next->next->next = child;
            child->prev = parent->lfdesc->next->next->next->next;
            break;
        default:
            for (r = parent->lfdesc; r->next != 0; r = r->next) {};
            r->next = child;
            child->prev = r;
    }
    parent->ndesc += 1;
    child->anc = parent;
}


static struct phy *phy_build(struct node *root, int nnode, int ntip)
{
    int i = 0, j = 0, k = 0;
    struct node *p, *q;
    struct phy *phy = malloc(sizeof(struct phy));
    if (!phy) {
        phy_errno = 1;
        return NULL;
    }
    phy->ntip = ntip;
    phy->nnode = nnode;
    phy->root = root;
    phy->nodes = malloc(nnode * sizeof(struct node *));
    if (!phy->nodes) {
        phy_errno = 1;
        free(phy);
        return NULL;
    }
    phy->vseq = malloc(nnode * sizeof(int));
    if (!phy->vseq) {
        phy_errno = 1;
        free(phy->nodes);
        free(phy);
        return NULL;
    }
    phy->inodes = malloc((nnode-ntip) * sizeof(struct node *));
    if (!phy->inodes) {
        phy_errno = 1;
        free(phy->vseq);
        free(phy->nodes);
        free(phy);
        return NULL;
    }
    /* Assign indices to each node and record their visitation sequence
    ** in a preorder traversal beginning from the root. Terminal nodes
    ** are numbered 0 ... ntip-1 and internal nodes are numbered from
    ** ntip ... nnode-1. The root node always has index ntip using this
    ** scheme. */

    p = phy->root;
    while (p) {
        phy->nodes[i] = p;
        if (p->ndesc) {
            p->index = ntip + j;
            phy->inodes[j++] = p;
        } else
            p->index = k++;
        phy->vseq[p->index] = i++;

        if (p->lfdesc)
            p = p->lfdesc;
        else if (p->next)
            p = p->next;
        else {
            // on entry p is a terminal node that marks a clade
            // boundary. this will be the last node visited in a
            // preorder traversal of the subtree rooted at each node
            // on the path back from this node up to and including the
            // first encountered node with a ->next sibling
            q = p;
            while (p->anc != 0 && p->next == 0) {
                p = p->anc;
                p->lastvisit = q;
            }
            p = p->next;
        }
    }

    /* Initialize the tree traversal structure state */
    phy->t.visit = 0;
    phy->t.order = 0;
    phy->t.cursor = 0;
    phy->t.begin = 0;
    phy->t.end = 0;
    phy->t.next = 0;
    return phy;
}


struct ReadCtx {
    unsigned int n;
    unsigned int nAlloc;
    unsigned int cursor;
    unsigned int ntip;
    unsigned int nnode;
    char *z;
    const char *newick;
    struct node *q;
    struct node *p;
    struct node *root;
};


static int read_label(struct ReadCtx *ctx)
{
    char c;
    int toread = 1;
    while (toread) {
        c = ctx->newick[ctx->cursor++];
        switch (c) {
            case ':':
            case ',':
            case ')':
            case ';':
            case '[':
                ctx->cursor--;
                toread = 0;
                break;
            case ' ':
            case '\n':
            case '\r':
            case '\v':
            case '\t':
            case '\f':
            case '(':
            case ']':
                phy_errno = 2;
                return 1;
            case '\0':
                phy_errno = 4;
                return 1;
            default:
                if (ctx->n == ctx->nAlloc) {
                    ctx->nAlloc += 100;
                    ctx->z = realloc(ctx->z, ctx->nAlloc);
                    if (!ctx->z) {
                        phy_errno = 1;
                        return 1;
                    }
                }
                ctx->z[ctx->n++] = c;
        }
    }
    ctx->z[ctx->n] = 0;
    if (ctx->n) {
        ctx->p->lab = calloc(ctx->n+1, 1);
        if (!ctx->p->lab) {
            phy_errno = 1;
            ctx->p->lab = 0;
            return 1;
        }
        strcpy(ctx->p->lab, ctx->z);
    }
    ctx->n = 0;
    return 0;
}


static int read_note(struct ReadCtx *ctx)
{
    char c;
    int opened;
    c = ctx->newick[ctx->cursor++];
    if (c == '[') {
        opened = 1;
        while ((c = ctx->newick[ctx->cursor++]), opened > 0) {
            if (c == '\0') {
                phy_errno = 4;
                return 1;
            }
            if (c == '[')
                ++opened;
            else if (c == ']')
                --opened;
            if (opened) {
                if (ctx->n == ctx->nAlloc) {
                    ctx->nAlloc += 100;
                    ctx->z = realloc(ctx->z, ctx->nAlloc);
                    if (!ctx->z) {
                        phy_errno = 1;
                        return 1;
                    }
                }
                ctx->z[ctx->n++] = c;
            }
        }
        ctx->z[ctx->n] = 0;
        if (ctx->n) {
            ctx->p->note = calloc(ctx->n+1, 1);
            if (!ctx->p->note) {
                phy_errno = 1;
                ctx->p->note = 0;
                return 1;
            }
            strcpy(ctx->p->note, ctx->z);
        }
        ctx->n = 0;
    } else
        ctx->cursor--;
    return 0;
}


static int read_brlen(struct ReadCtx *ctx)
{
    char c;
    int toread = 1;
    c = ctx->newick[ctx->cursor++];
    if (c == ':') {
        while (toread) {
            c = ctx->newick[ctx->cursor++];
            switch (c) {
                case 'e':
                case '-':
                case '+':
                case '0':
                case '1':
                case '2':
                case '3':
                case '4':
                case '5':
                case '6':
                case '7':
                case '8':
                case '9':
                case '.':
                    if (ctx->n == ctx->nAlloc) {
                        ctx->nAlloc += 100;
                        ctx->z = realloc(ctx->z, ctx->nAlloc);
                        if (!ctx->z) {
                            phy_errno = 1;
                            return 1;
                        }
                    }
                    ctx->z[ctx->n++] = c;
                    break;
                case ',':
                case ')':
                case ';':
                    ctx->cursor--;
                    toread = 0;
                    break;
                case '\0':
                    phy_errno = 4;
                    return 1;
                default:
                    phy_errno = 2;
                    return 1;
            }
        }
    } else
        ctx->cursor--;
    ctx->z[ctx->n] = 0;
    if (ctx->n)
        ctx->p->brlen = atof(ctx->z);
    ctx->n = 0;
    return 0;
}


static struct node *read_newick(struct ReadCtx *ctx)
{
    char c;

    if (ctx->newick[strlen(ctx->newick)-1] != ';') {
        phy_errno = 4;
        return 0;
    }

    ctx->nnode++;
    ctx->root = node_new();
    if (ctx->root == NULL) {
        phy_errno = 1;
        return NULL;
    }
    ctx->p = ctx->root;
    while ((c = ctx->newick[ctx->cursor++]) != ';') {
        switch (c) {
            case '(':
                if (ctx->cursor > 1 && !(
                    ctx->newick[ctx->cursor-2] == ','
                    || ctx->newick[ctx->cursor-2] == '(')) {
                    phy_errno = 4;
                    return NULL;
                }
                ctx->nnode++;
                ctx->q = node_new();
                if (ctx->q == NULL) {
                    phy_errno = 1;
                    return NULL;
                }
                node_add_child(ctx->p, ctx->q);
                ctx->p = ctx->q;
                break;
            case ',':
                ctx->nnode++;
                if (!ctx->p->anc && ctx->p != ctx->root) {
                    phy_errno = 4;
                    free(ctx->p->lab);
                    free(ctx->p);
                    return NULL;
                } else if (ctx->p == ctx->root) {
                    phy_errno = 4;
                    return NULL;
                }
                ctx->q = node_new();
                if (ctx->q == NULL) {
                    phy_errno = 1;
                    return NULL;
                }
                node_add_child(ctx->p->anc, ctx->q);
                if (!ctx->p->ndesc)
                    ctx->ntip++;
                ctx->p = ctx->q;
                break;
            case ')':
                if (!ctx->p->ndesc)
                    ctx->ntip++;
                if (!ctx->p->anc && ctx->p != ctx->root) {
                    phy_errno = 4;
                    free(ctx->p->lab);
                    free(ctx->p);
                    return NULL;
                } else if (ctx->p == ctx->root) {
                    phy_errno = 4;
                    return NULL;
                }
                ctx->p = ctx->p->anc;
                if (ctx->p->ndesc < 2) {
                    phy_errno = 3;
                    return NULL;
                }
                break;
            default:
                ctx->cursor--;
                if (read_label(ctx))
                    return NULL;
                if (read_note(ctx))
                    return NULL;
                if (read_brlen(ctx))
                    return NULL;
        }
    }
    if (ctx->p->ndesc < 2) {
        phy_errno = 3;
        return NULL;
    }
    return ctx->p;
}


struct WriteCtx {
    unsigned int n;
    unsigned int nAlloc;
    char *newick;
};


static void write_chars(char *z, struct WriteCtx *ctx)
{
    unsigned int len = strlen(z);
    if (len > 0) {
        int togo = ctx->nAlloc - ctx->n;
        while (togo <= len) {
            ctx->nAlloc += 100;
            ctx->newick = realloc(ctx->newick, ctx->nAlloc);
            if (!ctx->newick) {
                perror("Error in write_label");
                exit(0);
            }
            togo = ctx->nAlloc - ctx->n;
        }
        ctx->newick[ctx->n] = 0;
        strcat(ctx->newick, z);
        ctx->n += len;
    }
}


static void write_brlen(struct node *p, struct WriteCtx *ctx)
{
    if (p->brlen > 0) {
        char brlen[25] = "";
        sprintf(brlen, ":%f", p->brlen);
        write_chars(brlen, ctx);
    }
}


static void write_label(struct node *p, struct WriteCtx *ctx)
{
    if (p->lab != 0) {
        if (strlen(p->lab) > 0)
            write_chars(p->lab, ctx);
    }
}


static void write_newick(struct node *node, struct WriteCtx *ctx)
{
    struct node *d = 0;
    if (node->ndesc) {
        for (d = node->lfdesc; d != 0; d = d->next) {
            if (d == node->lfdesc)
                write_chars("(", ctx);
            write_newick(d, ctx);
            if (d->next)
                write_chars(",", ctx);
        }
        write_chars(")", ctx);
    }
    write_label(node, ctx);
    write_brlen(node, ctx);
}


char *phy_write(struct phy *phy)
{
    struct WriteCtx ctx = {0, 0, 0};
    write_newick(phy->root, &ctx);
    write_chars(";", &ctx);
    return ctx.newick;
}


void phy_node_add_data(struct node *node, void *data, void (*data_free)(void *))
{
    if (node->data && node->data_free)
        node->data_free(node->data);
    node->data = data;
    node->data_free = data_free;
}


void phy_traverse_prepare(struct phy *phy, struct node *node, int visit, int order)
{
    phy->t.visit = visit;
    phy->t.order = order;
    if (order == PREORDER) {
        phy->t.begin = node;
        phy->t.next = node;
        if (node->ndesc) {
            if (visit == ALL_NODES) {
                phy->t.end = phy->nodes[phy->vseq[node->lastvisit->index]];
            } else {
                struct node *p;
                p = node->lastvisit->anc->lfdesc;
                while (p->ndesc)
                    p = p->lastvisit->anc->lfdesc;
                phy->t.end = phy->inodes[p->anc->index - phy->ntip];
            }
        } else {
            phy->t.end = node;
        }
        if (node->ndesc) {
            if (visit == ALL_NODES)
                phy->t.cursor = phy->vseq[node->index] + 1;
            else
                phy->t.cursor = node->index - phy->ntip + 1;
        }
    } else { // POSTORDER
        phy->t.end = node;
        if (node->ndesc) {
            if (visit == ALL_NODES) {
                phy->t.begin = phy->nodes[phy->vseq[node->lastvisit->index]];
            } else {
                struct node *p;
                p = node->lastvisit->anc->lfdesc;
                while (p->ndesc)
                    p = p->lastvisit->anc->lfdesc;
                phy->t.begin = phy->inodes[p->anc->index - phy->ntip];
            }
        } else {
            phy->t.begin = node;
        }
        phy->t.next = phy->t.begin;
        if (node->ndesc) {
            if (visit == ALL_NODES)
                phy->t.cursor = phy->vseq[phy->t.begin->index] - 1;
            else
                phy->t.cursor = phy->t.begin->index - phy->ntip - 1;
        }
    }
}


struct node *phy_traverse_step(struct phy *phy)
{
    struct node *node = phy->t.next;
    if (node) {
        if (node != phy->t.end) {
            if (phy->t.visit == ALL_NODES) {
                if (phy->t.order == PREORDER)
                    phy->t.next = phy->nodes[phy->t.cursor++];
                else
                    phy->t.next = phy->nodes[phy->t.cursor--];
            } else {
                if (phy->t.order == PREORDER)
                    phy->t.next = phy->inodes[phy->t.cursor++];
                else
                    phy->t.next = phy->inodes[phy->t.cursor--];
            }
        } else {
            phy->t.visit = 0;
            phy->t.order = 0;
            phy->t.cursor = 0;
            phy->t.begin = 0;
            phy->t.end = 0;
            phy->t.next = 0;
        }
    }
    return node;
}


int phy_isbinary(struct phy *phy)
{
    return phy->nnode == (2 * phy->ntip - 1) ? 1 : 0;
}


void phy_free(struct phy *phy)
{
    if (phy) {
        int i;
        for (i = 0; i < phy->nnode; ++i)
            node_free(phy->nodes[i]);
        free(phy->nodes);
        free(phy->inodes);
        free(phy->vseq);
        free(phy);
    }
}


struct phy *phy_read(const char *newick)
{
    struct node *root = 0;
    struct phy *phy = 0;
    struct ReadCtx ctx = {0, 0, 0, 0, 0, 0, newick, 0, 0, 0};
    root = read_newick(&ctx);
    if (root) {
        phy = phy_build(root, ctx.nnode, ctx.ntip);
        if (!phy)
            cleanup(ctx.root);
    }
    else
        cleanup(ctx.root);
    free(ctx.z);
    return phy;
}


struct phy *phy_fread(FILE *in)
{
    struct phy *phy = 0;
    char *newick = 0;
    int nBytes = 0;
    fseek(in, 0, SEEK_END);
    nBytes = ftell(in);
    rewind(in);
    if (nBytes)
        newick = malloc(nBytes);
    if (newick) {
        fread(newick, nBytes, 1, in);
        phy = phy_read(newick);
        free(newick);
    }
    return phy;
}


void phy_fwrite(struct phy *phy, FILE *out)
{
    char *newick = phy_write(phy);
    fputs(newick, out);
    free(newick);
}


struct node *phy_getnode_with_index(struct phy *phy, int index)
{
    return phy->nodes[phy->vseq[index]];
}


struct node *phy_getnode_with_label(struct phy *phy, const char *label)
{
    struct node *node = 0;
    phy_traverse_prepare(phy, phy->root, ALL_NODES, POSTORDER);
    while ((node = phy_traverse_step(phy)) != 0) {
        if (node->lab) {
            if (strcmp(label, node->lab) == 0)
                break;
        }
    }
    return node;
}


struct phy *phy_extract_clade(struct node *node)
{
    struct WriteCtx ctx = {0, 0, 0};
    struct phy *phy = 0;
    write_newick(node, &ctx);
    write_chars(";", &ctx);
    phy = phy_read(ctx.newick);
    free(ctx.newick);
    if (!phy)
        return NULL;
    phy->root->brlen = 0;
    return phy;
}


struct phy *phy_extract_subtree(int ntip, struct node **tips, struct phy *phy)
{
    int i;
    int nnode = 0;
    int bitmask[phy->nnode];
    struct node *p;
    struct node *q;
    struct node *root;
    struct node *head;

    memset(bitmask, 0, phy->nnode * sizeof(int));
    for (i = 0; i < ntip; ++i) {
        p = tips[i];
        bitmask[p->index] = 1;
        while ((p = p->anc) != 0) {
            if (bitmask[p->index])
                break;
            bitmask[p->index] = 1;
        }
    }

    root = 0;
    head = 0;

    p = phy->root;
    while (p) {
        if (bitmask[p->index] && !root) {
            nnode++;
            root = node_new();
            if (!root) {
                phy_errno = 1;
                return NULL;
            }
            head = root;
        } else if (bitmask[p->index]) {
            nnode++;
            q = node_new();
            if (!q) {
                phy_errno = 1;
                cleanup(root);
                return NULL;
            }
            q->brlen = p->brlen;
            if (p->lab) {
                q->lab = malloc(strlen(p->lab)+1);
                if (!q->lab) {
                    phy_errno = 1;
                    cleanup(root);
                    return NULL;
                }
                strcpy(q->lab, p->lab);
            }
            node_add_child(head, q);
            head = q;
        }
        if (p->lfdesc)
            p = p->lfdesc;
        else if (p->next) {
            if (bitmask[p->index])
                head = head->anc;
            p = p->next;
        }
        else {
            while (p->anc != 0 && p->next == 0) {
                if (bitmask[p->index])
                    head = head->anc;
                p = p->anc;
            }
            if (bitmask[p->index])
                head = head->anc;
            p = p->next;
        }
    }

    root->brlen = 0;
    p = root;
    while (p) {
        while (p->ndesc == 1) {
            q = p->lfdesc;
            q->brlen += p->brlen;
            q->next = p->next;
            q->prev = p->prev;
            q->anc = p->anc;
            if (p->prev)
                p->prev->next = q;
            if (p->next)
                p->next->prev = q;
            if (p->anc && p->anc->lfdesc == p)
                p->anc->lfdesc = q;
            if (p == root) {
                root = q;
                root->brlen = 0;
            }
            node_free(p);
            p = q;
            nnode--;
        }
        if (p->lfdesc)
            p = p->lfdesc;
        else if (p->next)
            p = p->next;
        else {
            while (p->anc != 0 && p->next == 0)
                p = p->anc;
            p = p->next;
        }
    }

    return phy_build(root, nnode, ntip);
}


const char *phy_errmsg() {
    switch (phy_errno) {
        case 1:
            phy_errno = 0;
            return PHY_ERR1;
        case 2:
            phy_errno = 0;
            return PHY_ERR2;
        case 3:
            return PHY_ERR3;
        case 4:
            phy_errno = 0;
            return PHY_ERR4;
        default:;
    }
    return "no errors detected";
}

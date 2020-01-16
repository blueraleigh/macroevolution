#ifndef RCM_H
#define RCM_H

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Connections.h>
#include "assert.h"
#include "phy.h"

/*
** Functionality that is shared between the continuous and discrete versions
** of the Ultra Common Mechanism (rcm) Markov model of character evolution with
** unobserved character states.
*/


#define twotothe256 \
    115792089237316195423570985008687907853269984665640564039457584007913129639936.0

#define minlikelihood  (1.0/twotothe256)
#define minusminlikelihood -minlikelihood

/* Macros for accessing the downpass, stem, and uppass
** parsimony score arrays */
#define DCLK(state, node) (state)->dclk[(node)->index]
#define SCLK(state, node) (state)->sclk[(node)->index]
#define UCLK(state, node) (state)->uclk[(node)->index]


struct rcm_data;
struct rcm_stat;

/* States on the tree
**
** To iterate over the states on the tree we do
**
** struct statelist *sl = &model->sl;
**
** for (state = sl->head; state; state = state->next) {
**     // do something with state
** }
**
*/
struct rcm_statelist {
    /* Number of states that are currently on the tree plus 1 */
    int len;

    /* Number of states in the model */
    int r;

    /* Number of dimensions for each datum */
    int p;

    /* Number of tips in the tree */
    int ntip;

    /* Number of nodes in the tree */
    int nnode;

    /* Root node index */
    int root;

    struct rcm_state *head;

    /* tail state is always non-analog state (when len < len_max) */
    struct rcm_state *tail;

    /* Pointer to the end of the list of states malloc'd but
    ** then removed from the tree that are available for reuse
    ** i.e.,
    **
    ** NULL <-- state <-- state <-- state <-- fl
    */
    struct rcm_state *fl;

    /* Function for state stat allocation */
    struct rcm_stat *(*stat_alloc)(int);

    void (*stat_free)(struct rcm_stat *);

    /* Function for computing marginal likelihood of data
    ** in state */
    double (*stat_loglk)(struct rcm_stat *, struct rcm_data *);
};


struct rcm_state {

    /* Identifier */
    int index;

    /* Count of terminal nodes assigned to this state */
    int ntip;

    /* Downpass conditional likelihoods */
    double *dclk;

    /* Stem conditional likelihoods */
    double *sclk;

    /* Uppass conditional likelihoods */
    double *uclk;

    /* Marginal log likelihood of data generated from state */
    double loglk;

    /* Memory alloc'd */
    void *mem;

    /* Resource states stored as doubly linked list */
    struct rcm_state *next;

    struct rcm_state *prev;

    /* Running stat that tracks sufficient statistics for
    ** computing marginal likelihoods */
    struct rcm_stat *stat;
};


struct rcm {

    /* Number of states in the model. */
    int r;

    /* Number of dimensions for each datum */
    int p;

    /* Each item is the current resource class assignment for a node */
    struct rcm_state **stateid;

    int *stateid_r;

    /* Encodes state information about whether the conditional likelihoods
    ** of a node need to be updated to reflect new partition of terminals */
    int *needsupdate_up;

    int *needsupdate_dn;

    int *lzd;

    int *lzu;

    double rate;

    double rate_max;

    double dataloglk;

    double treeloglk;

    struct phy *phy;

    /* Observed data at the terminal nodes of phy */
    struct rcm_data *data;

    struct rcm_statelist sl;

    double (*push)(int, struct rcm_data *, struct rcm_stat *);
    void (*push0)(int, struct rcm_data *, struct rcm_stat *);

    void (*pop)(int, struct rcm_data *, struct rcm_stat *);
};


struct rcm_mcmc {
    int niter;

    int sample_freq;

    double update_node;

    double update_rate;

    double update_weight;

    double tune_rate;

    struct rcm *model;

    SEXP outputConn;
};


void rcm_clear(struct rcm *model);
struct rcm *rcm_init_start(int p, int r, double alpha, struct phy *phy);
double rcm_dataloglk(struct rcm *model);
double rcm_treeloglk(struct rcm *model);
void rcm_init_states(int *stateid, struct rcm *model);

#endif

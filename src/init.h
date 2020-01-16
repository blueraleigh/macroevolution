#ifndef MACROEVOLUTION_INIT_H
#define MACROEVOLUTION_INIT_H

#include <R.h>
#include <Rinternals.h>

/* Declarations for all native routine entry points */
SEXP phyr_read_newick(SEXP);
SEXP phyr_write_newick(SEXP);
SEXP phyr_tiplabels(SEXP);
SEXP phyr_node_brlens(SEXP);
SEXP phyr_node_ages(SEXP);
SEXP phyr_ancestors(SEXP, SEXP);
SEXP phyr_children(SEXP, SEXP);
SEXP phyr_descendants(SEXP, SEXP, SEXP, SEXP);
SEXP phyr_extract_clade(SEXP, SEXP);
SEXP phyr_extract_subtree(SEXP, SEXP, SEXP);
SEXP phyr_plot_ctree(SEXP, SEXP, SEXP);
SEXP phyr_plot_ptree(SEXP, SEXP);
SEXP phyr_ladderize(SEXP, SEXP);
SEXP phyr_rotate(SEXP, SEXP);

SEXP mkdmm_model_init(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP mkdmm_mcmc_run(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP mkdmm_marginal_asr(SEXP);
SEXP mkdmm_posterior_multinomial(SEXP);

SEXP rcm_dmm_model_init(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP rcm_model_mcmc_run(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP rcm_marginal_asr(SEXP);
SEXP rcm_stochastic_map_expected_counts(SEXP);
SEXP rcm_stochastic_map(SEXP, SEXP);
SEXP rcm_dmm_posterior_multinomial(SEXP);
SEXP rcm_expected_loss(SEXP, SEXP, SEXP);
SEXP rcm_posterior_coincidence(SEXP);

#endif

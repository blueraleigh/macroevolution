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

SEXP do_fitch_pscore(SEXP, SEXP, SEXP);
SEXP do_fitch_mpr(SEXP, SEXP, SEXP, SEXP);
SEXP do_fitch_history(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP do_fitch_count(SEXP, SEXP, SEXP);
SEXP do_fitch_count2(SEXP, SEXP, SEXP, SEXP);
SEXP do_fitch_phat(SEXP, SEXP, SEXP);

SEXP mk2_model_init(SEXP, SEXP);
SEXP mk2_loglk(SEXP, SEXP);
SEXP mk2_grad1(SEXP, SEXP);
SEXP mk2_grad2(SEXP, SEXP);
SEXP mk2_uclk(SEXP, SEXP);
SEXP mk2_marginal_asr(SEXP, SEXP);
SEXP mk2_perr(SEXP, SEXP);
SEXP mk2_mcmc_slice(SEXP, SEXP, SEXP, SEXP);
SEXP mk2_loglk_conditional(SEXP, SEXP, SEXP);
SEXP mk2_gradient_conditional(SEXP, SEXP, SEXP);

SEXP mkepoch_rate_init(SEXP, SEXP, SEXP, SEXP);
SEXP mkepoch_model_init(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP mkepoch_loglk(SEXP, SEXP);
SEXP mkepoch_marginal_asr(SEXP, SEXP);

SEXP rcm_dmm_model_init(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP rcm_kmeans_counts_model_init(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP rcm_model_mcmc_run(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP rcm_marginal_asr(SEXP);
SEXP rcm_stochastic_map_expected_counts(SEXP);
SEXP rcm_stochastic_map(SEXP, SEXP);
SEXP rcm_dmm_posterior_multinomial(SEXP);
SEXP rcm_expected_loss(SEXP, SEXP, SEXP);
SEXP rcm_posterior_coincidence(SEXP);
SEXP rcm_dmm_decomp(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP rcm_pij(SEXP, SEXP);


#endif

#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

#include "init.h"

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}


static const R_CallMethodDef CallEntries[] = {
    CALLDEF(phyr_read_newick, 1),
    CALLDEF(phyr_write_newick, 1),
    CALLDEF(phyr_tiplabels, 1),
    CALLDEF(phyr_node_brlens, 1),
    CALLDEF(phyr_node_ages, 1),
    CALLDEF(phyr_ancestors, 2),
    CALLDEF(phyr_children, 2),
    CALLDEF(phyr_descendants, 4),
    CALLDEF(phyr_extract_clade, 2),
    CALLDEF(phyr_extract_subtree, 3),
    CALLDEF(phyr_plot_ctree, 3),
    CALLDEF(phyr_plot_ptree, 2),
    CALLDEF(phyr_ladderize, 2),
    CALLDEF(phyr_rotate, 2),

    CALLDEF(mkdmm_model_init, 7),
    CALLDEF(mkdmm_mcmc_run, 9),
    CALLDEF(mkdmm_marginal_asr, 1),
    CALLDEF(mkdmm_posterior_multinomial, 1),

    CALLDEF(do_fitch_pscore, 3),
    CALLDEF(do_fitch_mpr, 4),
    CALLDEF(do_fitch_history, 5),
    CALLDEF(do_fitch_count, 3),
    CALLDEF(do_fitch_count2, 4),
    CALLDEF(do_fitch_phat, 3),

    CALLDEF(mk2_model_init, 2),
    CALLDEF(mk2_loglk, 2),
    CALLDEF(mk2_grad1, 2),
    CALLDEF(mk2_grad2, 2),
    CALLDEF(mk2_uclk, 2),
    CALLDEF(mk2_marginal_asr, 2),
    CALLDEF(mk2_perr, 2),
    CALLDEF(mk2_mcmc_slice, 4),
    CALLDEF(mk2_loglk_conditional, 3),
    CALLDEF(mk2_gradient_conditional, 3),

    CALLDEF(mkepoch_rate_init, 4),
    CALLDEF(mkepoch_model_init, 6),
    CALLDEF(mkepoch_loglk, 2),
    CALLDEF(mkepoch_marginal_asr, 2),

    CALLDEF(rcm_dmm_model_init, 7),
    CALLDEF(rcm_kmeans_counts_model_init, 6),
    CALLDEF(rcm_model_mcmc_run, 7),
    CALLDEF(rcm_marginal_asr, 1),
    CALLDEF(rcm_stochastic_map_expected_counts, 1),
    CALLDEF(rcm_stochastic_map, 2),
    CALLDEF(rcm_dmm_posterior_multinomial, 1),
    CALLDEF(rcm_expected_loss, 3),
    CALLDEF(rcm_posterior_coincidence, 1),

    {NULL, NULL, 0}
};


void attribute_visible R_init_macroevolution(DllInfo *info)
{
    R_registerRoutines(info, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    R_forceSymbols(info, TRUE);
}

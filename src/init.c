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

    CALLDEF(mkepoch_rate_init, 4),
    CALLDEF(mkepoch_model_init, 6),
    CALLDEF(mkepoch_loglk, 2),
    CALLDEF(mkepoch_marginal_asr, 2),

    {NULL, NULL, 0}
};


void attribute_visible R_init_macroevolution(DllInfo *info)
{
    R_registerRoutines(info, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    R_forceSymbols(info, TRUE);
}

make.mk2 = function(phy, x) {
    stopifnot(is.tree(phy))
    stopifnot(tree.isbinary(phy))
    stopifnot(!is.null(names(x)))
    stopifnot(storage.mode(x) == "integer")

    stopifnot(all(names(x) %in% tiplabels(phy))
        && all(tiplabels(phy) %in% names(x)))

    stopifnot(all(x[!is.na(x)] %in% c(0L, 1L)))

    x = x[tiplabels(phy)]

    clk = matrix(0, 2, Ntip(phy))

    for (i in 1L:Ntip(phy))
    {
        if (!is.na(x[i]))
            clk[x[i]+1L, i] = 1
        else
            clk[, i] = 1
    }

    model = .Call(mk2_model_init, clk, phy)

    loglk = function(par) {
        if (length(par) == 1L)
            par = c(1, par)
        .Call(mk2_loglk, par, model)
    }

    return (loglk)
}


mk2.grad = function(par, lik, intermediates=FALSE) {
    if (length(par) == 1L)
    {
        model = environment(lik)$model
        if (!intermediates)
            return (sum(.Call(mk2_grad1, c(1, par), model)))
        else
            return (.Call(mk2_grad1, c(1, par), model))
    }
    else
    {
        model = environment(lik)$model
        if (!intermediates)
            return (rowSums(.Call(mk2_grad2, par, model)))
        else
            return (.Call(mk2_grad2, par, model))
    }
}


mk2.asr = function(par, lik) {
    model = environment(lik)$model
    asr = .Call(mk2_marginal_asr, par, model)
    return (asr)
}

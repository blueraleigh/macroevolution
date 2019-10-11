#' Initialize epoch transition rate matrices
#'
#' @param layout An array with dim=c(nstate, nstate, nepoch). Each matrix in
#'  the array specifies the structure of the transition rate matrix for the
#'  corresponding epoch. A value of 0 means that the given transition rate is
#'  constrained to be 0. Any other value corresponds to the index of a transition
#'  rate parameter. Therefore, different transitions can be constrained
#'  to have the same rate by constraining them to share the same parameter.
#' @param covarion A boolean value indicating whether or not \code{layout}
#'  is setup for a covarion-like model.
#' @return An object that inherits class \code{ratematrix.mk.epoch} that can
#'  be used with \code{make.mk.epoch} to construct a likelihood function.
#' @author Michael C. Grundler
make.ratematrix.mk.epoch = function(layout, covarion=FALSE) {
    if (!inherits(layout, "array"))
        stop("Expected array argument")

    storage.mode(layout) = "integer"

    dims = dim(layout)

    if (dims[1L] != dims[2L])
        stop("Transition rate matrices must be square")

    pars = c(layout)
    pars = pars[pars > 0L]

    if (!length(pars))
        stop("Some transition rates must be non-zero")

    if (any(tabulate(pars) == 0L))
        stop("Transition rate parameter indices are not consecutive")

    nepoch = dims[3L]
    nstate = dims[1L]

    npar = length(unique(pars))

    # check to ensure that the layout is a valid covarion model
    if (covarion) {
        if (nstate %% 2L) {
            stop("Provided rate matrix layout is not a valid covarion model.
The number of states should be even.")
        }
        half = nstate / 2L
        for (i in 1L:nepoch) {
            bad = sum(layout[1L:half, 1L:half, i][upper.tri(layout[1L:half, 1L:half, i])])
            bad = bad + sum(layout[1L:half, 1L:half, i][lower.tri(layout[1L:half, 1L:half, i])])
            bad = bad + sum(layout[1L:half, (half+1L):nstate, i][upper.tri(layout[1L:half, (half+1L):nstate, i])])
            bad = bad + sum(layout[(half+1L):nstate, 1L:half, i][lower.tri(layout[(half+1L):nstate, 1L:half, i])])
            if (bad) {
                stop("Provided rate matrix layout is not a valid covarion model.
Ensure that all off-diagonal elements in upper left, upper right,
and lower left submatrices are zero.")
            }
        }
    }

    rate = .Call(mkepoch_rate_init, nepoch, npar, nstate, layout)

    if (!is.null(dimnames(layout)))
        dimnames(attr(rate, "layout")) = dimnames(layout)

    class(rate) = c("ratematrix.mk.epoch", "ratematrix.mk", "ratematrix")
    attr(rate, "covarion") = covarion

    return (rate)
}


print.ratematrix.mk.epoch = function(x, ...) {
    print(attr(x, "layout"))
}


#' Mk epoch model
#'
#' An Mk model that allows the evolution of a character in different
#' non-overlapping time periods (epochs) to evolve under different
#' transition rate matrices.
#'
#' @param phy An object of class \code{tree}.
#' @param x A vector of character states for tips in \code{phy}. May be a
#'  character vector or an integer vector.
#' @param rate An object of class \code{ratematrix.mk.epoch}. See documenatation
#'  for \code{make.ratematrix.mk.epoch}.
#' @param epochs An optional vector of epoch boundaries. Epoch boundaries are
#'  assumed to be specified as time units before present. If this argument is
#'  not missing it must be the case that the \code{rate} argument specifies a
#'  rate matrix for each epoch. The number of epochs is given by 1 plus the
#'  number of epoch boundaries.
#' @param levels. A vector of factor levels representing the different states
#'  for the character. If \code{x} is a character vector this argument must be
#'  supplied so that character states can be mapped to integers. In such a case,
#'  the integer value of a character state is its index position in this vector.
#'  If \code{x} is an integer vector this argument is unnecessary.
#' @param ambig. An optional list with names that match to character states
#'  in \code{x}. If present, each list item represents a mapping from the
#'  ambiguous coding state (represented by the name of the list item) to a
#'  set of non-ambiguous states (represented by the value of the list item).
#' @return A closure that inherits class \code{mk.epoch}, which is a function
#'  that computes the log likelihood of \code{x} given a vector of transtion
#'  rate parameters.
#' @author Michael C. Grundler
#' @examples
#'   library(macroevolution)
#'   data(snaketree)
#'
#'   phy = read.newick(text=snaketree)
#'
#'   # assume 2 character states, three epochs
#'   # assume a fully symmetric rate matrix for each
#'   # epoch, but allow for a different rate in each epoch
#'   layout = array(0, c(2, 2, 3))
#'   layout[, , 1] = 1
#'   layout[, , 1] = 2
#'   layout[, , 1] = 3
#'
#'   rate = make.ratematrix.mk.epoch(layout, covarion=FALSE)
#'
#'   # make up some tip states
#'   x = structure(c(rep(1, 10), rep(2, Ntip(phy)-10)), names=tiplabels(phy))
#'
#'   # epoch boundaries, which imply one epoch from 0 - 2 Ma, another from
#'   # 2 Ma - 8 Ma, and a final from 8 Ma - max(ages(phy)) Ma
#'   epochs = c(2, 8)
#'
#'   # quick visual of the "data" and epochs
#'   L = plot(phy)
#'   abline(v=max(ages(phy))-epochs, col=2, lty=2)
#'   points(L[[1]][1:10, 1], L[[1]][1:10, 3], pch=21, bg=8)
#'
#'   lik = make.mk.epoch(phy, x, rate, epochs)
#'
#'   # notice how param 1 is optimized to lower bound. this makes sense given
#'   # data can be explained by 2 transitions, one in each of the more ancient
#'   # epochs, and no transitions in the most recent epoch
#'   optim(runif(3, 0, .1), lik, method="L-BFGS-B", lower=1e-6, upper=1, control=list(fnscale=-1))
make.mk.epoch = function(phy, x, rate, epochs, levels, ambig) {
    stopifnot(is.tree(phy))
    stopifnot(tree.isbinary(phy))
    stopifnot(inherits(rate, "ratematrix.mk.epoch"))
    stopifnot(!is.null(names(x)))

    if (is.null(attr(rate, "covarion")))
        stop("Missing TRUE/FALSE value for rate matrix 'covarion' attribute")

    stopifnot(all(names(x) %in% tiplabels(phy))
        && all(tiplabels(phy) %in% names(x)))

    x = x[tiplabels(phy)]

    crown.age = max(ages(phy)) - ages(phy)
    stem.age = crown.age + brlens(phy)

    if (!missing(epochs)) {
        stopifnot(dim(rate)[3L] == (length(epochs)+1L))
        stopifnot(min(epochs) > min(ages(phy)))
        stopifnot(max(epochs) < max(ages(phy)))
        epochs = sort(epochs, decreasing=FALSE)
        epochs = c(min(ages(phy)), epochs, max(ages(phy)))
    } else {
        stopifnot(dim(rate)[3L] == 1L)
        epochs = c(min(ages(phy)), max(ages(phy)))
    }

    nstate = dim(attr(rate, "layout"))[1L]
    covarion = attr(rate, "covarion")
    half = nstate / 2L

    if (missing(levels) && inherits(x, "character"))
        stop("Must provide factor levels for tip state vector")

    if (inherits(x, "character")) {
        y = match(x, levels)
        if (anyNA(y)) {
            if (missing(ambig))
                stop("Cannot resolve ambiguous states")
        }
    } else {
        if (!covarion)
            y = match(x, 1L:nstate)
        else
            y = match(x, 1L:half)
        if (anyNA(y)) {
            if (missing(ambig))
                stop("Cannot resolve ambiguous states")
        }
    }

    storage.mode(y) = "integer"

    clk = matrix(0, nstate, Ntip(phy))

    for (i in 1L:Ntip(phy)) {
        if (is.na(y[i])) {
            if (is.na(match(y[i], names(ambig))))
                stop("Cannot resolve ambiguous states")
            z = ambig[[match(y[i], names(ambig))]]
            if (inherits(z, "character"))
                z = match(z, levels)
            if (anyNA(z))
                stop("Cannot resolve ambiguous states")
            clk[z, i] = 1
            if (covarion)
                clk[z+half, i] = 1
        } else {
            clk[y[i], i] = 1
            if (covarion)
                clk[y[i]+half, i] = 1
        }
    }

    model = .Call(mkepoch_model_init, epochs, clk, phy, rate,
        crown.age, stem.age)

    loglk = function(par) {
        .Call(mkepoch_loglk, par, model)
    }
    class(loglk) = c("mk.epoch", "mk")

    return (loglk)
}


print.mk.epoch = function(x, ...) {
    cat("Mk epoch model log likelihood function\n\n")
    cat("Epoch boundaries (time before present):\n")
    print(environment(x)$epochs)
    cat("Epoch rate matrices (parameter indices):\n")
    print(attr(environment(x)$rate, "layout"))
}


#' Marginal ancestral state reconstruction
#'
#' @param lik A likelihood function returned by \code{make.mk.epoch}
#' @return A function that computes marginal ancestral state estimates for
#'  the internal nodes of a phylogeny given a set of transition rate
#'  parameters.
#' @author Michael C. Grundler
make.asr.mk.epoch = function(lik) {
    stopifnot(inherits(lik, "mk.epoch"))

    model = environment(lik)$model

    asr = function(par) {
        .Call(mkepoch_marginal_asr, par, model)
    }

    return (asr)
}


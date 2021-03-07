#' Dirichlet-multinomial Markov process model (ultra-common mechanism)
#'
#' Bayesian MCMC implementation of the Dirichlet-multinomial Markov process
#' model described by Grundler and Rabosky (2020), https://doi.org/10.1093/sysbio/syaa031
#'
#' This function is deprecated in favor of \code{\link{make.rcm.dmm}}.
#'
#' @param phy An object of class \code{tree}.
#' @param x A three column \code{data.frame} containing count data. The first
#' column is expected to be the label of a terminal taxon; the second column, the
#' label of a resource category; the third column, the number of observations
#' for the particular combination of labels in columns one and two.
#' @param r The number of states in the prior model.
#' @param alpha.init An optional initial value for the Gamma prior hyperparameter on
#' branch lengths. If omitted defaults to a random value.
#' @param beta.init An optional initial value for the Dirichlet prior hyperparameter on
#' the multinomial distributions corresponding to each cluster. If omitted defaults
#' to the uniform Dirichlet distribution prior.
#' @param stateid.init An optional initial partition of terminal nodes into
#' states. If omitted all terminals are placed in a single cluster.
#' @return A function to perform MCMC inference that takes the following arguments
#' \describe{
#' \item{niter}{The number of iterations that the MCMC will run for. Default is 1000.}
#' \item{thin}{The number of iterations between each recorded sample. Default is 1.}
#' \item{tune.alpha}{Tuning parameter for the parameter of the Gamma prior. Should be between 0 and 1. Default is 0.1}
#' \item{tune.beta}{Tuning parameter for the parameter of the Dirichlet prior. Should be between 0 and 1. Default is 0.1}
#' \item{update.node}{Weight for proposing node state updates. Default is 1.}
#' \item{update.alpha}{Weight for proposing updates to the Gamma prior hyperparameter. Default is 1.}
#' \item{update.beta}{Weight for proposing updates to the Dirichlet prior hyperparameter. Default is 0.}
#' \item{output.file}{Name of file for recording the state of Markov chain. Default is "mcmc.out".}
#' \item{output.mode}{One of either "wb" or "ab", indicating whether the output file should be opened
#' for writing or appending. Default is "wb".}
#'}
#' The output from the MCMC function is written to a binary file that can be read
#' back into R with the \code{read.mk.dmm} function. For details on its structure
#' consult the documentation for \code{read.mk.dmm}.
#' @seealso \code{\link{read.mk.dmm}} for description of output file format,
#'   \code{\link{make.rcm.dmm}} for a version of the model that uses branch lengths.
#' @examples
#' \dontrun{
#'  data(snakediet)
#'  data(snaketree)
#'
#'  phy = read.newick(text=snaketree)
#'  mcmc = make.mk.dmm(phy, snakediet, 20L)
#'  mcmc(output.file="mcmc.out")
#'
#'  res = read.mk.dmm("mcmc.out")
#'
#'  # plot trace of data log likelihood
#'  plot(res$pars[, 1])
#'
#'  # plot trace of log prior probability
#'  plot(res$pars[, 2])
#'
#'  # plot trace of UCM hyperparameter
#'  plot(res$pars[, 3])
#'
#'  # look at the posterior distribution of the number of states among
#'  # terminal taxa
#'  barplot(table(apply(out$stateid, 1, function(p) length(unique(p)))))
#'
#'  # take a detailed look at a specific posterior sample, say the 250th
#'  obj = make.mk.dmm.from.sample(250, "mcmc.out")
#'
#'  str(obj)
#'
#'  # posterior mean estimate of each multinomial in this sample
#'  obj$dens.mean
#'
#'  # max a posteriori estimate of each multinomial in this sample
#'  obj$dens.map
#'
#'  # function to sample each multinomial in this sample from its posterior
#'  obj$dens.sampl
#'
#'  # marginal ancestral state probabilities. e.g. obj$asr[1, 1] is the marginal
#'  # posterior probability that the MRCA of all species in the phylogeny had a
#'  # pattern of resource use described by, e.g, the multinomial density in the
#'  # first row of obj$dens.map or obj$dens.mean.
#'  obj$asr
#'}
make.mk.dmm = function(phy, x, r, alpha.init, beta.init, stateid.init) {
    warning("make.mk.dmm is deprecated in favor of make.rcm.dmm")
    stopifnot(is.tree(phy))
    stopifnot(tree.isbinary(phy))
    stopifnot(is.data.frame(x))
    stopifnot(r >= 2L)

    x = aggregate(x[, 3L] ~ x[, 2L] + x[, 1L], FUN=sum)[, c(2L,1L,3L)]
    tips = match(x[, 1L], tiplabels(phy))
    x = x[!is.na(tips), ]
    rnames = as.character(unique(x[, 2L]))
    tips = match(x[, 1L], tiplabels(phy))
    cats = match(x[, 2L], rnames)
    cnts = as.integer(x[, 3L])
    x = cbind(tips, cats, cnts)
    x[, 1L] = x[, 1L] - 1L
    x[, 2L] = x[, 2L] - 1L
    storage.mode(x) = "integer"

    nrcls = as.integer(r)
    nrcat = length(rnames)

    f = (floor((Nnode(phy) - 1) / nrcls) + 1) / (Nnode(phy) - 1)
    alpha.max = -log((f*nrcls - 1) / (nrcls - 1)) / log((nrcls / (nrcls-1)) + 1)

    if (missing(alpha.init)) {
        alpha.init = runif(1, 0, alpha.max)
    } else {
        if (alpha.init > alpha.max)
            stop(sprintf("alpha.init should not exceed %f", alpha.max))
        stopifnot(alpha.init > 0)
    }

    if (missing(beta.init)) {
        beta.init = structure(rep(1, nrcat), names=rnames)
    } else {
        if (length(beta.init) == 1L)
            beta.init = structure(rep(beta.init, nrcat), names=rnames)
        else
            stopifnot(length(beta.init) == nrcat)

        if (is.null(names(beta.init))) {
            names(beta.init) = rnames
            warning("beta.init vector supplied without names")
        } else {
            stopifnot(all(names(beta.init) %in% rnames))
            stopifnot(all(rnames %in% names(beta.init)))
            beta.init = beta.init[rnames]
        }
        stopifnot(all(beta.init > 0))
    }

    if (missing(stateid.init)) {
        stateid.init = integer(Ntip(phy))
    } else {
        if (!is.null(names(stateid.init))) {
            stopifnot(all(names(stateid.init) %in% tiplabels(phy)))
            stopifnot(all(tiplabels(phy) %in% names(stateid.init)))
            stateid.init = stateid.init[tiplabels(phy)]
        }
        else {
            warning("ordering of initial state assignments assumed to match
                ordering of terminal nodes")
        }
        stopifnot(length(stateid.init) == Ntip(phy))
        if (length(unique(stateid.init)) > nrcls)
            stop("size of initial partition is too large")
    }

    model = .Call(
        mkdmm_model_init,
        phy,
        x,
        nrcat,
        nrcls,
        as.integer(stateid.init),
        as.double(alpha.init),
        as.double(beta.init))

    mcmc = function(niter=1000L, thin=1L, tune.alpha=0.1,
        tune.beta=0.1, update.node=1, update.alpha=1, update.beta=0,
        output.file="mcmc.out", output.mode=c("wb", "ab"))
    {
        if ((thin > 1) && ((thin %% 2) != 0))
            stop("Argument `thin` should be a power of 2")

        stopifnot(niter >= 0)
        stopifnot(tune.alpha > 0 && tune.alpha < 1)
        stopifnot(tune.beta > 0 && tune.beta < 1)
        stopifnot(update.node > 0)
        stopifnot(update.alpha >= 0)
        stopifnot(update.beta >= 0)

        output.mode = match.arg(output.mode)

        outputConn = file(output.file, output.mode)

        on.exit(close(outputConn))

        if (output.mode == "wb") {
            newick = write.newick(phy)
            writeBin(c(nrcat, nrcls, Ntip(phy), nrow(x), nchar(newick)), outputConn)
            writeBin(c("loglk", "prior", "alpha", rnames), outputConn)
            writeBin(tiplabels(phy), outputConn)
            writeBin(c(x), outputConn)
            # eos=NULL means newick string is not NUL terminated, in contrast
            # to other strings written by writeBin above
            writeChar(newick, eos=NULL, outputConn)
        }

        .Call(
            mkdmm_mcmc_run,
            model,
            as.integer(niter),
            as.integer(thin),
            as.double(tune.alpha),
            as.double(tune.beta),
            as.double(update.node),
            as.double(update.alpha),
            as.double(update.beta),
            outputConn)
    }

    return (mcmc)
}


read.header.mk.dmm = function(output.file) {
    con = file(output.file, "rb")
    on.exit(close(con))

    meta = readBin(con, integer(), 5L)
    nrcat = meta[1L]
    nrcls = meta[2L]
    ntip  = meta[3L]
    nnz = meta[4L]
    nc = meta[5L]

    ncols = 3 + nrcat
    header = readBin(con, character(), ncols)

    tip.names = readBin(con, character(), ntip)

    dataset = matrix(readBin(con, integer(), nnz * 3L), nnz, 3L)
    newick = readChar(con, nc)

    return (list(nrcat=nrcat, nrcls=nrcls, ntip=ntip,
        rnames=header[4L:ncols], dataset=dataset, newick=newick))
}


#' Read the output file from \code{make.mk.dmm} into R
#'
#'@details
#' The output file is a special binary file with the
#' following header format.
#'\tabular{rrl}{
#' Offset \tab Size  \tab  Description\cr
#' 0      \tab 4     \tab  Number of resource categories (= J)\cr
#' 4      \tab 4     \tab  Number of resource states (= K)\cr
#' 8      \tab 4     \tab  Number of terminal taxa (= N)\cr
#' 12     \tab 4     \tab  Number of rows in the input data matrix (= M)\cr
#' 16     \tab 4     \tab  Number of characters in the Newick string for the input phylogeny\cr
#' 20     \tab 6     \tab  The string "loglk"\cr
#' 26     \tab 6     \tab  The string "prior"\cr
#' 32     \tab 6     \tab  The string "alpha"\cr
#' 38     \tab X     \tab  An array of strings giving the names of the resource categories\cr
#' 38+X   \tab Y     \tab  An array of strings giving the names of the terminal taxa\cr
#' 38+X+Y \tab 12M   \tab  An array of integers giving the input dataset\cr
#' 38+X+Y+12M \tab Z \tab  The Newick string
#'}
#'
#' Note that byte offsets assume 1 byte chars, 4 byte integers, and 8 byte doubles.
#' After the header is the payload. The first entry in the payload is the initial state
#' configuration and hyperparameter values that were used to start the MCMC run.
#' The payload has the following repetitive structure, where W initially equals
#' 38+X+Y+12M+Z and is incremented by 24+8J+4N for each posterior sample.
#'
#'\tabular{rrl}{
#' Offset  \tab   Size \tab  Description\cr
#' W       \tab   8    \tab  The data log likelihood\cr
#' W+8     \tab   8    \tab  The log prior probability\cr
#' W+16    \tab   8    \tab  The hyperparameter of the UCM prior\cr
#' W+24    \tab   8J   \tab  The hyperparameters of the Dirichlet prior\cr
#' W+24+8J \tab   4N   \tab  The state assignments of the terminal taxa
#'}
#' @param output.file Name of the output file.
#' @param skip Number of rows in the file to skip over before reading. Defaults to 0.
#' @param n Number of rows to read. If -1 (default) all rows (after any skipped
#' rows) are read until the end of file is reached.
#' @return A list with the following structure
#' \describe{
#' \item{r}{The number of states in the UCM prior model.}
#' \item{pars}{A numeric matrix whose rows correspond to posterior samples. The
#' first column contains the log likelihood of the data conferred by the specific
#' configuration of state assignments in the sample. The second column contains
#' the prior probability of the state configuration under the UCM model. The
#' third column contains the hyperparameter of the UCM prior model. The subsequent
#' columns are the hyperparameters of the Dirichlet prior on the multinomial
#' distribution of each state.}
#' \item{stateid}{An integer matrix that records the index of the state sampled
#' for each terminal taxon in each posterior sample. Note that the indices of
#' states in one sample do not necessarily align with the indices in another.}}
#' \item{dataset}{The count dataset used for analysis. The first column holds the
#' terminal node indices, the second column holds the resource category indices,
#' and the third column holds the number of observations for the particular
#' combination of indices in columns one and two. Note that the indices in the
#' first two columns are 0-based, following C (not R) convention. Also, resource
#' category indices correspond to the column name ordering in the \code{pars}
#' list member. Thus, index 0 corresponds to the resource named in column 4 of
#' of \code{pars}; index 1, to column 5; and so on.}
#' \item{newick}{The Newick string of the phylogeny used for analysis.}
read.mk.dmm = function(output.file, skip=0, n=-1) {
    con = file(output.file, "rb")
    on.exit(close(con))

    meta = readBin(con, integer(), 5L)

    nrcat = meta[1L]
    nrcls = meta[2L]
    ntip  = meta[3L]
    nnz = meta[4L]
    nc = meta[5L]

    ncols = 3L + nrcat
    header = readBin(con, character(), ncols)

    tip.names = readBin(con, character(), ntip)

    dataset = matrix(readBin(con, integer(), nnz * 3L), nnz, 3L)
    newick = readChar(con, nc)

    header.sz = 20 + (3*nnz*4) + nc + sum(nchar(header, "bytes") + 1) + sum(nchar(tip.names, "bytes") + 1)
    file.sz = file.size(output.file) - header.sz
    line.sz = 8 * ncols + 4 * ntip

    header[4:ncols] = paste0("beta.", header[4:ncols])

    nlines = file.sz / line.sz

    if (skip > 0) {
        skip = min(skip, nlines)
        nlines = nlines - skip
        invisible(readChar(con, skip * line.sz, useBytes=TRUE))
    }

    if (n > 0)
        nlines = min(nlines, n)

    if (nlines) {
        pars = matrix(0, nlines, ncols, dimnames=list(NULL, header))
        stateid = matrix(0L, nlines, ntip, dimnames=list(NULL, tip.names))

        for (i in 1:nlines) {
            pars[i, ] = readBin(con, numeric(), ncols)
            stateid[i, ] = readBin(con, integer(), ntip)
        }

        return (list(r=nrcls, pars=pars, stateid=stateid,
            dataset=dataset, newick=newick))
    }

    return (NULL)
}


#' Create a summary object from a single posterior sample
#'
#' @param output.file The name of the output.file
#' @param i The index of the posterior sample to summarize
#' @return A list with following structure
#' \describe{
#'  \item{asr}{Marginal ancestral state probabilities.}
#'  \item{dens.mean}{Posterior mean estimate of each multinomial in the sample.}
#'  \item{dens.map}{Max a posteriori estimate of each multinomial in the sample.}
#'  \item{dens.sampl}{Function to sample from the posterior distribution of each
#'   multinomial in this sample.}
#'}
make.mk.dmm.from.sample = function(i, output.file) {

    rnames = read.header.mk.dmm(output.file)$rnames
    out = read.mk.dmm(output.file, skip=i, n=1)
    phy = read.newick(text=out$newick)

    model = .Call(
        mkdmm_model_init,
        phy,
        out$dataset,
        length(unique(out$dataset[, 2L])),
        out$r,
        as.integer(out$stateid[1L, ]),
        as.double(out$pars[1L, 3L]),
        as.double(out$pars[1L, -(1L:3L)]))

    # align tip state indices to rows of multinom.* objects
    tip.state = match(out$stateid[1L, ], unique(out$stateid[1L, ]))

    asr = .Call(mkdmm_marginal_asr, model)
    dirichlet.pars = .Call(mkdmm_posterior_multinomial, model)

    l = ncol(dirichlet.pars)
    k = nrow(dirichlet.pars)

    dens.mean = structure(
        sweep(dirichlet.pars, 1, rowSums(dirichlet.pars), "/")
        , dimnames=list(NULL, rnames))
    dens.map = structure(
        sweep(dirichlet.pars[-k, ] - 1, 1, rowSums(dirichlet.pars[-k, ] - 1), "/")
        , dimnames=list(NULL, rnames))

    dens.sampl = function(n, j) {
        if (j > k)
            stop("invalid state number")
        p = matrix(rgamma(n*l, dirichlet.pars[j, ]), ncol=l, byrow=TRUE,
            dimnames=list(NULL, rnames))
        return (sweep(p, 1, rowSums(p), "/"))
    }

    return (list(asr=asr, dens.map=dens.map, dens.mean=dens.mean
        , dens.sampl=dens.sampl))
}

#' Dirichlet-multinomial Markov process model (random cluster model)
#'
#' Bayesian MCMC implementation of the Dirichlet-multinomial Markov process
#' model described by Grundler and Rabosky (2020), https://doi.org/10.1093/sysbio/syaa031
#' that does not integrate out branch lengths.
#'
#' @param phy An object of class \code{tree}.
#' @param x A three column \code{data.frame} containing count data. The first
#' column is expected to be the label of a terminal taxon; the second column, the
#' label of a resource category; the third column, the number of observations
#' for the particular combination of labels in columns one and two.
#' @param r The number of states in the prior model.
#' @param stateid.init An optional initial partition of terminal nodes into
#' states. If omitted all terminals are placed in a single cluster.
#' @param rate.init An optional initial value for the rate of the Poisson process
#' prior on the evolution of states. If omitted defaults to a random value.
#' @param beta.init An optional initial value for the Dirichlet prior hyperparameter on
#' the multinomial distributions corresponding to each cluster. If omitted defaults
#' to the uniform Dirichlet distribution prior. This is never updated during the MCMC.
#' @return A function to perform MCMC inference that takes the following arguments
#' \describe{
#' \item{niter}{The number of iterations that the MCMC will run for. Default is 1000.}
#' \item{thin}{The number of iterations between each recorded sample. Default is 1.}
#' \item{tune.rate}{Tuning parameter for the rate. Should be between 0 and 1. Default is 0.1}
#' \item{update.node}{Weight for proposing node state updates. Default is 1.}
#' \item{update.rate}{Weight for proposing updates to the rate. Default is 1.}
#' \item{output.file}{Name of file for recording the state of Markov chain. Default is "mcmc.out".}
#' \item{output.mode}{One of either "wb" or "ab", indicating whether the output file should be opened
#' for writing or appending. Default is "wb".}
#'}
#' The output from the MCMC function is written to a binary file that can be read
#' back into R with the \code{read.rcm.dmm} function. For details on its structure
#' consult the documentation for \code{read.rcm.dmm}.
#' @details In the \code{make.mk.dmm} implementation, each branch is allowed to have its
#' own rate, and through a judicious choice of prior these are integrated out of
#' the model, resulting in a single set of transition probabilities that describe
#' the dynamics on all branches and that do not depend on the temporal durations
#' of branches. In this version, we assume instead that there is a single rate
#' across the entire phylogeny. In other words, the prior model for the evolution
#' of states is just the regular fully-symmetric, or equal-rates, model. In this
#' case, the temporal branch lengths do matter, and each branch has its own set
#' of transition probabilities.
#' @seealso \code{\link{read.rcm.dmm}} for description of output file format.
#' @examples
#' \dontrun{
#'  data(snakediet)
#'  data(snaketree)
#'
#'  phy = read.newick(text=snaketree)
#'  mcmc = make.rcm.dmm(phy, snakediet, 20L)
#'  mcmc(output.file="mcmc.out")
#'
#'  res = read.rcm.dmm("mcmc.out")
#'
#'  # plot trace of data log likelihood
#'  plot(res$pars[, 1])
#'
#'  # plot trace of log prior probability
#'  plot(res$pars[, 2])
#'
#'  # plot trace of rate hyperparameter
#'  plot(res$pars[, 3])
#'
#'  # look at the posterior distribution of the number of states among
#'  # terminal taxa
#'  barplot(table(apply(out$stateid, 1, function(p) length(unique(p)))))
#'
#'  # take a detailed look at a specific posterior sample, say the 250th
#'  obj = make.rcm.dmm.from.sample(250, "mcmc.out")
#'
#'  str(obj)
#'
#'  # posterior mean estimate of each multinomial in this sample
#'  obj$dens.mean
#'
#'  # max a posteriori estimate of each multinomial in this sample
#'  obj$dens.map
#'
#'  # marginal ancestral state probabilities. e.g. obj$asr[1, 1] is the marginal
#'  # posterior probability that the MRCA of all species in the phylogeny had a
#'  # pattern of resource use described by, e.g, the multinomial density in the
#'  # first row of obj$dens.map or obj$dens.mean.
#'  obj$asr
#'}
make.rcm.dmm = function(phy, x, r, stateid.init, rate.init, beta.init)
{
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

    r = as.integer(r)
    p = length(rnames)

    f = (floor((Nnode(phy) - 1) / r) + 1) / (Nnode(phy) - 1)
    rate.max = -log((f*r - 1) / (r - 1)) / (r * mean(brlens(phy)[-root(phy)]))

    if (missing(beta.init))
    {
        beta.init = structure(rep(1, p), names=rnames)
    }
    else {
        if (length(beta.init) == 1L)
            beta.init = structure(rep(beta.init, p), names=rnames)
        else
            stopifnot(length(beta.init) == p)

        if (is.null(names(beta.init)))
        {
            names(beta.init) = rnames
            warning("beta.init vector supplied without names")
        }
        else
        {
            stopifnot(all(names(beta.init) %in% rnames))
            stopifnot(all(rnames %in% names(beta.init)))
            beta.init = beta.init[rnames]
        }
        stopifnot(all(beta.init > 0))
    }

    if (missing(rate.init))
    {
        rate.init = runif(1, 0, rate.max)
    }
    else
    {
        storage.mode(rate.init) = "double"
        if (rate.init > rate.max)
            stop(sprintf("rate.init should not exceed %f", rate.max))
        stopifnot(rate.init > 0)
    }

    if (missing(stateid.init))
    {
        stateid.init = integer(Ntip(phy))
    }
    else
    {
        if (!is.null(names(stateid.init)))
        {
            stopifnot(all(names(stateid.init) %in% tiplabels(phy)))
            stopifnot(all(tiplabels(phy) %in% names(stateid.init)))
            stateid.init = stateid.init[tiplabels(phy)]
        }
        else
        {
            warning("ordering of initial state assignments assumed to match
                ordering of terminal nodes")
        }
        stopifnot(length(stateid.init) == Ntip(phy))
        if (length(unique(stateid.init)) > r)
            stop("size of initial partition is too large")
    }

    model = .Call(
        rcm_dmm_model_init,
        phy,
        x,
        p,
        r,
        rate.init,
        beta.init,
        as.integer(stateid.init))

    mcmc = function(niter=1000L, thin=1L, update.node=1, update.rate=1,
        tune.rate=0.1, output.file="mcmc.out", output.mode=c("wb", "ab"))
    {
        if ((thin > 1) && ((thin %% 2) != 0))
            stop("Argument `thin` should be a power of 2")

        stopifnot(niter >= 0)

        output.mode = match.arg(output.mode)

        outputConn = file(output.file, output.mode)

        on.exit(close(outputConn))

        if (output.mode == "wb")
        {
            newick = write.newick(phy)
            writeBin(c(p, r, nrow(x), Ntip(phy), nchar(newick)), outputConn)
            writeBin(tiplabels(phy), outputConn)
            writeBin(rnames, outputConn)
            writeBin(c(x), outputConn)
            # eos=NULL means newick string is not NUL terminated, in contrast
            # to other strings written by writeBin above
            writeChar(newick, eos=NULL, outputConn)
            writeBin(c(beta.init), outputConn)
        }

        .Call(
            rcm_model_mcmc_run,
            model,
            as.integer(niter),
            as.integer(thin),
            as.double(update.node),
            as.double(update.rate),
            as.double(tune.rate),
            outputConn)
    }

    return (mcmc)
}


#' Read the output file from \code{make.rcm.dmm} into R
#'
#'@details
#' The output file is a special binary file with the
#' following header format.
#'\tabular{rrl}{
#' Offset \tab Size  \tab  Description\cr
#' 0      \tab 4     \tab  Number of resource categories (= J)\cr
#' 4      \tab 4     \tab  Number of resource states (= K)\cr
#' 8      \tab 4     \tab  Number of rows in the input data matrix (= M)\cr
#' 12     \tab 4     \tab  Number of terminal taxa (= N)\cr
#' 16     \tab 4     \tab  Number of characters in the Newick string for the input phylogeny\cr
#' 20     \tab X     \tab  An array of strings giving the names of the terminal taxa\cr
#' 20+X   \tab Y     \tab  An array of strings giving the names of the resource categories\cr
#' 20+X+Y \tab 12M   \tab  An array of integers giving the input dataset\cr
#' 20+X+Y+12M \tab Z \tab  The Newick string\cr
#' 20+X+Y+12M+Z \tab 8J \tab The hyperparameters of the Dirichlet prior
#'}
#'
#' Note that byte offsets assume 1 byte chars, 4 byte integers, and 8 byte doubles.
#' After the header is the payload. The first entry in the payload is the initial state
#' configuration and hyperparameter values that were used to start the MCMC run.
#' The payload has the following repetitive structure, where W initially equals
#' 20+X+Y+12M+Z+8J and is incremented by 24+4N for each posterior sample.
#'
#'\tabular{rrl}{
#' Offset  \tab   Size \tab  Description\cr
#' W       \tab   8    \tab  The data log likelihood\cr
#' W+8     \tab   8    \tab  The log prior probability\cr
#' W+16    \tab   8    \tab  The rate of the Poisson prior\cr
#' W+24    \tab   4N   \tab  The state assignments of the terminal taxa
#'}
#' @param output.file Name of the output file.
#' @param skip Number of rows in the file to skip over before reading. Defaults to 0.
#' @param n Number of rows to read. If -1 (default) all rows (after any skipped
#' rows) are read until the end of file is reached.
#' @return A list with the following structure
#' \describe{
#' \item{r}{The number of states in the prior model.}
#' \item{pars}{A numeric matrix whose rows correspond to posterior samples. The
#' first column contains the log likelihood of the data conferred by the specific
#' configuration of state assignments in the sample. The second column contains
#' the prior probability of the state configuration under the Markov model. The
#' third column contains the rate of the Markov model.}
#' \item{stateid}{An integer matrix that records the index of the state sampled
#' for each terminal taxon in each posterior sample. Note that the indices of
#' states in one sample do not necessarily align with the indices in another.}}
#' \item{dirichlet.prior}{The hyperparameters of the Dirichlet prior.}
#' \item{dataset}{The count dataset used for analysis. The first column holds the
#' terminal node indices, the second column holds the resource category indices,
#' and the third column holds the number of observations for the particular
#' combination of indices in columns one and two. Note that the indices in the
#' first two columns are 0-based, following C (not R) convention. Also, resource
#' category indices correspond to the name ordering in the \code{dirichlet.prior}
#' list member. Thus, index 0 corresponds to the first named element of
#' of \code{dirichlet.prior}; index 1, to the second; and so on.}
#' \item{phy}{The phylogeny used for analysis.}
read.rcm.dmm = function(output.file, skip=0, n=-1)
{
    con = file(output.file, "rb")
    on.exit(close(con))

    meta = readBin(con, integer(), 5L)

    p = meta[1L]
    r = meta[2L]
    nr = meta[3L]
    ntip  = meta[4L]
    nc = meta[5L]

    tip.names = readBin(con, character(), ntip)
    rnames = readBin(con, character(), p)

    dataset = matrix(readBin(con, integer(), 3L * nr), nr, 3L)
    newick = readChar(con, nc)

    # Dirichlet hyperparameters for each state
    dprior = structure(readBin(con, double(), p), names=rnames)

    header.sz = 20 + (nr*3*4) + nc + sum(nchar(tip.names, "bytes") + 1) +
        p * 8 + sum(nchar(rnames, "bytes") + 1)
    file.sz = file.size(output.file) - header.sz
    line.sz = 8 * 3 + 4 * ntip

    nlines = file.sz / line.sz

    phy = read.newick(text=newick)

    #dataset = data.frame(predator=tiplabels(phy)[dataset[, 1] + 1],
    #    prey=rnames[dataset[, 2] + 1], count=dataset[, 3], stringsAsFactors=FALSE)

    if (skip > 0)
    {
        skip = min(skip, nlines)
        nlines = nlines - skip
        invisible(readChar(con, skip * line.sz, useBytes=TRUE))
    }

    if (n > 0)
        nlines = min(nlines, n)

    if (nlines)
    {
        pars = matrix(0, nlines, 3, dimnames=list(NULL, c("loglk", "logprior",
            "rate")))
        stateid = matrix(0L, nlines, ntip, dimnames=list(NULL, tip.names))

        for (i in 1:nlines)
        {
            pars[i, ] = readBin(con, numeric(), 3)
            stateid[i, ] = readBin(con, integer(), ntip)
        }

        return (list(r=r, pars=pars, stateid=stateid,
            dirichlet.prior=dprior, dataset=dataset, phy=phy))
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
#'  \item{smap}{A function to perform stochastic character mapping.}
#'  \item{dens.map}{Max a posteriori estimate of each multinomial in the sample.}
#'  \item{dens.mean}{Posterior mean estimate of each multinomial in the sample.}
#'  \item{tip.state}{Terminal node state assignments. These align with the row
#'  indices of \code{dens.map} and \code{dens.mean}.}
#'}
make.rcm.dmm.from.sample = function(i, output.file)
{
    out = read.rcm.dmm(output.file, skip=i, n=1)
    phy = out$phy
    dataset = out$dataset

    model = .Call(
        rcm_dmm_model_init,
        phy,
        dataset,
        length(out$dirichlet.prior),
        out$r,
        out$pars[1, 3],
        out$dirichlet.prior,
        out$stateid[1, ])

    # align tip state indices to rows of dens.* objects
    tip.state = match(out$stateid[1L, ], unique(out$stateid[1L, ]))

    asr = .Call(rcm_marginal_asr, model)
    dirichlet.pars = .Call(rcm_dmm_posterior_multinomial, model)

    dens.mean = structure(
        sweep(dirichlet.pars, 1, rowSums(dirichlet.pars), "/")
        , dimnames=list(NULL, names(out$dirichlet.prior)))
    dens.map = structure(
        sweep(dirichlet.pars - 1, 1, rowSums(dirichlet.pars - 1), "/")
        , dimnames=list(NULL, names(out$dirichlet.prior)))

    smap = function(n)
    {
        stopifnot(as.integer(n) > 0)
        return (.Call(rcm_stochastic_map, as.integer(n), model))
    }

    expected.counts = structure(t(.Call(rcm_stochastic_map_expected_counts, model))
        , dimnames=list(NULL, c("posterior", "prior")))

    return (list(asr=asr, smap=smap, dens.map=dens.map, dens.mean=dens.mean,
        tip.state=tip.state, expected.counts=expected.counts))
}

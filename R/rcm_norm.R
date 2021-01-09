#' Normal-inverse Wishart Markov process model (random cluster model)
#'
#'
#' @param phy An object of class \code{tree}.
#' @param x A \code{data.frame} containing continuous character data. Columns
#' correspond to characters and rows to species or specimens.
#' @param r The number of states in the prior model.
#' @param stateid.init An optional initial partition of terminal nodes into
#' states. If omitted all terminals are placed in a single cluster.
#' @param rate.init An optional initial value for the rate of the Poisson process
#' prior on the evolution of states. If omitted defaults to a random value.
#' @param kappa.init An optional initial value for the kappa hyperparameter of
#' the NIV prior on the normal distributions corresponding to each cluster.
#' If omitted defaults to 1. This is never updated during the MCMC.
#' @param nu.init An optional initial value for the nu hyperparameter of
#' the NIV prior on the normal distributions corresponding to each cluster.
#' If omitted defaults to 1. This is never updated during the MCMC.
#' @param mu.init An optional initial value for the mean vector hyperparameter of
#' the NIV prior on the normal distributions corresponding to each cluster.
#' If omitted defaults to the global mean vector. This is never updated during the MCMC.
#' @param lambda.init An optional initial value for the covariance matrix hyperparameter of
#' the NIV prior on the normal distributions corresponding to each cluster.
#' If omitted defaults to a diagonal matrix of sample variances.
#' This is never updated during the MCMC.
#' @param integrate.brlen An optional boolean value indicating whether or not to
#' assume branches have i.i.d. rates drawn from a Gamma prior. If TRUE, the branch
#' rates are integrated out of the likelihood function such that a single set of
#' transition probabilities applies to all branches. The default is FALSE, in
#' which case all branches share a single rate.
#' @param species.col The column index containing the species labels. The
#' default is 1.
#' @param specimen.col The column index containing the specimen labels. The
#' default is 0, meaning that character values correspond to species means.
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
#' @seealso \code{\link{read.rcm.dmm}} for description of output file format.
make.rcm.norm = function(phy, x, r, stateid.init, rate.init, kappa.init,
    nu.init, mu.init, lambda.init, integrate.brlen=FALSE, species.col=1L,
    specimen.col=0L)
{
    stopifnot(is.tree(phy))
    stopifnot(tree.isbinary(phy))
    stopifnot(is.data.frame(x))
    stopifnot(r >= 1L)

    species.col = species.col[1]
    specimen.col = specimen.col[1]

    y = x[, -c(species.col, specimen.col)]
    cnames = colnames(y)

    normalize.row = function(p) {
        sp = p[1, species.col]
        if (specimen.col > 0)
            specimen = p[1, specimen.col]
        else
            specimen = sprintf("%s1", sp)
        vals = as.numeric(p[-c(species.col, specimen.col)])
        idx = which(!is.na(vals))
        n = length(idx)
        data.frame(
            species=rep(sp, n),
            specimen=rep(specimen, n),
            trait=cnames[idx],
            value=vals[idx])
    }

    z = vector("list", nrow(x))
    for (i in 1:nrow(x))
        z[[i]] = normalize.row(x[i, ])

    x = do.call(rbind, z)

    tips = match(x[, 1L], tiplabels(phy))
    x = x[!is.na(tips), ]
    tips = as.numeric(match(x[, 1L], tiplabels(phy)))
    inds = as.numeric(factor(x[, 2L]))
    chars = as.numeric(match(x[, 3L], cnames))
    vals = as.numeric(x[, 4L])
    x = cbind(tips, inds, chars, vals)
    x[, 1L] = x[, 1L] - 1
    x[, 3L] = x[, 3L] - 1

    r = as.integer(r)
    p = length(cnames)

    if (missing(integrate.brlen))
    {
        integrate.brlen = 0L
    }
    else
    {
        stopifnot(inherits(integrate.brlen, "logical"))
        if (integrate.brlen[1])
            integrate.brlen = 1L
        else
            integrate.brlen = 0L
    }

    f = (floor((Nnode(phy) - 1) / r) + 1) / (Nnode(phy) - 1)
    if (!integrate.brlen)
        rate.max = -log((f*r - 1) / (r - 1)) / (r * mean(brlens(phy)[-root(phy)]))
    else
        rate.max = -log((f*r - 1) / (r - 1)) / log((r/(r - 1)) + 1)

    if (is.nan(rate.max))
        rate.max = 0

    if (missing(kappa.init)) {
        kappa.init = 1
    } else {
        kappa.init = kappa.init[1]
        storage.mode(kappa.init) = "double"
        stopifnot(kappa.init > 0)
    }
    if (missing(nu.init)) {
        nu.init = p+1
    } else {
        nu.init = nu.init[1]
        storage.mode(nu.init) = "double"
        stopifnot(nu.init >= p)
    }
    if (missing(mu.init)) {
        mu.init = colMeans(y, na.rm=TRUE)
    } else {
        storage.mode(mu.init) = "double"
        stopifnot(length(mu.init) == p)
    }
    if (missing(lambda.init)) {
        lambda.init = diag(apply(y, 2, var, na.rm=TRUE), p, p)
    } else {
        storage.mode(lambda.init) = "double"
        stopifnot(isSymmetric(lambda.init))
        stopifnot(nrow(lambda.init) == p)
        stopifnot(all(diag(lambda.init) > 0))
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
        rcm_norm_model_init,
        phy,
        x,
        p,
        r,
        integrate.brlen,
        rate.init,
        nu.init,
        kappa.init,
        mu.init,
        lambda.init,
        as.integer(stateid.init))

    mcmc = function(niter=1000L, thin=1L, update.node=1, update.rate=1,
        tune.rate=0.1, output.file="mcmc.out", output.mode=c("wb", "ab"))
    {
        if ((thin > 1) && ((thin %% 2) != 0))
            stop("Argument `thin` should be a power of 2")

        stopifnot(niter >= 0)

        if (r == 1)
            niter = 0

        output.mode = match.arg(output.mode)

        outputConn = file(output.file, output.mode)

        on.exit(close(outputConn))

        if (output.mode == "wb")
        {
            newick = write.newick(phy)
            writeBin(c(p, r, integrate.brlen, nrow(x), Ntip(phy),
                nchar(newick)), outputConn)
            writeBin(tiplabels(phy), outputConn)
            writeBin(cnames, outputConn)
            writeBin(c(x), outputConn)
            # eos=NULL means newick string is not NUL terminated, in contrast
            # to other strings written by writeBin above
            writeChar(newick, eos=NULL, outputConn)
            writeBin(c(kappa.init,nu.init,mu.init,lambda.init), outputConn)
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


#' Read the output file from \code{make.rcm.norm} into R
#'
#'@details
#' The output file is a special binary file with the
#' following header format.
#'\tabular{rrl}{
#' Offset \tab Size  \tab  Description\cr
#' 0      \tab 4     \tab  Number of quantitative characters (= J)\cr
#' 4      \tab 4     \tab  Number of resource states (= K)\cr
#' 8      \tab 4     \tab  Flag indicating if branch lengths are integrated out\cr
#' 12     \tab 4     \tab  Number of rows in the input data matrix (= M)\cr
#' 16     \tab 4     \tab  Number of terminal taxa (= N)\cr
#' 20     \tab 4     \tab  Number of characters in the Newick string for the input phylogeny\cr
#' 24     \tab X     \tab  An array of strings giving the names of the terminal taxa\cr
#' 24+X   \tab Y     \tab  An array of strings giving the names of the characters\cr
#' 24+X+Y \tab 32M   \tab  An array of doubles giving the input dataset\cr
#' 24+X+Y+32M \tab Z \tab  The Newick string\cr
#' 24+X+Y+32M+Z \tab 8(2+J+J*J) \tab The hyperparameters of the Normal-inverse Wishart prior
#'}
#'
#' Note that byte offsets assume 1 byte chars, 4 byte integers, and 8 byte doubles.
#' After the header is the payload. The first entry in the payload is the initial state
#' configuration and hyperparameter values that were used to start the MCMC run.
#' The payload has the following repetitive structure, where W initially equals
#' 24+X+Y+32M+Z+8(2+J+J*J) and is incremented by 24+4N for each posterior sample.
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
#' \item{integrate.brlen}{Were branch lengths integrated out.}
#' \item{pars}{A numeric matrix whose rows correspond to posterior samples. The
#' first column contains the log likelihood of the data conferred by the specific
#' configuration of state assignments in the sample. The second column contains
#' the prior probability of the state configuration under the Markov model. The
#' third column contains the rate of the Markov model.}
#' \item{stateid}{An integer matrix that records the index of the state sampled
#' for each terminal taxon in each posterior sample. Note that the indices of
#' states in one sample do not necessarily align with the indices in another.}}
#' \item{kappa.init}{The hyperparameters of the Normal-inverse Wishart prior.}
#' \item{nu.init}{The hyperparameters of the Normal-inverse Wishart prior.}
#' \item{mu.init}{The hyperparameters of the Normal-inverse Wishart prior.}
#' \item{lambda.init}{The hyperparameters of the Normal-inverse Wishart prior.}
#' \item{dataset}{The dataset used for analysis. The first column holds the
#' terminal node indices, the second colum, the specimen labels,
#' the third column holds the character indices,
#' and the fourth column holds the measurment for the particular
#' combination of indices in columns one, two, and three. Note that the indices in the
#' first and third columns are 0-based, following C (not R) convention. Also, resource
#' category indices correspond to the name ordering in the \code{mu.init}
#' list member. Thus, index 0 corresponds to the first named element of
#' of \code{mu.init}; index 1, to the second; and so on.}
#' \item{phy}{The phylogeny used for analysis.}
read.rcm.norm = function(output.file, skip=0, n=-1)
{
    con = file(output.file, "rb")
    on.exit(close(con))

    meta = readBin(con, integer(), 6L)

    p = meta[1L]
    r = meta[2L]
    integrate.brlen = meta[3L]
    nr = meta[4L]
    ntip  = meta[5L]
    nc = meta[6L]

    tip.names = readBin(con, character(), ntip)
    cnames = readBin(con, character(), p)

    dataset = matrix(readBin(con, double(), 4L * nr), nr, 4L)
    newick = readChar(con, nc)

    # NIV hyperparameters for each state
    kappa = readBin(con, double(), 1)
    nu = readBin(con, double(), 1)
    mu = structure(readBin(con, double(), p), names=cnames)
    lambda = matrix(readBin(con, double(), p*p), p, p, dimnames=list(cnames, cnames))

    header.sz = 24 + (nr*4*8) + nc + sum(nchar(tip.names, "bytes") + 1) +
        (2 + p + p*p) * 8 + sum(nchar(cnames, "bytes") + 1)
    file.sz = file.size(output.file) - header.sz
    line.sz = 8 * 3 + 4 * ntip

    nlines = file.sz / line.sz

    phy = read.newick(text=newick)

    if (skip > 0)
    {
        skip = min(skip, nlines)
        nlines = nlines - skip
        #invisible(readChar(con, skip * line.sz, useBytes=TRUE))
        invisible(readBin(con, "raw", skip * line.sz))
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

        return (list(r=r, integrate.brlen=integrate.brlen, pars=pars,
            stateid=stateid, kappa.init=kappa, nu.init=nu, mu.init=mu,
            lambda.init=lambda, dataset=dataset, phy=phy))
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
#'  \item{tip.state}{Terminal node state assignments. These align with the row
#'  indices of \code{dens.map} and \code{dens.mean}.}
#'  \item{niv.pars}{Parameters of the posterior distribution for each state.}
#'  \item{sample}{Function to sample from each state's posterior.}
#'}
make.rcm.norm.from.sample = function(i, output.file)
{
    out = read.rcm.norm(output.file, skip=i-1, n=1)
    if (!is.null(out)) {
        phy = out$phy
        dataset = out$dataset

        model = .Call(
            rcm_norm_model_init,
            phy,
            dataset,
            length(out$mu.init),
            out$r,
            out$integrate.brlen,
            out$pars[1, 3],
            out$nu.init,
            out$kappa.init,
            out$mu.init,
            out$lambda.init,
            out$stateid[1, ])

        # align tip state indices to rows of dens.* objects
        tip.state = match(out$stateid[1L, ], unique(out$stateid[1L, ]))

        asr = .Call(rcm_marginal_asr, model)
        niv.pars = .Call(rcm_norm_posterior_normal, model)

        smap = function(n)
        {
            stopifnot(as.integer(n) > 0)
            return (.Call(rcm_stochastic_map, as.integer(n), model))
        }

        sample = function(n, state) {
            p = length(out$mu.init)
            kappa = niv.pars[state, 1]
            nu = niv.pars[state, 2]
            mu = niv.pars[state, 3:(2+p)]
            lambda = matrix(niv.pars[state, -(1:(2+p))], p, p)
            Sigma = solve(rWishart(1, nu, solve(lambda))[,,1])
            xbar = MASS::mvrnorm(1, mu, Sigma/kappa)
            return (MASS::mvrnorm(n, xbar, Sigma))
        }

        return (list(asr=asr, smap=smap, tip.state=tip.state,
            niv.pars=niv.pars, sample=sample))
    }
    return (NULL)
}


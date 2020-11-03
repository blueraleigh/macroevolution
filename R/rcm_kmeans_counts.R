make.rcm.kmeans.counts = function(phy, x, r, stateid.init, rate.init,
    integrate.brlen)
{
    stopifnot(is.tree(phy))
    stopifnot(tree.isbinary(phy))
    stopifnot(is.data.frame(x))
    stopifnot(r >= 2L)

    x = aggregate(x[, 3L] ~ x[, 2L] + x[, 1L], FUN=sum)[, c(2L,1L,3L)]
    tips = match(x[, 1L], tiplabels(phy))
    x = x[!is.na(tips), ]
    x = do.call(function(...) rbind(..., make.row.names=FALSE),
            split(x, x[, 1L])[tiplabels(phy)])
    counts = xtabs(x[, 3] ~ x[, 1] + x[, 2])
    counts = counts[tiplabels(phy), ]
    rnames = colnames(counts)
    tips = match(x[, 1L], tiplabels(phy))
    cats = match(x[, 2L], rnames)
    cnts = as.integer(x[, 3L])
    x = cbind(tips, cats, cnts)
    storage.mode(x) = "integer"

    r = as.integer(r)
    p = length(rnames)

    f = (floor((Nnode(phy) - 1) / r) + 1) / (Nnode(phy) - 1)
    rate.max = -log((f*r - 1) / (r - 1)) / (r * mean(brlens(phy)[-root(phy)]))

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
        if (length(unique(stateid.init)) < r)
            stop("size of initial partition is too small")
    }

    if (missing(integrate.brlen))
    {
        integrate.brlen = 1L
    }
    else
    {
        if (integrate.brlen)
            integrate.brlen = 1L
        else
            integrate.brlen = 0L
    }

    model = .Call(
        rcm_kmeans_counts_model_init,
        phy,
        counts,
        p,
        r,
        integrate.brlen,
        rate.init,
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
            writeBin(c(p, r, integrate.brlen, nrow(x), Ntip(phy),
                nchar(newick)), outputConn)
            writeBin(tiplabels(phy), outputConn)
            writeBin(rnames, outputConn)
            writeBin(c(x), outputConn)
            # eos=NULL means newick string is not NUL terminated, in contrast
            # to other strings written by writeBin above
            writeChar(newick, eos=NULL, outputConn)
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


read.rcm.kmeans.counts = function(output.file, skip=0, n=-1)
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
    rnames = readBin(con, character(), p)

    dataset = matrix(readBin(con, integer(), 3L * nr), nr, 3L)
    newick = readChar(con, nc)

    header.sz = 24 + (nr*3*4) + nc + sum(nchar(tip.names, "bytes") + 1) +
        sum(nchar(rnames, "bytes") + 1)
    file.sz = file.size(output.file) - header.sz
    line.sz = 8 * 3 + 4 * ntip

    nlines = file.sz / line.sz

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

        return (list(r=r, integrate.brlen=integrate.brlen, pars=pars,
            stateid=stateid, dataset=dataset, newick=newick))
    }

    return (NULL)
}


make.rcm.kmeans.counts.from.sample = function(i, output.file)
{
    out = read.rcm.kmeans.counts(output.file, skip=i, n=1)
    phy = read.newick(text=out$newick)

    counts = xtabs(out$dataset[, 3] ~ out$dataset[, 1] + out$dataset[, 2])

    model = .Call(
        rcm_kmeans_counts_model_init,
        phy,
        counts,
        ncol(counts),
        out$r,
        out$integrate.brlen,
        out$pars[1, 3],
        out$stateid[1, ])

    asr = .Call(rcm_marginal_asr, model)

    smap = function(n, prior=FALSE)
    {
        stopifnot(as.integer(n) > 0)
        return (.Call(rcm_stochastic_map, as.integer(n), model))
    }

    expected.counts = structure(t(.Call(rcm_stochastic_map_expected_counts, model))
        , dimnames=list(NULL, c("posterior", "prior")))

    return (list(asr=asr, smap=smap, tip.state=out$stateid[1,],
        expected.counts=expected.counts))
}

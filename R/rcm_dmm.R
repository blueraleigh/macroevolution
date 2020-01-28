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

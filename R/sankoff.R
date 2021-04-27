prepare.sankoff = function(phy, data, levels, ambig) {
    I = rep(1L, Ntip(phy))
    stopifnot(!is.null(names(data)))
    if (missing(levels)) {
        stopifnot(class(data) == "integer")
        downpass = xtabs(I ~ names(data) + data)
    } else {
        if (missing(ambig)) {
            data = factor(data, levels=levels)
            downpass = xtabs(I ~ names(data) + data)
        } else {
            stopifnot(!is.null(names(ambig)))
            stopifnot(!any(levels %in% names(ambig)))
            poly = which(data %in% names(ambig))
            idata = factor(data[-poly], levels=levels)
            I = I[-poly]
            downpass = xtabs(I ~ names(idata) + idata)
            for (i in poly) {
                d = numeric(length(levels))
                d[ambig[[data[poly[i]]]]] = 1
                downpass = rbind(downpass, matrix(d, nrow=1, dimnames=list(
                    names(data)[i], NULL)))
            }
        }
    }
    downpass[downpass == 0] = Inf
    downpass[downpass == 1] = 0
    downpass = downpass[tiplabels(phy), , drop=FALSE]
    downpass = rbind(downpass, matrix(0, Nnode(phy)-Ntip(phy), ncol(downpass)))

    return (downpass)
}



# Perform a maximum parsimony reconstruction of ancestral
# states given an r-state character using Sankoff's algorithm
mpr.sankoff = function(phy, data, levels, ambig, cost) {
    stopifnot(is.tree(phy))

    g = prepare.sankoff(phy, data, levels, ambig)
    r = ncol(g)
    h = matrix(0, Nnode(phy), r)
    f = matrix(0, Nnode(phy), r)

    if (missing(cost)) {
        cost = matrix(1, r, r)
        diag(cost) = 0
    }
    storage.mode(cost) = "double"

    .Call(do_sankoff_downpass, phy, r, g, h, cost)
    .Call(do_sankoff_uppass, phy, r, g, h, f, cost)
    n = .Call(do_sankoff_count, phy, r, g, h, f, cost)

    sample = function(n, only.mpr=TRUE) {
        node.state = matrix(0L, Nnode(phy), n)
        .Call(do_sankoff_sample, phy, r, g, h, f, cost,
            as.integer(n), node.state, as.integer(only.mpr))
        return (node.state)
    }

    return (list(g=g,h=f,f=f,n=n, sample=sample))
}

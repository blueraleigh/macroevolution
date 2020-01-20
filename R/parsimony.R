prepare.downpass = function(phy, data, levels, ambig) {
    I = rep(1L, Ntip(phy))
    stopifnot(!is.null(names(data)))
    if (missing(levels))
    {
        stopifnot(class(data) == "integer")
        downpass = xtabs(I ~ data + names(data))
    }
    else
    {
        if (missing(ambig))
        {
            data = factor(data, levels=levels)
            downpass = xtabs(I ~ data + names(data))
        }
        else
        {
            poly = which(data %in% ambig)
            idata = factor(data[-poly], levels=levels)
            I = I[-poly]
            downpass = xtabs(I ~ idata + names(data))
            downpass = cbind(downpass, matrix(1L, length(levels), length(poly),
                dimnames=list(NULL, names(data)[poly])))
        }
    }

    downpass = downpass[, tiplabels(phy), drop=FALSE]
    downpass = cbind(downpass, matrix(0L, nrow(downpass), Nnode(phy)-Ntip(phy)))

    return (downpass)
}


# Calculate the Fitch parsimony score given an
# r-state character.
pscore.fitch = function(phy, data, levels, ambig, intermediates=FALSE) {
    stopifnot(is.tree(phy))
    stopifnot(tree.isbinary(phy))

    downpass = prepare.downpass(phy, data, levels, ambig)
    pscore = integer(Nnode(phy))

    .Call(do_fitch_pscore, phy, downpass, pscore)

    if (intermediates)
        return (pscore[-(1L:Ntip(phy))])

    return (pscore[root(phy)])
}


# Perform a maximum parsimony reconstruction of ancestral
# states given an r-state character using Fitch's algorithm
mpr.fitch = function(phy, data, levels, ambig) {
    stopifnot(is.tree(phy))
    stopifnot(tree.isbinary(phy))

    downpass = prepare.downpass(phy, data, levels, ambig)
    uppass = matrix(0L, nrow(downpass), Nnode(phy))

    pscore = integer(Nnode(phy))

    .Call(do_fitch_mpr, phy, uppass, downpass, pscore)

    pscore = pscore[root(phy)]

    # Sample a history of character evolution from a
    # maximum parsimony reconstruction of ancestral states
    sample.mpr = function(n) {
        storage.mode(n) = "integer"

        node.state = matrix(0L, Nnode(phy), n)
        node.state[1L:Ntip(phy), ] = data[tiplabels(phy)]

        .Call(do_fitch_history, phy, node.state, uppass, downpass, n)

        return (node.state)
    }

    return (sample.mpr)
}


# Count the number of MPR reconstructions
count.mpr.fitch = function(phy, data, levels, ambig) {
    stopifnot(is.tree(phy))
    stopifnot(tree.isbinary(phy))
    storage.mode(data) = "integer"
    stopifnot(!is.null(names(data)))

    downpass = prepare.downpass(phy, data, levels, ambig)
    uppass = matrix(0L, nrow(downpass), Nnode(phy))
    pscore = integer(Nnode(phy))

    .Call(do_fitch_mpr, phy, uppass, downpass, pscore)

    return (.Call(do_fitch_count, phy, uppass, downpass))
}

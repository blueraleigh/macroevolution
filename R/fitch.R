# Perform a maximum parsimony reconstruction of ancestral
# states given an r-state character using Fitch's algorithm
mpr.fitch2 = function(phy, data, levels, ambig) {
    stopifnot(is.tree(phy))
    stopifnot(tree.isbinary(phy))
    stopifnot(!is.null(names(data)))
    stopifnot(is.integer(data))
    stopifnot(all(data > 0L))
    stopifnot(max(data) < 32)

    r = length(unique(data))
    data = data[tiplabels(phy)]
    downpass = integer(Nnode(phy))
    uppass = integer(Nnode(phy))
    pscore = integer(Nnode(phy))

    downpass[1:Ntip(phy)] = bitwShiftL(1, data-1)

    .Call(do_fitch_mpr2, phy, r, uppass, downpass, pscore)
}

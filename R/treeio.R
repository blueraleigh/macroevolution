#' Phylogenetic tree input
#'
#' Parse a phylogenetic tree in Newick string format
#'
#' @param file A filename pointing to a file containing a Newick character string.
#' @param text A Newick character string. If not \code{NULL} any \code{file}
#' argument is ignored.
#' @return An object of class \code{tree}
#' @seealso \code{\link{write.newick}}
read.newick = function(file, text=NULL) {
    if (is.null(text)) {
        stopifnot(file.exists(file))
        newick = paste0(scan(file, what=character(), quiet=TRUE, sep=";"), ";")
    } else {
        newick = text
    }
    tree = .Call(phyr_read_newick, newick)
    class(tree) = "tree"
    return (tree)
}


#' Test whether an object is a valid \code{tree} object.
is.tree = function(x) {
    return (inherits(x, "tree"))
}


print.tree = function(x, ...) {
    str(x)
}


#' Phylogenetic tree output
#'
#' Serialize a \code{tree} object to a Newick string
#'
#' @param phy An object of class \code{tree}.
#' @param file A filename pointing to a file into which the Newick string will
#' be written. If empty the output is written to the R console.
#' @seealso \code{\link{read.newick}}
write.newick = function(phy, file="") {
    stopifnot(is.tree(phy))
    if (file != "")
        cat(.Call(phyr_write_newick, phy), "\n", sep="", file=file)
    .Call(phyr_write_newick, phy)
}

#' Root node index
#'
#' @param phy An object of class \code{tree}.
#' @return The index of the root node
root = function(phy) {
    attr(phy, "root")
}

#' Number of terminal taxa in a phylogeny
#'
#' @param phy An object of class \code{tree}.
#' @return The number of terminal taxa in the phylogeny
Ntip = function(phy) {
    attr(phy, "Ntip")
}


#' Number of nodes in a phylogeny
#'
#' @param phy An object of class \code{tree}.
#' @return The number of terminal and internal nodes in the phylogeny.
Nnode = function(phy) {
    attr(phy, "Nnode")
}


#' Terminal taxa labels
#'
#' @param phy An object of class \code{tree}.
#' @return The labels of the terminal taxa in the phylogeny.
tiplabels = function(phy) {
    if (is.null(tip.label <- attr(phy, "tip.label"))) {
        stopifnot(is.tree(phy))
        tip.label = .Call(phyr_tiplabels, phy)
        attr(phy, "tip.label") = tip.label
    }
    return (tip.label)
}


#' Notes for nodes embedded in a Newick strings
#'
#' @param phy An object of class \code{tree}.
#' @return The embedded notes
notes = function(phy) {
    if (is.null(notes <- attr(phy, "notes"))) {
        stopifnot(is.tree(phy))
        notes = .Call(phyr_node_notes, phy)
        attr(phy, "notes") = notes
    }
    return (notes)
}


#' Branch lengths for all edges
#'
#' @param phy An object of class \code{tree}.
#' @return The branch lengths of all edges in the phylogeny.
brlens = function(phy) {
    if (is.null(brlen <- attr(phy, "brlen"))) {
        stopifnot(is.tree(phy))
        brlen = .Call(phyr_node_brlens, phy)
        attr(phy, "brlen") = brlen
    }
    return (brlen)
}


#' Ages for all nodes
#'
#' @param phy An object of class \code{tree}.
#' @return The age (time since root divergence) of all nodes in the phylogeny.
ages = function(phy) {
    if (is.null(age <- attr(phy, "age"))) {
        stopifnot(is.tree(phy))
        age = .Call(phyr_node_ages, phy)
        attr(phy, "age") = age
    }
    return (age)
}


#' Ancestors for all nodes
#'
#' @param phy An object of class \code{tree}.
#' @return A list holding the node indices of the sequence of ancestors that lead
#' to each node.
ancestors = function(phy) {
    if (is.null(anc <- attr(phy, "ancestor"))) {
        stopifnot(is.tree(phy))
        anc = lapply(1L:Nnode(phy), function(n) .Call(phyr_ancestors, phy, n)[-1L])
        attr(phy, "ancestor") = anc
    }
    return (anc)
}


#' Immediate descendants of all nodes
#'
#' @param phy An object of class \code{tree}.
#' @return A list holding the node indices of the immediate descendants of
#' each node.
children = function(phy) {
    if (is.null(kids <- attr(phy, "children"))) {
        stopifnot(is.tree(phy))
        kids = lapply(1L:Nnode(phy), function(n) .Call(phyr_children, phy, n))
        attr(phy, "children") = kids
    }
    return (kids)
}


#' Return the descendants of a node
#'
#' @param node Index of node whose descendants are desired.
#' @param phy An object of class \code{tree}.
#' @param visit If \code{ALL_NODES} returned descendants include terminal and
#' internal nodes. If \code{INTERNAL_NODES_ONLY} terminal taxa are omitted.
#' @param order Descendants may be returned in \code{PREORDER} or \code{POSTORDER}
#' sequence.
#' @return A vector holding the node indices of the descendants of \code{node}
descendants = function(node, phy, visit=c("ALL_NODES", "INTERNAL_NODES_ONLY"), order=c("PREORDER", "POSTORDER")) {
    stopifnot(is.tree(phy))
    storage.mode(node) = "integer"
    if (node <= 0 || node > Nnode(phy))
        stop("Invalid node index")

    order = match.arg(order)
    visit = match.arg(visit)

    order = switch(order, PREORDER=0L, POSTORDER=1L)
    visit = switch(visit, ALL_NODES=0L, INTERNAL_NODES_ONLY=1L)

    descendants = .Call(phyr_descendants, phy, node, visit, order)

    if (order == 0L)
        return (descendants[-1L])
    else
        return (descendants[-length(descendants)])
}


#' Return the tips descended from a node
tips = function(node, phy) {
    stopifnot(is.tree(phy))
    storage.mode(node) = "integer"
    if (node <= 0 || node > Nnode(phy))
        stop("Invalid node index")

    descendants = .Call(phyr_descendants, phy, node, 0L, 0L)

    descendants = descendants[-1L]

    return (descendants[descendants <= Ntip(phy)])
}


#' Extract a clade
#'
#' @param phy An object of class \code{tree}.
#' @param node Index of node whose subtree is desired.
#' @return An object of class \code{tree} containing all nodes and edges
#' descended from \code{node}.
extract.clade = function(phy, node) {
    stopifnot(is.tree(phy))
    storage.mode(node) = "integer"
    stopifnot(node > Ntip(phy))
    stopifnot(node < Nnode(phy))
    clade = .Call(phyr_extract_clade, phy, node)
    class(clade) = "tree"
    return (clade)
}


#' Prune a phylogeny
#'
#' @param phy An object of class \code{tree}.
#' @param tip A vector of terminal taxa labels
#' @return An object of class \code{tree} containing only those terminal taxa
#' specified by \code{tip}.
keep.tip = function(phy, tip) {
    stopifnot(is.tree(phy))
    stopifnot(class(tip) == "character")
    tips = match(tip, tiplabels(phy))
    tips = tips[!is.na(tips)]
    tips = unique(tips)
    ntip = length(tips)
    subtree = .Call(phyr_extract_subtree, phy, ntip, tips)
    class(subtree) = "tree"
    return (subtree)
}


#' Prune a phylogeny
#'
#' @param phy An object of class \code{tree}.
#' @param tip A vector of terminal taxa labels
#' @return An object of class \code{tree} containing only those terminal taxa
#' not among the set specified by \code{tip}.
drop.tip = function(phy, tip) {
    stopifnot(is.tree(phy))
    stopifnot(class(tip) == "character")
    tips = which(!tiplabels(phy) %in% tips)
    ntip = length(tips)
    subtree = .Call(phyr_extract_subtree, phy, ntip, tips)
    class(subtree) = "tree"
    return (subtree)
}


#' Determine if a phylogeny is strictly bifurcating
#'
#' @param phy An object of class \code{tree}.
tree.isbinary = function(phy) {
    stopifnot(is.tree(phy))
    if (Ntip(phy) == ((Nnode(phy) + 1) / 2))
        return (TRUE)
    else
        return (FALSE)
}


#' Determine if a phylogeny is ultrametric
#'
#' @param phy An object of class \code{tree}.
tree.isultrametric = function(phy) {
    stopifnot(is.tree(phy))
    if (var(ages(phy)[1L:Ntip(phy)]) < sqrt(.Machine$double.eps))
        return (TRUE)
    return (FALSE)
}


#' Create a copy of a phylogeny
#'
#' Because the R phylogeny object is just an external pointer to a
#' C structure any changes to the underlying C structure are propagated
#' to all R objects that point to that object. To avoid side-effects
#' of this type make a copy of the phylogeny first.
#'
#' @param phy An object of class \code{tree}.
tree.duplicate = function(phy) {
    stopifnot(is.tree(phy))
    return (read.newick(text=write.newick(phy)))
}


#' Ladderize a phylogeny
#'
#' Rotates all nodes such that the descendant with the larger
#' subtree is placed on the left.
#'
#' @param phy An object of class \code{tree}.
tree.ladderize = function(phy) {
    stopifnot(is.tree(phy))
    stopifnot(tree.isbinary(phy))
    ndesc = sapply(1L:Nnode(phy), function(n) length(descendants(n, phy)))
    phy.dup = tree.duplicate(phy)
    perm = .Call(phyr_ladderize, phy.dup, ndesc)
    r = root(phy.dup)
    ntip = Ntip(phy.dup)
    nnode = Nnode(phy.dup)
    attributes(phy.dup) = NULL
    attr(phy.dup, "root") = r
    attr(phy.dup, "Ntip") = ntip
    attr(phy.dup, "Nnode") = nnode
    # the perm attribute stores the original node index:
    # e.g. the node with index 1 in the ladderized tree has index perm[1]
    # in the unladderized tree
    attr(phy.dup, "perm") = perm
    class(phy.dup) = "tree"
    return (phy.dup)
}


#' Rotate a node
#'
#' Swaps the left and right descendants
#'
#' @param node A vector of node indices to rotate.
#' @param phy An object of class \code{tree}.
rotate.node = function(node, phy) {
    stopifnot(is.tree(phy))
    stopifnot(tree.isbinary(phy))
    storage.mode(node) = "integer"
    stopifnot(all(node > Ntip(phy) & node <= Nnode(phy)))
    phy.dup = tree.duplicate(phy)
    .Call(phyr_rotate, phy.dup, node-1L)
    r = root(phy.dup)
    ntip = Ntip(phy.dup)
    nnode = Nnode(phy.dup)
    attributes(phy.dup) = NULL
    attr(phy.dup, "root") = r
    attr(phy.dup, "Ntip") = ntip
    attr(phy.dup, "Nnode") = nnode
    class(phy.dup) = "tree"
    return (phy.dup)
}

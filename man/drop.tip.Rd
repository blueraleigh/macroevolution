% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/treeio.R
\name{drop.tip}
\alias{drop.tip}
\title{Prune a phylogeny}
\usage{
drop.tip(phy, tip)
}
\arguments{
\item{phy}{An object of class \code{tree}.}

\item{tip}{A vector of terminal taxa labels}
}
\value{
An object of class \code{tree} containing only those terminal taxa
not among the set specified by \code{tip}.
}
\description{
Prune a phylogeny
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rcm_dmm.R
\name{make.rcm.dmm}
\alias{make.rcm.dmm}
\title{Dirichlet-multinomial Markov process model (random cluster model)}
\usage{
make.rcm.dmm(
  phy,
  x,
  r,
  stateid.init,
  rate.init,
  beta.init,
  integrate.brlen = FALSE
)
}
\arguments{
\item{phy}{An object of class \code{tree}.}

\item{x}{A three column \code{data.frame} containing count data. The first
column is expected to be the label of a terminal taxon; the second column, the
label of a resource category; the third column, the number of observations
for the particular combination of labels in columns one and two.}

\item{r}{The number of states in the prior model.}

\item{stateid.init}{An optional initial partition of terminal nodes into
states. If omitted all terminals are placed in a single cluster.}

\item{rate.init}{An optional initial value for the rate of the Poisson process
prior on the evolution of states. If omitted defaults to a random value.}

\item{beta.init}{An optional initial value for the Dirichlet prior hyperparameter on
the multinomial distributions corresponding to each cluster. If omitted defaults
to the uniform Dirichlet distribution prior. This is never updated during the MCMC.}

\item{integrate.brlen}{An optional boolean value indicating whether or not to
assume branches have i.i.d. rates drawn from a Gamma prior. If TRUE, the branch
rates are integrated out of the likelihood function such that a single set of
transition probabilities applies to all branches. The default is FALSE, in
which case all branches share a single rate.}
}
\value{
A function to perform MCMC inference that takes the following arguments
\describe{
\item{niter}{The number of iterations that the MCMC will run for. Default is 1000.}
\item{thin}{The number of iterations between each recorded sample. Default is 1.}
\item{tune.rate}{Tuning parameter for the rate. Should be between 0 and 1. Default is 0.1}
\item{update.node}{Weight for proposing node state updates. Default is 1.}
\item{update.rate}{Weight for proposing updates to the rate. Default is 1.}
\item{output.file}{Name of file for recording the state of Markov chain. Default is "mcmc.out".}
\item{output.mode}{One of either "wb" or "ab", indicating whether the output file should be opened
for writing or appending. Default is "wb".}
}
The output from the MCMC function is written to a binary file that can be read
back into R with the \code{read.rcm.dmm} function. For details on its structure
consult the documentation for \code{read.rcm.dmm}.
}
\description{
Bayesian MCMC implementation of the Dirichlet-multinomial Markov process
model described by Grundler and Rabosky (2020), https://doi.org/10.1093/sysbio/syaa031
}
\examples{
\dontrun{
 data(snakediet)
 data(snaketree)

 phy = read.newick(text=snaketree)
 mcmc = make.rcm.dmm(phy, snakediet, 20L)
 mcmc(output.file="mcmc.out")

 res = read.rcm.dmm("mcmc.out")

 # plot trace of data log likelihood
 plot(res$pars[, 1])

 # plot trace of log prior probability
 plot(res$pars[, 2])

 # plot trace of rate hyperparameter
 plot(res$pars[, 3])

 # look at the posterior distribution of the number of states among
 # terminal taxa
 barplot(table(apply(out$stateid, 1, function(p) length(unique(p)))))

 # take a detailed look at a specific posterior sample, say the 250th
 obj = make.rcm.dmm.from.sample(250, "mcmc.out")

 str(obj)

 # posterior mean estimate of each multinomial in this sample
 obj$dens.mean

 # max a posteriori estimate of each multinomial in this sample
 obj$dens.map

 # marginal ancestral state probabilities. e.g. obj$asr[1, 1] is the marginal
 # posterior probability that the MRCA of all species in the phylogeny had a
 # pattern of resource use described by, e.g, the multinomial density in the
 # first row of obj$dens.map or obj$dens.mean.
 obj$asr
}
}
\seealso{
\code{\link{read.rcm.dmm}} for description of output file format.
}

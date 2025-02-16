% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phybreak_analysis_functions.R
\name{label.transmissions}
\alias{label.transmissions}
\title{Find transmission rates between individuals with certain labels}
\usage{
label.transmissions(
  MCMCstate,
  labels,
  infector.posterior.probabilities = phybreak.infector.posts(MCMCstate)$post_support,
  permute_test = FALSE,
  nPermute = 10000,
  between.clust.bs = FALSE,
  label.levels = NULL,
  label.probs = NULL,
  quantiles = c(0.025, 0.975)
)
}
\arguments{
\item{MCMCstate}{The output from sample_phybreak}

\item{labels}{The labels for each individual}

\item{infector.posterior.probabilities}{The matrix of posterior probabilities of the infectors, 
as supplied by phybreak.infector.posts(MCMCstate)$post_support}

\item{permute_test}{Whether or not to run a permutation test of how the individuals are labeled 
to find a null distribution for the label transmission rates}

\item{nPermute}{Number of permutations to use in the permutation test}

\item{between.clust.bs}{Whether or not the bootstrapping should include demographic label randomization between clusters as well as within clusters}

\item{label.levels}{The labels of the demographic groups of interest.}

\item{label.probs}{The overall proportions of labels across the clusters of interest}

\item{quantiles}{Quantiles at which to find quantiles from the null distribution}
}
\value{
A matrix containing the the probabilities for each individual to be infected by a sampled individual of a certain label
and a matrix with the number of transmissions between and within individuals with each label
}
\description{
Function to find the expected number of transmission between and within individuals 
that belong to certain groups
}

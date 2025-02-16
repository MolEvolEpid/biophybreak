% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phybreak_analysis_functions.R
\name{transmission.rate.violins}
\alias{transmission.rate.violins}
\title{Plot group transmission rates}
\usage{
transmission.rate.violins(
  transmission.list,
  include.null = FALSE,
  bcbs = FALSE,
  title = NULL
)
}
\arguments{
\item{transmission.list}{A list output as from "label.transmission.wrapper"}

\item{include.null}{Whether or not to include the null distributions of the transmission rates given the within-cluster label demographics}

\item{bcbs}{Whether or not to include the null distributions of the transmission rates given only the overall label demographics}

\item{title}{A string for the title of the plot}
}
\description{
Function to make violin plots of transmission rates 
between individuals with certain labels over multiple clusters
}

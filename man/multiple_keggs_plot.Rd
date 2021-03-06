% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AnnotationHelper.R
\name{multiple_keggs_plot}
\alias{multiple_keggs_plot}
\title{Generates a plot representing significant kegg pathways from an enrichment result.}
\usage{
multiple_keggs_plot(kegg.results.list, filename, p.threshold = 0.05,
  n.threshold = 2)
}
\arguments{
\item{kegg.results.list}{A list of enrichment results returned by 
\code{\link{kegg_enrichment}}.}

\item{filename}{The name of the file where the plot should be saved.}

\item{p.threshold}{Minimum p-value for a category to be plotted.}

\item{n.threshold}{Minimum number times a pathway is reported for it to be plotted.}
}
\description{
Generates a plot representing significant kegg pathways from an enrichment result.
}


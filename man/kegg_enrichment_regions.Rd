% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AnnotationHelper.R
\name{kegg_enrichment_regions}
\alias{kegg_enrichment_regions}
\title{Perform KEGG enrichment on a set of regions.}
\usage{
kegg_enrichment_regions(regions, annotations.list, ...)
}
\arguments{
\item{regions}{A \linkS4class{GRanges} object with regions to enriched for KEGG pathways.}

\item{annotations.list}{A list of annotation databases returned by
\code{\link{select_annotations}}.}

\item{...}{Parameters to be passed to \code{\link{kegg_enrichment}}.}
}
\value{
A vector of Entrez gene ids containing a non-redundant list of the
  genes represented by the given regions.
}
\description{
Conveniance function. Regions are converted to genes with
\code{\link{gene_from_regions}} using default parameters.
}


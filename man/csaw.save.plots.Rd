% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CsawAnalysis.R
\name{csaw.save.plots}
\alias{csaw.save.plots}
\title{Output of csaw analysis}
\usage{
csaw.save.plots(list.reads.to.counts, list.reads.filter,
  list.normalization.factors, results, list.multiple.testing, annotations,
  output.dir)
}
\arguments{
\item{list.reads.to.counts}{List returned by \code{\link{reads.to.counts}}.}

\item{list.reads.filter}{List returned by \code{\link{filter.reads.to.counts}}.}

\item{list.normalization.factors}{List returned by \code{\link{normalize.factors}}.}

\item{results}{Results of \code{\link{test.diff.binding}}.}

\item{list.multiple.testing}{List returned by \code{\link{multiple.testing}}}

\item{annotations}{The annotations to use.}

\item{output.dir}{The name of the directory where to save the files and graphs.}
}
\description{
Creates and saves plots and tables for csaw analysis: \describe{
\item{Fragements lengths.pdf}{A plot of fragments (reads) lengths in all bam files.}
\item{Window size.pdf}{A cross-correlation plot of the efficiency of the ChIP-seq data.}
\item{Window size (zoom).pdf}{A zommed version of "Window size.pdf".}
\item{Effect of filtering.pdf}{A histogram of the counted reads, with lines representing background and filtering.}
\item{Normalization efforts.pdf}{MA plots of the effect of normalization between the first file and every others.}
\item{TopTags results.txt}{The best results of \code{\link{test.diff.binding}}.}
\item{Representative windows (differential binding).txt}{Table of the most representative windows, by "differential binding"
clustering method, after the correction for multiple testing.}
\item{LogFC of best windows in each cluster.txt}{For each window cluster, the LogFC of the best windows, by "differential binding"
clustering method, after the correction for multiple testing.}
\item{Representative windows (promoter-based).txt}{Table of the most representative windows, by "promoter-based" clustering method,
after the correction for multiple testing.}
\item{LogFC of best windows in each cluster (promoter-based).txt}{For each window cluster, the LogFC of the best windows, by
"promoter-based" clustering method, after the correction for multiple testing.}
\item{csaw_clusters.gz}{The zipped table of results.}}
}


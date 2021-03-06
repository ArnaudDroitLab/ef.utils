% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ImportHelper.R
\name{read_identical}
\alias{read_identical}
\title{Import a set of identically structured files.}
\usage{
read_identical(file.names, header.columns, data.columns,
  file.labels = basename(file.names), sep = "\\t", header = TRUE, ...)
}
\arguments{
\item{file.names}{The files to be read.}

\item{header.columns}{Indices or names of row-identifying columns which should
be repeated across all files. Those columns are added only once to the
output, as the very first columns.}

\item{data.columns}{Indices or names fo the columns containing unique data in
each file. The values from each file will be added to the output.}

\item{file.labels}{A vector of labels for the imported files. This must be of
the of same length as \code{file.names}. The label is prefixed to column names
in the resulting data frame.}
}
\value{
A \code{data-frame} with the concatenated information from all files.
}
\description{
Given a set of files with identical row names and column names, this
function reads all files and concatenate the requested columns from each.
}


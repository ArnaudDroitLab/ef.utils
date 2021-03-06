% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EncodeHelper.R
\name{get.chrom.state.url}
\alias{get.chrom.state.url}
\title{Given a mnemonic and a number of states, builds the download URL for a chromatin state map.}
\usage{
get.chrom.state.url(mnemonic, number.of.states)
}
\arguments{
\item{mnemonic}{The ENCODE mnemonic for the cell type (for example, "E123")}

\item{number.of.states}{Which chromatin state map to download, either 15 or 18.}
}
\value{
The URL to download the requested chromatin state map.
}
\description{
Given a mnemonic and a number of states, builds the download URL for a chromatin state map.
}


#' Import a bed file and discards its metadata.
#'
#' @param x The name of the BED file to import.
#' @param extraCols The names and types of the extra columns.
#' @return A \linkS4class{GRanges} object representing the regions of x, without the metadata.
#' @importFrom GenomicRanges mcols
import.and.discard.metadata <- function(x, extraCols) {
    regions = rtracklayer::import(x, format="BED", extraCols=extraCols)
    GenomicRanges::mcols(regions) = NULL
    return(regions)
}

#' Import a set of bed files in a directory into a \linkS4class{GRangesList}.
#'
#' @param input.dir The name of the directory from which BED files must be imported.
#' @param file.format The format of the files to be imported. Can be "bed",
#'   "broad" or "narrow".
#' @param file.ext The extension of the files to be imported.
#' @param discard.metadata If \code{TRUE}, metadata will be discarded after importing.
#' @param dir.type The architecture of the directory from which the files are to
#'   be imported. If "plain", files are sought directly within the directory. If
#'   "mugqic", the output structure of the MUGQIC pipeline is searched for.
#' @return A \linkS4class{GRanges} object representing the regions of the files in input.dir.
#' @importFrom GenomicRanges GRangesList
#' @importFrom rtracklayer import
#' @export
import.into.grl <- function(input.dir=".", file.format="bed", file.ext=NULL, discard.metadata=FALSE, dir.type="plain") {
    # Define certain parameters based on the file format.
    if(file.format=="bed") {
        if(is.null(file.ext)) {
            file.ext="bed"
        }
        extraCols = c()
    } else if(file.format == "broad") {
        if(is.null(file.ext)) {
            file.ext="_peaks.broadPeak"
        }
        extraCols <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric")
    } else {
        if(is.null(file.ext)) {
            file.ext="_peaks.narrowPeak"
        }
        extraCols <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")
    }

    if(dir.type=="mugqic") {
        peak.input.dir = file.path(input.dir, "peak_call")

        # Make a list of all directories in the peak_call directory.
        # Each directory will contain one peak file.
        peak.dirs = list.files(peak.input.dir, include.dirs=TRUE)

        # Generate a list of all peak files.
        all.files = file.path(peak.input.dir, peak.dirs, paste(peak.dirs, file.ext, sep=""))

        list.names = peak.dirs
    } else {
        # Grab all files in directory.
        all.filenames = list.files(input.dir, pattern=file.ext, include.dirs=TRUE)
        all.files = file.path(input.dir, all.filenames)
        list.names = gsub("\\.$", "", gsub(file.ext, "", all.filenames))
    }

    # Import regions and discard extra data if requested.
    if(discard.metadata) {
        grl <- GenomicRanges::GRangesList(lapply(all.files, import.and.discard.metadata, extraCols =extraCols))
    } else {
        grl <- GenomicRanges::GRangesList(lapply(all.files, rtracklayer::import, format="BED", extraCols =extraCols))
    }
    names(grl) <- gsub(file.ext, "", list.names)

    return(grl)
}


#' Import a set of identically structured files.
#'
#' Given a set of files with identical row names and column names, this
#' function reads all files and concatenate the requested columns from each.
#'
#' @param file.names The files to be read.
#' @param header.columns Indices or names of row-identifying columns which should
#'   be repeated across all files. Those columns are added only once to the
#'   output, as the very first columns.
#' @param data.columns Indices or names fo the columns containing unique data in
#'   each file. The values from each file will be added to the output.
#' @param file.labels A vector of labels for the imported files. This must be of
#'   the of same length as \code{file.names}. The label is prefixed to column names
#'   in the resulting data frame.
#' @return A \code{data-frame} with the concatenated information from all files.
#' @export
read.identical <- function(file.names, header.columns, data.columns, file.labels=basename(file.names)) {
    results=NULL
    for(i in 1:length(file.names)) {
        file.name = file.names[i]
        file.label = file.labels[i]

        file.data = read.table(file.name, sep="\t", header=TRUE)
        if(is.null(results)) {
            results = file.data[, header.columns]
        }

        colnames(file.data) <- paste(file.label, colnames(file.data), sep=".")

        results = cbind(results, file.data[,data.columns])
    }

    return(results)
}

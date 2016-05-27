library(GenomicRanges)
library(rtracklayer)

#' Import a bed file and discards its metadata.
#'
#' @param x The name of the BED file to import.
#' @param extraCols The names and types of the extra columns.
#' @return A GRanges object representing the regions of x, without the metadata.
import.and.discard.metadata <- function(x, extraCols) {
    regions = import(x, format="BED", extraCols=extraCols)
    mcols(regions) = NULL
    return(regions)
}

#' Import a set of bed files in a directory into a GRangesList.
#'
#' @param input.dir The name of the directory from which BED files must be imported.
#' @param file.format The format of the files to be imported. Can be "bed",
#'   "broad" or "narrow".
#' @param file.ext The extension of the files to be imported.
#' @param discard.metadata If TRUE, metadata will be discarded after importing.
#' @param dir.type The architecture of the directory from which the files are to
#'   be imported. If "plain", files are sought directly within the directory. If
#'   "mugqic", the output structure of the MUGQIC pipeline is searched for.
#' @return A GRanges object representing the regions of the files in input.dir.
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
        all.file.names = list.files(input.dir, pattern=file.ext, include.dirs=TRUE)
        all.files = file.path(input.dir, all.file.name)
        list.names = gsub(file.ext, "", all.file.names)
    }

    # Import regions and discard extra data if requested.
    if(discard.metadata) {
        grl <- GRangesList(lapply(all.files, import.and.discard.metadata, extraCols =extraCols))
    } else {
        grl <- GRangesList(lapply(all.files, import, format="BED", extraCols =extraCols))
    }
    names(grl) <- gsub(file.ext, "", list.namess)
    
    return(grl)    
}
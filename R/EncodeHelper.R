#' Given a biosample's name, tries to match it with an ENCODE mnemonic.
#'
#' @param biosample The biosample identifier from ENCODE. Valid examples are
#'   GM12878, K562 or HeLa.
#'
#' @return The ENCODE mnemonic associated with the biosample.
#' @export
get_encode_mnemonic <- function(biosample) {
    return(biosample.code$EID[grep(biosample, biosample.code$E.Mnemonic, ignore.case = TRUE)])
}

#' Given a mnemonic and a number of states, builds the download URL for a chromatin state map.
#'
#' @param mnemonic The ENCODE mnemonic for the cell type (for example, "E123")
#' @param number.of.states Which chromatin state map to download, either 15 or 18.
#'
#' @return The URL to download the requested chromatin state map.
get.chrom.state.url <- function(mnemonic, number.of.states) {
    base.marks = ifelse(number.of.states==18, "core_K27ac", "coreMarks")
    url = paste0("http://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/",
                 base.marks,
                 "/jointModel/final/",
                 mnemonic, "_", number.of.states, "_",
                 base.marks, "_mnemonics.bed.gz")
}

#' Download a single chromatin state map from ENCODE.
#'
#' @param mnemonic The ENCODE mnemonic for the cell type (for example, "E123")
#' @param number.of.states Which chromatin state map to download, either 15 or 18.
#' @param download.dir The folder where the downloaded files should be stored.
#'
#' @return The path of the downloaded file, if it was found.
#'
#' @importFrom utils download.file
download.one.chrom.state <- function(mnemonic, number.of.states, download.dir) {
    download.url = get.chrom.state.url(mnemonic, number.of.states)
    dest.file = file.path(download.dir, paste0(number.of.states, ".chrom.states.bed"))
    if(download.file(download.url, dest.file) == 0) {
        return(dest.file)
    } else {
        return(NULL)
    }
}

#' Attempts to download 15 and 18 chromatin state maps from ENCODE for a given cell type.
#'
#' @param mnemonic The ENCODE mnemonic for the cell type (for example, "E123")
#' @param download.dir The folder where the downloaded files should be stored.
#'
#' @return A list containing the names of the unzipped chromatin state files, if they were found.
internal.download_chrom_states <- function(mnemonic, download.dir=".") {
    results <- list()
    for(number.of.states in c(15, 18)) {
        gzipped.states = download.one.chrom.state(mnemonic, 18, download.dir)
        if(!is.null(gzipped.states)) {
            system(paste0("gzip -d -k ", gzipped.states))
            results[[as.character(number.of.states)]] <- gsub(".gz", "", gzipped.states)
        }
    }
    
    return(results)
}

#' Import HMM chromatin states for the right cell type.
#'
#' @param biosample The biosample identifier from ENCODE. Valid examples are
#'   GM12878, K562 or HeLa.
#' @param download.dir The folder where the downloaded files should be stored.
#'
#' @return The name and path of the downloaded and unzipped files.
#'
#' @importFrom utils download.file
#' @export
import_chrom_states <- function(biosample, download.dir="."){
    # Grab the ENCODE identifier for the cell line. If there are more than
    # one (IE the biosample is ambiguous), we'll use the first.
    mnemonic <- get_encode_mnemonic(biosample)
    if(length(mnemonic) > 1) {
        warning("biosample is ambiguous: the first matchign entry will be used.")
    }

    if(length(mnemonic) > 0) {
        downloaded.states = internal.download_chrom_states(mnemonic, download.dir)
        if(!is.null(downloaded.states[["18"]])) {
            return(downloaded.states[["18"]])
        } else if(!is.null(downloaded.states[["15"]])) {
            return(downloaded.states[["15"]])
        } else {
            warning("No chromatin state segmentation were found for the given biosample.")
            return(NULL)
        }
    } else {
        warning("No ENCODE tissue match the given biosample.")
        return(NULL)    
    }
}


#' Default filtering function for \code{\link{download_encode_chip}}.
#'
#' The filtering function does three things: \enumerate{
#'   \item It removes all files which do not have the correct genome assembly.
#'   \item It removes broad marks chips, such as histones and Pol II.
#'   \item If the provided results have been re-analyzed by the ENCODE Consortium,
#'     the ENCODE results are kept and the original ones discarded.}
#'
#' @param query.results A partial \code{data.frame} obtained from the \code{\link[ENCODExplorer]{queryEncode}}
#'   function.The biosample identifier from ENCODE. Valid examples are GM12878, K562.
#' @param genome.assembly Which genome assembly should the results come from?
#' @return A filtered \code{data frame}.
#' @importFrom plyr ddply
#' @export
default_download_filter_chip <- function(query.results, genome.assembly) {
    filtered.results = plyr::ddply(query.results, ~accession, function(x, genome.assembly) {
        x = subset(x, assembly==genome.assembly)

        if(grepl("^H\\d", x$target[1]) || grepl("^POL", x$target[1])) {
            return(NULL)
        }

        result_set=x
        
        if(sum(x$lab=="ENCODE Consortium Analysis Working Group") > 0) {
            result_set = subset(x, lab=="ENCODE Consortium Analysis Working Group")
        }
        
        if(sum(result_set$output_type=="optimal idr thresholded peaks") > 0) {
            result_set = subset(result_set, output_type=="optimal idr thresholded peaks")
        } else if(sum(result_set$output_type=="conservative idr thresholded peaks") > 0) {
            result_set = subset(result_set, output_type=="conservative idr thresholded peaks")
        }
        
        return(result_set)
        
    }, genome.assembly=genome.assembly)

    return(filtered.results)
}

#' Alternative filtering function for \code{\link{download_encode_chip}}.
#'
#' The filtering function does three things: \enumerate{
#'   \item It removes all files which do not have the correct genome assembly.
#'   \item It removes all marks chips, except histone marks chips.
#'   \item If the provided results have been re-analyzed by the ENCODE Consortium,
#'     the ENCODE results are kept and the original ones discarded.}
#'
#' @param query.results A partial \code{data.frame} obtained from the \code{\link[ENCODExplorer]{queryEncode}}
#'   function.The biosample identifier from ENCODE. Valid examples are GM12878, K562.
#' @param genome.assembly Which genome assembly should the results come from?
#' @return A filtered \code{data frame}.
#' @importFrom plyr ddply
#' @export
histone_download_filter_chip <- function(query.results, genome.assembly) {
  filtered.results = plyr::ddply(query.results, ~accession, function(x, genome.assembly) {
    x = subset(x, assembly==genome.assembly)

    if(!grepl("^H\\d", x$target[1])) {
      return(NULL)
    }

    if(sum(x$lab=="ENCODE Consortium Analysis Working Group") > 0) {
      return(subset(x, lab=="ENCODE Consortium Analysis Working Group"))
    } else {
      return(x)
    }
  }, genome.assembly=genome.assembly)

  return(filtered.results)
}

#' Alternative filtering function for \code{\link{download_encode_chip}}.
#'
#' The filtering function does three things: \enumerate{
#'   \item It removes all files which do not have the correct genome assembly.
#'   \item It removes all marks chips, except Pol II marks chips.
#'   \item If the provided results have been re-analyzed by the ENCODE Consortium,
#'     the ENCODE results are kept and the original ones discarded.}
#'
#' @param query.results A partial \code{data.frame} obtained from the \code{\link[ENCODExplorer]{queryEncode}}
#'   function.The biosample identifier from ENCODE. Valid examples are GM12878, K562.
#' @param genome.assembly Which genome assembly should the results come from?
#' @return A filtered \code{data frame}.
#' @importFrom plyr ddply
#' @export
pol2_download_filter_chip <- function(query.results, genome.assembly) {
  filtered.results = plyr::ddply(query.results, ~accession, function(x, genome.assembly) {
    x = subset(x, assembly==genome.assembly)

    if(!grepl("^POL", x$target[1])) {
      return(NULL)
    }

    if(sum(x$lab=="ENCODE Consortium Analysis Working Group") > 0) {
      return(subset(x, lab=="ENCODE Consortium Analysis Working Group"))
    } else {
      return(x)
    }
  }, genome.assembly=genome.assembly)

  return(filtered.results)

}

#' Default filtering function for \code{\link{download_encode_rna}}.
#'
#' The filtering function does three things: \enumerate{
#'   \item It removes all files which do not have the correct genome assembly.
#'   \item It only keeps files with counts for genes, not transcripts or isoforms.
#'   \item If removes any file where treatment is not NA.}
#'
#' @param query.results A partial \code{data.frame} obtained from the \code{\link[ENCODExplorer]{queryEncode}}
#'   function.The biosample identifier from ENCODE. Valid examples are GM12878, K562.
#' @param genome.assembly Which genome assembly should the results come from?
#' @return A filtered \code{data frame}.
#' @importFrom plyr ddply
#' @export
default_download_filter_rna <- function(query.results, genome.assembly) {
    filtered.results = plyr::ddply(query.results, ~accession, function(x, genome.assembly) {
        return(subset(x, grepl("genes", x$submitted_file_name) & assembly==genome.assembly & is.na(treatment)))
    }, genome.assembly=genome.assembly)

    return(filtered.results)
}

#' Alternative filtering function for \code{\link{download_encode_rna}}.
#'
#' The filtering function does three things: \enumerate{
#'   \item It removes all files which do not have the correct genome assembly.
#'   \item It only keeps files with counts for isoforms, not transcripts or genes.
#'   \item If removes any file where treatment is not NA.}
#'
#' @param query.results A partial \code{data.frame} obtained from the \code{\link[ENCODExplorer]{queryEncode}}
#'   function.The biosample identifier from ENCODE. Valid examples are GM12878, K562.
#' @param genome.assembly Which genome assembly should the results come from?
#' @return A filtered \code{data frame}.
#' @importFrom plyr ddply
#' @export
isoform_download_filter_rna <- function(query.results, genome.assembly) {
  filtered.results = plyr::ddply(query.results, ~accession, function(x, genome.assembly) {
    return(subset(x, grepl("isoform", x$submitted_file_name) & assembly==genome.assembly & is.na(treatment)))
  }, genome.assembly=genome.assembly)

  return(filtered.results)
}

consensus.and.signal.mean <- function(grl, keep.signal) {
    overlap.results = intersect_overlap(build_intersect(grl, keep.signal = keep.signal))
    mcols(overlap.results) <- rowMeans(as.data.frame(mcols(overlap.results)), na.rm = TRUE)
    names(mcols(overlap.results)) <- "signalValue"
    return(unlist(overlap.results))
}

#' @importFrom rtracklayer import
import.plus.consensus <- function(files.to.import, file.format, file.ext = ".bed", keep.signal = FALSE) {
    grl = import_files_into_grl(files.to.import, file.format, file.ext=file.ext, discard.metadata=!keep.signal)
    return(consensus.and.signal.mean(grl, keep.signal))
}

download.chip.and.import <- function(query.results, peak.type, out.dir=".", keep.signal = FALSE) {
    # Create a directory for downloaded files.
    download.dir = file.path(out.dir, peak.type)
    dir.create(download.dir, recursive = TRUE, showWarnings=FALSE)

    # Subset the query results to only keep files in the right format.
    query.results.subset = query.results
    query.results.subset$experiment = query.results$experiment[grepl(peak.type, query.results$experiment$file_format_type),]

    # Download the files.
    downloaded.files = ENCODExplorer::downloadEncode(resultSet=query.results.subset, resultOrigin="queryEncode", dir=download.dir, force=FALSE)

    # Write the metadata about the downloaded files.
    write.table(query.results$experiment, file=file.path(out.dir, paste0(peak.type, ".metadata.txt")))

    # Unzip the files. Use gzip -d since on windows, gunzip is not installed by default.
    if(dir.exists(download.dir)) {
      system(paste0("gzip -d -k ", download.dir, "/*.gz"))
      
      results.list = list()
      for(target in unique(query.results.subset$experiment$target)) {
        query.subset = query.results.subset$experiment[query.results.subset$experiment$target == target,]
        
        # Perform consensuss across replicates.
        accession.list = list()
        for(accession in unique(query.subset$accession)) {
          accession.files = file.path(download.dir, paste0(query.subset$file_accession, ".bed"))
          accession.list[[accession]] = import.plus.consensus(accession.files, peak.type, keep.signal=keep.signal)
        }
        
        # Perform consensus across experiments.
        results.list[[target]] = consensus.and.signal.mean(GRangesList(accession.list), keep.signal)
      }
      
      out.gr = GRangesList(results.list)
      out.gr@unlistData@elementMetadata@listData$peak <- NULL
    } else {
      out.gr = NULL
    }
    
    return(list(Metadata=query.results.subset$experiment,
                Downloaded=downloaded.files,
                Regions=out.gr))
}

#' Obtains and process transcription factor data from ENCODE.
#'
#' Given an encode biosample identifier (GM12878, MCF-7, etc.), this function
#' retrieved processed bed files giving the genomic location of binding peaks.
#'
#' @param biosample The biosample identifier from ENCODE. Valid examples are
#'   GM12878, K562 or MCF-7.
#' @param genome.assembly Which genome assembly should the results come from?
#' @param download.filter A filtering function to be applied to the experiment
#'   data-frame returned by \code{\link[ENCODExplorer]{queryEncode}}. The function should return a filtered
#'   \code{data-frame} containing only the files which should be downloaded.
#' @param download.dir The folder where the downloaded files should be stored.
#'   defaults to \code{file.path("input/ENCODE", biosample, "chip-seq")}.
#' @param keep.signal Should the signalValue from the ChIP-seq be kept?
#' @return A list containing three elements: \describe{
#' \item{Metadata}{The metadata returned by \code{\link[ENCODExplorer]{queryEncode}}, containing information
#'     about all files which matched the query.}
#' \item{Downloaded}{The list of files which were downloaded.}
#' \item{Regions}{The processed peak regions.}}
#' @importFrom ENCODExplorer queryEncode
#' @importFrom ENCODExplorer downloadEncode
#' @importFrom GenomicRanges GRangesList
#' @importFrom GenomicRanges mcols
#' @importFrom stats aggregate
#' @importMethodsFrom GenomicRanges findOverlaps
#' @export
download_encode_chip <- function(biosample, assembly, download.filter=default_download_filter_chip,
                                   download.dir=file.path("input/ENCODE", biosample, assembly, "chip-seq"), keep.signal = FALSE) {
    # Query ENCODE to obtain appropriate files.
    # queryEncode has a bug and will fail if encode_df is not loaded.
    data("encode_df", package="ENCODExplorer")
    query.results = ENCODExplorer::queryEncode(assay="ChIP-seq", biosample=biosample, file_format="bed", status="released")

    # Filter the ENCODE files using the supplied functions.  Only download relevant files.
    query.results$experiment = download.filter(query.results$experiment, assembly)
    dir.create(download.dir, recursive=TRUE, showWarnings=FALSE)

    results = list()
    for(peak.type in c("narrow", "broad")) {
        results[[peak.type]] = download.chip.and.import(query.results, peak.type, download.dir, keep.signal)
    }
    
    return(results)
}

#' Obtains and process expression data from ENCODE.
#'
#' Given an encode biosample identifier (GM12878, MCF-7, etc.), this function
#' retrieved processed expression levels.
#'
#' @param biosample The biosample identifier from ENCODE. Valid examples are
#'   GM12878, K562 or MCF-7.
#' @param genome.assembly Which genome assembly should the results come from?
#' @param download.filter A filtering function to be applied to the experiment
#'   data-frame returned by \code{\link[ENCODExplorer]{queryEncode}}. The function should return a filtered
#'   \code{data-frame} containing only the files which should be downloaded.
#' @param download.dir The folder where the downloaded files should be stored.
#'   defaults to \code{file.path("input/ENCODE", biosample, "chip-seq")}.
#' @return A list containing three elements: \describe{
#' \item{Metadata}{The metadata returned by queryEncode, containing information
#'     about all files which matched the query.}
#' \item{Downloaded}{The list of files which were downloaded.}
#' \item{Expression}{A data-frame with processed expression levels.}}
#' @importFrom ENCODExplorer queryEncode
#' @importFrom ENCODExplorer downloadEncode
#' @export
download_encode_rna <- function(biosample, assembly, download.filter=default_download_filter_rna, download.dir=file.path("input/ENCODE", biosample, assembly, "rna-seq")) {
    # Query ENCODE to obtain appropriate files.
    # queryEncode has a bug and will fail if encode_df is not loaded.
    data("encode_df", package="ENCODExplorer")
    query.results = ENCODExplorer::queryEncode(assay="RNA-seq", biosample=biosample, file_format="tsv", status="released")

    # Filter the ENCODE files using the supplied functions.  Only download relevant files.
    query.results$experiment = download.filter(query.results$experiment, assembly)
    dir.create(download.dir, recursive=TRUE, showWarnings=FALSE)
    downloaded.files = ENCODExplorer::downloadEncode(resultSet=query.results, resultOrigin="queryEncode", dir=download.dir, force=FALSE)

    # Read the files.
    if(!is.null(query.results)) {
        #rna.filenames = list.files(download.dir)
        rna.data = read_identical(downloaded.files, 1:5, 6:7, file.labels=gsub(".tsv", "", downloaded.files))
        
        # Calculate mean of metrics.
        for(metric in c("TPM", "FPKM")) {
            mean.metric = apply(rna.data[,grepl(metric, colnames(rna.data))], 1, mean, na.rm=TRUE)
            #sd.metric = apply(rna.data[,grepl(metric, colnames(rna.data))], 1, sd, na.rm=TRUE)
            rna.data = cbind(rna.data, mean.metric)
            colnames(rna.data)[ncol(rna.data)] <- paste0("Mean.", metric)
            #rna.data = cbind(rna.data, mean.metric)
            #colnames(rna.data)[ncol(rna.data)] <- paste0("SD.", metric)
        }
    } else {
        rna.data = NULL
    }
    
    # Return results
    return(list(Metadata=query.results$experiment,
                Downloaded=downloaded.files,
                Expression=rna.data))    
}

#' Helper function for obtaining transcription factor data through download_encode_chip.
#'
#' @param biosample The biosample identifier from ENCODE. Valid examples are
#'   GM12878, K562 or MCF-7.
#' @param genome.assembly Which genome assembly should the results come from?
#' @param download.dir The folder where the downloaded files should be stored.
#'   defaults to \code{file.path("input/ENCODE", biosample, assembly, "chip-seq", "tf")}.
#' @param ... Other parameters to be passed to download_encode_chip.
#' @return A list containing three elements: \describe{
#' \item{Metadata}{The metadata returned by \code{\link[ENCODExplorer]{queryEncode}}, containing information
#'     about all files which matched the query.}
#' \item{Downloaded}{The list of files which were downloaded.}
#' \item{Regions}{The processed peak regions.}}
#' @export
download_encode_tf <- function(biosample, assembly, 
                               download.dir=file.path("input/ENCODE", biosample, assembly, "chip-seq", "tf"),
                               ...) {
    return(download_encode_chip(biosample, assembly,
                                download.filter=default_download_filter_chip, 
                                download.dir=download.dir, ...))
}

#' Helper function for obtaining histone data through download_encode_chip.
#'
#' @param biosample The biosample identifier from ENCODE. Valid examples are
#'   GM12878, K562 or MCF-7.
#' @param genome.assembly Which genome assembly should the results come from?
#' @param download.dir The folder where the downloaded files should be stored.
#'   defaults to \code{file.path("input/ENCODE", biosample, assembly, "chip-seq", "tf")}.
#' @param ... Other parameters to be passed to download_encode_chip.
#' @return A list containing three elements: \describe{
#' \item{Metadata}{The metadata returned by \code{\link[ENCODExplorer]{queryEncode}}, containing information
#'     about all files which matched the query.}
#' \item{Downloaded}{The list of files which were downloaded.}
#' \item{Regions}{The processed peak regions.}}
#' @export
download_encode_histones <- function(biosample, assembly,
                                     download.dir=file.path("input/ENCODE", biosample, assembly, "chip-seq", "histones"),
                                     ...) {
    return(download_encode_chip(biosample, assembly, 
                                download.filter=histone_download_filter_chip, 
                                download.dir=download.dir, ...))
}

#' Helper function for obtaining polymerase data through download_encode_chip.
#'
#' @param biosample The biosample identifier from ENCODE. Valid examples are
#'   GM12878, K562 or MCF-7.
#' @param genome.assembly Which genome assembly should the results come from?
#' @param download.dir The folder where the downloaded files should be stored.
#'   defaults to \code{file.path("input/ENCODE", biosample, assembly, "chip-seq", "tf")}.
#' @param ... Other parameters to be passed to download_encode_chip.
#' @return A list containing three elements: \describe{
#' \item{Metadata}{The metadata returned by \code{\link[ENCODExplorer]{queryEncode}}, containing information
#'     about all files which matched the query.}
#' \item{Downloaded}{The list of files which were downloaded.}
#' \item{Regions}{The processed peak regions.}}
#' @export
download_encode_polymerases <- function(biosample, assembly,
                                        download.dir=file.path("input/ENCODE", biosample, assembly, "chip-seq", "polymerases"),
                                        ...) {
    return(download_encode_chip(biosample, assembly,
                                download.filter=pol2_download_filter_chip, 
                                download.dir=download.dir, ...))
}

#' Helper function for obtaining chromatin states through import_chrom_states.
#'
#' @param biosample The biosample identifier from ENCODE. Valid examples are
#'   GM12878, K562 or MCF-7.
#' @param genome.assembly Which genome assembly should the results come from?
#' @param download.dir The folder where the downloaded files should be stored.
#'   defaults to \code{file.path("input/ENCODE", biosample, assembly, "chip-seq", "tf")}.
#' @return A GRanges object with the loaded chromatin states, or NULL if the do not exist.
#' @export
download_chromatin_states <- function(biosample, assembly, download.dir=file.path("input/ENCODE", biosample, assembly)) {
    downloaded.file = import_chrom_states(biosample, download.dir)
    if(!is.null(downloaded.file)) {
        return(load_chrom_state(downloaded.file))
    } else {
        return(NULL)
    }
}
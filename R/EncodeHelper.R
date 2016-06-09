#' Default filtering function for download.encode.chip.
#' 
#' The filtering function does three things:
#'   - It removes all files which do not have the correct genome assembly.
#'   - It removes broad marks chips, such as histones and Pol II.
#'   - If the provided results have been re-analyzed by the ENCODE Consortium,
#'     the ENCODE results are kept and the original ones discarded.
#'
#' @param query.results A partial data.frame obtained from the queryEncode
#'   function.The biosample identifier from ENCODE. Valid examples are
#' @param genome.assembly Which genome assembly should the results come from?
#' @return A filtered data frame.
#' @export
default.download.filter.chip <- function(query.results, genome.assembly) {
    filtered.results = plyr::ddply(query.results, ~accession, function(x, genome.assembly) {
        x = subset(x, assembly==genome.assembly)
        
        if(grepl("^H\\d", x$target) | grepl("^POL", x$target)) {
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

#' Default filtering function for download.encode.rna.
#' 
#' The filtering function does three things:
#'   - It removes all files which do not have the correct genome assembly.
#'   - It only keeps files with counts for genes, not transcripts or isoforms.
#'   - If removes any file where treatment is not NA.
#'
#' @param query.results A partial data.frame obtained from the queryEncode
#'   function.The biosample identifier from ENCODE. Valid examples are
#' @param genome.assembly Which genome assembly should the results come from?
#' @return A filtered data frame.
#' @export
default.download.filter.rna <- function(query.results, genome.assembly) {
    filtered.results = plyr::ddply(query.results, ~accession, function(x, genome.assembly) {
        return(subset(x, grepl("genes", x$submitted_file_name) & assembly==genome.assembly & is.na(treatment)))
    }, genome.assembly=genome.assembly)
    
    return(filtered.results)
}

#' Obtains and process transcription factor data from ENCODE.
#' 
#' Given an encode biosample identifier (GM12878, MCF-7, etc.), this function
#' retrieved processed bed files giving the genomic location of binding peaks.
#'
#' @param biosample The biosample identifier from ENCODE. Valid examples are
#'   GM12878, K562 or MCF-7.
#' @param download.filter A filtering function to be applied to the experiment
#'   data-frame returned by queryEncode. The function should return a filtered
#'   data-frame containing only the files which should be downloaded.
#' @param download.dir The folder where the downloaded files should be stored.
#'   defaults to file.path("input/ENCODE", biosample, "chip-seq").
#' @return A list containing three elements:
#'   - Metadata: The metadata returned by queryEncode, containing information 
#'     about all files which matched the query.
#'   - Downloaded: The list of files which were downloaded.
#'   - Regions: The processed peak regions.
#' @export
download.encode.chip <- function(biosample, assembly, download.filter=default.download.filter, download.dir=file.path("input/ENCODE", biosample, "chip-seq")) {
    # Query ENCODE to obtain appropriate files.
    query.results = ENCODExplorer::queryEncode(assay="ChIP-seq", biosample=biosample, file_format="bed", assembly=assembly, status="released")
    
    # Filter the ENCODE files using the supplied functions.  Only download relevant files.
    query.results$experiment = download.filter(query.results$experiment, assembly)
    dir.create(download.dir, recursive=TRUE, showWarnings=FALSE)
    downloaded.files = ENCODExplorer::downloadEncode(resultSet=query.results, resultOrigin="queryEncode", dir=download.dir, force=FALSE)

    # Unzip the files. Use gzip -d since on windows, gunzip is not installed by default.
    # Also, use -k to keep the original file so that future calls to downloadEncode will not 
    # redownload the files.
    system(paste0("gzip -d -k ", download.dir, "/*.gz"))
    
    # Write the metadata about the downloaded files.
    write.table(query.results$experiment, file=file.path(download.dir, "metadata.txt"))
    
    # Import the downloaded files.
    all.gr = import.into.grl(download.dir, file.format="narrow", file.ext="bed", discard.metadata=TRUE, dir.type="plain")

    # Combine biological/technical replicates using consensus regions.
    accession.replicates = GenomicRanges::GRangesList(dlply(query.results$experiment, ~accession, function(x) {
        gr.subset = all.gr[names(all.gr) %in% x$file_accession]
        return(intersect.overlap(build.intersect(gr.subset)))
    }))
    
    # Combine all same-target replicates using consensus regions.
    target.replicates = GenomicRanges::GRangesList(dlply(query.results$experiment, ~target, function(x) {
        gr.subset = accession.replicates[names(accession.replicates) %in% x$accession]
        return(intersect.overlap(build.intersect(gr.subset)))
    }))
    
    return(list(Metadata=query.results$experiment,
                Downloaded=downloaded.files,
                Regions=target.replicates))
}

#' Obtains and process expression data from ENCODE.
#' 
#' Given an encode biosample identifier (GM12878, MCF-7, etc.), this function
#' retrieved processed expression levels.
#'
#' @param biosample The biosample identifier from ENCODE. Valid examples are
#'   GM12878, K562 or MCF-7.
#' @param download.filter A filtering function to be applied to the experiment
#'   data-frame returned by queryEncode. The function should return a filtered
#'   data-frame containing only the files which should be downloaded.
#' @param download.dir The folder where the downloaded files should be stored.
#'   defaults to file.path("input/ENCODE", biosample, "chip-seq").
#' @return A list containing three elements:
#'   - Metadata: The metadata returned by queryEncode, containing information 
#'     about all files which matched the query.
#'   - Downloaded: The list of files which were downloaded.
#'   - Expression: A data-frame with processed expression levels.
#' @export
download.encode.rna <- function(biosample, assembly, download.filter=default.download.filter.rna, download.dir=file.path("input/ENCODE", biosample, "rna-seq")) {
    # Query ENCODE to obtain appropriate files.
    query.results = ENCODExplorer::queryEncode(assay="RNA-seq", biosample=biosample, file_format="tsv", status="released")
    
    # Filter the ENCODE files using the supplied functions.  Only download relevant files.
    query.results$experiment = download.filter(query.results$experiment, assembly)
    dir.create(download.dir, recursive=TRUE, showWarnings=FALSE)
    downloaded.files = ENCODExplorer::downloadEncode(resultSet=query.results, resultOrigin="queryEncode", dir=download.dir, force=FALSE)
    
    # Read the files.
    rna.filenames = list.files(download.dir)
    rna.data = read.identical(file.path(download.dir, rna.filenames), 1:5, 6:7, file.labels=gsub(".tsv", "", rna.filenames))
    
    # Calculate mean of metrics.
    for(metric in c("TPM", "FPKM")) {
        mean.metric = apply(rna.data[,grepl(metric, colnames(rna.data))], 1, mean, na.rm=TRUE)
        #sd.metric = apply(rna.data[,grepl(metric, colnames(rna.data))], 1, sd, na.rm=TRUE)
        rna.data = cbind(rna.data, mean.metric)
        colnames(rna.data)[ncol(rna.data)] <- paste0("Mean.", metric)
        #rna.data = cbind(rna.data, mean.metric)
        #colnames(rna.data)[ncol(rna.data)] <- paste0("SD.", metric)
    }
    
    
    # Return results
    return(list(Metadata=query.results$experiment,
                Downloaded=downloaded.files,
                Expression=rna.data))
}
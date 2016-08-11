#' Create an object storing parameters for ChIA processing.
#'
#' @param input.chrom.state The name of the file containing the information about chromatin states.
#' @param biosample The biosample identifier from ENCODE. Valid examples are GM12878, K562 or MCF-7.
#' @param genome.build The name of the chosen annotation ("hg38", "hg19").
#' @param tf.regions A \linkS4class{GRangesList} object containing TF regions.
#' @param histone.regions A \linkS4class{GRangesList} object containing histone regions.
#' @param pol.regions A \linkS4class{GRangesList} object containing the pol2 regions.
#' @param expression.levels A data frame containing the levels of expression of genes, according to their EMSEMBL id.
#'
#' @return An environment that can be passed to process.chia.pet, annotate.chia.pet or analyze.chia.pet.
#' @export
build.chia.params <- function(input.chrom.state = NULL, biosample = NULL, genome.build = NULL, tf.regions = NULL,
                             histone.regions = NULL, pol.regions = NULL, expression.data = NULL, tssRegion = c(-3000, 3000),
                             centrality.measures=c("Degree"), weight.attr=NULL) {
    chia.params = new.env()
    
    chia.params$biosample = biosample
    chia.params$genome.build = genome.build

    chia.params$input.chrom.state = input.chrom.state
    chia.params$tf.regions = tf.regions
    chia.params$histone.regions = histone.regions
    chia.params$pol.regions = pol.regions
    chia.params$expression.data = expression.data

    chia.params$tssRegion = tssRegion
    chia.params$centrality.measures = centrality.measures
    chia.params$weight.attr = weight.attr
               
    return(chia.params)
}

#' Retrieve whatever data it can from ENCODE and add it to the ChIA parameter object.
#'
#' @param chia.params The ChIA parameters object for which additional annotations must be fetched from ENCODE.
#'
#' @return An environment that can be passed to process.chia.pet, annotate.chia.pet or analyze.chia.pet.
#' @export
add.encode.data <- function(chia.params) {
    biosample = chia.params$biosample
    genome.build = chia.params$genome.build
    
    # If biosample is provided, download missing annotations from ENCODE.
    if(!is.null(biosample) && !is.null(genome.build)) {
        # Download transcription factors
        if (is.null(chia.params$tf.regions)) {
            chia.params$tf.regions <- download.encode.chip(biosample, genome.build)$Regions
        }

        # Download histone marks
        if (is.null(chia.params$histone.regions)) {
            chia.params$histone.regions <- download.encode.chip(biosample, genome.build, download.filter=histone.download.filter.chip,
                                                  download.dir=file.path("input/ENCODE", biosample, "chip-seq", "histone"))$Regions
        }

        # Download PolII regions.
        if (is.null(chia.params$pol.regions)) {
            chia.params$pol.regions <- download.encode.chip(biosample, genome.build, download.filter = pol2.download.filter.chip,
                                              download.dir = file.path("input/ENCODE", biosample, "chip-seq", "pol2"))$Regions
        }

        # Download expression data
        if (is.null(chia.params$expression.data)) {
            chia.params$expression.data <- download.encode.rna(biosample, genome.build)$Expression
            chia.params$expression.data$ENSEMBL = gsub("\\.\\d+$", "", chia.params$expression.data$gene_id)
            chia.params$expression.data$FPKM = log2(chia.params$expression.data$Mean.FPKM + 1)
        }

        # Download chromatin states
        if (is.null(chia.params$input.chrom.state) && genome.build=="hg19") {
            chia.params$input.chrom.state <- import.chrom.states(biosample, file.path("input/chrom_states", biosample))
        }
    }
    
    return(chia.params)
}
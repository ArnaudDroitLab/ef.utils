#############################################################################################################
# Functions that should always run:
#############################################################################################################

#' Install required annotation packages.
#'
#' Install \pkg{BiocInstaller} and the required annotation packages for the given
#' genome assemblies.
#'
#' @param which.assembly Which assembly's annotations should be installed.
#'   Supported assemblies are hg19, hg28, mm10 and mm9.
#' @importFrom BiocInstaller biocLite
#' @export
install_annotations <- function(which.assembly) {
    if(!requireNamespace("BiocInstaller")) {
        source("https://bioconductor.org/biocLite.R")
        biocLite("BiocInstaller")
    }

    if("all" %in% which.assembly) {
        which.assembly = c("hg19", "hg38", "mm9", "mm10")
    }

    if("hg19" %in% which.assembly) {
        if(requireNamespace("BiocInstaller")) {
            BiocInstaller::biocLite(pkgs=c("BSgenome.Hsapiens.UCSC.hg19",
                 "TxDb.Hsapiens.UCSC.hg19.knownGene",
                 "PWMEnrich.Hsapiens.background",
                 "org.Hs.eg.db"))
        }
    }

    if("hg38" %in% which.assembly) {
        if(requireNamespace("BiocInstaller")) {
            BiocInstaller::biocLite(pkgs=c("BSgenome.Hsapiens.UCSC.hg38",
                 "TxDb.Hsapiens.UCSC.hg38.knownGene",
                 "PWMEnrich.Hsapiens.background",
                 "org.Hs.eg.db"))
        }
    }

    if("mm10" %in% which.assembly) {
        if(requireNamespace("BiocInstaller")) {
            BiocInstaller::biocLite(pkgs=c("BSgenome.Mmusculus.UCSC.mm10",
                 "TxDb.Mmusculus.UCSC.mm10.knownGene",
                 "PWMEnrich.Mmusculus.background",
                 "org.Mm.eg.db"))
        }
    }

    if("mm9" %in% which.assembly) {
        if(requireNamespace("BiocInstaller")) {
            BiocInstaller::biocLite(pkgs=c("BSgenome.Mmusculus.UCSC.mm9",
                 "TxDb.Mmusculus.UCSC.mm9.knownGene",
                 "PWMEnrich.Mmusculus.background",
                 "org.Mm.eg.db"))
        }
    }
}

#' Annotations selection helper.
#'
#' Given a genome build identifier, creates a list of annotation databases
#' which can be passed to the annotation functions.
#'
#' @param genome.build The genome build to use for annotation.
#' @return A list with the following elements: \describe{
#' \item{TxDb}{A TxDb of class \linkS4class{AnnotationDBI} providing informations about the genes
#'   of the selected build.}
#' \item{OrgDb}{An \linkS4class{OrgDb} object for the selected species.}
#' \item{OrgDbStr}{A character string representing the name of the \linkS4class{OrgDb} object.}
#' \item{BSGenome}{A \linkS4class{BSGenome} object to retrieve DNA sequences.}
#' \item{KEGG}{A cache of KEGG pathways, as returned by \code{\link[gage]{kegg.gsets}} from the \pkg{gage} library.}
#' \item{PWMBG}{A PWMLogn background for motif enrichment of promoter regions.}}
#'
#' @export
select_annotations <- function(genome.build) {
    if(genome.build=="hg38") {
        if(requireNamespace(c("BSgenome.Hsapiens.UCSC.hg38",
                              "TxDb.Hsapiens.UCSC.hg38.knownGene",
                              "PWMEnrich.Hsapiens.background"), quietly=TRUE)) {
            data(PWMLogn.hg19.MotifDb.Hsap, package="PWMEnrich.Hsapiens.background")

            return(list(TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
                        OrgDbStr="org.Hs.eg.db",
                        OrgDb=org.Hs.eg.db::org.Hs.eg.db,
                        PWMBG=PWMLogn.hg19.MotifDb.Hsap,
                        BSGenome=BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
                        KEGG=hs.keggs))
        } else {
            stop("The required annotation packages are not available.")
        }
    } else if(genome.build=="mm10") {
        if(requireNamespace(c("BSgenome.Mmusculus.UCSC.mm10",
                              "TxDb.Mmusculus.UCSC.mm10.knownGene",
                              "PWMEnrich.Mmusculus.background"), quietly=TRUE)) {
            data(PWMLogn.mm9.MotifDb.Mmus, package="PWMEnrich.Mmusculus.background")

            return(list(TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene,
                        OrgDbStr="org.Mm.eg.db",
                        OrgDb=org.Mm.eg.db::org.Mm.eg.db,
                        PWMBG=PWMLogn.mm9.MotifDb.Mmus,
                        BSGenome=BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10,
                        KEGG=mm.keggs))
        } else {
            stop("The required annotation packages are not available.")
        }
    } else if(genome.build=="mm9") {
        if(requireNamespace(c("BSgenome.Mmusculus.UCSC.mm9",
                              "TxDb.Mmusculus.UCSC.mm9.knownGene",
                              "PWMEnrich.Mmusculus.background"), quietly=TRUE)) {
            data(PWMLogn.mm9.MotifDb.Mmus, package="PWMEnrich.Mmusculus.background")

            return(list(TxDb=TxDb.Mmusculus.UCSC.mm9.knownGene::TxDb.Mmusculus.UCSC.mm9.knownGene,
                        OrgDbStr="org.Mm.eg.db",
                        OrgDb=org.Mm.eg.db::org.Mm.eg.db,
                        PWMBG=PWMLogn.mm9.MotifDb.Mmus,
                        BSGenome=BSgenome.Mmusculus.UCSC.mm9::BSgenome.Mmusculus.UCSC.mm9,
                        KEGG=mm.keggs))
        } else {
            stop("The required annotation packages are not available.")
        }
    } else if(genome.build=="hg19") {
        if(requireNamespace(c("BSgenome.Hsapiens.UCSC.hg19",
                              "TxDb.Hsapiens.UCSC.hg19.knownGene",
                              "PWMEnrich.Hsapiens.background"), quietly=TRUE)) {
            data(PWMLogn.hg19.MotifDb.Hsap, package="PWMEnrich.Hsapiens.background")

            return(list(TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
                        OrgDbStr="org.Hs.eg.db",
                        OrgDb=org.Hs.eg.db::org.Hs.eg.db,
                        PWMBG=PWMLogn.hg19.MotifDb.Hsap,
                        BSGenome=BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
                        KEGG=hs.keggs))
        } else {
            stop("The required annotation packages are not available.")
        }
    } else {
        stop("The selected genome build is not supported.")
    }

    # Cannot be reached.
    stop("ERROR: Cannot be reached!")
}

#' Annotate "chip.data", given as parameter.
#'
#' @param chip.data A \linkS4class{GRanges} object containing the ChIP-seq data.
#' @param input.chrom.state The name of the file containing the information about chromatin states.
#' @param tf.regions A \linkS4class{GRangesList} object containing the TF regions.
#' @param histone.regions A \linkS4class{GRangesList} object containing the histone regions.
#' @param pol.regions A \linkS4class{GRangesList} object containing the pol2 regions.
#' @param expression.levels A data frame containing the levels of expression of genes. according to their EMSEMBL id.
#' @param genome.build The name of the chosen annotation ("hg38", "mm9", "mm10", "hg19").
#' @param biosample The biosample identifier from ENCODE. Valid examples are
#'   GM12878, K562.
#' @param tssRegion A vector with the region range to TSS.
#' @param output.dir The name of the directory where to write the selected annotations.
#' @param label The name of the file containing all annotation.
#'
#' @return The annotated "\code{chip.data}".
#'
#' @importFrom Biobase cache
#'
#' @export
annotate_chip <- function(chip.data, input.chrom.state, tf.regions, histone.regions, pol.regions, expression.levels,
                          genome.build = c("hg19", "mm9", "mm10", "hg38"), biosample = "GM12878",
                          tssRegion = c(-3000, 3000), output.dir, label) {
  dir.create(output.dir, recursive = TRUE)
  genome.build <- match.arg(genome.build)

  # Add an ID to every region.
  chip.data$ID = 1:length(chip.data)

  chip.data <- associate_genomic_region(chip.data, genome.build, output.dir, tssRegion = tssRegion)

  if(!is.null(input.chrom.state)) {
    chip.data = associate_chrom_state(chip.data, input.chrom.state)
  }

  if(!is.null(tf.regions)) {
    chip.data = associate_tf(chip.data, tf.regions)
  }

  if(!is.null(histone.regions)) {
    chip.data = associate_histone_marks(chip.data, histone.regions)
  }

  if(!is.null(pol.regions)) {
    chia.obj$Regions = associate_histone_marks(chia.obj$Regions, pol.regions)
  }

  chip.data = associate_gene_chip(chip.data, expression.data=NULL, biosample, genome.build)

  if(genome.build=="hg19" || genome.build=="hg38") {
    chip.data = associate_tissue_specificity_human(chip.data)
    chip.data = associate_fitness_genes(chip.data)
  }


  chip.data = associate_is_gene_active(chip.data)


  write.table(chip.data, file.path(output.dir, label), sep = "\t", row.names = FALSE)

  return(chip.data)
}

#' Annotates a set of regions with genomic features.
#'
#' Given a genome build identifier, creates a list of annotation databases
#' which can be passed to the annotation functions.
#'
#' @param region The regions to annotate.
#' @param annotations.list A list of annotation databases returned by
#' \code{\link{select_annotations}}.
#' @param tssRegion A vector with the region range to TSS.
#' @param filename The name of the file where the results should be saved.
#' If \code{NULL}, results are not saved to disk.
#' @return An annotation object.
#' @importFrom ChIPseeker annotatePeak
#' @export
annotate_region <- function(region, annotations.list, tssRegion = c(-3000, 3000), filename=NULL) {
    tfAnnotation = NULL
    if(length(region) > 0) {
        tfAnnotation <- ChIPseeker::annotatePeak(region,
                                                 tssRegion=tssRegion,
                                                 TxDb=annotations.list$TxDb,
                                                 annoDb=annotations.list$OrgDbStr)

        if(!is.null(filename)) {
            write.table(tfAnnotation,
                        file=filename,
                        sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
        }
    }
    return(tfAnnotation)
}

#' Given a \linkS4class{GRanges} object, finds enriched motifs within those regions.
#'
#' Performs motif enrichment against \code{pwm.bg}. If \code{pwm.bg} is not provided,
#' \code{annotations.list$PWMBG} is used.
#'
#' @param regions The regions on which motif enrichment should be performed.
#' @param annotations.list A list of annotation databases returned by
#'   \code{\link{select_annotations}}.
#' @param file.label A label for generating file names to save the results.
#'   If \code{NULL}, results are not saved to disk.
#' @param pwm.bg A PWMLogn background object against which enrichment
#'   should be performed.
#' @param top.x Number of top motifs for which logos are generated.
#' @return A list with the following elements: \describe{
#' \item{Region}{The subset of regions where motifs were sought.}
#' \item{Enrichment}{The result of the \pkg{PWMEnrich} \code{\link[PWMEnrich]{motifEnrichment}} call.}
#' \item{Report}{The result of the \pkg{PWMEnrich} \code{\link[PWMEnrich]{groupReport}} call.}}
#' @importFrom Biostrings getSeq
#' @importFrom Biostrings DNAStringSet
#' @importFrom PWMEnrich motifEnrichment
#' @importFrom PWMEnrich groupReport
#' @importFrom PWMEnrich plotMultipleMotifs
#' @export
motif_enrichment <- function(regions, annotations.list, file.label=NULL,  pwm.bg=NULL, top.x=20) {
  # Make sure we have a valid background.
  if(is.null(pwm.bg)) {
    if(is.null(annotations.list$PWMBG)) {
      stop("The provided annotations.list does not have a valid PWMBG element.")
    } else {
      pwm.bg = annotations.list$PWMBG
    }
  }

  # Get sequences for the given regions.
  intersectSeq <- Biostrings::getSeq(annotations.list$BSGenome, regions)

  # Remove N prefix/suffixes. Will also deal with all N sequences, which would cause a crash.
  intersectSeq <- Biostrings::DNAStringSet(gsub("N*$", "", gsub("^N*", "", as.character(intersectSeq))))

  # Remove sequences which are smaller than the maximum PWM length.
  max.pwm.length = max(unlist(lapply(pwm.bg$pwms, length)))
  sequence.subset = width(intersectSeq) >= max.pwm.length
  intersectSeq <- intersectSeq[sequence.subset]

  if(length(intersectSeq) > 0) {
    res <-  PWMEnrich::motifEnrichment(intersectSeq, pwm.bg)
    report <- PWMEnrich::groupReport(res)

    if(!is.null(file.label)) {
      # Plot top X motifs
      ordered.motifs = order(res$group.bg)
      for(i in 1:top.x) {
        # Perform some name sanitation.
        # For HOCOMOCO motif names.
        motif.name = gsub("Hsapiens-HOCOMOCOv9_AD_PLAINTEXT_H_PWM_hg19-", "", names(res$group.bg)[ordered.motifs[i]])

        # For mouse names.
        # Mouse motif names contain forward slashes, which are obviously not valid file name characters.
        motif.name = gsub("/", "", motif.name)

        # Generate logo file.
        pdf(paste(file.label, " ", i, " - ", motif.name, ".pdf"), width=7/1.5, height=11/6)
        PWMEnrich::plotMultipleMotifs(res$pwms[ordered.motifs[i]], xaxis=FALSE, titles="")
        dev.off()
      }

      # Write results to disc.
      write.table(as.data.frame(report), file=paste0(file.label, " MotifEnrichment.txt"),
                  sep="\t", row.names=FALSE, col.names=TRUE)
    }

    return(list(Region=regions[sequence.subset], Enrichment=res, Report=report))
  } else {
    warning("No sequences were eligeible for motif enrichment.")
    return(NULL)
  }
}

#' Given a set of Entrez IDS, retrieve the coordinates of those genes' promoters.
#'
#' @param selected.genes A vector of ENTREZ gene ids whose promoters should be
#'   retrieved.
#' @param annotations.list A list of annotation databases returned by
#'   \code{\link{select_annotations}}.
#' @param flank.size How many base pairs upstream of the TSS should we retrieve?
#' @return A \linkS4class{GRanges} objects representing the promoters of the given genes.
#' @importFrom AnnotationDbi select
#' @importFrom GenomicRanges reduce
#' @export
get_promoters <- function(selected.genes, annotations.list, flank.size=1000) {
  # Get the transcription regions from the database.
  tx.regions = AnnotationDbi::select(annotations.list$TxDb, selected.genes, c("TXCHROM", "TXSTART", "TXEND", "TXSTRAND"), "GENEID")

  # Keep only the first record for each gene.
  tx.regions = tx.regions[match(selected.genes, tx.regions$GENEID),]

  # Keep promoter only.
  promoter.regions = GenomicRanges::reduce(GenomicRanges::flank(GRanges(tx.regions), flank.size))

  return(promoter.regions)
}

#' Perform motif enrichment on the promoters of a set of genes.
#'
#' Utility function which performs the same operations as motif_enrichment,
#' but accepts a list of genes instead of a list of regions.
#'
#' @param selected.genes A vector of ENTREZ gene ids whose promoters should be
#'   subjected to motif enrichment.
#' @param annotations.list A list of annotation databases returned by
#'   \code{\link{select_annotations}}.
#' @param flank.size How many base pairs upstream of the TSS should we retrieve?
#' @param ... Additional arguments for \code{\link{motif_enrichment}}.
#' @return A list with the following elements: \describe{
#'   \item{Region}{The subset of regions where motifs were sought.}
#'   \item{Enrichment}{The result of the \pkg{PWMEnrich} \code{\link[PWMEnrich]{motifEnrichment call.}}}
#'   \item{Report}{The result of the \pkg{PWMEnrich} \code{\link[PWMEnrich]{groupReport call.}}}}
#' @export
motif_enrichment_genes <- function(selected.genes, annotations.list, flank.size=1000, ...) {
  promoter.regions = get_promoters(selected.genes, annotations.list, flank.size)

  return(motif_enrichment(promoter.regions, annotations.list=annotations.list, ...))
}

#' Perform KEGG pathway enrichment on a set of genes.
#'
#' Utility function which performs the same operations as \code{\link{motif_enrichment}},
#' but accepts a list of genes instead of a list of regions.
#'
#' @param selected.genes A vector of ENTREZ gene ids to be subjected
#'   to motif enrichment.
#' @param annotations.list A list of annotation databases returned by
#'   \code{\link{select_annotations}}.
#' @param filename The name of the file where the results should be saved.
#'    If \code{NULL}, results are not saved to disk.
#' @param diseases If \code{TRUE}, the enrichment is performed against the disease pathways.
#' @param gene.background A list of Entrez gene ids of the genes to be used
#'   as the background of the enrichment. If \code{NULL}, all genes in \code{annotations.list$TxDb}
#'   are used.
#' @return A data-frame with the enrichment results.
#' @importFrom AnnotationDbi keys
#' @export
kegg_enrichment <- function(selected.genes, annotations.list, filename=NULL, diseases=FALSE, gene.background=NULL) {
  # Retrieve all KEGG pathways.
  keptSets <- annotations.list$KEGG$kg.sets[annotations.list$KEGG$dise.idx]
  if(!diseases) {
    keptSets <- annotations.list$KEGG$kg.sets[annotations.list$KEGG$sigmet.idx]
  }

  if(is.null(gene.background)) {
    gene.background <- AnnotationDbi::keys(annotations.list$TxDb)
  }

  # For all pathways, perform enrichment analysis.
  inUniverse <- as.numeric(gene.background)
  inDataset <- as.numeric(selected.genes)

  # list.enrichment <- function(all.drawn, all.category, all.universe) {
  #     chosen <- sum(unique(all.drawn) %in% unique(all.category))
  #     universe <- length(unique(all.universe))
  #     possible <- length(unique(all.category))
  #     drawn <- length(unique(all.drawn))
  #     expected <- possible*(drawn/universe)
  #
  #     return(data.frame(Chosen=chosen,
  #                       Possible=possible,
  #                       Universe=universe,
  #                       Drawn=drawn,
  #                       Expected=expected,
  #                       PVal=phyper(chosen, possible, universe - possible, drawn, lower.tail=FALSE)))
  # }
  #
  # ldply(keptSets, list.enrichment, function(x) { list.enrichment(all.drawn=inDataset, x, inUniverse) })

  results <- data.frame(Pathway=character(0),
                        Chosen=numeric(0),
                        Possible=numeric(0),
                        Universe=numeric(0),
                        Drawn=numeric(0),
                        Expected=numeric(0),
                        PVal=numeric(0),
                        AdjPVal=numeric(0))
  for(kegg.set in names(keptSets)) {
    inPathway <- as.numeric(keptSets[[kegg.set]])

    chosen <- sum(unique(inDataset) %in% unique(inPathway))
    universe <- length(unique(inUniverse))
    possible <- length(unique(inPathway))
    drawn <- length(unique(inDataset))
    expected <- possible*(drawn/universe)

    # Perform the hypergeometric test.
    results <- rbind(results,
                     data.frame(Pathway=kegg.set,
                                Chosen=chosen,
                                Possible=possible,
                                Universe=universe,
                                Drawn=drawn,
                                Expected=expected,
                                PVal=phyper(chosen, possible, universe - possible, drawn, lower.tail=FALSE),
                                AdjPVal=1))
  }

  results$AdjPVal <- p.adjust(results$PVal, method="fdr")
  if(!is.null(filename)) {
    write.table(results, file=filename, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
  }

  return(results)
}

# Given a set of annotations, convert it to a set of Entrez gene ids before calling kegg_enrichment.
# kegg_enrichment.annotation <- function(annotations, file.name, diseases=FALSE, gene.background=NULL) {
#   selected.genes = unique(as.data.frame(annotations)$geneId)
#
#   return(kegg_enrichment(selected.genes, file.name, diseases, gene.background))
# }

#' Given a set of regions, convert it to a set of Entrez gene ids.
#'
#' Utility function to get a gene list from a set of regions.
#'
#' @param regions A \linkS4class{GRanges} object with regions to be associated to genes.
#' @param annotations.list A list of annotation databases returned by
#'   \code{\link{select_annotations}}.
#' @param flank.size The extent around the TSS which is considered the
#'   promoter region.
#' @param region.types A character vector representing region types which will
#'   cause a gene association. These can be: \describe{
#' \item{Promoter}{Promoter region (upstream and downstream of the TSS).}
#' \item{Gene body}{Promoter, Exons, introns, 5' and 3' UTRs.}
#' \item{Downstream}{Region downstream of the TES.}
#' \item{Distal}{Regions which do not fit any of the above category.}
#' \item{All}{All of the above.}}
#' @return A data-frame with the enrichment results.
#' @importFrom ChIPseeker annotatePeak
#' @export
gene_from_regions <- function(regions, annotations.list, flank.size=c(-3000, 3000), region.types=c("Gene body", "Promoter", "Downstream", "Distal", "All")) {

  region.types <- match.arg(region.types, several.ok = TRUE)

  # Annotate regions to retrieve gene names.
  overlap.annotation <- ChIPseeker::annotatePeak(regions,
                                                 tssRegion=flank.size,
                                                 TxDb=annotations.list$TxDb,
                                                 annoDb=annotations.list$OrgDbStr)

  if("All" %in% region.types) {
    region.types=c("Promoter", "Gene body", "Downstream", "Distal")
  }

  to.keep = rep(FALSE, length(regions))
  if("Promoter" %in% region.types) {
    to.keep = to.keep | grepl("Promoter", overlap.annotation@anno$annotation)
  }

  if("Gene body" %in% region.types) {
    to.keep = to.keep | grepl("Promoter", overlap.annotation@anno$annotation) | grepl("3' UTR", overlap.annotation@anno$annotation) |
      grepl("5' UTR", overlap.annotation@anno$annotation) | grepl("Exon", overlap.annotation@anno$annotation) | grepl("Intron", overlap.annotation@anno$annotation)

  }

  if("Downstream" %in% region.types) {
    to.keep = to.keep | grepl("Downstream", overlap.annotation@anno$annotation)
  }

  if("Distal" %in% region.types) {
    to.keep = to.keep | grepl("Distal", overlap.annotation@anno$annotation)
  }

  return(unique(overlap.annotation@anno$geneId[to.keep]))
}

#' Perform KEGG enrichment on a set of regions.
#'
#' Conveniance function. Regions are converted to genes with
#' \code{\link{gene_from_regions}} using default parameters.
#'
#' @param regions A \linkS4class{GRanges} object with regions to enriched for KEGG pathways.
#' @param annotations.list A list of annotation databases returned by
#'   \code{\link{select_annotations}}.
#' @param ... Parameters to be passed to \code{\link{kegg_enrichment}}.
#' @return A vector of Entrez gene ids containing a non-redundant list of the
#'   genes represented by the given regions.
#' @export
kegg_enrichment_regions <- function(regions, annotations.list, ...) {
  selected.genes = gene_from_regions(regions, annotations.list)

  return(kegg_enrichment(selected.genes, annotations.list=annotations.list, ...))
}

#' Generates a plot representing significant kegg pathways from an enrichment result.
#'
#' @param kegg.results.list A list of enrichment results returned by 
#'   \code{\link{kegg_enrichment}}.
#' @param filename The name of the file where the plot should be saved.
#' @param p.threshold Minimum p-value for a category to be plotted.
#' @param n.threshold Minimum number times a pathway is reported for it to be plotted.
#'
#' @export
multiple_keggs_plot <- function(kegg.results.list, filename, p.threshold = 0.05, n.threshold = 2) {
    # Identify all significant pathways.
    sig = rep(FALSE, nrow(kegg.results.list[[1]]))
    for(item.name in names(kegg.results.list)) {
        sig = sig | (kegg.results.list[[item.name]]$AdjPVal <= p.threshold & kegg.results.list[[item.name]]$Chosen >= n.threshold)
        kegg.results.list[[item.name]]$ResultSet = item.name
    }
    
    # Concatenate results of all significant pathways.
    plot.df = kegg.results.list[[1]][sig,]
    for(i in 2:length(kegg.results.list)) {
        plot.df = rbind(plot.df, kegg.results.list[[i]][sig,])
    }
    
    # Plot it.
    ggplot(plot.df, aes(x=ResultSet, y=Pathway, fill=-log10(AdjPVal), label=Chosen)) +
        geom_tile() +
        geom_text() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
    ggsave(filename, width=7, height=7)
}

#' Given a set of regions, perform annotation, motif enrichment and KEGG enrichment.
#'
#' @param region A \linkS4class{GRanges} object with regions to be characterized.
#' @param annotations.list A list of annotation databases returned by
#'   \code{\link{select_annotations}}.
#' @param file.label A label for generating file names to save the results.
#'   If \code{NULL}, results are not saved to disk.
#' @param skip.motif If \code{TRUE}, motif enrichment is skipped.
#' @param skip.kegg If \code{TRUE}, KEGG pathway enrichment is skipped.
#' @param output.dir The directory where output should be stored.A directory for storing output.
#' @return A list containing the characterization results.
#' @export
characterize_region <- function(region, annotations.list, skip.motif=FALSE, skip.kegg=FALSE, output.dir="output/") {
  results = list()

  if(is.null(output.dir)){
    file.label.annotation <- NULL
    file.label.motif <- NULL
    file.label.KEGG.sig <- NULL
    file.label.KEGG.dis <- NULL
  } else {
    dir.create(output.dir, showWarnings = FALSE, recursive = TRUE)
    file.label.annotation <- file.path(output.dir, "Annotations.txt")
    file.label.motif <- file.path(output.dir, "Motifs")
    file.label.KEGG.sig <- file.path(output.dir, "KEGG signalisation and metabolism.txt")
    file.label.KEGG.dis <- file.path(output.dir, "KEGG diseases.txt")
  }

  #results[["Annotation"]] = annotate_region(region, file.path(base.dir, "Annotations.txt"))
  results[["Annotation"]] = annotate_region(region, annotations.list, filename = file.label.annotation)

  if(!skip.motif) {
    #results[["Motif"]] = motif_enrichment(region, file.path(base.dir, "Motifs"), use.HOCOMOCO=(!exists("GENOME_VERSION") || GENOME_VERSION=="hg38"))
    results[["Motif"]] = motif_enrichment(region, annotations.list, file.label = file.label.motif)
  }

  if(!skip.kegg) {
    #results[["KEGG.sig"]] = kegg_enrichment.annotation(results[["Annotation"]], file.path(base.dir, "KEGG signalisation and metabolism.txt"))
    #results[["KEGG.dis"]] = kegg_enrichment.annotation(results[["Annotation"]], file.path(base.dir, "KEGG diseases.txt"), diseases=TRUE)
    results[["KEGG.sig"]] = kegg_enrichment_regions(region, annotations.list, filename = file.label.KEGG.sig)
    results[["KEGG.dis"]] = kegg_enrichment_regions(region, annotations.list, filename = file.label.KEGG.dis, diseases = TRUE)
  }

  return(results)
}

#' Given a set of Entrez gene ids, perform  KEGG enrichment and motif enrichement at the promoters.
#'
#' @param gene.set A vector of Entrez IDs representing the genes to be characterized.
#' @param annotations.list A list of annotation databases returned by
#'   \code{\link{select_annotations}}.
#' @param file.label A label for generating file names to save the results.
#'   If \code{NULL}, results are not saved to disk.
#' @param skip.motif If \code{TRUE}, motif enrichment is skipped.
#' @param skip.kegg If \code{TRUE}, KEGG pathway enrichment is skipped.
#' @param output.dir The directory where output should be stored.A directory for storing output.
#' @return A list containing the characterization results.
#' @export
characterize_gene_set <- function(gene.set, annotations.list, skip.motif=FALSE, skip.kegg=FALSE, output.dir="output/", promoter.size=1000) {
  results = list()

  if(is.null(output.dir)){
    file.label.motif <- NULL
    file.label.KEGG.sig <- NULL
    file.label.KEGG.dis <- NULL
  } else {
    dir.create(output.dir, showWarnings = FALSE, recursive = TRUE)
    file.label.motif <- file.path(output.dir, "Motifs")
    file.label.KEGG.sig <- file.path(output.dir, "KEGG signalisation and metabolism.txt")
    file.label.KEGG.dis <- file.path(output.dir, "KEGG diseases.txt")
  }

  if(!skip.motif) {
    #results[["Motif"]] = motif_enrichment_genes(gene.set, file.path(base.dir, "Motifs"), use.HOCOMOCO=(!exists("GENOME_VERSION") || GENOME_VERSION=="hg38"), flank.size=promoter.size)
    results[["Motif"]] = motif_enrichment_genes(gene.set, annotations.list, flank.size=promoter.size, file.label = file.label.motif)

  }

  if(!skip.kegg) {
    #results[["KEGG.sig"]] = kegg_enrichment(gene.set, file.path(base.dir, "KEGG signalisation and metabolism.txt"))
    #results[["KEGG.dis"]] = kegg_enrichment(gene.set, file.path(base.dir, "KEGG diseases.txt"), diseases=TRUE)
    results[["KEGG.sig"]] = kegg_enrichment(gene.set, annotations.list, filename = file.label.KEGG.sig)
    results[["KEGG.dis"]] = kegg_enrichment(gene.set, annotations.list, filename = file.label.KEGG.dis, diseases=TRUE)

  }

  return(results)
}


################################################################################################################
# Functions that may not work:
# Possible reasons:
# associate_is_gene_active: The conditions of active genes may not be observed because
#   1) The annotation for CDK9 is not provided;
#   2) The annotation of the Pol2 phosphoS2 is not provided;
#   3) The FPKM is not provided.
# etc.
################################################################################################################

#' Associates a boolean in fonction of the activity of the gene
#'
#' Is.Gene.Active is \code{TRUE} if FPKM > 1

#' Associate gene activity marker to a \linkS4class{GRanges} object.
#'
#' @param regions A \linkS4class{GRanges} object to annotate.
#'
#' @return The \linkS4class{GRanges} object with associated histone overlap percentage.
#' @export
associate_is_gene_active <- function(regions){
  regions.df <- as.data.frame(regions)
  #colname.CDK9 <- colnames(regions.df)[grep("CDK9", colnames(regions.df))][1]
  #colname.serin <- colnames(regions.df)[grep("phosphoS2", colnames(regions.df))][1]
  active.genes <- regions.df$SYMBOL[regions.df$Expr.mean > 1]
  #if (!is.na(colname.CDK9)){
  #  active.genes <- active.genes[active.genes %in% regions.df$SYMBOL[regions.df[,colname.CDK9] > 0]]
  #}
  #if (!is.na(colname.serin)){
  #  active.genes <- active.genes[active.genes %in% regions.df$SYMBOL[regions.df[,colname.serin] > 0]]
  #}
  regions$Is.Gene.Active <- regions$SYMBOL %in% active.genes
  return(regions)
}

#' Associate histone marks to a \linkS4class{GRanges} object.
#'
#' @param regions A \linkS4class{GRanges} object to annotate.
#' @param genome.build The name of the chosen annotation ("hg38", "mm9", "mm10", "hg19").
#' @param biosample The biosample identifier from ENCODE. Valid examples are
#'   GM12878, K562 or MCF-7.
#'
#' @return The \linkS4class{GRanges} object with associated histone overlap percentage.
#' @importMethodsFrom GenomicRanges findOverlaps range
#' @importMethodsFrom BSgenome width
#' @importFrom stats aggregate
#' @importFrom GenomicRanges reduce
#' @importFrom GenomicRanges mcols
#' @export
associate_histone_marks <- function(regions, histone.regions){
  all.columns <- region_overlap(regions, histone.regions)
  
  colnames(all.columns) <- paste0(colnames(histone.regions), ".OverlapPercentage")
  if(dim(mcols(regions))[2] == 0){
    mcols(regions) <- all.columns
  } else {
    mcols(regions) <- cbind(as.data.frame(mcols(regions)), all.columns)
  }

  return(regions)
}

#' Generate a matrix detailing the percentage of overlap between query regions and
#' a list of target regions.
#'
#' @param query.regions The GRanges object detailing the regions against which 
#'   the overlaps must be calculated.
#' @param target.regions A GRangesList object detailing the regions whose overlaps
#'   against each region in query.regions must be calculated.
#' @param proportion If TRUE, proportions are returned instead of base pair counts.
#'
#' @return A matrix detailing the overlap between teh query and the targets.
#' @importMethodsFrom GenomicRanges findOverlaps reduce ranges
#' @importMethodsFrom BSgenome width
#' @importFrom stats aggregate
#' @export
region_overlap <- function(query.regions, target.regions, proportion=TRUE) {
  # Initialize the overlap matrix.
  overlap.matrix <- matrix(0, nrow = length(query.regions), ncol = length(target.regions),
                           dimnames=list(NULL, names(target.regions)))
  for (target in names(target.regions)) {
    # Flatten the target regions so we won't double dip.
    target.flat <- reduce(target.regions[[target]])
    
    # Find the overlaps between query and target regions.
    overlap.index <- findOverlaps(query.regions, target.flat)

    if (length(overlap.index) > 0){
      # Get the overlap length for each target, the aggregate the results oper query region.
      overlap.length <- width(ranges(overlap.index, ranges(query.regions), ranges(target.flat)))
      overlap.total <- aggregate(overlap.length, by=list(From=from(overlap.index)), FUN=sum)
      overlap.matrix[overlap.total$From, target] <- overlap.total$x
    }
  }
  
  if(proportion) {
    overlap.matrix = apply(overlap.matrix, 2, '/', width(query.regions))
  }
  
  return(overlap.matrix)
}

#' Associate genomic regions to a \linkS4class{GRanges} object.
#'
#' @param regions A \linkS4class{GRanges} object to annotate.
#' @param genome.build The name of the chosen annotation ("hg38", "mm9", "mm10", "hg19").
#' @param tssRegion A vector with the region range to TSS.
#' @param output.dir The name of the directory where to write the selected annotations.
#'
#' @return The \linkS4class{GRanges} object with associated genomic regions.
#' @importFrom GenomicRanges mcols
#' @export
associate_genomic_region <- function(regions, genome.build, output.dir, tssRegion = c(-3000, 3000)) {
  # Annotate with proximity to gene regions
  annotations.list = select_annotations(genome.build)
  annotations = annotate_region(regions, annotations.list, tssRegion = tssRegion, file.path(output.dir, "CHIA-PET annotation.txt"))
  annotations.df = as.data.frame(annotations)
  mcols(regions) = annotations.df[, 6:ncol(annotations.df)]

  # Write porportions of annotated region types.
  write.table(annotations@annoStat, file.path(output.dir, "Annotation summary.txt"), sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)

  # Simplify genomic region annotation
  simplified.annotation = mcols(regions)$annotation
  simplified.annotation[grepl("Promoter", simplified.annotation)] <- "Promoter"
  simplified.annotation[grepl("Exon", simplified.annotation)] <- "Exon"
  simplified.annotation[grepl("Intron", simplified.annotation)] <- "Intron"
  simplified.annotation[grepl("Downstream", simplified.annotation)] <- "Downstream"
  regions$Simple.annotation = factor(simplified.annotation, levels=c("Distal Intergenic", "Promoter", "Intron", "Exon", "Downstream", "3' UTR", "5' UTR"))

  return(regions)
}

#' Associate genes to a \linkS4class{GRanges} object from ChIP-seq data.
#'
#' @param regions A \linkS4class{GRanges} object to annotate.
#' @param expression.data A data frame containing the levels of expression of genes, according to their EMSEMBL id.
#' @param biosample The biosample identifier from ENCODE. Valid examples are
#'   GM12878, K562 or MCF-7.
#' @param genome.build Which genome assembly should the results come from?
#'
#' @return The \linkS4class{GRanges} object with associated genes.
#' @importFrom plyr ddply mutate
#' @export
associate_gene_chip <- function(regions, expression.data=NULL, biosample, genome.build){
  # Associate a gene to a contact only if it's in the promoter.
  promoter.subset = regions$Simple.annotation == "Promoter"

  # Load isoform data
  downloaded.rna <- download_encode_rna(biosample, genome.build, isoform_download_filter_rna)

  # Subset the rna isoforms to only keep the isoforms with the highest expression
  promoter.subset <- as.data.frame(regions[promoter.subset])
  promoter.subset$ENSEMBL.transcript <- ucsc.to.ensembl$ENSEMBL[match(promoter.subset$transcriptId, ucsc.to.ensembl$UCSC)]
  promoter.subset$Mean.FPKM <- downloaded.rna$Expression$Mean.FPKM[match(promoter.subset$ENSEMBL.transcript, gsub("\\..", "", downloaded.rna$Expression$transcript_id))]
  max.FPKM.df <- aggregate(Mean.FPKM~ENSEMBL, data = promoter.subset, max)
  promoter.subset$Max.FPKM <- max.FPKM.df$Mean.FPKM[match(promoter.subset$ENSEMBL, max.FPKM.df$ENSEMBL)]
  isoform.subset = subset(promoter.subset, Mean.FPKM == Max.FPKM)

  # Further subset to keep the one closest to the TSS
  distance.df <- aggregate(abs(distanceToTSS)~ENSEMBL, data = promoter.subset, min)
  isoform.subset$min.distance <- distance.df$`abs(distanceToTSS)`[match(isoform.subset$ENSEMBL, distance.df$ENSEMBL)]
  distance.subset = subset(isoform.subset, distanceToTSS == min.distance)

  # If there are still more than one row, just pick the first one.
  gene.subset = aggregate(min.distance~ENSEMBL, data = distance.subset, function(x) { return(x[1]) })

  # Add the gene marker to the annotations.
  regions$Gene.Representative = regions$ENSEMBL %in% gene.subset$ENSEMBL

  if(!is.null(expression.data)) {
    index.match = match(regions$ENSEMBL, expression.data$ENSEMBL)
    regions$Expr.mean = expression.data$Mean.FPKM[index.match]
  }

  return(regions)
}

#' Given the name of a bed file containing a chromatin state map, loads it into a GRanges object.
#'
#' @param input.chrom.state The name of the chromatin state bed file.
#' @return A Granges object with the chromatin states.
#' @importFrom rtracklayer import
#' @export
load_chrom_state <- function(input.chrom.state) {
  # Load and rename chromatin states (for easier lexical ordering)
  chrom.states = import(input.chrom.state)
  chrom.states$name = gsub("^(.)_", "0\\1_", chrom.states$name)

  # Sort according to name, which will cause the first match to also be the
  # most relevant states (01_TSS first, 18_Quies last)
  chrom.states = chrom.states[order(chrom.states$name)]
  return(chrom.states)
}

#' Given a set of regions and a chromatin state map, match the "best" chromatin state to each region.
#'
#' @param regions The regions to which chromatin states should be associated.
#' @param chrom.states.gr A GRanges containing the chromatin state map.
#' @return A vector of the same length as regions, indicating the "best" chromatin state mapping to that region.
#' @importFrom GenomicRanges findOverlaps
#' @export
chrom_state_match <- function(regions, chrom.states.gr) {
  unique.states = sort(unique(chrom.states.gr$name))

  # Add state annotation to regions
  state.overlap.indices = findOverlaps(regions, chrom.states.gr, select="first")
  matching.chrom.state = factor(chrom.states.gr$name[state.overlap.indices], levels=unique.states)
    
  return(matching.chrom.state)
}

#' Associate chromatin sates to a \linkS4class{GRanges} object.
#'
#' @param regions A \linkS4class{GRanges} object to annotate.
#' @param input.chrom.state The path to a bed file representing a genome
#'   segmentation by chromatin state.
#' @return The \linkS4class{GRanges} object with associated chromatin states.
#' @importFrom rtracklayer import
#' @importMethodsFrom GenomicRanges findOverlaps
#' @export
associate_chrom_state <- function(regions, input.chrom.state) {
  # Load chromatin states
  chrom.states = load_chrom_state(input.chrom.state)
  matching = chrom_state_match(regions, chrom.states)
  regions$Chrom.State = matching

  return(regions)
}

#' Simplifies chromatin state names.
#'
#' Given a vector of canonical 18-state HMM-derived chromatin states,
#' return a vector of the same length where each state is associated
#' with a simplified version of itself.
#'
#' @param chrom.state.names The chromatin states to be simplified.
#' @return The simplified chromatin states.
#' @export
simplify_chrom_state <- function(chrom.state.names) {
    simple.chrom <- simple_chrom_state_map()
    return(factor(simple.chrom[chrom.state.names], levels=unique(simple_chrom_state_map())))
}

#' Returns the mapping used to go from full chromatin states to simple ones.
#'
#' @return The mapping used to go from full chromatin states to simple ones.
#' @export
simple_chrom_state_map <- function() {
    mapping = c("01_TssA"     = "TSS",
                "02_TssFlnk"  = "TSS", 
                "03_TssFlnkU" = "TSS",
                "04_TssFlnkD" = "TSS",
                "05_Tx"       = "Transcribed",
                "06_TxWk"     = "Transcribed", 
                "07_EnhG1"    = "Enhancer",
                "08_EnhG2"    = "Enhancer", 
                "09_EnhA1"    = "Enhancer",
                "10_EnhA2"    = "Enhancer",
                "11_EnhWk"    = "Enhancer",
                "12_ZNF/Rpts" = "Repressed",
                "13_Het"      = "Repressed", 
                "14_TssBiv"   = "TSS",
                "15_EnhBiv"   = "Enhancer",
                "16_ReprPC"   = "Repressed",
                "17_ReprPCWk" = "Repressed",
                "18_Quies"    = "Quiescent")
                
    mapping = factor(mapping, levels=c("TSS", "Transcribed", "Enhancer", "Repressed", "Quiescent"))
    
    return(mapping)
}

#' Associate TF overlaps to a \linkS4class{GRanges} object.
#'
#' @param regions A \linkS4class{GRanges} object to annotate.
#' @param tf.regions A data frame containing the TF data.
#'
#' @return The \linkS4class{GRanges} object with associated TF overlaps.
#' @importMethodsFrom GenomicRanges countOverlaps
#' @importFrom GenomicRanges mcols
#' @export
associate_tf <- function(regions, tf.regions) {
  # Calculate TF overlap with contact regions
  overlap.matrix = matrix(0, nrow=length(regions), ncol=length(tf.regions))
  colnames(overlap.matrix) <- paste0("TF.overlap.", names(tf.regions))
  for(i in 1:length(tf.regions)) {
    overlap.matrix[,i] <- countOverlaps(regions, tf.regions[[i]])
  }

  mcols(regions) = cbind(mcols(regions), as.data.frame(overlap.matrix))

  return(regions)
}

#' Associate tissue specificity of genes to a \linkS4class{GRanges} object.
#'
#' @param regions A \linkS4class{GRanges} object to annotate.
#'
#' @return The \linkS4class{GRanges} object with associated tissue scpecificity.
#' @export
associate_tissue_specificity_human <- function(regions) {

  # Annotate regions with Tau, Expression category.
  calculate.tau <- function(x) {
    return(sum(1 - (x/max(x))) / (length(x) - 1))
  }
  tissue.expression$Tau = apply(tissue.expression[,c(-1, -ncol(tissue.expression))], 1, calculate.tau)

  tissue.match = match(regions$ENSEMBL, tissue.expression$Ensembl.gene.id)
  regions$Expression.Category = factor(tissue.expression$Category[tissue.match],
                                                levels=c("Not detected",
                                                         "Mixed low", "Mixed high",
                                                         "Moderately tissue enriched", "Highly tissue enriched", "Group enriched",
                                                         "Expressed in all low", "Expressed in all high"))

  regions$Expression.Tau = tissue.expression$Tau[tissue.match]

  return(regions)
}

#' Identify essential genes of a \linkS4class{GRanges} object.
#'
#' @param regions A \linkS4class{GRanges} object to annotate.
#'
#' @return The \linkS4class{GRanges} object with identified essential genes.
#' @export
associate_fitness_genes <- function(regions) {

  # Add the "essential ratio" to the data
  fitness.match <- match(regions$SYMBOL, essential.genes$Gene)
  regions$Fitness <- essential.genes$numTKOHits[fitness.match]
  regions$Fitness <- ifelse(is.na(regions$Fitness), 0, (regions$Fitness / 6))

  return(regions)
}

#' Generate a GRanges object partitioning the genome into Promoter, Exon, introns,
#' downstream and distal intergenic regions.
#'
#' @param TxDb A TxDb object giving the locations of all genes for the given genome. 
#' @param available.genome A Granges object to be used as a background of all 
#'   the available genome.
#' @param BSgenome A BSgenome object to be used for generating an available genome
#'   background. Ignored if available.genome is provided.
#' @param flank.width The distance from the TSS/TES which are considered Promoters/Downstream.
#'
#' @return A GRanges object partitioning the genome into functional genomic regions.
#' @export
#' @importFrom GenomicFeatures transcriptsBy
#' @importFrom GenomicFeatures exonsBy
#' @importFrom GenomicFeatures intronsByTranscript
partition_genomic_regions <- function(TxDb, available.genome=NULL, BSgenome=NULL, flank.width=1000) {
    stopifnot(!is.null(BSgenome) | !is.null(available.genome))
    
    # Genomic region enrichment of the whole network.            
    geneList = unlist(GenomicFeatures::transcriptsBy(TxDb))
    
    # Generate "available.genome" regions containing the entire genome.
    if(is.null(available.genome)) {
      all.chr = unique(seqnames(geneList))
      all.chr.lengths = seqlengths(BSgenome)[all.chr]
      available.genome = GRanges(data.frame(seqnames=all.chr, start=1, end=all.chr.lengths, strand="+"))
    }
    
    region.list = list(Promoter=flank(geneList, width=flank.width, start=TRUE, both=TRUE),
                       Exon=unlist(GenomicFeatures::exonsBy(TxDb)),
                       Intron=unlist(GenomicFeatures::intronsByTranscript(TxDb)),
                       Downstream=flank(geneList, width=flank.width, start=FALSE, both=FALSE),
                       "Distal intergenic"=available.genome)
    genomic.regions = collapse_regions(region.list)
}


#' Generate a summary of the distances separating a set of genes from a set of peaks.
#'
#' @param peaks A GRanges object describing the peaks whose distance from the genes must be summarized.
#' @param genes A vector of gene symbols.
#' @param annot An annotation list returned by select.annotation for putting th esymbols in context.
#'
#' @return A numeric vector containing the proportion of genes within each distance category.
#' @export
gene_peaks_distances <- function(peaks, genes, annot) {
  # Annotate the NIPBL peaks.
  peaks.annotation = as.data.frame(annotatePeak(peaks, TxDb=annot$TxDb))
  
  # Associate symbols.
  peaks.annotation$SYMBOL = mapIds(annot$OrgDb, as.data.frame(peaks.annotation)$transcriptId, c("SYMBOL"), "UCSCKG")
  
  # Reorder them by distance to TSS.
  peaks.annotation = peaks.annotation[order(abs(peaks.annotation$distanceToTSS)),]
  
  # Keep closest to TSS
  peaks.best = peaks.annotation[match(unique(peaks.annotation$SYMBOL), peaks.annotation$SYMBOL),]
  peaks.best$Simple.Annotation = gsub(" \\(.*\\)", "", peaks.best$annotation)
  
  # Cross reference with gene list
  gene.peaks = peaks.best[peaks.best$SYMBOL %in% genes,]

  nb.peaks = table(gene.peaks$Simple.Annotation)
  nb.peaks["None"] = length(genes) - sum(nb.peaks)
  
  return(nb.peaks / sum(nb.peaks))
}  
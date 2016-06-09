#' Install required annotation packages.
#' 
#' Install BiocInstaller and the required annotation packages for the given
#' genome assemblies.
#'
#' @param which.assembly Which assembly's annotations should be installed.
#'   Supported assemblies are hg19, hg28, mm10 and mm9.
#' @export
install.annotations <- function(which.assembly) {
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
#' @return A list with the following elements:
#'   - *TxDb*: A TxDb of class AnnotationDBI providing informations about the genes
#'     of the selected build.
#'   - *OrgDb*: An OrgDb object for the selected species.
#'   - *OrgDbStr*: A character string representing the name of the OrgDb object.
#'   - *BSGenome*: A BSGenome object to retrieve DNA sequences
#'   - *KEGG*: A cache of KEGG pathways, as returned by kegg.gsets from the gage library.
#'   - *PWMBG*: A PWMLogn background for motif enrichment of promoter regions.
#' @export
select.annotations <- function(genome.build) {
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
            
            return(list(TxDb=TxDxDb.Mmusculus.UCSC.mm10.knownGene::TxDxDb.Mmusculus.UCSC.mm10.knownGene,
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

#' Annotates a set of regions with genomic features.
#' 
#' Given a genome build identifier, creates a list of annotation databases
#' which can be passed to the annotation functions.
#'
#' @param region The regions to annotate.
#' @param annotations.list A list of annotation databases returned by 
#'   select.annotations.
#' @param filename The name of the file where the results should be saved.
#    If NULL, results are not saved to disk.
#' @return An annotation object.
#' @export
annotate.region <- function(region, annotations.list, filename=NULL) {
    tfAnnotation = NULL
    if(length(region) > 0) {
        tfAnnotation <- ChIPseeker::annotatePeak(region,
                                                 tssRegion=c(-3000, 3000), 
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

#' Given a GRanges object, finds enriched motifs within those regions. 
#' 
#' Performs motif enrichment against pwm.bg. If pwm.bg is not provided,
#' annotations.list$PWMBG is used.
#'
#' @param regions The regions on which motif enrichment should be performed.
#' @param annotations.list A list of annotation databases returned by 
#'   select.annotations.
#' @param file.label A label for generating file names to save the results. 
#'   If NULL, results are not saved to disk.
#' @param pwm.bg A PWMLogn background object against which enrichment
#'   should be performed.
#' @param top.x Number of top motifs for which logos are generated.
#' @return A list with the following elements:
#'   - *Region*: The subset of regions where motifs were sought.
#'   - *Enrichment*: The result of the PWMEnrich motifEnrichment call.
#'   - *Report*: The result of the PWMEnrich groupReport call.
#' @export
motif.enrichment <- function(regions, annotations.list, file.label=NULL,  pwm.bg=NULL, top.x=20) {
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
#'   select.annotations.
#' @param flank.size How many base pairs upstream of the TSS should we retrieve?
#' @return A GRanges objects representing the promoters of the given genes.
#' @export
get.promoters <- function(selected.genes, annotations.list, flank.size=1000) {
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
#' Utility function which performs the same operations as motif.enrichment,
#' but accepts a list of genes instead of a list of regions.
#'
#' @param selected.genes A vector of ENTREZ gene ids whose promoters should be
#'   subjected to motif enrichment.
#' @param annotations.list A list of annotation databases returned by 
#'   select.annotations.
#' @param flank.size How many base pairs upstream of the TSS should we retrieve?
#' @param ... Additional arguments for motif.enrichment.
#' @return A list with the following elements:
#'   - *Region*: The subset of regions where motifs were sought.
#'   - *Enrichment*: The result of the PWMEnrich motifEnrichment call.
#'   - *Report*: The result of the PWMEnrich groupReport call.
#' @export
motif.enrichment.genes <- function(selected.genes, annotations.list, flank.size=1000, ...) {
    promoter.regions = get.promoters(selected.genes, annotations.list, flank.size)
    
    return(motif.enrichment(promoter.regions, annotations.list=annotations.list, ...))
}

#' Perform KEGG pathway enrichment on a set of genes.
#' 
#' Utility function which performs the same operations as motif.enrichment,
#' but accepts a list of genes instead of a list of regions.
#'
#' @param selected.genes A vector of ENTREZ gene ids to be subjected 
#'   to motif enrichment.
#' @param annotations.list A list of annotation databases returned by 
#'   select.annotations.
#' @param filename The name of the file where the results should be saved.
#    If NULL, results are not saved to disk.
#' @param diseases If true, the enrichment is performed against the disease pathways.
#' @param gene.background A list of Entrez gene ids of the genes to be used
#'   as the background of the enrichment. If NULL, all genes in annotations.list$TxDb
#'   are used.
#' @return A data-frame with the enrichment results.
#' @export
kegg.enrichment <- function(selected.genes, annotations.list, filename=NULL, diseases=FALSE, gene.background=NULL) {
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

# Given a set of annotations, convert it to a set of Entrez gene ids before calling kegg.enrichment.
# kegg.enrichment.annotation <- function(annotations, file.name, diseases=FALSE, gene.background=NULL) {
#   selected.genes = unique(as.data.frame(annotations)$geneId)                                   
# 
#   return(kegg.enrichment(selected.genes, file.name, diseases, gene.background))
# }

#' Given a set of regions, convert it to a set of Entrez gene ids.
#' 
#' Utility function to get a gene list from a set of regions.
#'
#' @param regions A GRanges object with regions to be associated to genes.
#' @param annotations.list A list of annotation databases returned by 
#'   select.annotations.
#' @param flank.size The extent around the TSS which is considered the 
#'   promoter region.
#' @param region.types A character vector representing region types which will
#'   cause a gene association. These can be:
#'     - *Promoter*: Promoter region (upstream and downstream of the TSS)
#'     - *Gene body*: Exons, introns, 5' and 3' UTRs.
#'     - *Downstream*: Region downstream of the TES.
#'     - *Distal*: Regions which do not fit any of the above category.
#'     - *All*: All of the above.
#'   disease pathways.
#' @return A data-frame with the enrichment results.
#' @export
gene.from.regions <- function(regions, annotations.list, flank.size=c(-3000, 3000), region.types=c("Gene body", "Promoter", "Downstream", "Distal", "All")) {
  
    region.types <- match.arg(region.types, several.ok = TRUE)
  
    # Annotate regions to retrieve gene names.
    overlap.annotation <- ChIPseeker::annotatePeak(regions,
                                                   tssRegion=flank.size, 
                                                   TxDb=annotations.list$TxDb, 
                                                   annoDb=annotations.list$OrgDbStr)
                                       
    if("All" %in% regions.types) {
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
#' gene.from.regions using default parameters.
#'
#' @param regions A GRanges object with regions to enriched for KEGG pathways.
#' @param annotations.list A list of annotation databases returned by 
#'   select.annotations.
#' @param ... Parameters to be passed to kegg.enrichment.
#' @return A vector of Entrez gene ids containing a non-redundant list of the 
#'   genes represented by the given regions.
#' @export
kegg.enrichment.regions <- function(regions, annotations.list, ...) {
    selected.genes = gene.from.regions(regions, annotations.list)
    
    return(kegg.enrichment(selected.genes, annotations.list=annotations.list, ...))
}

#' Given a set of regions, perform annotation, motif enrichment and KEGG enrichment.
#'  
#' @param region A GRanges object with regions to be characterized.
#' @param annotations.list A list of annotation databases returned by 
#'   select.annotations.
#' @param file.label A label for generating file names to save the results. 
#'   If NULL, results are not saved to disk.
#' @param skip.motif If TRUE, motif enrichment is skipped.
#' @param skip.kegg If TRUE, KEGG pathway enrichment is skipped.
#' @param output.dir The directory where output should be stored.A directory for storing output.
#' @return A list containing the characterization results.
#' @export
characterize.region <- function(region, annotations.list, skip.motif=FALSE, skip.kegg=FALSE, output.dir="output/") {
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

    #results[["Annotation"]] = annotate.region(region, file.path(base.dir, "Annotations.txt"))
    results[["Annotation"]] = annotate.region(region, annotations.list, filename = file.label.annotation)
    
    if(!skip.motif) {
        #results[["Motif"]] = motif.enrichment(region, file.path(base.dir, "Motifs"), use.HOCOMOCO=(!exists("GENOME_VERSION") || GENOME_VERSION=="hg38"))
      results[["Motif"]] = motif.enrichment(region, annotations.list, file.label = file.label.motif)
    }
    
    if(!skip.kegg) {
        #results[["KEGG.sig"]] = kegg.enrichment.annotation(results[["Annotation"]], file.path(base.dir, "KEGG signalisation and metabolism.txt"))
        #results[["KEGG.dis"]] = kegg.enrichment.annotation(results[["Annotation"]], file.path(base.dir, "KEGG diseases.txt"), diseases=TRUE)
        results[["KEGG.sig"]] = kegg.enrichment.regions(region, annotations.list, filename = file.label.KEGG.sig)
        results[["KEGG.dis"]] = kegg.enrichment.regions(region, annotations.list, filename = file.label.KEGG.dis, diseases = TRUE)
    }
    
    return(results)
}

#' Given a set of Entrez gene ids, perform  KEGG enrichment and motif enrichement at the promoters.
#'  
#' @param gene.set A vector of Entrez IDs representing the genes to be characterized.
#' @param annotations.list A list of annotation databases returned by 
#'   select.annotations.
#' @param file.label A label for generating file names to save the results. 
#'   If NULL, results are not saved to disk.
#' @param skip.motif If TRUE, motif enrichment is skipped.
#' @param skip.kegg If TRUE, KEGG pathway enrichment is skipped.
#' @param output.dir The directory where output should be stored.A directory for storing output.
#' @return A list containing the characterization results.
#' @export
characterize.gene.set <- function(gene.set, annotations.list, skip.motif=FALSE, skip.kegg=FALSE, output.dir="output/", promoter.size=1000) {
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
        #results[["Motif"]] = motif.enrichment.genes(gene.set, file.path(base.dir, "Motifs"), use.HOCOMOCO=(!exists("GENOME_VERSION") || GENOME_VERSION=="hg38"), flank.size=promoter.size)
        results[["Motif"]] = motif.enrichment.genes(gene.set, annotations.list, flank.size=promoter.size, file.label = file.label.motif)
        
    }
    
    if(!skip.kegg) {
        #results[["KEGG.sig"]] = kegg.enrichment(gene.set, file.path(base.dir, "KEGG signalisation and metabolism.txt"))
        #results[["KEGG.dis"]] = kegg.enrichment(gene.set, file.path(base.dir, "KEGG diseases.txt"), diseases=TRUE)
        results[["KEGG.sig"]] = kegg.enrichment(gene.set, annotations.list, filename = file.label.KEGG.sig)
        results[["KEGG.dis"]] = kegg.enrichment(gene.set, annotations.list, filename = file.label.KEGG.dis, diseases=TRUE)
        
    }
    
    return(results)
}

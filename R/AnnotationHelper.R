library(GenomicRanges)
library(VennDiagram)

# Perform KEGG enrichment
library(gage)

# Libraries to do motif search.
library(PWMEnrich)
library(MotifDb)

# If we're on guillimin, increase the core number.
#if(ON_SERVER) {
#    registerCoresPWMEnrich(8)
#}

# Load whichever libraries are necessary to annotate the genome of interest.
#if(exists("GENOME_VERSION") && GENOME_VERSION=="hg38") {
#    library(BSgenome.Hsapiens.UCSC.hg38)
#    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#    library(PWMEnrich.Hsapiens.background)
#    
#    txDb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#    orgDb <- "org.Hs.eg.db"
#    pwmBG = PWMLogn.hg19.MotifDb.Hsap
#    bsGenome = BSgenome.Hsapiens.UCSC.hg38
#} else if(exists("GENOME_VERSION") && GENOME_VERSION=="mm10") {
#    library(BSgenome.Mmusculus.UCSC.mm10)
#    library(TxDb.Mmusculus.UCSC.mm10.knownGene)    
#    library(PWMEnrich.Mmusculus.background)
#    data(PWMLogn.mm9.MotifDb.Mmus)
#    
#    txDb <- TxDb.Mmusculus.UCSC.mm10.knownGene
#    orgDb <- "org.Mm.eg.db"    
#    pwmBG = PWMLogn.mm9.MotifDb.Mmus
#    bsGenome = BSgenome.Mmusculus.UCSC.mm10
#} else if(exists("GENOME_VERSION") && GENOME_VERSION=="mm9") {
#    library(BSgenome.Mmusculus.UCSC.mm9)
#    library(TxDb.Mmusculus.UCSC.mm9.knownGene)    
#    library(PWMEnrich.Mmusculus.background)
#    data(PWMLogn.mm9.MotifDb.Mmus)
#    
#    txDb <- TxDb.Mmusculus.UCSC.mm9.knownGene
#    orgDb <- "org.Mm.eg.db"
#    pwmBG = PWMLogn.mm9.MotifDb.Mmus
#    bsGenome = BSgenome.Mmusculus.UCSC.mm9
#    
#} else {
#    library(BSgenome.Hsapiens.UCSC.hg19)
#    library(PWMEnrich.Hsapiens.background)
#    data(PWMLogn.hg19.MotifDb.Hsap)
#    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#
#    txDb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#    orgDb <- "org.Hs.eg.db"    
#    pwmBG = PWMLogn.hg19.MotifDb.Hsap
#    bsGenome = BSgenome.Hsapiens.UCSC.hg19
#}

library(ChIPseeker)

# Convert the annotation to a list for easier manipulation.
#list.txDb = as.list(txDb)

# # Load HOCOMOCO background.
# load(file.path(cache.dir, "bg.denovo.RData"))
# 
# # Load KEGG backgrounds. We can't rely on fetching them from the web when executing on a computing
# # node on guillimin/colosse.
# if(exists("GENOME_VERSION") && (GENOME_VERSION=="mm10" || GENOME_VERSION=="mm9")) {
#     load(file.path(cache.dir, "cached.mm.keggs.RData"))
# } else {
#     load(file.path(cache.dir, "cached.hs.keggs.RData"))
# }


# Given a GRanges object, annotates it using default parameters.
annotate.region <- function(region, filename) {
    tfAnnotation = NULL
    if(length(region) > 0) {
        tfAnnotation <- annotatePeak2(region,
                                      tssRegion=c(-3000, 3000), 
                                      TxDb=txDb, 
                                      annoDb=orgDb)
        
        write.table(tfAnnotation,
                    file=filename, 
                    sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
    }
    return(tfAnnotation)
}

# Given a GRanges object, finds enriched motifs within those regions.
motif.enrichment <- function(regions, file.label, use.HOCOMOCO=FALSE, top.x=20) {
    # Remove regions smaller than 30nt.
    regions.subset = regions[as.data.frame(regions)$width >= 30]
    
    intersectSeq <- getSeq(bsGenome, regions.subset)
    if(!use.HOCOMOCO) {
        res <-  motifEnrichment(intersectSeq, pwmBG)
    } else {
        res <-  motifEnrichment(intersectSeq, bg.denovo)
    }
    report <- groupReport(res)
    
    # Plot top X motifs
    ordered.motifs = order(res$group.bg)
    for(i in 1:top.x) {
        motif.name = gsub("Hsapiens-HOCOMOCOv9_AD_PLAINTEXT_H_PWM_hg19-", "", names(res$group.bg)[ordered.motifs[i]])
        # Mouse motif names contain forward slashes, which are obviously not valid file name characters.
        motif.name = gsub("/", "", motif.name)
        pdf(paste(file.label, " ", i, " - ", motif.name, ".pdf"), width=7/1.5, height=11/6)
        plotMultipleMotifs(res$pwms[ordered.motifs[i]], xaxis=FALSE, titles="")
        dev.off()
    }
    
    
    write.table(as.data.frame(report), file=paste0(file.label, " MotifEnrichment.txt"),
                sep="\t", row.names=FALSE, col.names=TRUE)
    
    return(list(Region=regions.subset, Enrichment=res, Report=report))
}


motif.enrichment.genes <- function(selected.genes, file.label, use.HOCOMOCO=FALSE, top.x=20, flank.size=1000) {
    which.tx.id =  list.txDb$genes$tx_id[list.txDb$genes$gene_id %in% selected.genes]
    ranges.df = list.txDb$transcripts[which.tx.id,c("tx_chrom", "tx_start", "tx_end", "tx_strand")]
    promoter.regions = reduce(flank(GRanges(ranges.df), flank.size))

    return(motif.enrichment(promoter.regions, file.label, use.HOCOMOCO, top.x))
}

kegg.enrichment <- function(selected.genes, file.name, diseases=FALSE, gene.background=NULL) {
  # Retrieve all KEGG pathways.
  # allKeggs <- kegg.gsets()
  keptSets <- allKeggs$kg.sets[allKeggs$dise.idx]
  if(!diseases) {
      keptSets <- allKeggs$kg.sets[allKeggs$sigmet.idx]
  }
  
  if(is.null(gene.background)) {
      gene.background <- unique(as.list(txDb)$genes$gene_id)
  }
  
  # For all pathways, perform enrichment analysis.
  inUniverse <- as.numeric(gene.background)
  inDataset <- as.numeric(selected.genes)

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
  write.table(results, file=file.name, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
  
  return(results)
}

# Given a set of annotations, convert it to a set of Entrez gene ids before calling kegg.enrichment.
kegg.enrichment.annotation <- function(annotations, file.name, diseases=FALSE, gene.background=NULL) {
  selected.genes = unique(as.data.frame(annotations)$geneId)                                   

  return(kegg.enrichment(selected.genes, file.name, diseases, gene.background))
}


# Given a set of regions, convert it to a set of annotations before calling kegg.enrichment.annotation.
kegg.enrichment.regions <- function(regions, file.name, diseases=FALSE, gene.background=NULL) {
  # Annotate regions to retrieve gene names.
  overlap.annotation <- annotatePeak(regions,
                                     tssRegion=c(-3000, 3000), 
                                     TxDb=txDb, 
                                     annoDb=orgDb)
                                     
  return(kegg.enrichment.annotation(overlap.annotation, file.name, diseases, gene.background))
}

# Given a GRanges object, annotate it  and find enriched motifs/KEGG pathways.
characterize.region <- function(region, label, skip.motif=FALSE, skip.kegg=FALSE, output.dir="output/") {
    results = list()

    base.dir = file.path(output.dir, label)
    dir.create(base.dir, showWarnings = FALSE, recursive = TRUE)

    results[["Annotation"]] = annotate.region(region, file.path(base.dir, "Annotations.txt"))
    
    if(!skip.motif) {
        results[["Motif"]] = motif.enrichment(region, file.path(base.dir, "Motifs"), use.HOCOMOCO=(!exists("GENOME_VERSION") || GENOME_VERSION=="hg38"))
    }
    
    if(!skip.kegg) {
        results[["KEGG.sig"]] = kegg.enrichment.annotation(results[["Annotation"]], file.path(base.dir, "KEGG signalisation and metabolism.txt"))
        results[["KEGG.dis"]] = kegg.enrichment.annotation(results[["Annotation"]], file.path(base.dir, "KEGG diseases.txt"), diseases=TRUE)
    }
    
    return(results)
}

# Given a set of Entrez gene ids, perform  KEGG enrichment and motif enrichement at the promoters.
characterize.gene.set <- function(gene.set, label, skip.motif=FALSE, skip.kegg=FALSE, output.dir="output/", promoter.size=1000) {
    results = list()

    base.dir = file.path(output.dir, label)
    dir.create(base.dir, showWarnings = FALSE, recursive = TRUE)

    if(!skip.motif) {
        results[["Motif"]] = motif.enrichment.genes(gene.set, file.path(base.dir, "Motifs"), use.HOCOMOCO=(!exists("GENOME_VERSION") || GENOME_VERSION=="hg38"), flank.size=promoter.size)
    }
    
    if(!skip.kegg) {
        results[["KEGG.sig"]] = kegg.enrichment(gene.set, file.path(base.dir, "KEGG signalisation and metabolism.txt"))
        results[["KEGG.dis"]] = kegg.enrichment(gene.set, file.path(base.dir, "KEGG diseases.txt"), diseases=TRUE)
    }
    
    return(results)
}

# Fix bug in annotatePeak whereas the id type is not found because the TxDb metadata
# is accessed by numerical index rather than name.
annotatePeak2 <- function (peak, tssRegion = c(-3000, 3000), TxDb = NULL, level = "transcript",
    assignGenomicAnnotation = TRUE, genomicAnnotationPriority = c("Promoter",
        "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"),
    annoDb = NULL, addFlankGeneInfo = FALSE, flankDistance = 5000,
    verbose = TRUE)
{
    level <- match.arg(level, c("transcript", "gene"))
    if (all(genomicAnnotationPriority %in% c("Promoter", "5UTR",
        "3UTR", "Exon", "Intron", "Downstream", "Intergenic")) ==
        FALSE) {
        stop("genomicAnnotationPriority should be any order of c(\"Promoter\", \"5UTR\", \"3UTR\", \"Exon\", \"Intron\", \"Downstream\", \"Intergenic\")")
    }
    if (is(peak, "GRanges")) {
        input <- "gr"
        peak.gr <- peak
    }
    else {
        input <- "file"
        peak.gr <- loadPeak(peak, verbose)
    }
    peakNum <- length(peak.gr)
    if (verbose)
        cat(">> preparing features information...\t\t", format(Sys.time(),
            "%Y-%m-%d %X"), "\n")
    TxDb <- loadTxDb(TxDb)
    if (level == "transcript") {
        features <- getGene(TxDb, by = "transcript")
    }
    else {
        features <- getGene(TxDb, by = "gene")
    }
    if (verbose)
        cat(">> identifying nearest features...\t\t", format(Sys.time(),
            "%Y-%m-%d %X"), "\n")
    idx.dist <- getNearestFeatureIndicesAndDistances(peak.gr,
        features)
    nearestFeatures <- features[idx.dist$index]
    if (verbose)
        cat(">> calculating distance from peak to TSS...\t",
            format(Sys.time(), "%Y-%m-%d %X"), "\n")
    distance <- idx.dist$distance
    peak.gr <- idx.dist$peak
    if (verbose)
        cat(">> assigning genomic annotation...\t\t", format(Sys.time(),
            "%Y-%m-%d %X"), "\n")
    if (assignGenomicAnnotation == TRUE) {
        anno <- getGenomicAnnotation(peak.gr, distance, tssRegion,
            TxDb, level, genomicAnnotationPriority)
        annotation <- anno[["annotation"]]
        detailGenomicAnnotation <- anno[["detailGenomicAnnotation"]]
    }
    else {
        annotation <- NULL
        detailGenomicAnnotation <- NULL
    }
    names(nearestFeatures) <- NULL
    nearestFeatures.df <- as.data.frame(nearestFeatures)
    if (level == "transcript") {
        colnames(nearestFeatures.df) <- c("geneChr", "geneStart",
            "geneEnd", "geneLength", "geneStrand", "geneId",
            "transcriptId")
        nearestFeatures.df$geneId <- TXID2EG(as.character(nearestFeatures.df$geneId),
            geneIdOnly = TRUE)
    }
    else {
        colnames(nearestFeatures.df) <- c("geneChr", "geneStart",
            "geneEnd", "geneLength", "geneStrand", "geneId")
    }
    if (!is.null(annotation))
        mcols(peak.gr)[["annotation"]] <- annotation
    for (cn in colnames(nearestFeatures.df)) {
        mcols(peak.gr)[[cn]] <- unlist(nearestFeatures.df[, cn])
    }
    mcols(peak.gr)[["distanceToTSS"]] <- distance
    if (!is.null(annoDb)) {
        if (verbose)
            cat(">> adding gene annotation...\t\t\t", format(Sys.time(),
                "%Y-%m-%d %X"), "\n")
        # Won't work. The Type of Gene ID key is now at row index 9!
        # IDType <- metadata(TxDb)[8, 2]
        IDType <- metadata(TxDb)$value[metadata(TxDb)$name == "Type of Gene ID"]
        geneAnno <- addGeneAnno(annoDb, peak.gr$geneId, type = IDType)
        if (!all(is.na(geneAnno))) {
            for (cn in colnames(geneAnno)[-1]) {
                mcols(peak.gr)[[cn]] <- geneAnno[, cn]
            }
        }
    }
    if (addFlankGeneInfo == TRUE) {
        if (verbose)
            cat(">> adding flank feature information from peaks...\t",
                format(Sys.time(), "%Y-%m-%d %X"), "\n")
        flankInfo <- getAllFlankingGene(peak.gr, features, flankDistance)
        mcols(peak.gr)[["flank_txIds"]] <- NA
        mcols(peak.gr)[["flank_geneIds"]] <- NA
        mcols(peak.gr)[["flank_gene_distances"]] <- NA
        mcols(peak.gr)[["flank_txIds"]][flankInfo$peakIdx] <- flankInfo$flank_txIds
        mcols(peak.gr)[["flank_geneIds"]][flankInfo$peakIdx] <- flankInfo$flank_geneIds
        mcols(peak.gr)[["flank_gene_distances"]][flankInfo$peakIdx] <- flankInfo$flank_gene_distances
    }
    if (verbose)
        cat(">> assigning chromosome lengths\t\t\t", format(Sys.time(),
            "%Y-%m-%d %X"), "\n")
    peak.gr@seqinfo <- seqinfo(TxDb)[names(seqlengths(peak.gr))]
    if (verbose)
        cat(">> done...\t\t\t\t\t", format(Sys.time(), "%Y-%m-%d %X"),
            "\n")
    if (assignGenomicAnnotation) {
        res <- new("csAnno", anno = peak.gr, tssRegion = tssRegion,
            level = level, hasGenomicAnnotation = TRUE, detailGenomicAnnotation = detailGenomicAnnotation,
            annoStat = getGenomicAnnoStat(peak.gr), peakNum = peakNum)
    }
    else {
        res <- new("csAnno", anno = peak.gr, tssRegion = tssRegion,
            level = level, hasGenomicAnnotation = FALSE, peakNum = peakNum)
    }
    return(res)
}
environment(annotatePeak2) <- asNamespace('ChIPseeker')
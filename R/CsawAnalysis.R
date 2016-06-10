library(csaw)
library(ggplot2)
library(edgeR)
library(GenomicRanges)

#' reads.to.counts
#' 
#' Converts reads from bam files to counts
#' 
#' @param bam.files, vector, path to bam files
#' @param param, vector, list of parameters
#' 
#' @return A list with the following elements :
#' [[1]] counted.reads, RangedSummarizedExperiment, output of WindowCounts
#' [[2]] plot.df, data frame, necessary information for the construction of a plot of window sizes
#' [[3]] correlate.reads, vector, carrelated reads
#' 
#' @export
reads.to.counts <- function(bam.files, param){
  # Setting up parameters
  no.dup <- reform(param, dedup=TRUE)
  
  # Finding average fragment length
  correlate.reads <- correlateReads(bam.files, max.dist = 500, param = no.dup)
  frag.length <- maximizeCcf(correlate.reads)
  
  # Assess window size
  collected <- list()
  window.widths = c()
  for(curbam in bam.files) {
    windowed <- windowCounts(curbam, spacing=50, width=150, ext = frag.length, param=no.dup, filter=10)
    rwsms <- rowSums(assay(windowed))
    maxed <- findMaxima(rowRanges(windowed), range=1000, metric=rwsms)
    collected[[curbam]] <- profileSites(curbam, rowRanges(windowed)[maxed],
                                        param=no.dup, weight=1/rwsms[maxed])
    
    # Returned widths are 1. Obviously a problem. Won't use downstream for now.
    window.widths[curbam] = wwhm(collected[[curbam]], rowRanges(windowed)[maxed], ext = frag.length)
  }
  
  plot.df = data.frame(Val=unlist(collected), 
                       Position=as.integer(unlist(lapply(collected, names))), 
                       File=rep(bam.files, times=lapply(collected, length)[1]))
  
  # Converting reads to counts
  counted.reads <- windowCounts(bam.files, spacing = 50, width = 150, ext = frag.length, filter = 10, param = param)
  
  # Returning counted.reads
  return (list(counted.reads, correlate.reads, plot.df))
}

#' filter
#' 
#' Filters out unintersing windows, by global enrichment
#' 
#' @param bam.files, vector, path to bam files
#' @param param, vector, list of parameters
#' @param counted.reads, RangedSummarizedExperiment, output of WindowCounts
#' @param threshold, numeric, threshold to consider in order to cut out windows
#' 
#' @return A list with the following elements :
#' [[1]] filtered.counted.reads, RangedSummarizedExperiment, filtered output of WindowCounts
#' [[2]] filter.stat, list, statistics of windows aboundances and their filter
#' 
#' @export
filter <- function(bam.files, param, counted.reads, threshold){
  # Filtering by global enrichment
  binned.filter <- windowCounts(bam.files, bin = TRUE, width = 2000L, param = param)
  filter.stat <- filterWindows(counted.reads, background = binned.filter, type = "global")
  filtered.counted.reads <- counted.reads[filter.stat$filter > log2(threshold),]
  
  # Returning filtered.counted.reads
  return(list(filtered.counted.reads, filter.stat))
}

#' normalization.factors
#' 
#' Calculates normalization factors, by using TMM method
#' 
#' @param bam.files, vector, path to bam files
#' @param param, vector, list of parameters
#' 
#' @return A list with the following elements :
#' [[1]] norm.factors, vector, normalization factors
#' [[2]] binned.normalize, RangedSummarizedExperiment, output of WindowCounts considered as "bin"
#' 
#' @export
normalization.factors <- function(bam.files, param){
  # Removing composition biases by using TMM method
  binned.normalize <- windowCounts(bam.files, bin = TRUE, width = 10000, param = param)
  norm.factors <- normOffsets(binned.normalize)
  
  # Returning norm.factors
  return(list(norm.factors, binned.normalize))
}

#' test.diff.binding
#' 
#' Tests differential binding
#' 
#' @param bam.files, vector, path to bam files
#' @param filtered.counted.reads, RangedSummarizedExperiment, filtered output of WindowCounts
#' @param norm.factors, vector, normalization factors
#' 
#' @return results, DGELRT, results of glmQLFTest
#' 
#' @export
test.diff.binding <- function(bam.files, filtered.counted.reads, norm.factors, correspondances, reference){
  # Setting up data
  setting.up <- asDGEList(filtered.counted.reads, norm.factors)
  
  # Creating a design matrix
  condition = relevel(correspondances$Factors, ref=reference)
  design <- model.matrix(~condition)
  
  # Estimating dispersions
  estimated.disp <- estimateDisp(setting.up, design)
  
  # If there are no replicates
  if (length(unique(correspondances$Factors)) == length(correspondances$Factors)){
    if (length(correspondances$Factors) == 2) contrast <- c(0, 1) else contrast <- NULL
    fit <- glmFit(estimated.disp, design, dispersion=0.05)
    results <- glmLRT(fit, contrast=contrast) 
    
  } else {
    fit <- glmQLFit(estimated.disp, design, robust = TRUE)
    
    # Testing for DB windows
    results <- glmQLFTest(fit, coef=2:ncol(design))
  }
  
  # Returning results
  return(results)
}

#' multiple.testing
#' 
#' Makes corrections for multiple testing
#' @param filtered.counted.reads, RangedSummarizedExperiment, filtered output of WindowCounts
#' @param results, DGELRT, results of glmQLFTest
#' 
#' @return A list with the following elements :
#' [[1]] merged, GRanges, result of mergeWindows
#' [[2]] tabcom, data frame, results of combined tests
#' 
#' @export
multiple.testing <- function(filtered.counted.reads, results, txdb){
  # Clustering windows into regions
  broads <- genes(txdb)
  broads <- resize(broads, width(broads) + 3000, fix = "end")
  overlaps <- findOverlaps(broads, rowRanges(filtered.counted.reads))
  table.broads <- combineOverlaps(overlaps, results$table)
  
  # Quick clustering
  merged <- mergeWindows(rowRanges(filtered.counted.reads), tol=1000L)
  tabcom <- combineTests(merged$id, results$table)
  
  # Summarizing the direction of differential binding per cluster
  is.sig <- tabcom$FDR <= 0.05
  has.up <- tabcom$logFC.up > 0
  has.down <- tabcom$logFC.down > 0
  cluster.summary <- data.frame(Total.DB=sum(is.sig), DB.up=sum(is.sig & has.up & !has.down),
                                DB.down=sum(is.sig & !has.up & has.down), DB.both=sum(is.sig & has.up & has.down))
  
  # Summarizing the direction of differential binding per cluster, for promoter-based clusters
  is.sig <- table.broads$FDR <= 0.05 & !is.na(table.broads$FDR)
  has.up <- table.broads$logFC.up > 0
  has.down <- table.broads$logFC.down > 0
  promoter.cluster.summary <- data.frame(Total.DB=sum(is.sig), DB.up=sum(is.sig & has.up & !has.down),
                                         DB.down=sum(is.sig & !has.up & has.down), DB.both=sum(is.sig & has.up & has.down))
  
  # Choosing representative windows, based on differential binding
  table.best <- getBestTest(merged$id, results$table)
  #Log-fold change of the best window in each cluster
  tabcom$best.start <- start(rowRanges(filtered.counted.reads))[table.best$best]
  #Same principle, with overlaps
  table.best.broad <- getBestOverlaps(overlaps, results$table)
  table.broads$best.start <- start(rowRanges(filtered.counted.reads))[table.best.broad$best]
  
  if (length(results$comparison) == 1){
    tabcom$best.logFC <- table.best$logFC
    table.broads$best.logFC <- table.best.broad$logFC
    logFC.best.windows <- tabcom[,c("best.logFC", "best.start")]
    logFC.best.windows.broads <- (table.broads[!is.na(table.broads$PValue),c("best.logFC", "best.start")])
    
  } else {
    comparison <- paste0("logFC.", results$comparison)
    for (i in comparison){
      new.col.name <- paste0("best.", i)
      tabcom[,new.col.name] <- table.best[,i]
      table.broads[,new.col.name] <- table.best.broad[,i]
    }
    logFC.best.windows <- tabcom[,c(paste0("best.", comparison), "best.start")]
    logFC.best.windows.broads <- (table.broads[!is.na(table.broads$PValue),c(paste0("best.", comparison), "best.start")])
    
  }
  
  # Returning merged and tabcom
  return(list(merged, tabcom, table.best, logFC.best.windows, table.best.broad, logFC.best.windows.broads))
}

#' post.processing
#' 
#' Saves results of MulipleTesting() with annotations
#' @param merged, GRanges, result of mergeWindows
#' @param tabcom, data frame, results of combined tests
#' 
#' @return annotations, list, annotation of considered regions
#' 
#' @export
post.processing <- function(merged, tabcom, txdb, orgdb){
  # Adding gene-based annotation
  annotations <- detailRanges(merged$region, txdb = txdb, orgdb = orgdb, promoter = c(3000, 1000), dist = 5000)
  
  return(annotations)
}

#' csaw.analyze
#' 
#' Analyzes bam files with csaw and produces corresponding graphs.
#' 
#' @param correspondances A data frame with two columns :
#'      $Name Name (with the path) of the bam files.
#'      $Factors The factors to considerate during the analysis. Files with the same factor are replicas.
#' @param reference The name of the refenrece factor.
#' @param genome.build The name of the chosen annotation ("hg38", "hg19", "mm9", "mm10").
#' @param output.dir The name of the directory where to save the files and graphs.
#' @param threshold The threshold to consider while filtering data.
#' 
#' @export
csaw.analyze <- function(correspondances, reference, genome.build, output.dir = "csaw_output/", threshold = 3) {
  correspondances <- read.table(correspondances, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  bam.files <- correspondances[,2]
  correspondances$Factors <- factor(correspondances$Factors)
  if (!reference %in% correspondances$Factors){
    stop("The reference is not part of the given factors.")
  }
  
  param <- readParam(pe="none", minq = 50)
  dir.create(output, recursive = TRUE)
  dir.create("csaw_cache")
  
  #Load appropriate libraries
  install.annotations(genome.build)
  annotations <- select.annotations(genome.build)
  txdb <- annotations$TxDb
  orgdb <- annotations$OrgDb
  
  if (!file.exists("csaw_cache/reads.to.counts.RData")){
    list.reads.to.counts <- reads.to.counts(bam.files, param)
    counted.reads <- list.reads.to.counts[[1]]
    correlate.reads <- list.reads.to.counts[[2]]
    plot.df <- list.reads.to.counts[[3]]
    # Saving data
    save(counted.reads, file = "csaw_cache/reads.to.counts.RData")
    
    pdf(paste0(output, "/Fragements lengths.pdf"))
    plot(0:500, correlate.reads, type="l", ylab="CCF", xlab="Delay (bp)")
    dev.off()
    
    ggplot(plot.df, aes(x=Position, y=Val, color=File)) +
      geom_line() +
      theme(legend.position="none")
    ggsave(paste0(output, "/Window size.pdf"))
    
    subset <- plot.df[plot.df$Position >= -50 & plot.df$Position <= 200,]
    ggplot(subset, aes(x=Position, y=Val, color=File)) +
      geom_line() +
      theme(legend.position="none")
    ggsave(paste0(output, "/Window size (zoom).pdf"))
    
  } else {
    load("csaw_cache/reads.to.counts.RData")
  }
  
  if (!file.exists("csaw_cache/filter.RData")){
    list.filter <- filter(bam.files, param, counted.reads, threshold)
    filtered.counted.reads <- list.filter[[1]]
    filter.stat <- list.filter[[2]]
    # Saving data
    save(filtered.counted.reads, file = "csaw_cache/filter.RData")
    
    # Illustrate effect of filtering in historgam
    plot.df=data.frame(log.CPM=filter.stat$back.abundances)
    global.bg <- filter.stat$abundances - filter.stat$filter
    annot.df = data.frame(Type=c("background", "threshold"),
                          x=c(global.bg[1], global.bg[1]+log2(threshold)))
    ggplot(plot.df, aes(x=log.CPM)) +
      geom_histogram(color="black", fill="lightblue") +
      geom_vline(data=annot.df, mapping=aes(xintercept=x, color=Type)) +
      geom_vline(data=annot.df, mapping=aes(xintercept=x, color=Type)) + 
      scale_colour_manual(values=c(background="red", threshold="blue"))
    ggsave(paste0(output, "/Effect of filtering.pdf"))
    
  } else {
    load("csaw_cache/filter.RData")
  }
  
  if (!file.exists("csaw_cache/normalization.factors.RData")){
    list.normalization.factors <- normalization.factors(bam.files, param)
    norm.factors <- list.normalization.factors[[1]]
    binned.normalize <- list.normalization.factors[[2]]
    # Saving data
    save(norm.factors, file = "csaw_cache/normalization.factors.RData")
    
    # Save normalization efforts with MA plots
    pdf(paste0(output, "/Normalization efforts.pdf"))
    par(mfrow=c(1, 3), mar=c(5, 4, 2, 1.5))
    adj.counts <- cpm(asDGEList(binned.normalize), log=TRUE)
    for (i in 1:(length(bam.files)-1)) {
      cur.x <- adj.counts[,1]
      cur.y <- adj.counts[,1+i]
      smoothScatter(x=(cur.x+cur.y)/2+6*log2(10), y=cur.x-cur.y, xlab="A", ylab="M", main=paste("1 vs", i+1))
      all.dist <- diff(log2(norm.factors[c(i+1, 1)]))
      abline(h=all.dist, col="red")
    }
    dev.off()
    
  } else {
    load("csaw_cache/normalization.factors.RData")
  }
  
  if (!file.exists("csaw_cache/test.diff.binding.RData")){
    results <- test.diff.binding(bam.files, filtered.counted.reads, norm.factors, correspondances, reference)
    # Saving data
    save(results, file="csaw_cache/test.diff.binding.RData")
    write.table(topTags(results)[[1]], file = paste0(output,"/TopTags results.txt"))
    
  } else {
    load("csaw_cache/test.diff.binding.RData")
  }
  
  if (!file.exists("csaw_cache/multiple.testing.RData")){
    list.multiple.testing <- multiple.testing(filtered.counted.reads, results, txdb)
    merged <- list.multiple.testing[[1]]
    tabcom <- list.multiple.testing[[2]]
    table.best <- list.multiple.testing[[3]]
    logFC.best.windows <- list.multiple.testing[[4]]
    table.best.broad <- list.multiple.testing[[5]]
    logFC.best.windows.broads <- list.multiple.testing[[6]]
    # Saving data
    save(merged, tabcom, file = "csaw_cache/multiple.testing.RData")
    
    write.table(table.best, file = paste0(output, "/Representative windows (differential binding).txt"), sep = "\t")
    write.table(logFC.best.windows, file = paste0(output, "/LogFC of best windows in each cluster.txt"), sep = "\t")
    write.table(table.best.broad, file = paste0(output, "/Representative windows (promoter-based).txt"), sep = "\t")
    write.table(logFC.best.windows.broads, file = paste0(output, "/LogFC of best windows in each cluster (promoter-based).txt"), sep = "\t")
    
  } else {
    load("csaw_cache/multiple.testing.RData")
  }
  
  if (!file.exists(file = paste0(output, "/csaw_clusters.gz.RData"))){
    annotations <- post.processing(merged, tabcom, txdb, orgdb)
    
    # Saving the results
    csaw.output <- gzfile(paste0(output, "/csaw_clusters.gz"), open="w")
    write.table(data.frame(as.data.frame(merged$region)[,1:3], tabcom, annotations), file = csaw.output,
                row.names = FALSE, quote = FALSE, sep = "\t")
    close(csaw.output)
  }
}
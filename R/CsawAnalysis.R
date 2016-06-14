#' Converts reads from bam files to counts
#'
#' @param bam.files A vector with the path to bam files.
#' @param param A vector with the list of parameters to give to the functions.
#'
#' @return A list with the following elements: \describe{
#'   \item{[[1]]} {A \linkS4class{RangedSummarizedExperiment} object returned by \code{\link[csaw]{WindowCounts}}.}
#'   \item{[[2]]} {A data frame with the necessary information for the construction of a plot of window sizes.}
#'   \item{[[3]]} {A vector with correlated reads.}}
#'
#' @import csaw
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

#' Filters out unintersing windows, by global enrichment method.
#'
#' @param bam.files A vector with the path to bam files.
#' @param param A vector with the list of parameters to give to the functions.
#' @param counted.reads A \linkS4class{RangedSummarizedExperiment} object returned by \code{\link[csaw]{WindowCounts}}.
#' @param threshold A numeric, the threshold to consider in order to cut out windows.
#'
#' @return A list with the following elements: \describe{
#'   \item{[[1]]} {A filtered \linkS4class{RangedSummarizedExperiment} object returned by \code{\link[csaw]{WindowCounts}}.}
#'   \item{[[2]]} {filter.stat, list, statistics of windows aboundances and their filter}}
#'
#' @importFrom csaw windowCounts
#' @importFrom csaw filterWindows
filter.reads.to.counts <- function(bam.files, param, counted.reads, threshold){
  # Filtering by global enrichment
  binned.filter <- windowCounts(bam.files, bin = TRUE, width = 2000L, param = param)
  filter.stat <- filterWindows(counted.reads, background = binned.filter, type = "global")
  filtered.counted.reads <- counted.reads[filter.stat$filter > log2(threshold),]

  # Returning filtered.counted.reads
  return(list(filtered.counted.reads, filter.stat))
}

#' Calculates normalization factors, by using TMM method.
#'
#' @param bam.files A vector with the path to bam files.
#' @param param A vector with the list of parameters to give to the functions.
#'
#' @return A list with the following elements: \describe{
#'   \item{[[1]]} {A vector with the normalization factors to use to normalize.}
#'   \item{[[2]]} {A \linkS4class{RangedSummarizedExperiment} object returned by \code{\link[csaw]{WindowCounts}} to consider as a "bin".}}
#'
#' @importFrom csaw windowCounts
#' @importFrom csaw normOffsets
normalization.factors <- function(bam.files, param){
  # Removing composition biases by using TMM method
  binned.normalize <- windowCounts(bam.files, bin = TRUE, width = 10000, param = param)
  norm.factors <- normOffsets(binned.normalize)

  # Returning norm.factors
  return(list(norm.factors, binned.normalize))
}

#' Tests the differential binding.
#'
#' @param bam.files A vector with the path to bam files.
#' @param filtered.counted.reads A filtered \linkS4class{RangedSummarizedExperiment} object returned by \code{\link[csaw]{WindowCounts}}.
#' @param norm.factors A vector with the normalization factors to use to normalize.
#'
#' @return A \linkS4class{DGELRT} object returned as results by \code{\link[edgeR]{glmQLFTest}}.
#'
#' @importFrom csaw asDGEList
#' @importFrom stats relevel
#' @importFrom stats model.matrix
#' @import edgeR
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

#' Makes corrections for multiple testing
#' @param filtered.counted.reads A filtered \linkS4class{RangedSummarizedExperiment} object returned by \code{\link[csaw]{WindowCounts}}.
#' @param results A \linkS4class{DGELRT} object returned as results by \code{\link[edgeR]{glmQLFTest}}.
#'
#' @return A list with the following elements: \describe{
#'   \item{[[1]]} {A \linkS4class{GRanges} object returned by \code{\link[csaw]{mergeWindows}}}
#'   \item{[[2]]} {A data frame returned by \code{\link[csaw]{combineTests}}}}
#'
#' @import GenomicRanges
#' @import csaw
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

#' Saves results of \code{\link{multiple.testing}} with annotations
#' @param merged A \linkS4class{GRanges} object returned by \code{\link[csaw]{mergeWindows}}
#' @param tabcom A data frame returned by \code{\link[csaw]{combineTests}}
#'
#' @return annotations, list, annotation of considered regions
#'
#' @export
post.processing <- function(merged, tabcom, txdb, orgdb){
  # Adding gene-based annotation
  annotations <- detailRanges(merged$region, txdb = txdb, orgdb = orgdb, promoter = c(3000, 1000), dist = 5000)

  return(annotations)
}

#' Analyzes bam files with csaw and produces corresponding graphs.
#'
#' @param correspondances A data frame with two columns: \describe{
#'   \item{$Name} {Name (with the path) of the bam files.}
#'   \item{$Factors} {The factors to considerate during the analysis. Files with the same factor are replicas.}}
#' @param reference The name of the reference factor.
#' @param genome.build The name of the chosen annotation ("hg38", "hg19", "mm9", "mm10").
#' @param output.dir The name of the directory where to save the files and graphs.
#' @param threshold The threshold to consider while filtering data.
#' @param save.plots Should plots and files created by \code{\link{csaw.save.plots}} be saved?
#'
#' @export
csaw.analyze <- function(correspondances, reference, genome.build, output.dir = "csaw_output/", threshold = 3, save.plots = TRUE) {
  correspondances <- read.table(correspondances, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  bam.files <- correspondances[,2]
  correspondances$Factors <- factor(correspondances$Factors)
  if (!reference %in% correspondances$Factors){
    stop("The reference is not part of the given factors.")
  }

  param <- readParam(pe="none", minq = 50)
  dir.create(output.dir, recursive = TRUE)
  dir.create("csaw_cache")

  #Load appropriate libraries
  install.annotations(genome.build)
  annotations <- select.annotations(genome.build)
  txdb <- annotations$TxDb
  orgdb <- annotations$OrgDb

  # reads.to.counts
  if (!file.exists("csaw_cache/reads.to.counts.RData")){
    list.reads.to.counts <- reads.to.counts(bam.files, param)
    # Saving data
    save(list.reads.to.counts, file = "csaw_cache/reads.to.counts.RData")
  } else {
    load("csaw_cache/reads.to.counts.RData")
  }
  counted.reads <- list.reads.to.counts[[1]]

  # filter.reads.to.counts
  if (!file.exists("csaw_cache/filter.RData")){
    list.filter <- filter.reads.to.counts(bam.files, param, counted.reads, threshold)
    # Saving data
    save(list.filter, file = "csaw_cache/filter.RData")
  } else {
    load("csaw_cache/filter.RData")
  }
  filtered.counted.reads <- list.filter[[1]]

  # normalization.factors
  if (!file.exists("csaw_cache/normalization.factors.RData")){
    list.normalization.factors <- normalization.factors(bam.files, param)
    # Saving data
    save(list.normalization.factors, file = "csaw_cache/normalization.factors.RData")
  } else {
    load("csaw_cache/normalization.factors.RData")
  }
  norm.factors <- list.normalization.factors[[1]]

  # test.diff.binding
  if (!file.exists("csaw_cache/test.diff.binding.RData")){
    results <- test.diff.binding(bam.files, filtered.counted.reads, norm.factors, correspondances, reference)
    # Saving data
    save(results, file="csaw_cache/test.diff.binding.RData")
  } else {
    load("csaw_cache/test.diff.binding.RData")
  }

  # multiple.testing
  if (!file.exists("csaw_cache/multiple.testing.RData")){
    list.multiple.testing <- multiple.testing(filtered.counted.reads, results, txdb)
    # Saving data
    save(list.multiple.testing, file = "csaw_cache/multiple.testing.RData")
  } else {
    load("csaw_cache/multiple.testing.RData")
  }
  merged <- list.multiple.testing[[1]]
  tabcom <- list.multiple.testing[[2]]

  # post.processing
  annotations <- post.processing(merged, tabcom, txdb, orgdb)

  # Save plots and files
  if (save.plots){
    csaw.save.plots(list.reads.to.counts, list.reads.filter, list.normalization.factors,
                    results, list.multiple.testing, annotations, output.dir)
  }
}

#' Creates and saves plots and tables for csaw analysis: \describe{
#'   \item{Fragements lengths.pdf} {A plot of fragments (reads) lengths in all bam files.}
#'   \item{Window size.pdf} {A cross-correlation plot of the efficiency of the ChIP-seq data.}
#'   \item{Window size (zoom).pdf} {A zommed version of "Window size.pdf".}
#'   \item{Effect of filtering.pdf} {A histogram of the counted reads, with lines representing background and filtering.}
#'   \item{Normalization efforts.pdf} {MA plots of the effect of normalization between the first file and every others.}
#'   \item{TopTags results.txt} {The best results of \code{\link{test.diff.binding}}.}
#'   \item{Representative windows (differential binding).txt}
#' {Table of the most representative windows, by "differential binding" clustering method, after the correction for multiple testing.}
#'   \item{LogFC of best windows in each cluster.txt}
#' {For each window cluster, the LogFC of the best windows, by "differential binding" clustering method, after the correction for multiple testing.}
#'   \item{Representative windows (promoter-based).txt}
#' {Table of the most representative windows, by "promoter-based" clustering method, after the correction for multiple testing.}
#'   \item{LogFC of best windows in each cluster (promoter-based).txt}
#' {For each window cluster, the LogFC of the best windows, by "promoter-based" clustering method, after the correction for multiple testing.}
#'   \item{csaw_clusters.gz} {The zipped table of results.}}
#'
#' @param list.reads.to.counts List returned by \code{\link{reads.to.counts}}.
#' @param list.reads.filter List returned by \code{\link{filter.reads.to.counts}}.
#' @param list.normalization.factors List returned by \code{\link{normalize.factors}}.
#' @param results Results of \code{\link{test.diff.binding}}.
#' @param list.multiple.testing List returned by \code{\link{multiple.testing}}
#' @param annotations The annotations to use.
#' @param output.dir The name of the directory where to save the files and graphs.
#'
#' @importFrom grDevices pdf
#' @importFrom graphics plot
#' @importFrom grDevices dev.off
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 ggsave
#' @importFrom utils write.table
csaw.save.plots <- function(list.reads.to.counts, list.reads.filter, list.normalization.factors,
                            results, list.multiple.testing, annotations, output.dir){
  correlate.reads <- list.reads.to.counts[[2]]
  plot.df <- list.reads.to.counts[[3]]
  filter.stat <- list.filter[[2]]
  binned.normalize <- list.normalization.factors[[2]]
  merged <- list.multiple.testing[[1]]
  tabcom <- list.multiple.testing[[2]]
  table.best <- list.multiple.testing[[3]]
  logFC.best.windows <- list.multiple.testing[[4]]
  table.best.broad <- list.multiple.testing[[5]]
  logFC.best.windows.broads <- list.multiple.testing[[6]]

  # Save results of reads.to.counts
  pdf(file.path(output.dir, "/Fragements lengths.pdf"))
  plot(0:500, correlate.reads, type="l", ylab="CCF", xlab="Delay (bp)")
  dev.off()

  ggplot(plot.df, aes(x=Position, y=Val, color=File)) +
    geom_line() +
    theme(legend.position="none")
  ggsave(file.path(output.dir, "Window size.pdf"))

  subset <- plot.df[plot.df$Position >= -50 & plot.df$Position <= 200,]
  ggplot(subset, aes(x=Position, y=Val, color=File)) +
    geom_line() +
    theme(legend.position="none")
  ggsave(file.path(output.dir, "Window size (zoom).pdf"))

  # Save results of filter.reads.to.counts
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
  ggsave(file.path(output.dir, "Effect of filtering.pdf"))

  # Save results of normalization.factors
  # Save normalization efforts with MA plots
  pdf(file.path(output.dir, "Normalization efforts.pdf"))
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

  # Save results of test.diff.binding
  write.table(topTags(results)[[1]], file = file.path(output.dir,"TopTags results.txt"))

  # Save results of multiple.testing
  write.table(table.best, file = file.path(output.dir, "Representative windows (differential binding).txt"), sep = "\t")
  write.table(logFC.best.windows, file = file.path(output.dir, "LogFC of best windows in each cluster.txt"), sep = "\t")
  write.table(table.best.broad, file = file.path(output.dir, "Representative windows (promoter-based).txt"), sep = "\t")
  write.table(logFC.best.windows.broads, file = file.path(output.dir, "LogFC of best windows in each cluster (promoter-based).txt"), sep = "\t")

  # Save results of post.processing
  # Saving the results
  csaw.output <- gzfile(file.path(output.dir, "csaw_clusters.gz"), open="w")
  write.table(data.frame(as.data.frame(merged$region)[,1:3], tabcom, annotations), file = csaw.output,
              row.names = FALSE, quote = FALSE, sep = "\t")
  close(csaw.output)
}

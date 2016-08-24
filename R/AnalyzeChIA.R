#' Create a heatmap by connectivity group.
#'
#' @param chia.obj A list containing the annotated ChIA-PET data, as returned by \code{\link{annotate.chia}}.
#' @param variable.name The name of the variable according to which the heatmap should be computed.
#' @param label The Name to give te the variable name in the resulting heatmap.
#' @param output.dir The name of the directory where to save the heatmaps.
#'
#' @return The matrix used to create the heatmap.
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 ggsave
#' @importFrom reshape2 melt
contact.heatmap <- function(chia.obj, variable.name, label, output.dir) {
  type.df = data.frame(Left=mcols(chia.left(chia.obj))[,variable.name],
                       Right=mcols(chia.right(chia.obj))[,variable.name],
                       stringsAsFactors=FALSE)

  var.levels = levels(mcols(chia.left(chia.obj))[,variable.name])
  results.matrix = matrix(NA, nrow=length(var.levels), ncol=length(var.levels), dimnames=list(var.levels, var.levels))

  for(i in 1:(length(var.levels))) {
    for(j in i:length(var.levels)) {
      count = sum(type.df$Left==var.levels[i] & type.df$Right==var.levels[j] |
                    type.df$Left==var.levels[j] & type.df$Right==var.levels[i],
                  na.rm=TRUE)
      results.matrix[i, j] = count
    }
  }

  results.df = melt(results.matrix)
  results.df = results.df[!is.na(results.df$value),]

  results.df$Var1 = factor(results.df$Var1, levels = rev(var.levels))
  ggplot(results.df, aes(y=Var1, x=Var2, fill=value)) +
    geom_tile() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.key.size = unit(1, "cm"))
  ggsave(file.path(output.dir, paste0("Contact heatmap for ", label, ".pdf")))

  return(results.matrix)
}

#' Produces boxplots for every transcription factor, in relation to its 3D connectivity
#'
#' For every transcription factor in the ChIP data, creates a boxplot of the force of the ChIP-seq signal
#' in function of the contact frequence of the region.
#'
#' @param chip.data A \linkS4class{GRangesList} containing the regions of the ChIp-seq data, with signal values.
#' @param biosample The biosample identifier from ENCODE. Valid examples are GM12878, K562 or MCF-7.
#' @param genome.build The name of the chosen annotation ("hg38", "hg19").
#' @param chia.obj Annotated ChIA-PET data, as returned by \link{analyze.chia.pet} or \link{annotate.chia}.
#' @param output.dir The directory where to write the boxplots.
#' @param TSS Should only the TSS regions be kept?
#' @param tssRegion tssRegion A vector with the region range to TSS.
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_boxplot
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 ggsave
#' @importFrom cowplot plot_grid
#' @importFrom GenomicRanges flank
#' @importFrom GenomicFeatures genes
#' @export
boxplot.per.tf <- function(chip.data, biosample, genome.build, chia.obj, output.dir, TSS = TRUE, tssRegion = c(-3000, 3000)) {

  # Extract ChIA-PET regions
  chia.data <- chia.obj$Regions

  # Exctract all TF
  if (genome.build == "hg19"){
    TxDb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
  } else if (genome.build == "hg38"){
    TXDB <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  } else if (genome.build == "mm9"){
    TxDb <- TxDxDb.Mmusculus.UCSC.mm9.knownGene::TxDxDb.Mmusculus.UCSC.mm9.knownGene
  } else if (genome.build == "mm10"){
    TxDb <- TxDxDb.Mmusculus.UCSC.mm10.knownGene::TxDxDb.Mmusculus.UCSC.mm10.knownGene
  } else {
    stop("Error : genome.build is invalid.")
  }

  tss.regions <- genes(TxDb)
  tss.regions <- promoters(tss.regions, abs(tssRegion[1]), tssRegion[2])


  # Function to create boxplot with histogram
  create.boxplot <- function(chip.data, chia.data, label, label.x, label.y, output.dir, tss.regions, TSS=TRUE){
    if (TSS){
      chip.data <- chip.data[chip.data$distanceToTSS == 0]
      chia.data <- chia.data[chia.data$distanceToTSS == 0]
    }
    indices <- GenomicRanges::findOverlaps(chip.data, chia.data)
    if (length(indices) != 0) {

      signal.degree.df <- data.frame(Signal = log2(chip.data$signalValue[indices@from]),
                                     Degree = chia.data$Degree[indices@to])
      signal.degree.df$CutDegree <- cut(signal.degree.df$Degree, breaks = c(1, 5, 10, 20, 40, Inf), right = FALSE)
      if (length(chip.data$signalValue[-indices@from]) != 0){
        signal.degree.df <- rbind(data.frame(Signal = log2(chip.data$signalValue[-indices@from]),
                                             Degree = 0, CutDegree = "0"), signal.degree.df)
      }
      box <- ggplot(signal.degree.df) + geom_boxplot(aes(CutDegree, Signal)) + ylab(label.y) + xlab(label.x) + ggtitle(label)


      if (TSS) {
        tss.indices <- findOverlaps(tss.regions, chia.data)
        tss.degree.df <- data.frame(Region = tss.indices@from, Degree = chia.data$Degree[tss.indices@to])
        tss.degree.df$CutDegree <- cut(tss.degree.df$Degree, breaks = c(1, 5, 10, 20, 40, Inf), right = FALSE)
        if (length(chip.data$signalValue[-indices@from]) != 0){
          tss.degree.df <- rbind(data.frame(Region = c(1:length(tss.regions))[-tss.indices@from],
                                            Degree = 0, CutDegree = "0"), tss.degree.df)
        }

        mapping.df <- data.frame(x.pos = levels(tss.degree.df$CutDegree), y.pos = as.vector(table(signal.degree.df$CutDegree)),
                                 label = paste0(round(as.vector(table(signal.degree.df$CutDegree) / table(tss.degree.df$CutDegree))*100, digits = 1), "%"))
        hist <- ggplot() +
          geom_bar(aes(CutDegree), data = tss.degree.df, fill = "light blue") +
          geom_bar(aes(CutDegree), data = signal.degree.df) +
          geom_text(data = mapping.df, aes(x = x.pos, y = y.pos, label = label), vjust = -1) +
          xlab(label.x) +
          ggtitle(paste("Histogram of ", label.x))


      } else {
        hist <- ggplot(signal.degree.df) + geom_bar(aes(CutDegree)) + xlab(label.x) + ggtitle(paste("Histogram of ", label.x))
      }

      plot_grid(box, hist, nrow = 2, align = "v")
      ggsave(file.path(output.dir, paste0(label, ".pdf")), height = 14, width = 7)
    }
  }

  # Create plot for every TF
  dir.create(file.path(output.dir, biosample), recursive = TRUE, showWarnings=FALSE)
  for (tf in names(chip.data)){
    cat("Factor : ", tf, "\n")
    chip.subset <- chip.data[tf][[1]]
    chip.subset <- annotate.chip(chip.subset, input.chrom.state = NULL, tf.regions = NULL, histone.regions = NULL,
                                 pol.regions = NULL, expression.levels = NULL, genome.build = genome.build,
                                 biosample = biosample, tssRegion = tssRegion, output.dir = "output/annotations",
                                 label = paste0(tf, " annotated"))

    if (TSS){
      create.boxplot(chip.subset, chia.data,
                     paste0("Boxplot of log2(Signal) in fct of Degree of ", tf ," at TSS"),
                     "Contact frequency at TSS", "log2(Signal)", file.path(output.dir, biosample), tss.regions, TSS = TSS)
    } else {
      create.boxplot(chip.subset, chia.data,
                     paste0("Boxplot of log2(Signal) in fct of Degree of ", tf),
                     "Contact frequency", "log2(Signal)", file.path(output.dir, biosample), tss.regions, TSS = TSS)
    }

  }
}

#' Boxplot by connectivity
#'
#' Creates a boxplot of a specific value (given as parameter) according to the connectivity
#'
#' @param chia.obj A list containing the annotated ChIA-PET data, as returned by \code{\link{annotate.chia}}.
#' @param variable.name The name of the variable according to which the boxplot should be computed.
#' @param label The Name to give te the variable name in the resulting heatmap.
#' @param output.dir The name of the directory where to save the heatmaps.
#'
#' @import ggplot2
#'
#' @export
boxplot.by.connectivity <- function(chia.obj, variable.name, label, output.dir){
  data.for.boxplot <- as.data.frame(mcols(chia.obj$Regions)[, c(variable.name, "Degree")])
  data.for.boxplot$CutDegree <- cut(data.for.boxplot$Degree, breaks = c(1, 2, 6, 21, Inf), right = FALSE,
                                    labels = c("Singles", "Low", "Intermediate", "High"))
  data.for.boxplot <- data.for.boxplot[!is.na(data.for.boxplot[,variable.name]),]
  ggplot(data.for.boxplot) +
    geom_boxplot(aes(CutDegree, get(variable.name))) +
    ylab(variable.name) +
    xlab("Connectivity") +
    scale_y_log10() +
    ggtitle(label)
  ggsave(file.path(output.dir, paste0(label, ".pdf")))
}

#' Histogram of the proportion of "essential genes", by connectivity categories.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{annotate.chia}}.
#' @param essential.genes The ID's of the genes to anlayze. The ID's refer the the ChIA-PET ID's of the nodes.
#' @param label.x The title of the histogram.
#' @param output.dir The name of the directory where to save the plot.
#' @import ggplot2
#' @importFrom GenomicRanges mcols
#' @export
histogram.essential.genes <- function(chia.obj, essential.genes, label.x, output.dir){
  # Convert to data.frame
  chia.annotated.df <- as.data.frame(mcols(chia.obj$Regions))
  # Cut the size into categories
  chia.annotated.df$CutDegree <- cut(chia.annotated.df$Degree, breaks = c(1, 2, 6, 21, Inf), right = FALSE,
                                     labels = c("Singles", "Low", "Intermediate", "High"))
  # Keep only the ids that represent "essential genes"
  chia.essential <- chia.annotated.df[chia.annotated.df$ID %in% essential.genes,]
  # Prepare the mapping to add the percentage label over the bars of the histogram
  mapping.df <- data.frame(x.pos = levels(chia.essential$CutDegree), y.pos = as.vector(table(chia.essential$CutDegree)),
                           label = paste0(round(as.vector(table(chia.essential$CutDegree) / table(chia.annotated.df$CutDegree[chia.annotated.df$Gene.Representative]))*100, digits = 1), "%"))
  # Create the histogram
  ggplot(chia.essential) +
    geom_bar(aes(CutDegree)) +
    geom_text(data = mapping.df, aes(x = x.pos, y = y.pos, label = label), vjust = -1) +
    xlab("connectivity") +
    ggtitle(paste("Histogram of ", label.x))
  ggsave(file.path(output.dir, paste0("Histogram of ", label.x, ".pdf")))
}

#' Plot heatmaps in fonction of a ChIA-PET annotation
#'
#' Creates heatmaps of the proportion of a variable from the ChIA-PET annotation, given as parameter, in fonction of the networks.
#' To use the transcription factors as variable, write "TF". This will create the same heatmaps, but with presence or absence of each factor instead of proportions.
#'
#' @param chia.obj A list containing the annotated ChIA-PET data, as returned by \code{\link{annotate.chia}}.
#' @param size.limit The networks must have a size over size.limit.
#' @param variable.name Name of the column to use, or "TF" for transcription factors.
#' @param label The label to add to the heatmap title (name of the variable).
#' @param output.dir The name of the directory where output should be saved.
#'
#' @importFrom plyr ddply
#' @importFrom reshape2 dcast
#' @importFrom NMF aheatmap
#'
#' @export
chia.plot.network.heatmap <- function(chia.obj, size.limit, variable.name, label, output.dir) {
  dir.create(output.dir, recursive = TRUE, showWarnings=FALSE)

  presence.by.tf <- function(chia.df){
    tf.col <- grep("TF.overlap", colnames(chia.df))
    tf.presence <- apply(chia.df[,tf.col], 2, sum)
    tf.presence <- ifelse(tf.presence > 0, 1, 0)
    names(tf.presence) <- sub("TF.overlap.", "", colnames(chia.df)[tf.col])
    return(tf.presence)
  }

  if (variable.name == "TF") {
    per.component.cs = ddply(as.data.frame(chia.obj$Regions[chia.obj$Regions$Component.size>size.limit]), ~Component.Id,
                             presence.by.tf)

  } else {
    per.component.cs = ddply(as.data.frame(chia.obj$Regions[chia.obj$Regions$Component.size>size.limit]), ~Component.Id,
                             function(x) { as.data.frame(table(x[, variable.name])/nrow(x)) })
    # Un-melt it.
    per.component.cs = dcast(per.component.cs, Component.Id~Var1)
  }

  component.ids = per.component.cs$Component.Id

  per.component.cs.matrix = per.component.cs[,-1]
  component.sizes = chia.obj$Regions$Component.size[match(component.ids, chia.obj$Regions$Component.Id)]
  pdf(file.path(output.dir, paste0("Heatmap of the ", label, " - networks with more than ", size.limit, " nodes.pdf")), width=7, height=7)
  aheatmap(t(per.component.cs.matrix), annCol=list("log2(Network size)"=log2(component.sizes)), color="-Blues")
  dev.off()


  pdf(file.path(output.dir, paste0("Heatmap of the ", label, " - networks with more than ", size.limit, " nodes no clustering.pdf")), width=7, height=7)
  aheatmap(t(per.component.cs.matrix[order(component.sizes),]), annCol=list("log2(Network size)"=log2(component.sizes)[order(component.sizes)]), Colv=NA, color="-Blues")
  dev.off()

  per.component.cs.matrix[per.component.cs.matrix > 0.4] = 0.3
  pdf(file.path(output.dir, paste0("Heatmap of the ", label, " - networks with more than ", size.limit, " nodes limit 0.3.pdf")), width=7, height=7)
  aheatmap(t(per.component.cs.matrix), annCol=list("log2(Network size)"=log2(component.sizes)), color="-Blues")
  dev.off()


  pdf(file.path(output.dir, paste0("Heatmap of the ", label, " - networks with more than ", size.limit, " nodes no clustering limit 0.3.pdf")), width=7, height=7)
  aheatmap(t(per.component.cs.matrix[order(component.sizes),]), annCol=list("log2(Network size)"=log2(component.sizes)[order(component.sizes)]), Colv=NA, color="-Blues")
  dev.off()
}
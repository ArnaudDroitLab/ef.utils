#' Perform enrichments by connectivity group.
#'
#' @param chia.obj A list containing the annotated ChIA-PET data, as returned by annotate.chia.
#' @param variable.name The name of the variable according to which the enrichment should be performed.
#' @param label The Name to give te the variable name in the resulting graph.
#' @param wrap.row The number of rows in the final file.
#' @param wrap.col The number of columns in the final file.
#' @param output.dir The name of the directory where to save the graphs.
#'
#' @return The plots of enrichments by connectivity group.
#'
#' @importFrom GenomicRanges mcols
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 ggsave
connectivity.enrichment <- function(chia.obj, variable.name, label, wrap.row=3, wrap.col=6, output.dir="output") {
  # Define thresholds for different connectivity categories.
  thresholds.list = list(Singles=c(0, 1), Low=c(1, 5), Intermediate=c(5, 20), High=c(20, 1000))

  all.values = mcols(chia.obj$Regions)[,variable.name]
  results = matrix(0, nrow=length(levels(all.values)), ncol=length(thresholds.list))
  rownames(results) = levels(all.values)
  colnames(results) = names(thresholds.list)

  # Loop over all thresholds, calculating the frequency of the given variable
  # for each category.
  for(i in 1:length(thresholds.list)) {
    thresholds = thresholds.list[[i]]
    indices = chia.obj$Regions$Degree > thresholds[1] & chia.obj$Regions$Degree <= thresholds[2]
    proportions = table(all.values[indices]) / sum(indices)
    results[,i] <- proportions[rownames(results)]
  }

  melt.df = melt(results, factorsAsStrings=FALSE)
  colnames(melt.df) <- c("Var.of.interest", "Connectivity", "Proportion")
  melt.df$Connectivity = factor(melt.df$Connectivity, levels = names(thresholds.list))

  ggplot(melt.df, aes(x=Connectivity, y=Proportion)) +
    geom_line(group=1) +
    facet_wrap(~Var.of.interest, nrow=wrap.row, ncol=wrap.col) +
    ylab("Proportion of nodes in category") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(file.path(output.dir, paste0("Proportion of ", label, " as a function of connectivity category.pdf")))

  return(results)
}

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

#############################################################################################################
# Functions that may not work:
# Possible reasons:
# analyze.chromatin.states: The column "Chrom.State" is not in chia.obj$Regions.
# analyze.annotation: The column "Simple.annotation" is not in chia.obj$Regions.
# analyze.expression: The column "Expr.mean" is not in chia.obj$Regions.
# etc.
#############################################################################################################

#' Analyze the chromatin state of ChIA-PET data
#'
#' Analyzes ChIA-PET data and produces graphs according to its chromatin states: \describe{
#' \item{Chromatin states summary.txt}{A file with the summary of chromatin states.}
#' \item{log Degree histogram per chromatin state.pdf}{A histogram of the degree of each node according to the chromatin state.}
#' \item{Proportion of chromatin state as a function of connectivity category.pdf}{A plot of the proportion of chromatin state as a fonction of connectivity category.}
#' \item{Contact heatmap for chromatin states.pdf}{A contact heatmap of chromatin states.}}
#'
#' @param chia.obj A list containing the annotated ChIA-PET data, as returned by \code{\link{annotate.chia}}
#' @param output.dir The name of the directory where to save the graphs.
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 ggsave
#' @importFrom utils write.table
analyze.chromatin.states <- function(chia.obj, output.dir="output") {
    if(!has.chrom.state(chia.obj)) {
        warning("No chromatin states to analyze!")
    } else {
        # Write proportion of chromatin states/annotation types
        state.proportions = table(chia.obj$Regions$Chrom.State)/length(chia.obj$Regions)
        write.table(data.frame(State=names(state.proportions), Proportion=as.vector(state.proportions)),
                    file.path(output.dir, "Chromatin states summary.txt"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
        
        # Plot histogram of degrees for each chromatin state
        degree.per.state.df = data.frame(Chrom.State=chia.obj$Regions$Chrom.State, Degree=chia.obj$Regions$Degree)
        ggplot(degree.per.state.df, aes(x=Degree)) +
            geom_histogram() +
            scale_y_log10() +
            facet_wrap(~Chrom.State)
        ggsave(file.path(output.dir, "log Degree histogram per chromatin state.pdf"))
        
        connectivity.enrichment(chia.obj, "Chrom.State", "chromatin state", 3, 6, output.dir=output.dir)
        contact.heatmap(chia.obj, "Chrom.State", "chromatin states", output.dir=output.dir)
    }
}

#' Analyze the annotation of ChIA-PET data
#'
#' Analyzes ChIA-PET data and produces graphs according to the annotation: \describe{
#' \item{Proportion of genomic location as a function of connectivity category.pdf}{A plot of the proportion of the genomic location of the regions as a fonction of connectivity category.}
#' \item{Contact heatmap for genomic location.pdf}{A contact heatmap of the genomic location of the regions.}}
#'
#' @param chia.obj A list containing the annotated ChIA-PET data, as returned by \code{\link{annotate.chia}}
#' @param output.dir The name of the directory where to save the graphs.
analyze.annotation <- function(chia.obj, output.dir="output") {
    if(!has.gene.annotation(chia.obj)) {
        warning("No gene annotation to analyze!")
    } else {
        connectivity.enrichment(chia.obj, "Simple.annotation", "genomic location", 3, 3, output.dir)
        contact.heatmap(chia.obj, "Simple.annotation", "genomic location", output.dir)
    }
}

#' Analyze the expression of ChIA-PET data
#'
#' Produce a plot to the expression of a node according to its degree.
#'
#' @param chia.obj A list containing the annotated ChIA-PET data, as returned by \code{\link{annotate.chia}}.
#' @param output.dir The name of the directory where to save the graphs.
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 ggsave
analyze.expression <- function(chia.obj, output.dir="output") {
    if(!(has.degree(chia.obj) && has.expression.levels(chia.obj))) {
        warning("No expression levels to analyze!")
    } else {
        gene.reps = chia.obj$Regions[chia.obj$Regions$Gene.Representative==TRUE]
        degree.exp.df <- data.frame(Degree=gene.reps$Degree,
                                    Exp.Mean=gene.reps$Expr.mean)
        ggplot(degree.exp.df, aes(x=log2(Degree), y=log2(Exp.Mean))) +
            geom_point() +
            geom_smooth(method='lm')
        ggsave(file.path(output.dir, "Expression vs Degree at promoter.pdf"))
        
        boxplot.by.connectivity(chia.obj, "Expr.mean", "Boxplot of the expression in fct of Connectivity", output.dir)
    }
}

#' Analyze the gene specificity of ChIA-PET data
#'
#' Analyzes ChIA-PET data and produces graphs according to gene specificity: \describe{
#' \item{Tau vs degree.pdf}{A plot of the Tau of the nodes according to their degree.}
#' \item{Boxplot of degrees by expression category.pdf}{A boxplot the degrees of the nodes in each expression category.}}
#'
#' @param chia.obj A list containing the annotated ChIA-PET data, as returned by \code{\link{annotate.chia}}.
#' @param output.dir The name of the directory where to save the graphs.
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 ggsave
analyze.gene.specificity <- function(chia.obj, output.dir="output") {
    if(!has.gene.specificity(chia.obj)) {
        warning("No gene specificity to analyze!")
    } else {
        # Plot Tau and category vs degree
        tissue.specificity.df = with(chia.obj$Regions[chia.obj$Regions$Gene.Representative],
                                    data.frame(Degree=degree(chia.obj$Graph),
                                               Tau=chia.obj$Regions$Expression.Tau,
                                               Category=chia.obj$Regions$Expression.Category))
        tissue.specificity.df$Category = factor(tissue.specificity.df$Category,
                                                levels=c("Not detected",
                                                        "Mixed low", "Mixed high",
                                                        "Moderately tissue enriched", "Highly tissue enriched", "Group enriched",
                                                        "Expressed in all low", "Expressed in all high"))
        
        ggplot(tissue.specificity.df, aes(x=Tau, y=log2(Degree))) + geom_point()
        ggsave(file.path(output.dir, "Tau vs degree.pdf"))
        
        # Plot category vs degree
        ggplot(tissue.specificity.df, aes(x=Category, y=log2(Degree))) +
            geom_boxplot() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))
        ggsave(file.path(output.dir, "Boxplot of degrees by expression category.pdf"))
    }
}

#' Analyze the TF of ChIA-PET data
#' Produce a plot of the presence of TF at the connection points of ChIA-PET data, according to the
#' connectivity of the nodes.
#'
#' @param chia.obj A list containing the annotated ChIA-PET data, as returned by \code{\link{annotate.chia}}.
#' @param output.dir The name of the directory where to save the graphs.
#'
#' @importFrom GenomicRanges mcols
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 ggsave
#' @importFrom reshape2 melt
analyze.tf <- function(chia.obj, output.dir="output") {
    tf.columns = grepl("^TF\\.", colnames(mcols(chia.obj$Regions)))

    if(sum(tf.columns) > 0) {
        # Extract the overlap matrix from the region annotations.
        overlap.matrix <- as.matrix(mcols(chia.obj$Regions)[,tf.columns])

        # Look at TF presence curves as a function of connectivity
        boundaries.list = list(Singles=c(0, 1), Low=c(1, 5), Intermediate=c(5, 20), High=c(20, 1000))
        results = matrix(0, nrow=ncol(overlap.matrix), ncol=length(boundaries.list),
                         dimnames=list(Rows=colnames(overlap.matrix), Columns=names(boundaries.list)))

        # Loop over boundaries and TFs, calculating percentages of overlap.
        for(i in 1:length(boundaries.list)) {
            for(j in 1:ncol(overlap.matrix)) {
                boundaries = boundaries.list[[i]]
                indices = chia.obj$Regions$Degree > boundaries[1] & chia.obj$Regions$Degree <= boundaries[2]

                results[j,i] <- sum(overlap.matrix[indices,j] > 0) / sum(indices)
            }
        }

        results.df =  melt(results)
        colnames(results.df) = c("TF", "Connectivity", "Proportion")
        results.df$Connectivity = factor(results.df$Connectivity, levels = names(boundaries.list))

        # Reorganize TF by slope
        results.df$TF = factor(results.df$TF, levels=rownames(results)[order(results[,4]-results[,1])])

        ggplot(data=results.df, aes(x=Connectivity, y=Proportion)) +
            geom_line(group=1) +
            facet_wrap(~TF, ncol=10) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))
        ggsave(file.path(output.dir, "TF presence on contact point by connectivity.pdf"), width=14, height=14)
    } else {
        warning("No transcription factor to analyze!")
    }
}

#' Analyze the topology of ChIA-PET data
#'
#' Analyzes ChIA-PET data and produces graphs according to the topology of the nodes in the whole set
#' of regions and in the separated components: \describe{
#' \item{Histogram of number of edges.pdf}{Histogram of the number of edges in the whole set of data.}
#' \item{Scatter plot of degrees of the left node versus the right node.pdf}{Scatter plot showing the relation between the left and right vertices' degree.}
#' \item{Degree vs cluster width.pdf}{A plot of the degree of the nodes in one cluster according to the size of its cluster.}
#' \item{Log2(size of component) vs Log2(Number of components).pdf}{A plot of the number of components according to the size of the components.}
#' \item{Component table.txt}{A file with the information about TSS for each component.}
#' \item{Proportion of TSS in low connectivity nodes.pdf}{A plot of the proportion of TSS in nodes with low connectivity (component with 5 nodes or less).}
#' \item{Proportion of TSS in high connectivity nodes.pdf}{A plot of the proportion of TSS in nodes with high connectivity (component with more than 5 nodes).}}
#'
#' @param chia.obj A list containing the annotated ChIA-PET data, as returned by \code{\link{annotate.chia}}.
#' @param output.dir The name of the directory where to save the graphs.
#'
#' @import ggplot2
#' @importFrom igraph components
#' @importFrom igraph get.data.frame
#' @importFrom utils write.table
analyze.generic.topology <- function(chia.obj, output.dir="output") {
  # Plot an histogram of the number of edges.
  #hist(log2(degree(chia.obj$Graph)), breaks=seq(0, 300, by=5))
  ggplot(as.data.frame(chia.obj$Regions)) + geom_histogram(aes(Degree)) + scale_y_log10() + scale_x_log10()
  ggsave(file.path(output.dir, "Histogram of number of edges.pdf"))

  # Plot a scatter plot showing the relation between the left and right vertices' degree.
  chia.df = get.data.frame(chia.obj$Graph)
  count.per.index = table(c(chia.df$from, chia.df$to))
  left.count = count.per.index[as.character(chia.df$from)]
  right.count = count.per.index[as.character(chia.df$to)]

  ggplot(data.frame(Left=as.vector(left.count), Right=as.vector(right.count)), aes(x=Left, y=Right)) +
    geom_point(alpha=0.01)
  ggsave(file.path(output.dir, "Scatter plot of degrees of the left node versus the right node.pdf"))

  # Degree vs cluster size
  cluster.size.df = data.frame(Degree=degree(chia.obj$Graph), Width=width(chia.obj$Regions))
  ggplot(cluster.size.df, aes(x=Width, y=Degree)) + geom_point() + scale_x_log10() + scale_y_log10()
  ggsave(file.path(output.dir, "Degree vs cluster width.pdf"))

  # Analyze components
  if(has.components(chia.obj)) {
    component.table <- as.data.frame(mcols(chia.obj$Regions)[,c("Component.Id", "Component.size")])
    component.df <- data.frame(Size = unique(component.table)$Component.size, Number = 1)
    component.df <- aggregate(Number~Size, data = component.df, FUN = sum)
    
    ggplot(component.df) + geom_point(mapping=aes(x=log2(Size), y=log2(Number)))
    ggsave(file.path(output.dir, "Log2(size of component) vs Log2(Number of components).pdf"))
    
    annotate.component <- function(x) {
      result.df = data.frame(NumberOfTSS=sum(x$distanceToTSS==0),
                             Size=nrow(x),
                             data.frame(as.list((table(x$Simple.annotation)/nrow(x))), check.names=FALSE))
    
      if(!is.null(x$Chrom.State)) {
        result.df = cbind(result.df, data.frame(as.list((table(x$Chrom.State)/nrow(x))), check.names=FALSE))
      }
    
      return(result.df)
    }
    
    component.table = ddply(as.data.frame(chia.obj$Regions), "Component.Id", annotate.component)
    write.table(component.table, file=file.path(output.dir, "Component table.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
    
    component.table = ddply(component.table, "Size", function(x) { return(cbind(x, NumberOfComponentsOfThisSize=nrow(x)))})
    proportion.per.size = ddply(component.table, ~NumberOfTSS * Size, function(x) { return(nrow(x)/x$NumberOfComponentsOfThisSize[1])})
    proportion.per.size[order(proportion.per.size$Size, proportion.per.size$NumberOfTSS),]
    
    ggplot(subset(proportion.per.size, Size <= 5), aes(x=NumberOfTSS, y=V1)) + geom_bar(stat="identity") + facet_wrap(~Size)
    ggsave(file.path(output.dir, "Proportion of TSS in low connectivity nodes.pdf"))
    
    ggplot(subset(component.table, Size >5), aes(x=NumberOfTSS/Size)) + geom_histogram()
    ggsave(file.path(output.dir, "Proportion of TSS in high connectivity nodes.pdf"))
  }
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
plot.network.heatmap <- function(chia.obj, size.limit, variable.name, label, output.dir) {
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

#' Performs all possible analysis steps on a ChIA-PET object.
#'
#' @param chia.obj A list containing the annotated ChIA-PET data, as returned by \code{\link{annotate.chia}}.
#' @param output.dir The directory where output should be saved.
#' @param verbose Set to TRUE to get progress output to the console.
#' @param label The label to add to the heatmap title (name of the variable).
#'
#' @export
analyze.chia.pet <- function(chia.obj, output.dir=".", verbose=TRUE) {
    # If verbose output is turned off, redirect output to a NULL stream.
    cat.sink = ifelse(verbose, "", textConnection(NULL, w))
    
    # Output the results of the annotation.
    output.annotated.chia(chia.obj, output.dir)
    
    # Perform further in-depth analysis of the networks.
    cat(date(), " : Analyzing network topologies...\n",cat.sink)
    analyze.generic.topology(chia.obj, output.dir)
    
    cat(date(), " : Analyzing genomic annotations...\n",cat.sink)
    analyze.annotation(chia.obj, output.dir)
    
    cat(date(), " : Analyzing chromatin states...\n",cat.sink)
    analyze.chromatin.states(chia.obj, output.dir)
    
    cat(date(), " : Analyzing gene expression...\n",cat.sink)
    analyze.expression(chia.obj, output.dir)
    
    cat(date(), " : Analyzing transcription factor overlaps...\n",cat.sink)
    analyze.tf(chia.obj, output.dir)

    cat(date(), " : Analyzing gene specificity...\n",cat.sink)
    analyze.gene.specificity(chia.obj, output.dir)
    
    # Close dummy verbose stream.
    if(!verbose) {
        close(cat.sink)
    }
}
#############################################################################################################
# Functions That should always run:
#############################################################################################################

#' Return the left part of the ChIA-PET data.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load.chia}}.
#'
#' @return A \linkS4class{GRanges} object with the \code{Regions} from the "\code{chia.obj}" parameter corresponding the the "left side".
#' of the original data.
#' @importFrom igraph as_edgelist
chia.left <- function(chia.obj) {
    return(chia.obj$Regions[as_edgelist(chia.obj$Graph)[,1]])
}

#' Return the right part of the ChIA-PET data.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load.chia}}.
#'
#' @return A \linkS4class{GRanges} object with the \code{Regions} from the "\code{chia.obj}" parameter corresponding the the "right side".
#' of the original data.
#' @importFrom igraph as_edgelist
chia.right <- function(chia.obj) {
    return(chia.obj$Regions[as_edgelist(chia.obj$Graph)[,2]])
}

#' Read and load a ChIA-PET output file.
#'
#' @param input.chia The path of the file containing the ChIA-PET data.
#'
#' @return A list of 4 elements: \describe{
#' \item{$Left}{A \linkS4class{GRanges} object containing the information about the "left side" of the ChIA-PET data.}
#' \item{$Right}{A \linkS4class{GRanges} object containing the information about the "right side" of the ChIA-PET data.}
#' \item{$Regions}{A \linkS4class{GRanges} object containing the reduced information of both sides.}
#' \item{$Graph}{A directed \linkS4class{igraph} object picturing every interaction between the left side and the right
#' side of the ChIA-PET data.}}
#'
#' @importFrom GenomicRanges GRangesList
#' @importFrom GenomicRanges GRanges
#' @importMethodsFrom GenomicRanges findOverlaps
#' @importFrom igraph make_graph
#' @importFrom igraph set_edge_attr
#'
#' @export
load.chia <- function(input.chia) {
    chia.raw = read.table(input.chia, sep="\t")

    # Separate and extend on both sides
    split.raw.chia <- function(chia.raw, columns, flank.size=0) {
       result = chia.raw[,columns]
       colnames(result) = c("chr", "start", "end")
       result$start = pmax(0, result$start - flank.size)
       result$end = result$end + flank.size

       return(GRanges(result))
    }

    chia.left.ranges = split.raw.chia(chia.raw, 1:3)
    chia.right.ranges = split.raw.chia(chia.raw, 4:6)

    # Reduce to a single set of coordinates
    single.set = reduce(unlist(GRangesList(chia.left.ranges, chia.right.ranges)))

    # Build a graph.
    # Map back to original contact points
    chia.left.indices = findOverlaps(chia.left.ranges, single.set, select="first")
    chia.right.indices = findOverlaps(chia.right.ranges, single.set, select="first")

    # Create iGraph object.
    chia.graph = make_graph(c(rbind(chia.left.indices, chia.right.indices)))
    chia.graph = set_edge_attr(chia.graph, "Reads", value=chia.raw[,7])

    return(list(Left=chia.left.ranges, Right=chia.right.ranges, Regions=single.set, Graph=chia.graph))
}

#' Annotate "chia.obj", given as parameter.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load.chia}}.
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
#' @param split Should the data be divided into communities?
#' @param oneByOne Sould the netwoks be treated one by one or as a whole?
#' @param method What method sould be used to split data (ignored if split = \code{FALSE}).
#' @param verbose If TRUE, console output is suppressed.
#'
#' @return The annotated "\code{chia.obj}".
#'
#' @importFrom igraph degree
#'
#' @export
annotate.chia <- function(chia.obj, input.chrom.state, tf.regions, histone.regions, pol.regions, expression.levels,
                          genome.build = c("hg19", "mm9", "mm10", "hg38"), biosample = "GM12878",
                          tssRegion = c(-3000, 3000), output.dir, split = TRUE, oneByOne = FALSE,
                          method = igraph::cluster_fast_greedy, verbose=TRUE) {
  # Make sure the genome build is one of the supported ones.
  genome.build <- match.arg(genome.build)

  # Create the output directory.
  dir.create(output.dir, recursive = TRUE, showWarnings=FALSE)

  # If verbose output is turned off, redirect output to a NULL stream.
  cat.sink = ifelse(verbose, "", textConnection(NULL, w))

  # Add an ID to every region.
  chia.obj$Regions$ID = 1:length(chia.obj$Regions)

  # Add degree count directly to chia.obj$Regions
  chia.obj$Regions$Degree = degree(chia.obj$Graph)

  # Associate various
  cat(date(), " : Associating genomic regions...\n",cat.sink)
  chia.obj$Regions <- associate.genomic.region(chia.obj$Regions, genome.build, output.dir, tssRegion = tssRegion)

  if(!is.null(input.chrom.state)) {
    cat(date(), " : Associating chromatin states...\n",cat.sink)
    chia.obj$Regions = associate.chrom.state(chia.obj$Regions, input.chrom.state)
  }

  if(!is.null(tf.regions)) {
    cat(date(), " : Associating transcription factors...\n",cat.sink)
    chia.obj$Regions = associate.tf(chia.obj$Regions, tf.regions)
  }

  if(!is.null(histone.regions)) {
    cat(date(), " : Associating histone marks...\n",cat.sink)
    chia.obj$Regions = associate.histone.marks(chia.obj$Regions, histone.regions)
  }

  if(!is.null(pol.regions)) {
    cat(date(), " : Associating polymerase II regions...\n",cat.sink)
    chia.obj$Regions = associate.histone.marks(chia.obj$Regions, pol.regions)
  }

  cat(date(), " : Associating genes...\n",cat.sink)
  chia.obj$Regions = associate.gene(chia.obj$Regions, expression.levels)

  if(genome.build=="hg19" || genome.build=="hg38") {
    cat(date(), " : Associating tissue specificity...\n",cat.sink)
    chia.obj$Regions = associate.tissue.specificity.human(chia.obj$Regions)
    cat(date(), " : Associating fitness score...\n",cat.sink)
    chia.obj$Regions = associate.fitness.genes(chia.obj$Regions)
  }


  cat(date(), " : Associating components...\n",cat.sink)
  chia.obj = associate.components(chia.obj, split = split, oneByOne = oneByOne, method = method)

  cat(date(), " : Associating centrality scores...\n",cat.sink)
  chia.obj = associate.centralities(chia.obj)
  chia.obj$Regions = associate.is.in.factory(chia.obj$Regions)

  chia.obj$Regions = associate.is.gene.active(chia.obj$Regions)

  # If verbose output was turned off, close the dummy stream.
  if(!verbose) {
    close(cat.sink)
  }
  return(chia.obj)
}

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

#' Save annotated ChIA-PET data.
#'
#' Writes files: \itemize{
#'   \item The first one separates the "left" and "right" sides of the ChIA-PET data with annotations.
#'   \item The second set of files is made of components files. They group all interactions in a single component, with
#'        an extra column containing the number of reads supporting the data. Their format is supported by Cytoscape.}
#'
#' @param chia.obj A list containing the annotated ChIA-PET data, as returned by \code{\link{annotate.chia}}.
#' @param output.dir The name of the directory where to save the files.
#'
#' @importFrom utils write.table
#' @importFrom igraph edge_attr
#'
#' @export
output.annotated.chia <- function(chia.obj, output.dir="output") {
  # Create output directory if it does not exist.
  dir.create(output.dir, recursive = TRUE)

  # Export region annotation, putting the ID in the first column.
  chia.data <- as.data.frame(chia.obj$Regions)
  chia.data <- cbind(chia.data$ID, chia.data[,-which(colnames(chia.data) == "ID")])
  write.table(chia.data, file = file.path(output.dir, "Annotated CHIA-PET regions.txt"), row.names = FALSE, sep = "\t")

  # Write out annotated interactions by concatening left and right annotations.
  left.df = as.data.frame(chia.left(chia.obj))
  colnames(left.df) <- paste("Left", colnames(left.df), sep=".")
  right.df = as.data.frame(chia.right(chia.obj))
  colnames(right.df) <- paste("Right", colnames(right.df), sep=".")

  write.table(data.frame(left.df, right.df), file=file.path(output.dir, "Interactions.txt"), sep="\t", col.names=TRUE)


  # Output Cytoscape components

  # Create interactions table
  ids <- data.frame(left.df$Left.ID, right.df$Right.ID, edge_attr(chia.obj$Graph)$Reads, left.df$Left.Component.Id)
  colnames(ids) <- c("Source", "Target", "Reads", "Component")


  # Export networks in csv files
  dir.create(file.path(output.dir, "Size between 3 and 5 nodes (incl)"), recursive = TRUE)
  dir.create(file.path(output.dir, "Size between 6 and 20 nodes (incl)"), recursive = TRUE)
  dir.create(file.path(output.dir, "Size between 21 and 50 nodes (incl)"), recursive = TRUE)
  dir.create(file.path(output.dir, "Size between 51 and 100 nodes (incl)"), recursive = TRUE)
  dir.create(file.path(output.dir,"Size over 100 nodes"), recursive = TRUE)

  for (i in unique(ids$Component)){
    network <- ids[ids$Component == i,]
    size <- length(unique(c(network$Source, network$Target)))
    if (size > 2){
      if (size < 6){
        dir <- "Size between 3 and 5 nodes (incl)"
      } else if (size < 21){
        dir <- "Size between 6 and 20 nodes (incl)"
      } else if (size < 51){
        dir <- "Size between 21 and 50 nodes (incl)"
      } else if (size < 101){
        dir <- "Size between 51 and 100 nodes (incl)"
      } else {
        dir <- "Size over 100 nodes"
      }
      write.table(network[1:3],
                  file = file.path(output.dir, dir, paste0("Create network for components ", i, "(", size, " nodes)", ".csv")),
                  sep = ",", row.names = FALSE)
    }
  }
  chia.data <- as.data.frame(chia.obj$Regions)
  chia.data <- cbind(chia.data$ID, chia.data[,-which(colnames(chia.data) == "ID")])
  colnames(chia.data)[1] <- "ID"
  write.table(chia.data, file = file.path(output.dir, "Annotated CHIA-PET regions.txt"), row.names = FALSE, sep = "\t")
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
  dir.create(file.path(output.dir, biosample), recursive = TRUE)
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
  data.for.boxplot <- data.for.boxplot[!is.na(data.for.boxplot$Expr.mean),]
  ggplot(data.for.boxplot) +
    geom_boxplot(aes(CutDegree, get(variable.name))) +
    ylab(variable.name) +
    xlab("Connectivity") +
    scale_y_log10() +
    ggtitle(label)
  ggsave(file.path(output.dir, paste0(label, ".pdf")))
}

#' Separates a large networks into communities and saves the sub-netwoks in cytoscape-friendly format
#'
#' @param network.input The csv file containing the network to divide.
#' @param method The algorithm to use to divide the graph.
#' @return A list with two elements: the annotated chia.obj, with changed components ids ans sizes and the new networks.
#' @importFrom igraph make_graph
separate.into.communities <- function(network.input, method = igraph::cluster_fast_greedy){

  # Separate the network into smaller ones
  network.input <- unique(network.input)
  whole.graph <- make_graph(c(rbind(network.input$Source, network.input$Target)), directed = FALSE)

  communities <- method(whole.graph, weights = NULL)
  communities.df <- data.frame(cbind(unique(c(network.input$Source, network.input$Target)),
                                     communities$membership[unique(c(network.input$Source, network.input$Target))]))
  colnames(communities.df) <- c("Node.Id", "Group")
  simplify.group.id <- data.frame(Old = unique(communities.df$Group), New = 1:length(unique(communities.df$Group)))
  communities.df$Group <- simplify.group.id$New[match(communities.df$Group, simplify.group.id$Old)]
  communities.df$Size <- 1
  communities.df$Size <- aggregate(Size~Group, data = communities.df, FUN = sum)$Size[communities.df$Group]


  return(communities.df)
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

#' Analyze the annotation of ChIA-PET data
#'
#' Analyzes ChIA-PET data and produces graphs according to the annotation: \describe{
#' \item{Proportion of genomic location as a function of connectivity category.pdf}{A plot of the proportion of the genomic location of the regions as a fonction of connectivity category.}
#' \item{Contact heatmap for genomic location.pdf}{A contact heatmap of the genomic location of the regions.}}
#'
#' @param chia.obj A list containing the annotated ChIA-PET data, as returned by \code{\link{annotate.chia}}
#' @param output.dir The name of the directory where to save the graphs.
analyze.annotation <- function(chia.obj, output.dir="output") {
    connectivity.enrichment(chia.obj, "Simple.annotation", "genomic location", 3, 3, output.dir)
    contact.heatmap(chia.obj, "Simple.annotation", "genomic location", output.dir)
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
    gene.reps = chia.obj$Regions[chia.obj$Regions$Gene.Representative==TRUE]
    degree.exp.df <- data.frame(Degree=gene.reps$Degree,
                                Exp.Mean=gene.reps$Expr.mean)
    ggplot(degree.exp.df, aes(x=log2(Degree), y=log2(Exp.Mean))) +
        geom_point() +
        geom_smooth(method='lm')
    ggsave(file.path(output.dir, "Expression vs Degree at promoter.pdf"))

    boxplot.by.connectivity(chia.obj, "Expr.mean", "Boxplot of the expression in fct of Connectivity", output.dir)
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
    # Plot Tau and category vs degree
    tissue.specificity.df = with(chia.obj$Regions[chia.obj$Regions$Gene.Representative],
                                data.frame(Degree=chia.obj$Regions$Degree,
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

#' Analyze ChIA-PET data and produce graphs.
#'
#' @param input.chia The file containing processed ChIA-PET data.
#' @param input.chrom.state The name of the file containing the information about chromatin states.
#' @param biosample The biosample identifier from ENCODE. Valid examples are GM12878, K562 or MCF-7.
#' @param genome.build The name of the chosen annotation ("hg38", "hg19").
#' @param tf.regions A \linkS4class{GRangesList} object containing TF regions.
#' @param histone.regions A \linkS4class{GRangesList} object containing histone regions.
#' @param pol.regions A \linkS4class{GRangesList} object containing the pol2 regions.
#' @param expression.levels A data frame containing the levels of expression of genes, according to their EMSEMBL id.
#' @param output.dir The name of the directory where output should be saved.
#' @param ... Additional parameters to be passed to \code{\link{annotate.chia}}
#' @return The annotated chia.obj.
#' @importFrom Biobase cache
#' @export
analyze.chia.pet <- function(input.chia, input.chrom.state = NULL, biosample = NULL, genome.build = NULL, tf.regions = NULL,
                             histone.regions = NULL, pol.regions = NULL, expression.data = NULL, output.dir="output", verbose=TRUE,
                             ...) {
    # Create output directory.
    dir.create(file.path(output.dir), recursive=TRUE, showWarnings=FALSE)

    # Load interaction data.
    chia.obj = load.chia(input.chia)

    # If biosample is provided, download missing annotations from ENCODE.
    if(!is.null(biosample) && !is.null(genome.build)) {
        # Download transcription factors
        if (is.null(tf.regions)){
            cache(tf.regions <- download.encode.chip(biosample, genome.build)$Regions, output.dir)
        }

        # Download histone marks
        if (is.null(histone.regions)){
            cache(histone.regions <- download.encode.chip(biosample, genome.build, download.filter=histone.download.filter.chip,
                                                  download.dir=file.path("input/ENCODE", biosample, "chip-seq", "histone"))$Regions,
                output.dir)
        }

        # Download PolII regions.
        if (is.null(pol.regions)){
            cache(pol.regions <- download.encode.chip(biosample, genome.build, download.filter = pol2.download.filter.chip,
                                              download.dir = file.path("input/ENCODE", biosample, "chip-seq", "pol2"))$Regions,
                output.dir)
        }

        # Download expression data
        if (is.null(expression.data)){
            cache(expression.data <- download.encode.rna(biosample, genome.build)$Expression, output.dir)
            expression.data$ENSEMBL = gsub("\\.\\d+$", "", expression.data$gene_id)
            expression.data$FPKM = log2(expression.data$Mean.FPKM + 1)
        }

        # Download chromatin states
        if (is.null(input.chrom.state) && genome.build=="hg19") {
            input.chrom.state <- import.chrom.states(biosample, file.path("input/chrom_states", biosample))
        }
    }

    # Annotate the ChIA object.
    cache(chia.obj <- annotate.chia(chia.obj, input.chrom.state = input.chrom.state, tf.regions = tf.regions,
                              histone.regions=histone.regions, pol.regions = pol.regions,
                              expression.levels=expression.data, genome.build = genome.build, biosample=biosample,
                              output.dir = output.dir, verbose=verbose, ...), output.dir)


    # If verbose output is turned off, redirect output to a NULL stream.
    cat.sink = ifelse(verbose, "", textConnection(NULL, w))

    # Output the results of the annotation.
    output.annotated.chia(chia.obj, output.dir)

    # Perform further in-depth analysis of the networks.
    cat(date(), " : Analyzing network topologies...\n",cat.sink)
    analyze.generic.topology(chia.obj, output.dir)

    cat(date(), " : Analyzing genomic annotations...\n",cat.sink)
	analyze.annotation(chia.obj, output.dir)

	if(!is.null(input.chrom.state)) {
        cat(date(), " : Analyzing chromatin states...\n",cat.sink)
        analyze.chromatin.states(chia.obj, output.dir)
    }

	if(!is.null(expression.data)) {
        cat(date(), " : Analyzing gene expression...\n",cat.sink)
        analyze.expression(chia.obj, output.dir)
    }

    if(!is.null(tf.regions)) {
        cat(date(), " : Analyzing transcription factor overlaps...\n",cat.sink)
        analyze.tf(chia.obj, output.dir)
    }

	if(genome.build %in% c("hg19", "hg38")) {
        cat(date(), " : Analyzing gene specificity...\n",cat.sink)
        analyze.gene.specificity(chia.obj, output.dir)
    }

    # Close dummy verbose stream.
    if(!verbose) {
        close(cat.sink)
    }

	return(chia.obj)
}

#' Associates boolean to regions in focntion of their centrality
#'
#' @param chia.obj chia.obj ChIA-PET data, as returned by \code{\link{annotate.chia}}.
#'
#' @return The annotated chia.obj.
#'
#' @importFrom igraph make_graph
#' @importFrom igraph degree
#' @importFrom igraph estimate_closeness
#' @importFrom igraph betweenness
#' @importFrom igraph eigen_centrality
associate.centralities <- function(chia.obj){
  left.df = as.data.frame(chia.left(chia.obj))
  colnames(left.df) <- paste("Left", colnames(left.df), sep=".")
  right.df = as.data.frame(chia.right(chia.obj))
  colnames(right.df) <- paste("Right", colnames(right.df), sep=".")

  ids <- data.frame(left.df$Left.ID, right.df$Right.ID, left.df$Left.Component.Id, right.df$Right.Component.Id)
  ids <- ids[ids[,4] == ids[,3],1:3]
  colnames(ids) <- c("Source", "Target", "Component")

  # Creation of the new columns to fill
  chia.obj$Regions$Centrality.score <- 0
  chia.obj$Regions$Is.central <- 0

  for (id in unique(chia.obj$Regions$Component.Id)){
    network <- ids[ids$Component == id,]
    graph <- make_graph(c(rbind(network$Source, network$Target)))

    data <- data.frame(Nodes = unique(c(network$Source, network$Target)))
    data$Degree <- degree(graph)[data$Nodes]
    data$Betweenness <- betweenness(graph, directed = FALSE)[data$Nodes]
    data$Eigen <- eigen_centrality(graph, directed = FALSE)$vector[data$Nodes]

    data$Degree <- data$Degree / max(data$Degree)
    data$Betweenness <- data$Betweenness / max(data$Betweenness)
    data$Eigen <- data$Eigen / max(data$Eigen)

    data$Centrality.score <- with(data, (Degree + Betweenness + Eigen)/3)
    data$Centrality.score[is.nan(data$Centrality.score)] <- 0
    #data$Is.central <- data$Centrality.score > quantile(data$Centrality.score, probs = 0.95)
    chia.obj$Regions$Centrality.score[match(data$Nodes, chia.obj$Regions$ID)] <- data$Centrality.score

    chia.obj$Regions$Is.central[match(data$Nodes, chia.obj$Regions$ID)] <- data$Centrality.score > quantile(data$Centrality.score, probs = 0.95)
  }

  return(chia.obj)
}

#' Associates boolean to regions in fonction of their presence in factories
#'
#' Is.In.Factory is \code{TRUE} if the region is in a network with 3 genes or more.
#'
#' @param regions A \linkS4class{GRanges} object to annotate.
#' @return The annotated regions.
associate.is.in.factory <- function(regions){
  factories <- as.data.frame(regions)
  factories <- aggregate(Gene.Representative~Component.Id, data = regions, FUN = sum)
  factories <- factories$Component.Id[factories$Gene.Representative > 2]
  regions$Is.In.Factory <- (regions$Component.Id %in% factories) & (regions$Component.size > 4)
  return(regions)

}

#' Associates components ids and sizes to chia data, as returned by \code{\link{load.chia}}.
#'
#' @param chia.obj ChIA-PET data, as returned by \code{\link{annotate.chia}}.
#' @param split Should the data be divided into communities?
#' @param oneByOne Sould the netwoks be treated one by one or as a whole?
#' @param method What method sould be used to split data (ignored if split = \code{FALSE}).
#' @return The annotated chia.obj.
#' @importFrom igraph components
#' @importFrom igraph as.undirected
#' @importFrom igraph delete_edges
#' @importFrom igraph E
associate.components <- function(chia.obj, split = TRUE, oneByOne = FALSE, method = igraph::cluster_fast_greedy){

  left.df = as.data.frame(chia.left(chia.obj))
  right.df = as.data.frame(chia.right(chia.obj))

  ids <- data.frame(left.df$ID, right.df$ID)
  colnames(ids) <- c("Source", "Target")

  # Generate components
  components.out <- components(chia.obj$Graph)
  reorder.components <- components.out$membership[ids$Source]
  ids <- cbind(ids, Component = reorder.components)

  # Add data about components
  chia.obj$Regions$Component.Id <- components.out$membership
  chia.obj$Regions$Component.size <- components.out$csize[components.out$membership]


  if (split){
    if (oneByOne){
      for (i in 1:components.out$no){
        communities.df <- separate.into.communities(ids[ids$Component == i, 1:2], method = method)
        # If the network is divided...
        if (length(unique(communities.df$Group)) > 1){
          # ...Change the annotation
          chia.obj$Regions$Component.size[chia.obj$Regions$Component.Id == i] <- communities.df[order(communities.df$Node.Id),"Size"]
          chia.obj$Regions$Component.Id[chia.obj$Regions$Component.Id == i] <-
            paste0(network, ".", communities.df[order(communities.df$Node.Id),"Group"])
        }
      }
    } else {
      communities.df <- separate.into.communities(ids[, 1:2], method = method)

      chia.obj$Regions$Component.size <- communities.df[order(communities.df$Node.Id),"Size"]
      chia.obj$Regions$Component.Id <- communities.df[order(communities.df$Node.Id),"Group"]
    }
    # Update Left, Right and Graph
    # Reload Left and Right info, with new component ids ans sizes
    left.df = as.data.frame(chia.left(chia.obj))
    right.df = as.data.frame(chia.right(chia.obj))
    # Identify deleted edges
    to.delete <- which(left.df$Component.Id != right.df$Component.Id)

    chia.obj$Left <- chia.obj$Left[-to.delete]
    chia.obj$Right <- chia.obj$Right[-to.delete]

    chia.obj$Graph <- delete_edges(chia.obj$Graph, E(chia.obj$Graph)[to.delete])
  } else {
    left.df = as.data.frame(chia.left(chia.obj))
  }

  # Add a column with the number of edges in each community
  left.df$Edges <- 1
  number.edges.df <- aggregate(Edges~Component.Id, data = left.df, FUN = sum)
  chia.obj$Regions$Component.edges <- number.edges.df$Edges[match(chia.obj$Regions$Component.Id, number.edges.df$Component.Id)]

  return (chia.obj)
}

#' Associate genes to a \linkS4class{GRanges} object from ChIA-PET data.
#'
#' @param regions A \linkS4class{GRanges} object to annotate.
#' @param expression.data A data frame containing the levels of expression of genes, according to their EMSEMBL id.
#'
#' @return The \linkS4class{GRanges} object with associated genes.
#' @importFrom plyr ddply mutate
associate.gene <- function(regions, expression.data=NULL) {
  # Associate a gene to a contact only if it's in the promoter.
  promoter.subset = regions$Simple.annotation == "Promoter"

  # Subset the promoter contacts to only keep the highest degrees
  degree.info = ddply(as.data.frame(regions[promoter.subset]), "ENSEMBL", plyr::mutate, max.degree=max(Degree))
  degree.subset = subset(degree.info, Degree == max.degree)

  # Further subset to keep the one closest to the TSS
  distance.info = ddply(degree.subset, "ENSEMBL", plyr::mutate, min.distance=min(abs(distanceToTSS)))
  distance.subset = subset(distance.info, distanceToTSS == min.distance)

  # If there are still more than one row, just pick the first one.
  gene.subset = ddply(distance.subset, "ENSEMBL", function(x) { return(x[1,]) })

  # Add the gene marker to the annotations.
  regions$Gene.Representative = regions$ID %in% gene.subset$ID

  if(!is.null(expression.data)) {
    index.match = match(regions$ENSEMBL, expression.data$ENSEMBL)
    regions$Expr.mean = expression.data$Mean.FPKM[index.match]
  }

  return(regions)
}

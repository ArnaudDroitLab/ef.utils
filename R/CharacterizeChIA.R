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
#'   \item{$Left} {A \linkS4class{GRanges} object containing the information about the "left side" of the ChIA-PET data.}
#'   \item{$Right} {A \linkS4class{GRanges} object containing the information about the "right side" of the ChIA-PET data.}
#'   \item{$Regions} {A \linkS4class{GRanges} object containing the reduced information of both sides.}
#'   \item{$Graph} {A directed \linkS4class{igraph} object picturing every interaction between the left side and the right
#' side of the ChIA-PET data.}}
#'
#' @importFrom GenomicRanges GRangesList
#' @importFrom GenomicRanges GRanges
#' @importMethodsFrom GenomicRanges findOverlaps
#' @importFrom igraph make_graph
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
    #chia.left.merged = single.set[chia.left.indices]
    chia.right.indices = findOverlaps(chia.right.ranges, single.set, select="first")
    #chia.right.merged = single.set[chia.right.indices]

    # Create iGraph object.
    chia.graph = make_graph(c(rbind(chia.left.indices, chia.right.indices)))

    return(list(Left=chia.left.ranges, Right=chia.right.ranges, Regions=single.set, Graph=chia.graph))
}


#' Annotate "chia.obj", given as parameter.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load.chia}}.
#' @param input.chrom.state The name of the file containing the information about chromatin states.
#' @param tf.regions A data frame containing the TF data.
#' @param expression.levels A data frame containing the levels of expression of genes. according to their EMSEMBL id.
#' @param genome.build The name of the chosen annotation ("hg38", "mm9", "mm10", "hg19").
#' @param biosample The biosample identifier from ENCODE. Valid examples are
#'   GM12878, K562.
#' @param histone Should the overlap percentage of histone marks be added?
#' @param output.dir The name of the directory where to write the selected annotations.
#'
#' @return The annotated "\code{chia.obj}".
#'
#' @importFrom igraph degree
#'
#' @export
annotate.chia <- function(chia.obj, input.chrom.state, tf.regions, histone.regions, expression.levels, genome.build = c("hg19", "mm9", "mm10", "hg38"),
                          biosample = "GM12878", histone = FALSE, output.dir) {
    single.set = chia.obj$Regions
    genome.build <- match.arg(genome.build)

    # Add an ID to every region.
    chia.obj$Regions$ID = 1:length(chia.obj$Regions)

    # Add degree count to chia.obj$Regions
    chia.obj$Regions$Degree = degree(chia.obj$Graph)

    chia.obj <- associate.genomic.region(chia.obj, genome.build, output.dir)

    if(!is.null(input.chrom.state)) {
        chia.obj = associate.chrom.state(chia.obj, input.chrom.state)
    }

    if(!is.null(tf.regions)) {
        chia.obj = associate.tf(chia.obj, tf.regions)
    }

    if(!is.null(histone.regions)) {
        chia.obj = associate.histone.marks(chia.obj, histone.regions)
    }

    chia.obj = associate.gene(chia.obj, expression.levels)

    if(genome.build=="hg19" || genome.build=="hg38") {
        chia.obj = associate.tissue.specificity.human(chia.obj)
        chia.obj = associate.fitness.genes(chia.obj)
    }

    return(chia.obj)
}

#' Associate histone marks with the \code{Regions} of "chia.obj".
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by load.chia.
#' @param genome.build The name of the chosen annotation ("hg38", "mm9", "mm10", "hg19").
#' @param biosample The biosample identifier from ENCODE. Valid examples are
#'   GM12878, K562 or MCF-7.
#'
#' @return "\code{chia.obj}" with associated histone overlap percentage.
#' @importMethodsFrom GenomicRanges findOverlaps range
#' @importMethodsFrom BSgenome width
#' @importFrom stats aggregate
associate.histone.marks <- function(chia.obj, histone.regions){
  # find overlaps
  overlap.index <- findOverlaps(chia.obj$Regions, unlist(histone.regions))

  # to get the overlap length
  overlap.length <- width(ranges(overlap.index, ranges(chia.obj$Regions), ranges(unlist(histone.regions))))

  # correspondances between idexes from chia.obj$Regions and the lengths
  overlap.length.df <- data.frame(index=overlap.index@from, length=overlap.length)
  # addition of the rows with repeated indices
  overlap.length.df <- aggregate(overlap.length.df$length, list(overlap.length.df$index), sum)
  colnames(overlap.length.df) <- c("index", "length")

  # create the new vector
  histone.overlap.percentage <- vector(length = length(chia.obj$Regions@seqnames))
  histone.overlap.percentage[overlap.length.df$index] <- overlap.length.df$length
  chia.obj$Regions$Histone.overlap.percentage <- histone.overlap.percentage / width(chia.obj$Regions)

  return(chia.obj)
}


#' Associate genomic regions with the \code{Regions} of "chia.obj".
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by load.chia.
#' @param genome.build The name of the chosen annotation ("hg38", "mm9", "mm10", "hg19").
#' @param output.dir The name of the directory where to write the selected annotations.
#'
#' @return "\code{chia.obj}" with associated genomic regions.
#' @importFrom GenomicRanges mcols
associate.genomic.region <- function(chia.obj, genome.build, output.dir) {
    # Annotate with proximity to gene regions
    annotations.list = select.annotations(genome.build)
    annotations = annotate.region(chia.obj$Regions, annotations.list, file.path(output.dir, "CHIA-PET annotation.txt"))
    annotations.df = as.data.frame(annotations)
    mcols(chia.obj$Regions) = annotations.df[, 6:ncol(annotations.df)]

    # Write porportions of annotated region types.
    write.table(annotations@annoStat, file.path(output.dir, "Annotation summary.txt"), sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)

    # Simplify genomic region annotation
    simplified.annotation = mcols(chia.obj$Regions)$annotation
    simplified.annotation[grepl("Promoter", simplified.annotation)] <- "Promoter"
    simplified.annotation[grepl("Exon", simplified.annotation)] <- "Exon"
    simplified.annotation[grepl("Intron", simplified.annotation)] <- "Intron"
    simplified.annotation[grepl("Downstream", simplified.annotation)] <- "Downstream"
    chia.obj$Regions$Simple.annotation = factor(simplified.annotation, levels=c("Distal Intergenic", "Promoter", "Intron", "Exon", "Downstream", "3' UTR", "5' UTR"))

    return(chia.obj)
}


#' Associate genes with the \code{Regions} of "chia.obj".
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load.chia}}.
#' @param expression.data A data frame containing the levels of expression of genes, according to their EMSEMBL id.
#'
#' @return "\code{chia.obj}" with associated genes.
#' @importFrom plyr ddply mutate
associate.gene <- function(chia.obj, expression.data=NULL) {
    # Associate a gene to a contact only if it's in the promoter.
    promoter.subset = chia.obj$Regions$Simple.annotation == "Promoter"

    # Subset the promoter contacts to only keep the highest degrees
    degree.info = ddply(as.data.frame(chia.obj$Regions[promoter.subset]), "ENSEMBL", plyr::mutate, max.degree=max(Degree))
    degree.subset = subset(degree.info, Degree == max.degree)

    # Further subset to keep the one closest to the TSS
    distance.info = ddply(degree.subset, "ENSEMBL", plyr::mutate, min.distance=min(abs(distanceToTSS)))
    distance.subset = subset(distance.info, distanceToTSS == min.distance)

    # If there are still more than one row, just pick the first one.
    gene.subset = ddply(distance.subset, "ENSEMBL", function(x) { return(x[1,]) })

    # Add the gene marker to the annotations.
    chia.obj$Regions$Gene.Representative = chia.obj$Regions$ID %in% gene.subset$ID

    if(!is.null(expression.data)) {
        index.match = match(chia.obj$Regions$ENSEMBL, expression.data$ENSEMBL)
        chia.obj$Regions$Expr.mean = expression.data$Mean.FPKM[index.match]
    }

    return(chia.obj)
}


#' Associate chromatin sates with the \code{Regions} of "chia.obj".
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load.chia}}.
#' @param input.chrom.state The path to a bed file representing a genome
#'   segmentation by chromatin state.
#' @return "\code{chia.obj}" with associated chromatin states.
#' @importFrom rtracklayer import
#' @importMethodsFrom GenomicRanges findOverlaps
associate.chrom.state <- function(chia.obj, input.chrom.state) {
    # Annotate with chromatin states
    # Load and rename chromatin states (for easier lexical ordering)
    chrom.states = rtracklayer::import(input.chrom.state)
    chrom.states$name = gsub("^(.)_", "0\\1_", chrom.states$name)

    # Sort according to name, which will cause the first match to also be the
    # most relevant states (01_TSS first, 18_Quies last)
    chrom.states = chrom.states[order(chrom.states$name)]
    unique.states = sort(unique(chrom.states$name))

    # Add state annotation to chia.obj$Regions
    state.overlap.indices = findOverlaps(chia.obj$Regions, chrom.states, select="first")
    chia.obj$Regions$Chrom.State = factor(chrom.states$name[state.overlap.indices], levels=unique.states)

    return(chia.obj)
}


#' Associate TF overlaps with the \code{Regions} of "chia.obj".
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load.chia}}.
#' @param tf.regions A data frame containing the TF data.
#'
#' @return "\code{chia.obj}" with associated TF overlaps.
#' @importMethodsFrom GenomicRanges countOverlaps
#' @importFrom GenomicRanges mcols
associate.tf <- function(chia.obj, tf.regions) {
    # Calculate TF overlap with contact regions
    overlap.matrix = matrix(0, nrow=length(chia.obj$Regions), ncol=length(tf.regions))
    colnames(overlap.matrix) <- paste0("TF.overlap.", names(tf.regions))
    for(i in 1:length(tf.regions)) {
        overlap.matrix[,i] <- countOverlaps(chia.obj$Regions, tf.regions[[i]])
    }

    mcols(chia.obj$Regions) = cbind(mcols(chia.obj$Regions), as.data.frame(overlap.matrix))

    return(chia.obj)
}


#' Associate tissue specificity of genes with the \code{Regions} of "chia.obj".
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load.chia}}.
#'
#' @return "\code{chia.obj}" with associated tissue scpecificity.
associate.tissue.specificity.human <- function(chia.obj) {

    # Annotate chia.obj$Regions with Tau, Expression category.
    calculate.tau <- function(x) {
        return(sum(1 - (x/max(x))) / (length(x) - 1))
    }
    tissue.expression$Tau = apply(tissue.expression[,c(-1, -ncol(tissue.expression))], 1, calculate.tau)

    tissue.match = match(chia.obj$Regions$ENSEMBL, tissue.expression$Ensembl.gene.id)
    chia.obj$Regions$Expression.Category = factor(tissue.expression$Category[tissue.match],
                                            levels=c("Not detected",
                                                    "Mixed low", "Mixed high",
                                                    "Moderately tissue enriched", "Highly tissue enriched", "Group enriched",
                                                    "Expressed in all low", "Expressed in all high"))

    chia.obj$Regions$Expression.Tau = tissue.expression$Tau[tissue.match]

    return(chia.obj)
}

#' Identify essential genes of the \code{Regions} of "chia.obj".
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load.chia}}.
#'
#' @return "\code{chia.obj}" with identified essential genes.
associate.fitness.genes <- function(chia.obj){

  # Add the "essential ratio" to the data
  fitness.match <- match(chia.obj$Regions$SYMBOL, essential.genes$Gene)
  chia.obj$Regions$Fitness <- essential.genes$numTKOHits[fitness.match]
  chia.obj$Regions$Fitness <- ifelse(is.na(chia.obj$Regions$Fitness), 0, (chia.obj$Regions$Fitness / 6))

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
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
    ggsave(file.path(output.dir, paste0("Contact heatmap for ", label, ".pdf")))

    return(results.matrix)
}

#' Save annotated ChIA-PET data.
#'
#' Writes files: \itemize{
#'   \item The first one separates the "left" and "right" sides of the ChIA-PET data with annotations.
#'   \item The second set of files is made of components files. They group all interactions in a single component, with
#'        an extra column containing the number of reads supporting the data. Their format is supported by Cytoscape.
#'   \item The final file contains the annotated Regions of the ChIA-PET data, with two extra columns: \describe{
#'        \item{$Component.Id} {Id (number) of the hub in which the node is found.}
#'        \item{$Componend.Size} {Number of nodes of the component in which the node is found.}}}
#'
#' @param chia.obj A list containing the annotated ChIA-PET data, as returned by \code{\link{annotate.chia}}.
#' @param output.dir The name of the directory where to save the files.
#'
#' @importFrom utils write.table
#' @importFrom igraph components
#'
#' @export
output.annotated.chia <- function(chia.obj, output.dir="output") {

    # Write out annotated interactionsleft.df = as.data.frame(chia.left.merged)
    left.df = as.data.frame(chia.left(chia.obj))
    colnames(left.df) <- paste("Left", colnames(left.df), sep=".")
    right.df = as.data.frame(chia.right(chia.obj))
    colnames(right.df) <- paste("Right", colnames(right.df), sep=".")

    write.table(data.frame(left.df, right.df), file=file.path(output.dir, "Interactions.txt"), sep="\t", col.names=TRUE)


    # Output Cytoscape components

    # Generate components
    component.out <- components(chia.obj$Graph)

    # Create interactions table
    ids <- data.frame(left.df$ID, right.df$ID, chia.raw[,7])
    colnames(ids) <- c("Source", "Target", "Reads")
    dir.create(output, recursive = TRUE)
    write.table(ids, file = paste0(file.path, "all.csv"), sep=",", row.names = FALSE)

    # Export networks in csv files
    reorder.components <- components.out$membership[ids$Source]
    ids.components <- cbind(ids, reorder.components)
    colnames(ids.components) <- c(colnames(ids), "Component")
    dir.create(file.path(output, "Size of 5 nodes and less"), recursive = TRUE)
    dir.create(file.path(output, "Size between 6 and 20 nodes (incl)"), recursive = TRUE)
    dir.create(file.path(output, "Size between 21 and 50 nodes (incl)"), recursive = TRUE)
    dir.create(file.path(output, "Size between 51 and 100 nodes (incl)"), recursive = TRUE)
    dir.create(file.path(output,"Size over 100 nodes"), recursive = TRUE)
    for (i in 1:components.out$no){
      network <- ids[ids.components$Component == i,]
      if (components.out$csize[i] < 6){
        dir <- "Size of 5 nodes and less"
      } else if (components.out$csize[i] < 21){
        dir <- "Size between 6 and 20 nodes (incl)"
      } else if (components.out$csize[i] < 51){
        dir <- "Size between 21 and 50 nodes (incl)"
      } else if (components.out$csize[i] < 101){
        dir <- "Size between 51 and 100 nodes (incl)"
      } else {
        dir <- "Size over 100 nodes"
      }
      write.table(network,
                  file = file.path(output.dir, dir, paste0("Create network for components ", i, "(", components.out$csize[i], " nodes)", ".csv")),
                  sep = ",", row.names = FALSE)
    }

    # Export data on nodes
    chia.obj$Regions$Component.Id <- components.out$membership
    chia.obj$Regions$Component.size <- components.out$csize[components.out$membership]
    # Output region annotation.
    write.table(as.data.frame(chia.obj$Regions), file = file.path(output.dir, "Annotated CHIA-PET regions.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

}

#' Analyze the topology of ChIA-PET data
#'
#' Analyzes ChIA-PET data and produces graphs according to the topology of the nodes in the whole set
#' of regions and in the separated components: \describe{
#'   \item{Histogram of number of edges.pdf} {Histogram of the number of edges in the whole set of data.}
#'   \item{Scatter plot of degrees of the left node versus the right node.pdf} {Scatter plot showing the relation between the left and right vertices' degree.}
#'   \item{Degree vs cluster width.pdf} {A plot of the degree of the nodes in one cluster according to the size of its cluster.}
#'   \item{Log2(size of component) vs Log2(Number of components).pdf} {A plot of the number of components according to the size of the components.}
#'   \item{Component table.txt} {A file with the information about TSS for each component.}
#'   \item{Proportion of TSS in low connectivity nodes.pdf} {A plot of the proportion of TSS in nodes with low connectivity (component with 5 nodes or less).}
#'   \item{Proportion of TSS in high connectivity nodes.pdf} {A plot of the proportion of TSS in nodes with high connectivity (component with more than 5 nodes).}}
#'
#' @param chia.obj A list containing the annotated ChIA-PET data, as returned by \code{\link{annotate.chia}}.
#' @param output.dir The name of the directory where to save the graphs.
#'
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @importFrom graphics hist
#' @import ggplot2
#' @importFrom igraph components
#' @importFrom igraph get.data.frame
#' @importFrom utils write.table
analyze.generic.topology <- function(chia.obj, output.dir="output") {
    # Plot an histogram of the number of edges.
    pdf(file.path(output.dir, "Histogram of number of edges.pdf"))
    hist(degree(chia.obj$Graph), breaks=seq(0, 300, by=5))
    dev.off()

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
    chia.components = components(chia.obj$Graph)
    chia.obj$Regions$Component = chia.components$membership
    component.table <- table(chia.components$csize)

    component.df = data.frame(Size=as.integer(names( component.table)), Number=as.vector(component.table))
    ggplot(component.df) + geom_point(mapping=aes(x=log2(Size), y=log2(Number)))
    ggsave("Log2(size of component) vs Log2(Number of components).pdf")

    annotate.component <- function(x) {
        result.df = data.frame(NumberOfTSS=sum(x$distanceToTSS==0),
                          Size=nrow(x),
                          data.frame(as.list((table(x$Simple.annotation)/nrow(x))), check.names=FALSE))

        if(!is.null(x$Chrom.State)) {
            result.df = cbind(result.df, data.frame(as.list((table(x$Chrom.State)/nrow(x))), check.names=FALSE))
        }

        return(result.df)
    }

    component.table = ddply(as.data.frame(chia.obj$Regions), "Component", annotate.component)
    write.table(component.table, file=file.path(output.dir, "Component table.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

    component.table = ddply(component.table, "Size", function(x) { return(cbind(x, NumberOfComponentsOfThisSize=nrow(x)))})
    proportion.per.size = ddply(component.table, ~NumberOfTSS * Size, function(x) { return(nrow(x)/x$NumberOfComponentsOfThisSize[1])})
    proportion.per.size[order(proportion.per.size$Size, proportion.per.size$NumberOfTSS),]

    ggplot(subset(proportion.per.size, Size <= 5), aes(x=NumberOfTSS, y=V1)) + geom_bar(stat="identity") + facet_wrap(~Size)
    ggsave(file.path(output.dir, "Proportion of TSS in low connectivity nodes.pdf"))

    ggplot(subset(component.table, Size >5), aes(x=NumberOfTSS/Size)) + geom_histogram()
    ggsave(file.path(output.dir, "Proportion of TSS in high connectivity nodes.pdf"))
}

#' Analyze the chromatin state of ChIA-PET data
#'
#' Analyzes ChIA-PET data and produces graphs according to its chromatin states: \describe{
#'   \item{Chromatin states summary.txt} {A file with the summary of chromatin states.}
#'   \item{log Degree histogram per chromatin state.pdf} {A histogram of the degree of each node according to the chromatin state.}
#'   \item{Proportion of chromatin state as a function of connectivity category.pdf} {A plot of the proportion of chromatin state as a fonction of connectivity category.}
#'   \item{Contact heatmap for chromatin states.pdf} {A contact heatmap of chromatin states.}}
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
#'   \item{Proportion of genomic location as a function of connectivity category.pdf} {A plot of the proportion of the genomic location of the regions as a fonction of connectivity category.}
#'   \item{Contact heatmap for genomic location.pdf} {A contact heatmap of the genomic location of the regions.}}
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
    ggplot(degree.exp.df, aes(x=log2(Degree), y=Exp.Mean)) +
        geom_point() +
        geom_smooth(method='lm')
    ggsave(file.path(output.dir, "Expression vs Degree at promoter.pdf"))
}

#' Analyze the gene specificity of ChIA-PET data
#'
#' Analyzes ChIA-PET data and produces graphs according to gene specificity: \describe{
#'   \item{Tau vs degree.pdf} {A plot of the Tau of the nodes according to their degree.}
#'   \item{Boxplot of degrees by expression category.pdf} {A boxplot the degrees of the nodes in each expression category.}}
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
analyze.tf <- function(chia.obj, tf.regions, output.dir="output") {
    # Look at TF presence curves as a function of connectivity
    results = matrix(0, nrow=length(tf.regions), ncol=4)
    rownames(results) <- names(tf.regions)
    boundaries.list = list(Singles=c(0, 1), Low=c(1, 5), Intermediate=c(5, 20), High=c(20, 1000))
    colnames(results) <- names(boundaries.list)

    # Extract the overlap matrix from the region annotations.
    overlap.matrix <- as.matrix(mcols(chia.obj$Regions)[,grepl("^TF\\.", colnames(mcols(chia.obj$Regions)))])
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


#' Analyze ChIA-PET data and produce graphs.
#'
#' @param chia.obj A list containing the annotated ChIA-PET data, as returned by \code{\link{annotate.chia}}
#' @param input.chrom.state The name of the file containing the information about chromatin states.
#' @param biosample The biosample identifier from ENCODE. Valid examples are GM12878, K562 or MCF-7.
#' @param genome.build The name of the chosen annotation ("hg38", "hg19").
#' @param output.dir The name of the directory where to save the graphs.
#' @importFrom Biobase cache
#' @export
analyze.chia.pet <- function(input.chia, input.chrom.state = NULL, biosample = NULL, genome.build = NULL, output.dir="output/") {
    dir.create(file.path(output.dir), recursive=TRUE, showWarnings=FALSE)

    chia.obj = load.chia(input.chia)

    tf.regions = NULL
    histone.regions = NULL
    expression.data = NULL
    if(!is.null(biosample) && !is.null(genome.build)) {
        tf.regions = download.encode.chip(biosample, genome.build)$Regions
        histone.regions <- download.encode.chip(biosample, genome.build, download.filter=histone.download.filter.chip)$Regions

        expression.data = download.encode.rna(biosample, genome.build)$Expression
        expression.data$ENSEMBL = gsub("\\.\\d+$", "", expression.data$gene_id)
        expression.data$FPKM = log2(expression.data$Mean.FPKM + 1)
    }

    # Download chromatin states
    if(!is.null(biosample) && is.null(input.chrom.state)) {
        input.chrom.state <- import.chrom.states(biosample, file.path("input/chrom_states", biosample))
    }

    Biobase::cache(chia.obj <- annotate.chia(chia.obj,
                                             input.chrom.state = input.chrom.state,
                                             tf.regions = tf.regions,
                                             histone.regions=histone.regions,
                                             expression.levels=expression.data,
                                             genome.build = genome.build,
                                             biosample=biosample,
                                             histone=TRUE,
                                             output.dir = output.dir), dir=output.dir, prefix="cached_objects")



    analyze.generic.topology(chia.obj, output.dir)
	analyze.annotation(chia.obj, output.dir)

	if(!is.null(input.chrom.state)) {
        analyze.chromatin.states(chia.obj, output.dir)
    }

	if(!is.null(expression.data)) {
        analyze.expression(chia.obj, output.dir)
    }

  if(!is.null(tf.regions)) {
      analyze.tf(chia.obj, tf.regions, output.dir)
  }

	if(genome.build %in% c("hg19", "hg38")) {
        analyze.gene.specificity(chia.obj, output.dir)
    }

}

#' Separetes a large networks into communities and saves the sub-netwoks in cytoscape-friendly format
#'
#' @param network.input The csv file containing the network to divide.
#' @param network The id of the network.
#' @param annotated.chia The txt file containing all information about the networks.
#' @param output.dir The directory where to save the files.
#' @param method The algorithm to use to divide the graph.
#' @importFrom igraph make_graph
#' @export
separate.into.communities <- function(network.input, network, annotated.chia, output.dir, method = igraph::cluster_fast_greedy){
  # Separate the network into smaller ones
  network.df <- read.table(network.input, header = TRUE, sep = ",")
  network.df <- aggregate(Reads~(Source + Target), data = network.df, sum)
  whole.graph <- make_graph(c(rbind(network.df$Source, network.df$Target)), directed = FALSE)

  communities <- method(whole.graph, weights = NULL)
  communities.df <- data.frame(cbind(unique(c(network.df$Source, network.df$Target)),
                                     communities$membership[unique(c(network.df$Source, network.df$Target))]))
  colnames(communities.df) <- c("Node.Id", "Group")
  simplify.group.id <- data.frame(Old = unique(communities.df$Group), New = 1:length(unique(communities.df$Group)))
  communities.df$Group <- simplify.group.id$New[match(communities.df$Group, simplify.group.id$Old)]
  communities.df$Size <- 1
  communities.df$Size <- aggregate(Size~Group, data = communities.df, FUN = sum)$Size[communities.df$Group]

  # Change the annotation
  annotated.chia <- read.csv(annotated.chia, header = TRUE, sep = "\t")
  annotated.chia$Component.size[annotated.chia$Component.Id == network] <- communities.df[order(communities.df$Node.Id),"Size"]
  annotated.chia$Component.Id[annotated.chia$Component.Id == network] <-
    paste0(network, ".", communities.df[order(communities.df$Node.Id),"Group"])

  # Write the new annotation and the new network tables
  dir.create(output.dir, recursive = TRUE)
  write.table(communities.df, file = file.path(output.dir, "Annotated CHIA-PET regions.txt"), sep = "\t", row.names = FALSE)
  for (group in simplify.group.id$New) {
    nodes.to.keep <- communities.df$Node.Id[communities.df$Group == group]
    sub.network.df <- network.df[network.df$Source %in% nodes.to.keep & network.df$Target %in% nodes.to.keep,]
    file.name <- file.path(output.dir, paste0("Create network for component", network, ".", group, "(",
                                              communities.df$Size[communities.df$Group == group][1], " nodes).csv"))
    write.table(sub.network.df, file = file.name, sep = ",", row.names = FALSE)
  }
  return(paste0(group, " sub-networks were found."))
}

# analyze.chia.pet(input.chia="input/ChIA-PET/GSM1872887_GM12878_RNAPII_PET_clusters.txt",
#                  input.chrom.state="input/E116_18_core_K27ac_mnemonics.bed",
#                  biosample="GM12878", genome.build="hg19",
#                  output.dir="output/GM12878 PolII Ruan")
#
# analyze.chia.pet(input.chia="input/ChIA-PET/MCF7 Ruan 2012.txt",
#                  input.chrom.state=NULL,
#                  biosample="MCF-7", genome.build="hg19",
#                  output.dir="output/MCF-7 PolII Ruan")



# For testing in console
# input.chia="input/ChIA-PET/GSM1872887_GM12878_RNAPII_PET_clusters.txt"
# input.chrom.state="input/E116_18_core_K27ac_mnemonics.bed"
# biosample="GM12878"
# genome.build="hg19"
# output.dir="output/GM12878 PolII Ruan"

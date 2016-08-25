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

        connectivity.df <- categorize.by.connectivity(chia.obj)
        chia.plot.metrics(chia.obj, level.counts, connectivity.df, "Connectivity", "Proportion of nodes in category",
                     graph.type = "line", facet.rows = 3, facet.cols = 6,
                     file.out = file.path(output.dir, "Proportion of chromatin state as a function of connectivity category.pdf"),
                     variable.name = "Chrom.State")
        contact.heatmap(chia.obj, "Chrom.State", "chromatin states", output.dir=output.dir)
        
        if(has.components(chia.obj)) {
            size.categories <- categorize.by.components.size(chia.obj)
            chia.plot.metrics(chia.obj, level.counts, size.categories, 
              x.lab = "Size category", y.lab = "Proportion",
              graph.type = "line", facet.rows = 3,
              file.out = file.path(output.dir, "Proportion of chromatin state as a function of size category.pdf"),
              variable.name = "Chrom.State", proportion = TRUE)
        }
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
        connectivity.df <- categorize.by.connectivity(chia.obj)
        chia.plot.metrics(chia.obj, level.counts, connectivity.df, "Connectivity", "Proportion of nodes in category",
                   graph.type = "line", facet.rows = 3, facet.cols = 3,
                   file.out = file.path(output.dir, "Proportion of genomic location as a function of connectivity category.pdf"),
                   variable.name = "Simple.annotation")
        contact.heatmap(chia.obj, "Simple.annotation", "genomic location", output.dir)
        
        if(has.components(chia.obj)) {
            size.categories <- categorize.by.components.size(chia.obj)
            chia.plot.metrics(chia.obj, level.counts, size.categories,
                 x.lab = "Size category", y.lab = "Proportion", graph.type = "line", facet.rows = 3,
                 file.out = file.path(output.dir, "Proportion of genomic location as a function of size category.pdf"),
                 variable.name = "Simple.annotation", proportion = TRUE)
        }
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
    if(has.transcription.factors(chia.obj)) {
        # Extract the overlap matrix from the region annotations.
        overlap.matrix <- as.matrix(get.tf(chia.obj))

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

#' Analyze the components of ChIA-PET data
#'
#' Analyzes ChIA-PET data and produces graphs according to the topology of the nodes in the whole set
#' of regions and in the separated components: \describe{
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
analyze.components <- function(chia.obj, output.dir="output") {
  # Analyze components
  if(has.components(chia.obj)) {
    component.table <- chia.obj$Regions[,c("Component.Id", "Component.size")]
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
    
    if(has.gene.annotation(chia.obj)) {
        tss.metric.function <- functor.constructor(boolean.count, "Is.TSS", proportion=TRUE)
        chia.plot.metrics(chia.obj, apply.single.metric.by.component(tss.metric.function, "TSS proportion"),
                          categorize.by.components.size(chia.obj), graph.type="boxplot",
                          file.out=file.path(output.dir, "Boxplot of the proportion of TSS in fonction of the size of the networks.pdf"))
    }
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
}

#' Performs all possible analysis regarding the central nodes.
#'
#' @param chia.obj A list containing the annotated ChIA-PET data, as returned by \code{\link{annotate.chia}}.
#' @param output.dir The directory where output should be saved.
#'
#' @export
analyze.central.nodes <- function(chia.obj, output.dir=".") {
  centrality.categories <- NULL
  if(has.centrality(chia.obj)) {
    centrality.categories <- categorize.by.centrality(chia.obj)
  }
  
  if(has.chrom.state(chia.obj)) {
    chia.plot.metrics(chia.obj, level.counts, centrality.categories, graph.type = "heatmap",
                 file.out = file.path(output.dir, "Heatmap of centrality vs chromatin states.pdf"), variable.name = "Chrom.State")
  }
  
  if(has.gene.annotation(chia.obj)) {
    chia.plot.metrics(chia.obj, level.counts, centrality.categories, graph.type = "heatmap",
                 file.out = file.path(output.dir, "Heatmap of centrality vs genomic locations.pdf"), variable.name = "Simple.annotation")
  }
  
  if(has.transcription.factors(chia.obj)) {
    chia.plot.metrics(chia.obj, calculate.tf.presence, centrality.categories, graph.type = "heatmap",
     x.lab = "Centrality category", y.lab = "Proportion", 
     file.out = file.path(output.dir, "Proportion of TF as a function of size centrality.pdf"),
     proportion = TRUE)
  }
}

#' Performs analyses which do not fit in any of the other sub-analyses.
#'
#' @param chia.obj A list containing the annotated ChIA-PET data, as returned by \code{\link{annotate.chia}}.
#' @param output.dir The directory where output should be saved.
#'
#' @export
analyze.misc <- function(chia.obj, output.dir=".") {
    if(has.fitness(chia.obj)) {
        chia.plot.metrics(chia.obj, subset.counts, categorize.by.connectivity(chia.obj), x.lab = "Connectivity",
                          graph.type = "histogram", file.out = file.path(output.dir, "Histogram of the number of essential genes by connectivity category.pdf"),
                          all.conditions = "(Gene.Representative & (Fitness > 0.5))", proportion = FALSE)
    }
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

    cat(date(), " : Analyzing network components...\n",cat.sink)
    analyze.components(chia.obj, output.dir)
    
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

    cat(date(), " : Analyzing central nodes...\n",cat.sink)
    analyze.central.nodes(chia.obj, output.dir)
    
    cat(date(), " : Analyzing miscellaneous...\n",cat.sink)
    analyze.misc(chia.obj, output.dir)
    
    # Close dummy verbose stream.
    if(!verbose) {
        close(cat.sink)
    }
}

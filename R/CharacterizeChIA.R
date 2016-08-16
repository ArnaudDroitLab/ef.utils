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

has.chrom.state <- function(chia.obj) {
    return(!is.null(chia.obj$Regions$Chrom.State))
}

has.components <- function(chia.obj) {
    return(!is.null(chia.obj$Regions$Component.Id) && !is.null(chia.obj$Regions$Component.size))
}

has.gene.specificity <- function(chia.obj) {
    return(!is.null(chia.obj$Regions$Gene.Representative) &&
           !is.null(chia.obj$Regions$Expression.Tau) &&
           !is.null(chia.obj$Regions$Expression.Category))
}


has.degree <- function(chia.obj) {
    return(!is.null(chia.obj$Regions$Degree))
}

has.gene.representative <- function(chia.obj) {
    return(!is.null(chia.obj$Regions$Gene.Representative))
}

has.expression.levels <- function(chia.obj) {
    return(!is.null(chia.obj$Regions$Gene.Representative) &&
           !is.null(chia.obj$Regions$Expr.mean))
}

has.gene.annotation <- function(chia.obj) {
    return(!is.null(chia.obj$Regions$Simple.annotation))
}


#' Return the number of nodes in a CHIA object.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load.chia}}.
#'
#' @return The number of nodes in the chia object.
#' @importFrom igraph vcount
#' @export
number.of.nodes <- function(chia.obj) {
    node.count = vcount(chia.obj$Graph)
    stopifnot(nrow(chia.obj$Regions) == node.count)

    return(vcount(chia.obj$Graph))
}

#' Return the number of contacts in a CHIA object.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load.chia}}.
#'
#' @return The number of contacts in the chia object.
#' @importFrom igraph ecount
#' @export
number.of.contacts <- function(chia.obj) {
    return(ecount(chia.obj$Graph))
}

#' Return the number of components in a CHIA object.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load.chia}}.
#'
#' @return The number of components in the chia object.
#' @importFrom igraph components
#' @export
number.of.components <- function(chia.obj) {
    return(components(chia.obj$Graph)$no)
}

#' Return the mean component size of a CHIA object.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load.chia}}.
#'
#' @return The mean component size of the chia object.
#' @importFrom igraph components
#' @export
mean.component.size <- function(chia.obj) {
    return(mean(components(chia.obj$Graph)$csize))
}

#' Return the number of genes represented in a CHIA object.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load.chia}}.
#'
#' @return The number of genes represented in the chia object.
#' @export
number.of.genes <- function(chia.obj) {
    return(sum(chia.obj$Regions$Gene.Representative))
}

#' Return the number of active genes represented in a CHIA object.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load.chia}}.
#'
#' @return The number of active genes represented in the chia object.
#' @export
number.active.genes <- function(chia.obj) {
    return(sum(chia.obj$Regions$Gene.Representative & chia.obj$Regions$Is.Gene.Active))
}

#' Return the number genes per component in the CHIA object.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load.chia}}.
#'
#' @return The number of active genes represented in the chia object.
#' @export
genes.by.component <- function(chia.obj) {
    stopifnot(has.gene.representative(chia.obj) && has.components(chia.obj))
    return(aggregate(Gene.Representative~Component.Id, as.data.frame(chia.obj$Regions), sum))
}

#' Return the proportion of nodes representing genes in a CHIA object.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load.chia}}.
#'
#' @return The number of active genes represented in the chia object.
proportion.genes <- function(chia.obj) {
    return(number.of.genes(chia.obj) / number.of.nodes(chia.obj))
}

#' Return the proportion of active genes among all genes in a CHIA object.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load.chia}}.
#'
#' @return The proportion of active genes among all genes in the chia object.
proportion.active.genes <- function(chia.obj) {
    return(number.active.genes(chia.obj) / number.of.genes(chia.obj))
}

regions.to.vertex.attr <- function(chia.obj) {
  vertex_attr(chia.obj$Graph) <- as.data.frame(chia.obj$Regions, stringsAsFactors = FALSE)
  return(chia.obj)
}

vertex.attr.to.regions <- function(graph.obj) {
  region.df = as.data.frame(vertex_attr(graph.obj), stringsAsFactors = FALSE)
  region.df$strand = '*'
  return(GRanges(region.df))
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
#' @importFrom igraph edge_attr
#' @importFrom igraph set_edge_attr
#' @importFrom igraph edge_attr<-
#' @importFrom plyr ddply
#' @importFrom plyr summarize
#'
#' @export
load.chia <- function(input.chia) {
    chia.raw = read.table(input.chia, sep="\t",
                          col.names=c("L.chr", "L.start", "L.end", "R.chr", "R.start", "R.end", "Reads"))

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

    # Find and remove self loops.
    mapped.df = cbind(chia.raw, Left=chia.left.indices, Right=chia.right.indices)
    mapped.df = mapped.df[mapped.df$Left != mapped.df$Right,]

    # Summarize multiple edges.
    max.df = ddply(mapped.df, ~Left*Right, summarize, L.chr=head(L.chr, n=1), L.start=min(L.start), L.end=max(L.end),
                                                      R.chr=head(R.chr, n=1), R.start=min(R.start), R.end=max(R.end), Reads=sum(Reads))

    # Create iGraph object and set the original coordinates and the number of supporting reads as edge attributes.
    chia.graph = make_graph(c(rbind(max.df$Left, max.df$Right)), directed=FALSE)
    edge_attr(chia.graph) <- max.df

    return(list(Regions=single.set, Graph=chia.graph))
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
  dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

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
  dir.create(file.path(output.dir, "Size between 3 and 5 nodes (incl)"), recursive = TRUE, showWarnings=FALSE)
  dir.create(file.path(output.dir, "Size between 6 and 20 nodes (incl)"), recursive = TRUE, showWarnings=FALSE)
  dir.create(file.path(output.dir, "Size between 21 and 50 nodes (incl)"), recursive = TRUE, showWarnings=FALSE)
  dir.create(file.path(output.dir, "Size between 51 and 100 nodes (incl)"), recursive = TRUE, showWarnings=FALSE)
  dir.create(file.path(output.dir,"Size over 100 nodes"), recursive = TRUE, showWarnings=FALSE)

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


#' Identifies edges crossing community borders.
#'
#' @param input.graph The igraph whose community-crossing edges must be identified,
#' @param method A method returning a \code{community} object, such as igraph::cluster_fast_greedy.
#' @param weight.attr The name fo the edge attribute to be used as edge weight.
#'
#' @return The ids of the edges to be removed.
#' @importFrom igraph crossing
#' @importFrom igraph as.undirected
#' @importFrom igraph edge_attr
#' @export
identify.crossing.edges <- function(input.graph, method = igraph::cluster_fast_greedy, weight.attr=NULL){
  communities <- method(as.undirected(input.graph), weights = weight.attr)
  to.delete = crossing(communities, input.graph)
  return(edge_attr(input.graph)$original.id[to.delete])
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
#' @importFrom igraph induced_subgraph
#' @importFrom igraph set_edge_attr
#' @export
split.by.community <- function(chia.obj, oneByOne = FALSE, method = igraph::cluster_fast_greedy, weight.attr=NULL) {
  edge_attr(chia.obj$Graph)$original.id = 1:ecount(chia.obj$Graph)
  if (oneByOne){
    # Keep a record of edges marked for deletion, so that we can delete them all
    # at once. This prevents issues with edges being relabeled.
    marked.for.deletion = c()

    # Loop over components one by one.
    components.out = components(chia.obj$Graph)
    for (i in 1:components.out$no){
      # Get the component subgraph.
      component.subgraph = induced_subgraph(chia.obj$Graph, components.out$membership==i)

      # Split it into communities and record teh deleted edges.
      crossing.edges = identify.crossing.edges(component.subgraph, method = method, weight.attr=weight.attr)
      marked.for.deletion = c(marked.for.deletion, crossing.edges)
    }

    # Delete removed edges in the original chia object.
    chia.obj$Graph = delete_edges(chia.obj$Graph, marked.for.deletion)
  } else {
    # Split it into communities and delete the necessary edges.
    crossing.edges = identify.crossing.edges(chia.obj$Graph, method = method, weight.attr=weight.attr)
    chia.obj$Graph = delete_edges(chia.obj$Graph, crossing.edges)
  }

  # Remove the original.id edge attribute, since it won't be needed anymore.
  edge_attr(chia.obj$Graph)$original.id = NULL

  # Update the degree attribute of regions if it is present.
  chia.obj$Regions$Degree = degree(chia.obj$Graph)
  return(chia.obj)
}

#' Subset a CHIA object.
#'
#' @param chia.obj The chia object to eb subset.
#' @param indices The indices of the vertices to be kept.
#'
#' @return A chia object containing only the selected vertices.
#' @importFrom igraph induced_subgraph
#' @export
chia.vertex.subset <- function(chia.obj, indices) {
    return(list(Regions = chia.obj$Regions[indices],
                Graph = induced_subgraph(chia.obj$Graph, which(indices))))
}

#' Analyze ChIA-PET data and produce graphs.
#'
#' @param input.chia The file containing processed ChIA-PET data.
#' @param output.dir The name of the directory where output should be saved.
#' @return The annotated chia.obj.
#' @importFrom Biobase cache
#' @export
process.chia.pet <- function(input.chia, chia.param, output.dir="output", verbose=TRUE) {
    # Create output directory.
    dir.create(file.path(output.dir), recursive=TRUE, showWarnings=FALSE)

    # Load interaction data.
    chia.obj = load.chia(input.chia)

    # Annotate the ChIA object.
    chia.obj <- annotate.chia(chia.obj, chia.param, output.dir=output.dir, verbose=verbose)

    # Analyze the ChIA network.
    analyze.chia.pet(chia.obj, output.dir)

    # Return teh created object.
	return(chia.obj)
}

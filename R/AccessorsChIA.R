#' Return the left part of the ChIA-PET data.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load.chia}}.
#'
#' @return A \linkS4class{GRanges} object with the \code{Regions} from the "\code{chia.obj}" parameter corresponding the the "left side".
#' of the original data.
#' @importFrom igraph as_edgelist
chia.left <- function(chia.obj) {
    return(chia.obj$Regions[as_edgelist(chia.obj$Graph)[,1],])
}

#' Return the right part of the ChIA-PET data.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load.chia}}.
#'
#' @return A \linkS4class{GRanges} object with the \code{Regions} from the "\code{chia.obj}" parameter corresponding the the "right side".
#' of the original data.
#' @importFrom igraph as_edgelist
chia.right <- function(chia.obj) {
    return(chia.obj$Regions[as_edgelist(chia.obj$Graph)[,2],])
}

#' Return the right part of the ChIA-PET data.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load.chia}}.
#'
#' @return A \linkS4class{GRanges} object with the \code{Regions} from the "\code{chia.obj}" parameter corresponding the the "right side".
#' of the original data.
#' @export
get.granges <- function(chia.obj) {
    return(GRanges(chia.obj$Regions))
}

#' Determines if the given chia.obj has associated chromatin states.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load.chia}}.
#'
#' @return True if the object has associated chromatin states.
#' @export
has.chrom.state <- function(chia.obj) {
    return(!is.null(chia.obj$Regions$Chrom.State))
}

#' Determines if the given chia.obj has associated components.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load.chia}}.
#'
#' @return True if the object has associated components.
#' @export
has.components <- function(chia.obj) {
    return(!is.null(chia.obj$Regions$Component.Id) && !is.null(chia.obj$Regions$Component.size))
}

#' Determines if the given chia.obj has associated gene specificities.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load.chia}}.
#'
#' @return True if the object has associated gene specificities.
#' @export
has.gene.specificity <- function(chia.obj) {
    return(!is.null(chia.obj$Regions$Gene.Representative) &&
           !is.null(chia.obj$Regions$Expression.Tau) &&
           !is.null(chia.obj$Regions$Expression.Category))
}

#' Determines if the given chia.obj has associated node degrees.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load.chia}}.
#'
#' @return True if the object has associated node degrees.
#' @export
has.degree <- function(chia.obj) {
    return(!is.null(chia.obj$Regions$Degree))
}

#' Determines if the given chia.obj has associated gene representatives.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load.chia}}.
#'
#' @return True if the object has associated gene representatives.
#' @export
has.gene.representative <- function(chia.obj) {
    return(!is.null(chia.obj$Regions$Gene.Representative))
}

#' Determines if the given chia.obj has associated expression levels.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load.chia}}.
#'
#' @return True if the object has associated expression levels.
#' @export
has.expression.levels <- function(chia.obj) {
    return(!is.null(chia.obj$Regions$Gene.Representative) &&
           !is.null(chia.obj$Regions$Expr.mean))
}

#' Determines if the given chia.obj has associated gene annotations.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load.chia}}.
#'
#' @return True if the object has associated gene annotations.
#' @export
has.gene.annotation <- function(chia.obj) {
    return(!is.null(chia.obj$Regions$Simple.annotation))
}

#' Determines if the given chia.obj has associated node centralities.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load.chia}}.
#'
#' @return True if the object has associated node centralities.
#' @export
has.centrality <- function(chia.obj) {
    return(!is.null(chia.obj$Regions$Is.central))
}

#' Determines if the given chia.obj has TF binding information.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load.chia}}.
#'
#' @return True if the object has TF binding information.
#' @export
has.transcription.factors <- function(chia.obj) {
    return(sum(grepl("^TF", colnames(chia.obj$Regions))) > 0)
}

#' Determines if the given chia.obj has TF binding information.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load.chia}}.
#'
#' @return True if the object has TF binding information.
#' @export
has.fitness <- function(chia.obj) {
    return(!is.null(chia.obj$Regions$Fitness))
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
average.component.size <- function(chia.obj) {
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
#' @export
proportion.genes <- function(chia.obj) {
    return(number.of.genes(chia.obj) / number.of.nodes(chia.obj))
}

#' Return the proportion of active genes among all genes in a CHIA object.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load.chia}}.
#'
#' @return The proportion of active genes among all genes in the chia object.
#' @export
proportion.active.genes <- function(chia.obj) {
    return(number.active.genes(chia.obj) / number.of.genes(chia.obj))
}

#' Obtain the matrix of TF binding.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load.chia}}.
#'
#' @return A matrix indicating which TF binds to which regions.
#' @export
get.tf <- function(chia.obj) {
    stopifnot(has.transcription.factors(chia.obj))
    return(chia.obj$Regions[,grepl("TF", colnames(chia.obj$Regions))])
}

#' Returns a list of all statistics for a given ChIA object.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load.chia}}.
#'
#' @return A named list of calculated properties for the given ChIA object.
#' @export
get.all.statistics <- function(chia.obj) {
    results = lapply(list("Number of nodes"            = number.of.nodes,
                          "Number of contacts"         = number.of.contacts,
                          "Number of components"       = number.of.components,
                          "Average component size"     = average.component.size,
                          "Number of genes"            = number.of.genes,
                          "Number of active genes"     = number.active.genes,
                          "Proportion of genes"        = proportion.genes,
                          "Proportion of active genes" = proportion.active.genes),
                    function(f) f(chia.obj))

    return(results)
}
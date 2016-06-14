#' Internal exclusive overlap.
#'
#' Returns the set of regions within code{all.regions} where all
#' conditions/proteins at indices which.factors are present,
#' and no other factors are. Presence/absence at a given locus
#' is obtained from the \code{overlap.matrix} parameter.
#'
#' @param all.regions A \linkS4class{GRanges} object with the regions represented
#    by the rows of \code{overlap.matrix}.
#' @param overlap.matrix A matrix, with rows corresponding to regions
#'   in \code{all.rgions}, and columns corresponding to proteins/factors. A
#'   non-zero value in the matrix indicates the protein/factor of interest
#'   is present at this region.
#' @param which.factors Indices of columns where the factor should be present.
#' @return A \linkS4class{GRanges} object with the regions matching the given criteria.
exclusive.overlap.internal <- function(all.regions, overlap.matrix, which.factors) {
    has.factor = rep(TRUE, nrow(overlap.matrix))
    if(sum(which.factors) != 0) {
        has.factor = apply(overlap.matrix[,  which.factors, drop=FALSE] >= 1, 1, all)
    }

    no.others = rep(TRUE, nrow(overlap.matrix))
    if(sum(!which.factors) != 0) {
        no.others  = apply(overlap.matrix[, !which.factors, drop=FALSE] == 0, 1, all)
    }

    return(all.regions[has.factor & no.others])
}

#' Internal inclusive overlap.
#'
#' Returns the set of regions within \code{all.regions} where all
#' conditions/proteins at indices which.factors are present,
#' regardless of whether or not other factors are. Presence/absence
#' at a given locus is obtained from the overlap.matrix parameter.
#'
#' @param all.regions A \linkS4class{GRanges} object with the regions represented
#    by the rows of \code{overlap.matrix}.
#' @param overlap.matrix A matrix, with rows corresponding to regions
#'   in \code{all.rgions}, and columns corresponding to proteins/factors. A
#'   non-zero value in the matrix indicates the protein/factor of interest
#'   is present at this region.
#' @param which.factors Indices of columns where the factor should be present.
#' @return A \linkS4class{GRanges} object with the regions matching the given criteria.
inclusive.overlap.internal <- function(all.regions, overlap.matrix, which.factors) {
    has.factor = apply(overlap.matrix[,  which.factors, drop=FALSE] >= 1, 1, all)
    return(all.regions[has.factor])
}

#' Calculate an overlap of certain factors within an intersect.object.
#'
#' Given an \code{intersect.object} , finds all regions where all
#' factors described by either indices or names are present.
#' If both indices and names are \code{NULL}, the inner intersection
#' is calculated.
#'
#' @param intersect.object An \code{intersect.object} returned by \code{\link{build.intersect}}.
#'   by the rows of \code{overlap.matrix}.
#' @param indices A vector of indices into the factors of \code{intersect.object}.
#' @param names A vector of factor names in \code{intersect.object}.
#' @param exclusive If \code{TRUE}, a region will be returned if the factors in
#    indices or names are the ONLY the factors present at that region.
#' @return A \linkS4class{GRanges} object with the regions matching the given criteria.
#' @export
intersect.overlap <- function(intersect.object, indices=NULL, names=NULL, exclusive=FALSE) {
    if(is.null(indices) && is.null(names)) {
       which.factors = rep(TRUE, intersect.object$Length)
    } else if (!is.null(names)) {
       which.factors = intersect.object$Names %in% names
    } else {
        which.factors = indices
    }
    if(exclusive) {
        return(exclusive.overlap.internal(intersect.object$Region, intersect.object$Matrix, which.factors))
    } else {
        return(inclusive.overlap.internal(intersect.object$Region, intersect.object$Matrix, which.factors))
    }
}


#' Given a \linkS4class{GRangesList} object, determine which items overlap each others.
#'
#' @param grl The \linkS4class{GRangesList} object whose elements need to be overlapped with
#' each others.
#' @return A list with the following elements: \describe{
#'   \item{Regions} {A \linkS4class{GRanges} object with all genomic ranges occupied by at least one item.
#'   All ranges are "flattened", so if two of the initial ranges overlapped each other
#'   imperfectly, they are combined into a single region spanning both of them.}
#'   \item{Matrix} {A matrix, with \code{ncol=} the number of items in the initial \linkS4class{GRangesList} and \code{nrow=}
#'   the total number of loci, which is equal to the length of \code{Regions}. A value of 1 or more
#'   within the matrix indicates that the regions described by the column overlapped
#'   the region defined by the row.}
#'   \item{List} {A list of \code{length(grl)} numeric vectors indicating which indices of \code{Regions} overlap
#'   with the given condition/protein. Useful to translate the regions into unique names
#'   for drawing venn diagrams.}
#'   \item{Names} {The names of the initial grl items, corresponding to the column names of \code{Matrix} and the names
#'   of the element of \code{List}.}
#'   \item{Length} {The number of items in the initial \linkS4class{GRangesList}, corresponding to the number of columns in \code{Matrix}
#'   and the number of elements in \code{List}}.}
#' @importFrom GenomicRanges reduce
#' @importFrom GenomicRanges countOverlaps
#' @export
build.intersect <- function(grl) {
    # Flatten the GRangesList so we can get a list of all possible regions.
    #all.regions = biovizBase::flatGrl(GenomicRanges::reduce(unlist(grl)))
    all.regions = GenomicRanges::reduce(unlist(grl))

    # Build a matrix to hold the results.
    overlap.matrix <- matrix(0, nrow=length(all.regions), ncol=length(grl))
    overlap.list = list()

    # Loop over all ranges, intersecting them with the flattened list of all possible regions.
    for(i in 1:length(grl)) {
        overlap.matrix[,i] <- GenomicRanges::countOverlaps(all.regions, grl[[i]], type="any")
        overlap.list[[ names(grl)[i] ]] <- which(overlap.matrix[,i] != 0)
    }
    colnames(overlap.matrix) <-  names(grl)

    return(list(Regions = all.regions, Matrix=overlap.matrix, List=overlap.list, Names=colnames(overlap.matrix), Length=ncol(overlap.matrix)))
}

#' Generates a venn diagram from an \code{intersect.object}.
#'
#' @param intersect.object An intersect object returned by \code{\link{build.intersect}}.
#' @param filename A filename for the resulting venn diagram. Pass \code{NULL}
#'   to skip saving to a file.
#' @return The grid object representing the venn.diagram.
#' @importFrom VennDiagram venn.diagram
#' @export
plot.intersect.venn <- function(intersect.object, filename=NULL) {
    if(intersect.object$Length > 5) {
        stop("Cannot plot venn diagram of more than 5 groups!")
    }

    return(VennDiagram::venn.diagram(intersect.object$List,
                                     fill=c("red", "yellow", "green", "blue", "orange")[1:intersect.object$Length],
                                     filename=filename,
                                     print.mode="raw"))
}

#' Annotate the inner group of an \code{intersect.object}.
#'
#' Given an \code{intersect.object}, generate annotations for the regions
#' where all factors are present.
#'
#' @param intersect.object An intersect object returned by \code{\link{build.intersect}}.
#' @param annotations.list A list of annotation objects returned by \code{\link{select.annotations}}
#' @param filename A filename for the resulting annotations. Pass \code{NULL}
#'   to skip saving to a file.
#' @return The annotations for the intersect's inner regions.
#' @export
annotate.venn.center <- function(intersect.object, annotations.list, filename=NULL) {
    #overlap.regions = exclusive.overlap(intersect.object, rep(TRUE, intersect.object$Length))
    overlap.regions <- intersect.overlap(intersect.object, exclusive = TRUE)

    return(annotate.region(overlap.regions, annotations.list, filename))
}

#' Annotate the outer groups of an \code{intersect.object}.
#'
#' Given an \code{intersect.object}, generate annotations for the regions
#' where only one factor is present. One annotation per factor is generated.
#'
#' @param intersect.object An intersect object returned by \code{\link{build.intersect}}.
#' @param annotations.list A list of annotation objects returned by \code{\link{select.annotations}}
#' @param file.prefix A prefix for the names of the output files where the
#'   annotations will be written. Pass \code{NULL} to skip saving to a file.
#' @return A list with the generated annotations.
#' @export
annotate.venn.exclusive <- function(intersect.object, annotations.list, file.prefix=NULL) {
    results = list()
    for(i in 1:intersect.object$Length) {
        which.factors = 1:intersect.object$Length == i
        subset.regions = intersect.overlap(intersect.object, which.factors, exclusive = TRUE)

        if(length(subset.regions) > 0) {
            factor.name = intersect.object$Names[i]
            if(is.null(file.prefix)) {
                output.file = NULL
            } else {
                output.file = paste0(file.prefix, "Annotation for ", factor.name, " specific.txt")
            }

            results[[factor.name]] = annotate.region(subset.regions, annotations.list, output.file)
        }
    }

    return(results)
}

#' Calculate the pairwise overlaps of all factors within an \code{intersect.object}.
#'
#' @param intersect.object An intersect object returned by \code{\link{build.intersect}}.
#' @param filename A name for prefix for the names of the output files where the
#'   annotations will be written. Pass \code{NULL} to skip saving to a file.
#' @return A matrix containing the pairwise overlaps of all factors in
#' \code{intersect.object}. The row's factor is used as a denominator. Therefore, the
#    matrix is not symmetric.
#' @export
pairwise.overlap <- function(intersect.object, filename=NULL) {
    overlap.percentage <- matrix(0, nrow=intersect.object$Length, ncol=intersect.object$Length,
                                 dimnames=list(intersect.object$Names, intersect.object$Names))

    # Compare factors two by two.
    for(i in 1:intersect.object$Length) {
        for(j in 1:intersect.object$Length) {
            i.vector = intersect.object$Matrix[,i] >= 1
            j.vector = intersect.object$Matrix[,j] >= 1
            overlap.percentage[i,j] = sum(i.vector & j.vector) / sum(i.vector)
        }
    }

    if(!is.null(filename)) {
        write.table(overlap.percentage, file=filename, sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)
    }

    return(overlap.percentage)
}

#' Builds an overlap of the given regions and output all default annotations.
#'
#' After building the overlap, this function will generate: \enumerate{
#'   \item A venn diagram.
#'   \item An annotation of the intersect's inner regions.
#'   \item An annotation of the intersect's outer regions.
#'   \item The pairwise overlap of all factors.}
#'
#' @param regions A \linkS4class{GRangesList} of the regions to be intersected and analyzed.
#' @param annotations.list A list of annotation objects returned by \code{\link{select.annotations}}
#' @param label A label to use when generating file names.
#' @return The generated intersect object.
#' @export
build.intersect.all <- function(regions, annotations.list, label) {
    base.dir = file.path("output/", label)
    dir.create(base.dir, showWarnings = FALSE, recursive = TRUE)

    intersect.object = build.intersect(regions)

    plot.intersect.venn(intersect.object, file.path(base.dir, "Venn diagram.tiff"))
    annotate.venn.center(intersect.object, annotations.list, file.path(base.dir, "Venn intersection annotation.txt"))
    annotate.venn.exclusive(intersect.object, annotations.list, file.path(base.dir, "/"))
    pairwise.overlap(intersect.object, file.path(base.dir, "Pairwise overlap.txt"))

    return(intersect.object)
}

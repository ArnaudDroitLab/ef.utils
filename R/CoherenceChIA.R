#' Determines if certain attributes are coherent within certain categories of nodes.
#'
#' @param chia.obj The base ChIA object on which to assess intra-category coherence.
#' @param coherence.function A function which takes in a ChIA object and returns a 
#'   vector of a 2-level factor of the same length as the number fo regions in the 
#'   ChIA object.
#' @param node.categories A data-frame with a number of rows equal to the number of regions
#'   in chia.obj, where each column represents a set of indices indicating which nodes belong
#'   in a givenc ategory.
#' @param output.file The name of the file where the resulting graph should be saved.
#'
#' @return A list containing the computed data and the generated plot object.
#'
#' @import ggplot2
#' @export
coherence.test <- function(chia.obj, coherence.function, node.categories, output.file=NULL) {
  # Apply the coherence function
  coherence.list = category.apply(chia.obj, coherence.function, node.categories)
  
  # Count the success/failures from each category.
  coherence.counts = lapply(coherence.list, function(x) {
    tmp = c(sum(x==levels(x)[1], na.rm=TRUE), sum(x==levels(x)[2], na.rm=TRUE))
    names(tmp) = levels(x)
    return(tmp)
  })
  
  # Turn the count list into a data-frame.
  coherence.df = data.frame(matrix(unlist(coherence.counts), ncol=2, byrow=TRUE))
  colnames(coherence.df) <- names(coherence.counts[[1]])
  rownames(coherence.df) <- colnames(node.categories)
  
  # Calculate additional statistics
  coherence.df$Ratio = coherence.df[,2] / coherence.df[,1]
  coherence.df$Diff = coherence.df[,2] - coherence.df[,1]
  coherence.df$Category = factor(rownames(coherence.df), rownames(coherence.df)[order(coherence.df$Diff)])

  # Remove categories which do not have at least two elements.
  coherence.df.subset = coherence.df[coherence.df[,1] + coherence.df[,2] > 1,]
  
  # Generate a plot representing the data.
  plot.obj = ggplot(coherence.df.subset, aes(x=Category)) +
    geom_bar(mapping=aes(y=eval(parse(text=colnames(coherence.df)[2]))), stat="identity", fill="Blue", color="black") +
    geom_bar(mapping=aes(y=eval(parse(text=paste0("-",colnames(coherence.df)[1])))), stat="identity", fill="Red", color="black") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  if(!is.null(output.file)) {
    ggsave(output.file, width=14, height=7, plot=plot.obj)
  }
  
  return(list(Plot=plot.obj, Data=coherence.df))
}

#' Given a ChIA object, determine which genes are active and which are not.
#'
#' @param chia.obj The base ChIA object on which to assess gene activity.
#'
#' @return A 2-level vector with values "Inactive", "Active" or NA for each node.
#'
#' @importFrom plyr revalue
#' @export
active.gene.coherence <- function(chia.obj) {
    tmp = ifelse(chia.obj$Regions$Gene.Representative, chia.obj$Regions$Is.Gene.Active, NA)
    tmp = revalue(factor(tmp, levels=c(FALSE, TRUE)), c("FALSE"="Inactive", "TRUE"="Active"))
    return(tmp)
}

#' Generates a fucntion which takes a chia object, and determine which of its nodes
#' have an positive or negative fold-change based on the given column name.
#'
#' @param fc.column The column the new function should check for fold-change values.
#'
#' @return A function taking a ChIA object, and returning a 2-level vector with 
#'   "Upregulated", "Downregulated" or NA depending on the value of the given 
#'   fold-change column.
#'
#' @importFrom GenomicRanges mcols
#' @export
fold.change.coherence <- function(fc.column) {
    force(fc.column)
    function(chia.obj) {
        values = chia.obj$Regions[,fc.column]
        direction = factor(ifelse(values > 0, "Upregulated", "Downregulated"), levels=c("Downregulated", "Upregulated"))
        
        return(direction)
    }
}
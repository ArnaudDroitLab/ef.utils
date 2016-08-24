#' Converts a list of metrics into a data-frame for plotting.
#'
#' @param The list of metrics returned by mulitple call to a metric.function.
#' @return A data-frame representing the list values in a structure appropriate for ggplot.
metric.list.to.df <- function(metric.list) {
  # There are two posibilities: either we have multidimensional metrics (Multiple points per metric)
  # or unidimensional metrics (one point by metric). Infer which case we are in by input type.
  if(is.list(metric.list[[1]])) {
    # Get the list of all values into a single vector.
    all.values = unlist(lapply(metric.list, unlist))
    
    # Now, figure out which metric goes with each value. First, find out how
    # many times each metric is repeated within each category.
    per.category.metric.repeats = unlist(lapply(metric.list, function(x) { lapply(x, length) }))
  
    # Now just repeat the metric names this exact number of times.
    all.metrics = rep(unlist(lapply(metric.list, names)), unlist(per.category.metric.repeats))
    
    # Finally, figure out which category each value belong to.
    # Figure out the total number of metrics for each categoryé
    category.repeats = lapply(lapply(metric.list, unlist), length)
    
    # Now repeat each category name that exact number of times.
    all.categories = rep(names(metric.list), category.repeats)
    
    # Finally, get it all into the data.frame
    metric.df = data.frame(value    = all.values,
                           Metric   = factor(all.metrics, levels=names(metric.list[[1]])),
                           Category = factor(all.categories, levels=names(metric.list)))
  } else {
    # Get the metrics into a matrix
    metric.matrix = matrix(unlist(metric.list), ncol=length(metric.list),
                           dimnames = list(names(metric.list[[1]]), names(metric.list)))
    
    # Convert it to data frame for plotting
    metric.df = melt(metric.matrix, varnames=c("Metric", "Category"))  
  }
  
  return(metric.df)
}

#' Plots a subset of ChIA-PET nodes
#'
#' @param chia.obj A list containing the annotated ChIA-PET data, as returned by \code{\link{annotate.chia}}.
#' @param metric.function A function taking a chia.obj parameter and returning metrics in order to construct the plot.
#' @param node.categories Indices of nodes representing subsets of these nodes.
#' @param x.lab The label to add to the x axis of the gaph.
#' @param y.lab The label to add to the y axis of the graph.
#' @param title The title to add to the graph.
#' @param graph.type Either "line", "histogram", "heatmap".
#' @param facet.rows Optionnal graphic parameter.
#' @param facet.cols Optionnal graphic parameter.
#' @param file.out The name of the file where to save the plot, or NULL if none should be saved.
#' @param \dots Paramters to pass to metric.function.
#' @return The plot.
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export
chia.plot.metrics <- function(chia.obj, metric.function, node.categories, x.lab = NULL, y.lab = NULL,
                         title = NULL, graph.type = "line", facet.rows = NULL, facet.cols = NULL,
                         file.out=NULL, ...) {
  # Verify the valididy of arguments
  graph.type <- match.arg(graph.type, c("line", "histogram", "heatmap", "boxplot"))

  # Calculate metrics for all categories
  category.apply(chia.obj, metric.function, node.categories, ...)

  # Convert the list into a data.frame
  metric.df = metric.list.to.df(metric.list)

  # Plot metrics
  plot.obj = ggplot(data=metric.df) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

  # Type 1: line graph, comme nos graphs connectivity vs TF binding.
  # Chaque métrique devient une facette, avec une valeur par sous-catégorie de noeud.
  if(graph.type=="line") {
    plot.obj = plot.obj + geom_line(mapping=aes(x=Category, y=value), group = 1)
  }

  # Type 2: histogramme
  # Même chose que les graphes en ligne, mais en histogramme.
  if(graph.type=="histogram") {
    plot.obj = plot.obj + geom_bar(mapping = aes(x = Category, y = value), stat = "identity", colour = "black")
  }

  if(graph.type=="boxplot") {
    plot.obj = plot.obj + geom_boxplot(mapping = aes(x=Category, y=value))
  }
  
  if((graph.type=="line" || graph.type == "histogram" || graph.type == "boxplot") && (length(unique(metric.df$Metric)) > 1)) {
    plot.obj = plot.obj + facet_wrap(~Metric, nrow=facet.rows, ncol = facet.cols)
  }

  # Type 3: heatmap
  # Chaque sous-catégorie de noeuds est une colonne. Chaque métrique est une rangée.
  if(graph.type=="heatmap") {
    plot.obj = plot.obj + geom_tile(mapping=aes(x=Category, y=Metric, fill = value))
  }

  # Add labels and title to the graph
  plot.obj = plot.obj + xlab(x.lab) + ylab(y.lab) + ggtitle(title)

  if(!is.null(file.out)) {
    ggsave(plot.obj, file=file.out)
  }

  return(plot.obj)
}

#' Counts the occurences of each factor
#'
#' @param chia.subset A list containing a graph and ChIA-PET regions, as returned by \code{\link{chia.vertex.subset}}.
#' @param variable.name The name of the column containing a factor.
#' @param proportion Should the number of occurences be converted to proportions?
#' @return A named vector (table) with the number of occurences (or proportions) of each factor in the data.
#' @importFrom GenomicRanges mcols
#' @export
level.counts <- function(chia.subset, variable.name, proportion = TRUE){
  # Convert data into data frame
  variable.data <- mcols(chia.subset$Regions)[,variable.name]
  # Convert into factor if necessary
  if(class(variable.data) != "factor"){
    variable.data <- as.factor(variable.data)
    warning("The column given as parameter is not a factor. It has been converted as one.")
  }
  # Count the occurence of each factor
  table <- table(variable.data)
  # If proportion are needed, divide by the total number of nodes
  if (proportion){
    table <- table / length(chia.subset$Regions)
  }
  return(table)
}

#' Counts the occurences of categories made with a variable
#'
#' @param chia.subset A list containing a graph and ChIA-PET regions, as returned by \code{\link{chia.vertex.subset}}.
#' @param variable.name The name of the column containing a factor.
#' @param proportion Should the number of occurences be converted to proportions?
#' @return A named vector (table) with the number of occurences (or proportions) of each category.
#' @importFrom GenomicRanges mcols
#' @export
count.cut <- function(chia.subset, variable.name, proportion = TRUE, ...){
  # Convert data into data frame
  variable.data <- mcols(chia.subset$Regions)[,variable.name]
  # Cut variable in categories
  cut <- cut(variable.data, ...)
  # Count the occurence of each category
  table <- table(cut)
  # If proportion are needed, divide by the total number of nodes
  if (proportion){
    table <- table / length(chia.subset$Regions)
  }
  return(table)
}

#' Counts the number of nodes respecting the given condition(s)
#'
#' @param chia.subset A list containing a graph and ChIA-PET regions, as returned by \code{\link{chia.vertex.subset}}.
#' @param conditions a string with the condition to respect
#' @param proportion Should the number of occurences be converted to proportions?
#' @return A named vector (table) with the number of occurences (or proportions) of each category.
#' @importFrom GenomicRanges mcols
#' @export
subset.counts <- function(chia.subset, all.conditions, proportion = TRUE) {
  # Convert data into data frame
  variable.data <- mcols(chia.subset$Regions)
  # Count the number of nodes respecting the conditions
  count <- nrow(subset(variable.data, subset = eval(parse(text = all.conditions))))
  # If proportion are needed, divide by the total number of nodes
  if (proportion){
    count <- count / nrow(variable.data)
  }
  return(count)
}

#' Finds which transcription factors are present
#'
#' @param chia.subset A list containing a graph and ChIA-PET regions, as returned by \code{\link{chia.vertex.subset}}.
#' @param number Shoul the number of overlap be counted?
#' @return A named vector with the presence or absence of each transcription factor.
#' @importFrom GenomicRanges mcols
#' @export
TF.presence <- function(chia.subset, number = FALSE) {
  # Convert data into data frame
  variable.data <- mcols(chia.subset$Regions)
  # Extract the columns woth TF overlaps
  columns <- grep("TF.overlap.", colnames(variable.data))
  # Extract the names of these columns
  names <- sub("TF.overlap.", "", colnames(variable.data)[columns])
  # For each factor, counts the number of nodes overlapping
  presence <- vector()
  for (col in columns) {
    presence <- c(presence, sum(variable.data[,col] > 0))
  }
  # If only the presence is to consider, attribute the value "1" if the TF is present, "0" if absent
  if (!number) {
    presence <- ifelse(presence > 0, 1, 0)
  }
  # Add names
  names(presence) <- names
  return(presence)
}

calculate.tf.presence <- function(chia.obj, proportion=TRUE) {
    # Get the count of regions where the TF is presence, and rename the output vector
    # to remove the TF.overlap. prefix.
    results = apply(as.matrix(get.tf(chia.obj)) > 0, 2, sum)
    names(results) = gsub("TF.overlap.", "", names(results))
    
    # If we're calculating proportions, divide by the total number of regions.
    if(proportion) {
        results = results / sum(results)
    }
    
    # Return the results.
    return(results)
}

#' Finds which transcription factors are present
#'
#' @param chia.subset A list containing a graph and ChIA-PET regions, as returned by \code{\link{chia.vertex.subset}}.
#' @param variable.name The name of the boolean variable whose proportion should be calculated.
#' @return The count/proportion of nodes whose attribute named "variable.name" is TRUE.
#' @importFrom GenomicRanges mcols
#' @export
boolean.count <- function(chia.obj, variable.name, proportion=FALSE) {
  # Convert data into data frame
  variable.data <- as.logical(mcols(chia.obj$Regions)[[variable.name]])
  
  # Count the number of nodes respecting the conditions
  results <- sum(variable.data)
  
  # If proportion are needed, divide by the total number of nodes
  if (proportion){
    results <- results / length(variable.data)
  }
  return(results)
}

#' Builds a function object for extracting metrics using a variable name.
#'
#' @param method The method to be used by the resulting function object. Valid values could be
#'   levels.count, count.cut, boolean.count, etc.
#' @param variable.name The name of the variable whose proportion/counts should be calculated.
#' @param proportion If true, proportions are returned instead of counts.
#' @return A function to calculate the desired metric.
#' @export
functor.constructor <- function(method, variable.name, proportion=FALSE) {
    return(function(x) {
        method(x, variable.name=variable.name, proportion=proportion)
    })
}


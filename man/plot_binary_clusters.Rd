% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_binary_clusters.R
\name{plot_binary_clusters}
\alias{plot_binary_clusters}
\title{Plot binary clusters data (tile).}
\usage{
plot_binary_clusters(x, ...)
}
\arguments{
\item{x}{An \code{mtree} tree.}

\item{...}{Extra parameters, not used.}

\item{patient}{A patient id.}
}
\value{
A \code{ggplot} plot.
}
\description{
This function creates a \code{ggplot}-style
heatmap of the input CCF cluster of each clone in the data.
The heatmap is annotated for the drivers status of each
clone (with/ without driver). The CCF values are used to
colour the heatmap (`NA` values are in light gray).
}
\examples{
data(mtree_input)

x = mtrees(
mtree_input$binary_clusters, 
mtree_input$drivers,
mtree_input$samples,
mtree_input$patient,
mtree_input$sspace.cutoff,
mtree_input$n.sampling,
mtree_input$store.max
)

plot_binary_clusters(x[[1]])
}

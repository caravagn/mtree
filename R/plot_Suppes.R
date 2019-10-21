#' Plot Suppes' conditions for this tree
#' 
#' @description This function creates a joint plot assembled with
#' \code{cowplot} integrating two tile plots, one reporting which
#' conditions by Suppes' are satisfied, and the other reporting the
#' numbers used in the evaluation.
#'
#' @param x An \code{mtree} tree.
#' @param patient A patient id.
#' @param ... Extra parameters, not used.
#'
#' @return A \code{ggplot} plot.
#' 
#' @import reshape2
#' @import cowplot
#' @export
#'
#' @examples
#' data(mtree_input)
#' 
#' x = mtrees(
#' mtree_input$binary_clusters, 
#' mtree_input$drivers,
#' mtree_input$samples,
#' mtree_input$patient,
#' mtree_input$sspace.cutoff,
#' mtree_input$n.sampling,
#' mtree_input$store.max
#' )
#' 
#' plot_Suppes(x[[1]])
plot_Suppes = function(x, ...)
{
  p = x$patient
  sm = x$samples
  cl = x$CCF

  # Values 
  values = x$Suppes$Suppes %>%
    mutate(
      edge = paste0(from, ' \u2192 ', to)
    ) %>%
    select(-from, -to)
  
  colnames(values)[1:4] = c(
    'p(i)',
    'p(j)',
    'p(i, j)',
    'p(i) p(j)'
  )
  
  values = values %>% arrange(desc(Suppes))
  
  bool = reshape2::melt(
    values %>% select(edge, TP, PR, Suppes),
    id = 'edge'
  )
  
  bplot = ggplot(
    bool,
    aes(
      y = edge,
      x = variable,
      z = value,
      fill = value)
  ) +
    geom_tile(aes(width = .8, height = .8), size = 1) +
    my_ggplot_theme() +
    scale_fill_manual(values = c(`TRUE` = 'forestgreen', `FALSE` = 'gray')) +
    labs(
      title = paste("Suppes' conditions for", p),
      y = 'Tree edges',
      x = 'Conditions',
      subtitle = paste("Clonal cluster = ", cl$cluster[cl$is.clonal], 'evaluation as', x$Suppes$evaluation)
    ) +
    guides(
      fill = guide_legend("Binary data"),
      color = guide_legend("With driver")
    ) +
    scale_y_discrete(limits = values %>% pull(edge) %>% rev)
  
  
  nums = reshape2::melt(
    values %>% select(edge, starts_with('p(')),
    id = 'edge'
  )
  
  nplot = ggplot(
    nums,
    aes(
      y = edge,
      x = variable,
      z = value,
      fill = value)
  ) +
    geom_tile(aes(width = .8, height = .8), size = 1) +
    geom_text(aes(label = value), color = 'white') +
    my_ggplot_theme() +
    scale_fill_gradient(low = 'steelblue', high = 'darkred') +
    labs(
      title = paste("Marginal and joint probabilities for", p),
      y = 'Tree edges',
      x = 'Conditions',
      subtitle = paste("Clonal cluster = ", cl$cluster[cl$is.clonal], 'evaluation as', x$Suppes$evaluation)
    ) +
    guides(
      fill = guide_colorbar("Probability", barwidth = 6)
    ) +
    scale_y_discrete(limits = values %>% pull(edge) %>% rev)
  
  figure = cowplot::plot_grid(
    bplot, 
    nplot,
    ncol = 2,
    nrow = 1,
    align = 'h'
  )
  
  figure
}


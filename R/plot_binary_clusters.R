#' Plot binary clusters data (tile).
#' 
#' @description This function creates a \code{ggplot}-style
#' heatmap of the input CCF cluster of each clone in the data.
#' The heatmap is annotated for the drivers status of each
#' clone (with/ without driver). The CCF values are used to
#' colour the heatmap (`NA` values are in light gray).
#'
#' @param x An \code{mtree} tree.
#' @param patient A patient id.
#' @param ... Extra parameters, not used.
#'
#' @return A \code{ggplot} plot.
#' 
#' @import reshape2
#' 
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
#' plot_binary_clusters(x[[1]])
plot_binary_clusters = function(x, ...)
{
  p = x$patient
  sm = x$samples
  cl = x$CCF

  # Values CCF
  cl_tab = cl %>%
    select(!!sm, cluster) %>%
    reshape2::melt(id = 'cluster') %>%
    rename(
      region = variable,
      CCF = value) %>%
    as_tibble() %>%
    mutate(CCF = ifelse(CCF == 0, NA, CCF))
  
  # Annotations
  cl_tab_anno = cl %>%
     select(cluster, nMuts, is.driver, is.clonal)
  
  # Cluster ordering by sum of CCF
  cluster_ordering = cl_tab %>%
    group_by(cluster) %>%
    summarise(tot = sum(CCF, na.rm = TRUE)) %>%
    arrange(tot) %>%
    pull(cluster)
  
  # Combined
  cl_tab = cl_tab %>%
    left_join(cl_tab_anno %>% select(cluster, is.driver), by = 'cluster') %>%
    mutate(is.driver = ifelse(is.na(CCF), FALSE, is.driver))
  
  # Factor to sort
  cl_tab$cluster = factor(cl_tab$cluster, levels = cluster_ordering)
  
  ggplot(
    cl_tab %>% mutate(CCF = ifelse(is.na(CCF), '0', CCF)),
    aes(
      x = region,
      y = cluster,
      z = CCF,
      fill = paste(CCF),
      color = is.driver)
  ) +
    geom_tile(aes(width = .8, height = .8), size = 1) +
    scale_fill_manual(
      values = c(`1` = 'steelblue', `0` = 'gainsboro')
    ) +
    scale_color_manual(
      values = c(`TRUE` = 'red', `FALSE` = NA)
    ) +
    theme(
      legend.position = 'bottom',
      legend.key.size = unit(3, 'mm')
    ) +
    labs(
      title = paste("Binary data for", p),
      y = 'Binary cluster',
      x = 'Region',
      subtitle = paste("Clonal cluster = ", cl$cluster[cl$is.clonal])
    ) +
    guides(
      fill = guide_legend("Binary data"),
      color = guide_legend("With driver")
      )
}


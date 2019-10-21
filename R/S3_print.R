#' Print a \code{mtree} tree.
#' 
#' @details 
#' 
#' Print a summary for a \code{mtree} object, which includes a
#' shell-frinedly layout and other information.
#'
#' @param x A \code{mtree} tree.
#' @param ... Extra S3 parameters
#'
#' @return Nothing
#'
#' @export
#'
#' @import crayon
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
#' print(x[[1]])
print.mtree <- function(x, ...)
{
  stopifnot(inherits(x, "mtree"))
  
  M = x$adj_mat
  tb = x$tb_adj_mat
  
  printPretty = function(node, indent, last)
  {
    cat(indent)
    if (last) {
      cat("\\-")
      indent = paste(indent, " ", sep = '')
    }
    else {
      cat("|-")
      indent = paste(indent, "| ", sep = '')
    }
    cat(node)
    
    A = tb %>%
      activate(nodes) %>%
      filter(cluster == !!node) %>%
      pull(attachment)
    
    if (!is.na(A))
      cat(paste(' [', A, ']', sep = ''))
    
    D = tb %>%
      activate(nodes) %>%
      filter(cluster == !!node) %>%
      pull(driver)
    
    if (!is.na(D))
      cat(sprintf(' :: %s', D))
    
    cat('\n')
    
    cl = children(M, node)
    
    for (c in cl)
      printPretty(c, indent, c == cl[length(cl)])
  }
  
  pio::pioHdr(paste0('mtree - ', x$annotation))
  
  cat('\n')
  print(x$CCF)
  
  pio::pioStr("Tree shape (drivers annotated)",
              "",
              prefix = '\n',
              suffix = '\n\n')
  
  printPretty(node = "GL",
              indent = "  ",
              last = TRUE)
  
  pio::pioStr("Information transfer",
              "",
              prefix = '\n',
              suffix = '\n\n')
  
  apply(x$transfer$drivers, 1, function(w)
    cat('  ', w[1], '--->', w[2], '\n'))
  
  pio::pioStr("Tree score", x$score, suffix = '\n', prefix = '\n')
}




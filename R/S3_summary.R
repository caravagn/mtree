



#' Summary of a \code{"mtree"} object.
#' @description 
#' 
#' Reports some summary statistics for a mutation tree.
#'
#' @param x An \code{mtree} tree.
#' @param ... Extra parameters
#'
#' @return None.
#' @import crayon
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
#' summary(x[[1]])
summary.mtree <- function(x, ...) {
  print.ctree(x, ...)
  
  pio::pioStr("Binary clusters:",
              "",
              prefix = '\n',
              suffix = '\n\n')
  print(x$CCF)
  
  pio::pioStr("Drivers:",
              "",
              prefix = '\n',
              suffix = '\n\n')  
  print(x$drivers)
  
  stat = stats.rev_phylo(x)
  
  pio::pioStr("Pigeonhole principle:",
              green(
                nrow(stat$CCF.pigeonhole) * ncol(stat$CCF.pigeonhole) - stat$violations['pp']
              ),
              red(stat$violations['pp']),
              prefix = '\n',
              suffix = '\n\n')  
  print.noquote(stat$CCF.pigeonhole)
  
  pio::pioStr("Goodness-of-fit:",
              stat$gofit,
              "",
              prefix = '\n',
              suffix = '\n\n')  
}




stats.rev_phylo = function(x)
{
  
  M = MatrixToDataFrame(x$adj_mat)
  M = M[M$from != 'GL', , drop = FALSE]
  
  if (nrow(M) == 0) {
    return(
      list(
        probs = NULL,
        violations = NULL,
        gofit = 1,
        Suppes = NULL,
        CCF.crossing.rule = NULL,
        CCF.pigeonhole = NULL
      )
    )
    
  }
  
  # CCF data for this patient
  CCF = x$CCF[, x$samples, drop = FALSE]
  
  ## Branching penalty is the Pigeonhole Principle
  Mmatrix = DataFrameToMatrix(M)
  
  branching.penalty = sapply(x$CCF %>% pull(cluster),
                             function(n) {
                               # print(n)
                               cl = children(Mmatrix, n)
                               if (length(cl) == 0)
                                 return(NULL)
                               CCF[n, , drop = FALSE] >= colSums(CCF[cl, , drop = FALSE])
                             })
  branching.penalty = branching.penalty[!unlist(lapply(branching.penalty, is.null))]
  branchp = Reduce(rbind, branching.penalty)
  rownames(branchp) = names(branching.penalty)
  
  # cr.viol = sum(as.integer(!direction.penalty))
  ph.viol = sum(as.integer(!branchp))
  
  # violations = c(tp = tp.viol, pr = pr.viol, cr = cr.viol, pp = ph.viol)
  violations = c(pp = ph.viol)
  
  # Some alternative Goodness-Of-Fit measures. In the end, we take the
  # exact definition of the penalty for phylogenetic trees.
  # gofit = 2 * (x$numNodes - 1) * x$numRegions + x$numNodes  * x$numRegions
  # gofit = 1 - sum(violations)/gofit
  # gofit = nrow(branchp) * x$numRegions
  # gofit = 1 - ph.viol/gofit
  capture.output({
    gofit = node.penalty.for.branching(list(M), CCF)
  })
  
  
  return(
    list(
      violations = violations,
      gofit = gofit,
      CCF.pigeonhole = branchp
    )
  )
  
}
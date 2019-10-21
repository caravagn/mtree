#' Construct a `mtree` mutation tree computing its structure.
#'
#' @description
#'
#' This constructor creates a list of objects of class `'mtree'`, after using a
#' sampling strategy to determine possible trees that fit the data. The strategy
#' to sample trees can be controlled, a maximum number of trees can be sampled
#' with a Monte Carlo procedure and the actual process can be exhausted if there
#' are less than a number of available trees to fit the data.
#'
#' Note that the parameters of this function includes the same parmeters of
#' function \code{\link{mtree}}, plus the parameters of the sampler. See
#'  \code{\link{mtree}} for an explanation of the parameters.
#'
#' @param binary_clusters Clusters of binary annotations in the data of
#' this patient. See the package vignette to see the format in which this should
#' be specified.
#' @param drivers A list of driver events that should be annotated to each one
#' of the input clusters contained in the `CCF_clusters` parameter. See the package
#' vignette to see the format in which this should be specified.
#' @param samples A vector of samples names (e.g., the biopsies sequenced for
#' this patient).
#' @param patient A string id that represent this patient.
#' @param M The adjacency matrix defined to connect all the nodes of this tree.
#' @param score A scalar score that can be associated to this tree.
#' @param annotation Any string annotation that one wants to add to this `ctree`.
#' This will be used by some of the plotting functions that display `ctree` objects.
#' @param sspace.cutoff If there are less than this number of tree available, all the
#' structures are examined in an exhaustive fashion. Otherwise, if there are more than
#' this, a Monte Carlo sampler is used.
#' @param n.sampling If a Monte Carlo sampler is used, \code{n.sampling} distinct
#' trees are sampled and scored.
#' @param store.max When a number of trees are generated, scored and ranked, a maximum
#' of \code{store.max} are returned to the user (these are selected following the
#' ranking).
#' @param evaluation How Suppes conditions should be evaluated (`>=` or `>`).
#'
#' @return An list of objects of class \code{"mtree"} that represent the trees that
#' can be fit to the data of this patient.
#'
#' @export
#'
#' @import tidyverse
#' @import tidygraph
#' @import crayon
#' @import clisymbols
#' @import entropy
#' @import matrixcalc
#' @import reshape2
#' @import clisymbols
#' @import easypar
#' @import entropy
#'
#' @examples
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
#' plot(x[[1]])
mtrees = function(binary_clusters,
                  drivers,
                  samples,
                  patient,
                  sspace.cutoff = 10000,
                  n.sampling = 5000,
                  store.max = 100,
                  evaluation = '>=')
{
  # TODO - check input formats
  
  pio::pioHdr(paste("mtree ~ generate mutation trees for", patient))
  
  pioStr(
    'Sampler : ',
    sspace.cutoff,
    '(cutoff), ',
    n.sampling,
    '(sampling), ',
    store.max,
    '(max store)',
    suffix = '\n'
  )
  
  pioStr(
    'Suppes\' conditions : ', evaluation,
    suffix = '\n'
  )
  
  # Sample structure for all trees
  structures = trees_sampler(binary_clusters,
                             drivers,
                             samples,
                             patient,
                             sspace.cutoff,
                             n.sampling,
                             store.max)
  
  TREES = structures[[1]]
  SCORES = structures[[2]]
  
  # Trees assembly
  
  pio::pioStr(
    " Trees with non-zero sscore",
    length(TREES),
    'storing',
    min(length(TREES), store.max),
    prefix = crayon::green(clisymbols::symbol$tick),
    suffix = '\n'
  )
  
  if (length(TREES) > store.max)
  {
    TREES = TREES[1:store.max]
    SCORES = SCORES[1:store.max]
  }
  
  LSCORES = as.data.frame(SCORES)
  LSCORES = split(LSCORES, f = LSCORES[, 1])
  
  LSCORES = lapply(LSCORES, function(w)
    w[sample(1:nrow(w), replace = FALSE), , drop = FALSE])
  
  # Shuffle indexes of equal-scoring fits to avoid a sistematic bias in the way I coded this
  permuted.indexes = as.integer(rev(unlist(lapply(LSCORES, rownames))))
  names(permuted.indexes) = NULL
  
  TREES = TREES[permuted.indexes]
  SCORES = SCORES[permuted.indexes]
  
  # print(TREES)
  # print(SCORES)
  
  #####################################################################
  
  easypar::run(
    FUN = function(i)
    {
      mtree(
        binary_clusters,
        drivers,
        samples,
        patient,
        M = TREES[[i]],
        score = SCORES[i],
        annotation = paste('mtree rank ', i, '/', length(TREES), ' for ', patient, sep =
                             '')
      )
    },
    PARAMS = lapply(seq_along(TREES), list),
    parallel = FALSE
  )
  
  
}

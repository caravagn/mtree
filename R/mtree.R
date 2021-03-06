#' Construct a `mtree` mutation tree with known structure.
#'
#' @description 
#'
#' This constructor creates an object of class `'mtree'`, which represents a mutation tree. 
#' The tree is created from a set of binary clusters computed for a patient, here a cluster
#' is defined as the set of alterations (e.g., mutations) that are detected as present or
#' absent in the same set of sequenced biopsies.
#' 
#' To create a tree a list of drivers can be provided to be annotated to an input 
#' set of binary clusters. There are a minimum amount of information and formatting
#' fields that are required for tree construction to operate successfully. Please
#' refer to the package vignette and the provided input datasets for more instructions.
#' 
#' @seealso This function requires the input tree to be specified in the
#' format of an adjacency matrix; plese see function \code{\link{mtrees}} if you
#' need to create de novo also the adjacency matrices that fit your data.
#' 
#' @param binary_clusters Clusters of Cancer Cell Fractions available in the data of
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
#' @param evaluation How Suppes conditions should be evaluated (`>=` or `>`).
#'
#' @return An object of class \code{"mtree"} that represents this tree.
#' 
#' @export
#'
#' @import tidyverse
#' @import tidygraph
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
#' x = x[[1]]    
#'    
#'    
#' # Adj matrix inside of the objects, we remove the GL
#' # entry that is added as fake root by ctree
#' M = x$adj_mat
#' M = M[rownames(M) != 'GL', colnames(M) != 'GL']
#' 
#' print(M)
#' 
#' # Manual construction
#' y = mtree(
#' mtree_input$binary_clusters, 
#' mtree_input$drivers,
#' mtree_input$samples,
#' mtree_input$patient,
#' M,
#' score = 123456,
#' annotation = paste0("Some mutation tree")
#' )
#' 
#' # The same
#' print(x)
#' print(y)
mtree = 
  function(
    binary_clusters,
    drivers,
    samples,
    patient,
    M,
    score,
    annotation = paste0("Mutation tree for patient ", patient),
    evaluation = '>='
  )
  {
    # This function will create this output object
    obj <-
      structure(
        list(
          adj_mat = NA,
          # Adjacency matrix for the tree
          tb_adj_mat = NA,
          # Tidygraph object for this tree
          score = NA,
          # Tree score
          patient = NA,
          # Patient ID
          samples = NA,
          # Samples name
          drivers = NA,
          # Driver mapping to clusters
          CCF = NA,
          binary = NA,
          # CCF (aggregated per cluster) + binary (same)
          transfer = NA,
          # Information Transfer
          annotation = NA,
          # Some custom annotation
          tree_type = "Mutation tree",  # Mutation tree 
          evaluation = evaluation
        ),
        class = "mtree",
        call = match.call()
      )
    
    # =-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-
    # The information that we need for each tree are
    # binary clusters, data and driver events (mapped to clusters)
    # =-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-
    # For historical reasons and consistency with other packages (e.g., ctree from evoverse)
    # the binary_clusters are stored inside a CCF field.
    obj$CCF = obj$binary_data = binary_clusters %>% mutate(cluster = paste(cluster))

    obj$samples = samples
    obj$drivers = drivers
    
    # =-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-
    # We begin to create a representation of the tree
    # =-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-
    
    # Check that is must have drivers
    has_drivers = nrow(drivers) > 0
    if(!has_drivers) stop("No drivers are annotated in this clone tree, aborting.")
    
    # Handle the special case of empty matrix
    is_monoclonal = (binary_clusters$cluster %>% unique %>% length) == 1
    
    if(is_monoclonal == 1)
    {
      message('\nThis tree has 1 node, creating a monoclonal model disregarding the input matrix.')
      
      M = matrix(0, ncol = 1, nrow = 1)
      colnames(M) = rownames(M) = binary_clusters$cluster %>% unique
    }
    
    # Input can should be an adjacency matrix (which works also for monoclonal tumours)
    if (!inherits(M, 'matrix'))
      stop("Input `M` should be an adjacency matrix, aborting.")
    
    adj_mat =  M
    df_mat = MatrixToDataFrame(M)
    
    # We check for that to be a tree - empty is OK in this case if monoclonal
    if (!is_tree(adj_mat, allow.empty = nrow(obj$CCF) == 1))
      stop("The input adjacency matrix is not a valid tree, aborting.")
    
    # Add GL node, beware of the special case of empty adj_mat (might happen for a monoclonal tumour)
    M_root = ifelse(sum(adj_mat) == 0, rownames(adj_mat), root(adj_mat))
    df_mat = rbind(df_mat,
                   data.frame(
                     from = 'GL',
                     to = M_root,
                     stringsAsFactors = FALSE
                   ))
    
    # Update the adj matrix with GL, which can now go into obj
    obj$adj_mat = DataFrameToMatrix(df_mat)
    
    # =-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-
    # Create a tidygraph object for this tree
    # =-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-
    tb_adj_mat = as_tbl_graph(df_mat) %>%
      activate(nodes) %>%
      rename(cluster = name)
    
    # We can add information specific for this tree to the tidygraph
    tb_adj_mat = tb_adj_mat %>%
      left_join(obj$CCF, by = 'cluster')
    
    # Sample attachment for input data
    attachment = obj$CCF %>%
      select(cluster, obj$samples) %>%
      reshape2::melt(id = 'cluster') %>%
      mutate(value = ifelse(value > 0, 1, 0)) %>%
      group_by(variable) %>%
      filter(sum(value) == 1) %>%
      ungroup() %>%
      filter(value == 1) %>%
      select(-value) %>%
      rename(attachment = variable) %>%
      mutate(attachment = paste(attachment)) %>%
      group_by(cluster) %>%
      summarise(attachment = paste(attachment, collapse = ', '))
    
    tb_adj_mat = tb_adj_mat %>%
      left_join(attachment, by = 'cluster')
    
    # Drivers per node
    tb_adj_mat = tb_adj_mat %>%
      left_join(obj$drivers %>%
                  group_by(cluster) %>%
                  summarise(driver = paste(variantID, collapse = ', ')),
                by = 'cluster')
    
    # Store it in obj
    obj$tb_adj_mat = tb_adj_mat
    
    # Compute the information transfer
    obj$transfer = information_transfer(obj)
    
    # Extras
    obj$score = score
    obj$patient = patient
    
    obj$annotation = annotation
    
    obj$tree_type = "Mutation tree from binary data."
    
    # Add Suppes to this guy
    obj$Suppes = Suppes_poset(obj$binary_data, obj$samples, evaluation = evaluation)
    
    return(obj)
  }




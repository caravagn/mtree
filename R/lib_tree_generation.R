all.possible.trees = function(
  G,
  W,
  sspace.cutoff = 10000,
  n.sampling = 1000
  )
{
  M = DataFrameToMatrix(G)
  r = root(M)
  parents = sapply(colnames(M), pi, model = M)

  parents = parents[!unlist(lapply(parents, is.null))]
  # cat("* Nodes: ", colnames(M), '\n')
  # cat("* Root: ", r, '\n')
  # cat("* Parent set \n\t")
  #
  # nothing = lapply(parents, function(x) cat(paste('{', paste(x, collapse =', ', sep =''), '}', collapse = '')))

  # singletons -- template
  singletons = parents[unlist(lapply(parents, function(x) length(x) == 1))]
  nsingletons = length(singletons)


  # cat('\n* Singletons:', names(singletons))
  # cat(' [ n =', nsingletons, ']\n')

  sglt = expand.grid(singletons, stringsAsFactors = FALSE)
  if(ncol(sglt) == 0) {
    sglt = NULL
  }

  alternatives = parents[unlist(lapply(parents, function(x) length(x) > 1))]
  nalternatives = length(alternatives)

  altn = NULL
  if(nalternatives > 0)
  {
    combalternatives = prod(unlist(lapply(alternatives, length)))
    # cat('* Alternatives:', names(alternatives), '-- num.', combalternatives, '\n')
    
    # print(combalternatives)
    # print(sspace.cutoff)
    
    ex_search = ifelse(
      combalternatives < sspace.cutoff,
      "exahustive",
      paste0('Montecarlo for ', n.sampling, 'distinct trees')
    )
    
    # print(ex_search)
    
    pio::pioStr(
      " Structures", combalternatives,  '- search is', ex_search,
      prefix = crayon::green(clisymbols::symbol$tick),
      suffix = '\n')
    
    

    if(combalternatives > sspace.cutoff)
    {
      return(
        weighted.sampling(
          DataFrameToMatrix(G),
          W,
          n.sampling
        )
      )
    }
    altn = expand.grid(alternatives, stringsAsFactors = FALSE)
  }
  else cat(red('There are no alternatives!\n'))

  # all combinations
  if(is.null(altn) && is.null(sglt)) stop('Error -- no trees?')

  if(is.null(sglt) && !is.null(altn)) comb = altn
  if(!is.null(sglt) && is.null(altn)) comb = sglt
  if(!is.null(sglt) && !is.null(altn)) comb =  cbind(altn, sglt)


  
  pb =  dplyr::progress_estimated(n = nrow(comb), min_time = 2)
  progress_bar = getOption('ctree.progressBar', default = TRUE)

  models = NULL
  for(i in 1:nrow(comb))
  {
    if (progress_bar)
      pb$tick()$print()
    
    
    tree = data.frame(from = unlist(comb[i, ]), to = colnames(comb), stringsAsFactors = FALSE)
    test.tree = DataFrameToMatrix(tree)

    if(length(root(test.tree)) > 1 ) {
      # cat('Solution', DataFrameToEdges(tree), 'has multiple roots, removed\n')
      next;
    }

    if(!igraph::is_dag(igraph::graph_from_adjacency_matrix(test.tree))){
      # cat('Solution', DataFrameToEdges(tree), 'has loops, removed\n')
      next;
    }

    # # revert mapping
    # tree = DataFrameToMatrix(tree)
    # tree = reverse.mapping(tree, clusterIdsMapping)
    # tree = MatrixToDataFrame(tree)

    models = append(models, list(tree))
  }


  # cat(length(models), ' ')
  # cat('\r')

  return(models)
}

rankTrees = function(TREES, MI.table, structural.score)
{
  pb = dplyr::progress_estimated(n = length(TREES), min_time = 2)
  progress_bar = getOption('ctree.progressBar', default = TRUE)

  MI.TREES = NULL
  for(i in 1:length(TREES))
  {
    if (progress_bar)
      pb$tick()$print()
    
        M = DataFrameToMatrix(TREES[[i]])
        M = M[colnames(MI.table), colnames(MI.table)]

        M.entries = MI.table[which(M > 0, arr.ind = TRUE)]

        val = NA
        if(all(is.null(structural.score))) val = prod(M.entries)
        else val = prod(M.entries) * structural.score[i]

        # print(paste(prod(M.entries), sum(log(M.entries))))

        if(any(M.entries == 0))
        {
          n = sum(M.entries == 0)
          M.entries[M.entries == 0] = 1e-9

          # cat("\nMI correction for", n, "entries equal 0; set equal to 1e-9.", prod(M.entries))

          warning("Used MI correction for", n, "entries equal 0; set equal to 1e-9.")
        }
#
#         print(TREES[[i]])
#         print(M)
#         print(M.entries)
#         readline("")

        MI.TREES = c(MI.TREES, val)
  }

  # print(head(sort(table(MI.TREES), decreasing = T)))

  o = order(MI.TREES, decreasing = TRUE)
  MI.TREES = MI.TREES[o]
  TREES = TREES[o]
  structural.score = structural.score[o]

  return(list(TREES = TREES, SCORES = MI.TREES, PENALTIES = structural.score))
}

computeMI.table = function(binary.data, MI.Bayesian.prior = 0, add.control = FALSE)
{
  if(add.control) binary.data = rbind(binary.data, wt = 0)

  # • a=0:maximum likelihood estimator (see entropy.empirical)
  # • a=1/2:Jeffreys’ prior; Krichevsky-Trovimov (1991) entropy estimator
  # • a=1:Laplace’s prior
  # • a=1/length(y):Schurmann-Grassberger (1996) entropy estimator
  # • a=sqrt(sum(y))/length(y):minimax prior

  MI.table = matrix(
    apply(
      expand.grid(colnames(binary.data), colnames(binary.data)),
      1,
      function(x) {
        # Counting process.. with a Bayesian prior
        i = x[1]
        j = x[2]
        jo11 = (binary.data[, i] %*% binary.data[, j])/nrow(binary.data)
        jo10 = (binary.data[, i] %*% (1-binary.data[, j]))/nrow(binary.data)
        jo01 = ((1-binary.data[, i]) %*% binary.data[, j])/nrow(binary.data)
        jo00 = 1 - (jo10 + jo01 + jo11)
        entropy::mi.Dirichlet(matrix(c(jo11, jo10, jo01, jo00), nrow = 2), a = MI.Bayesian.prior)
      }),
    byrow = TRUE, ncol = ncol(binary.data))
  colnames(MI.table) = rownames(MI.table) = colnames(binary.data)

  return(MI.table)
}


weightMI.byMultinomial = function(MI.table, W)
{
  Coeff = matrix(0, ncol = ncol(MI.table), nrow = nrow(MI.table))
  rownames(Coeff) = rownames(MI.table)
  colnames(Coeff) = colnames(MI.table)

  for(j in names(W))
    Coeff[ names(W[[j]]) , j] = unlist(W[[j]])

  return(matrixcalc::hadamard.prod(MI.table, Coeff))
}

weighted.sampling = function(G, W, n)
{
  S = S.hashcodes = NULL

  sampleT = function()
  {
    r = root(G)

    E = setdiff(colnames(G), r)
    E = sample(E, length(E))

    tree = NULL
    repeat {
      pi = sapply(E, function(node){
        parents = W[[node]]

        draw = sample(names(parents), prob = unlist(parents), size = 1)
      })

      tree = data.frame(from = unlist(pi), to = names(pi), stringsAsFactors = FALSE)
      R = c(r, reach(tree, r))

      if(all(colnames(G) %in% R)) break
    }

    return(tree)
  }

  c = 0
  repeat{
    Tree = sampleT()
    hash = paste(sort(DataFrameToEdges(Tree)), collapse = ':')

    if(!(hash %in% S.hashcodes))
    {
      S = append(S, list(Tree))
      S.hashcodes = c(S.hashcodes, hash)
      c = c + 1
    }
    if(c == n) break;
  }
  
  return(S)

}

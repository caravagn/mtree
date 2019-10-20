trees_sampler = function(binary_clusters,
                         drivers,
                         samples,
                         patient,
                         sspace.cutoff = 10000,
                         n.sampling = 5000,
                         store.max = 100)
{
  TREES = SCORES = NULL
  
  # We will need to use the CCF clusters for this patient
  clusters = binary_clusters
  nclusters = nrow(clusters)
  
  clonal.cluster = clusters %>% filter(is.clonal) %>% pull(cluster)
  
  df_clusters = data.frame(clusters[, samples],
                           row.names = clusters$cluster)
  
  print(clusters)
  cat('\n')
  
  if (nclusters == 1)
  {
    cat(red('Sampler: this model has 1 node, it has trivial models.\n'))
    
    M = matrix(0, ncol = 1, nrow = 1)
    colnames(M) = rownames(M) = rownames(clusters)
    
    TREES = append(TREES, list(M))
    SCORES = c(SCORES, 1)
  }
  else
  {
    
    # ################## Generate Suppes poset
    SUPPES = Suppes_poset(binary_clusters, samples)
    
    CONSENSUS.TREE = SUPPES$Suppes %>% 
      filter(Suppes) %>%
      select(from, to)
    
    WEIGHTS.CONSENSUS.TREE = CONSENSUS.TREE %>% group_split(to)
    names(WEIGHTS.CONSENSUS.TREE) = sapply(WEIGHTS.CONSENSUS.TREE, function(x) x$to[1]) %>% paste
    
    WEIGHTS.CONSENSUS.TREE = lapply(WEIGHTS.CONSENSUS.TREE, function(w)
    {
      # 1/K uniform
      pio:::nmfy(
        w$from, 
        rep(
          1/length(w$from), length(w$from)
        )
      )
    })
    
    ################## Build all possible mutation trees

    # print(data.frame(CONSENSUS.TREE, stringsAsFactors = FALSE))
    # print(WEIGHTS.CONSENSUS.TREE)
    
    # # Sampling is carried out if there are more than 'sspace.cutoff' trees, in that case we
    # # sample 'n.sampling' possible trees. Otherwise all possible trees are generated.
    TREES = all.possible.trees(
      G = data.frame(CONSENSUS.TREE, stringsAsFactors = FALSE),
      W = WEIGHTS.CONSENSUS.TREE,
      sspace.cutoff,
      n.sampling
    )
    
    # print(TREES)
  
    # ################## Ranking trees. A tree is good according to the following factors:
    # # 1) the MI among the variables x and y, if they are connected by an edge x --> y 
    
    # 1) MI from binarized data -- different options, with a control sample which avoids 0log(0)
    # • a=0:maximum likelihood estimator (see entropy.empirical)
    # • a=1/2:Jeffreys’ prior; Krichevsky-Trovimov (1991) entropy estimator
    # • a=1:Laplace’s prior
    # • a=1/length(y):Schurmann-Grassberger (1996) entropy estimator
    # • a=sqrt(sum(y))/length(y):minimax prior
    binary.data = binary_clusters %>% select(samples) %>% data.frame
    
    # print(binary.data)
    rownames(binary.data) = binary_clusters$cluster
    
    MI.table = computeMI.table(binary.data %>% t,
                               MI.Bayesian.prior = 0,
                               add.control = TRUE)

    RANKED = rankTrees(TREES, MI.table, structural.score = NULL)
    TREES = RANKED$TREES
    SCORES = RANKED$SCORES
    
    TREES = TREES[SCORES > 0]
    SCORES = SCORES[SCORES > 0]
    
    TREES = lapply(TREES, DataFrameToMatrix)
  }
  
  return(list(adj_mat = TREES, scores = SCORES))
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



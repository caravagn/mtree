Suppes_poset = function(x, regions, evaluation = '>=')
{
  samples = x[, regions, drop = FALSE] %>% data.frame
  rownames(samples) = x$cluster
  
  variables = x$cluster
  
  # Marginal p(i)
  marginals = rowSums(samples)/length(regions)
  names(marginals) = variables
  
  # Joint p(i,j)
  joints = outer(
    1:nrow(samples), 1:nrow(samples),
    FUN = Vectorize( function(i,j) sum(samples[i, ] * samples[j,])/length(regions) ))
  colnames(joints) = rownames(joints) = rownames(samples)
  
  joints = reshape2::melt(joints)
  colnames(joints) = c('from', 'to', 'joint')
  
  joints$from = paste(joints$from)
  joints$to = paste(joints$to)
  
  # Table with Suppes inequalities - pairwise across all possible variables annotated
  pwise = expand.grid(from = variables, to = variables, stringsAsFactors = FALSE) %>%
    filter(from != to) %>% 
    as_tibble()
  
  # Add marginals
  pwise = pwise %>%
    mutate(
      m_f = marginals[paste(from)],
      m_t = marginals[paste(to)]
    )
  
  # print(pwise)
  # print(joints)
  
  # Add joints
  pwise = pwise %>% 
    left_join(joints, by = c('from', 'to')) %>%
    mutate(
      m_f_t = m_f * m_t
    )
  
  if(evaluation == '>=')
  {
    pwise = pwise %>%
      mutate(TP = m_f >= m_t,
             PR = joint >= m_f_t)
  }
  
  if(evaluation == '>')
  {
    pwise = pwise %>%
      mutate(TP = m_f > m_t,
             PR = joint > m_f_t)
  }
  
  pwise = pwise %>%
    mutate(Suppes = TP & PR) 
  
  incoming = pwise %>% filter(Suppes) %>% pull(to)
  multi_root = length(setdiff(variables, incoming)) > 1
  
  if(multi_root) {
    warning("Suppes conditions with ", evaluation, " leave >1 disconnected components: ", 
            paste(setdiff(variables, incoming), collapse = ', '),
            ". Will use just Temporal Priority.")
    
    pwise = pwise %>% mutate(Suppes = TP)
  }
  
  list(Suppes = pwise, evaluation = evaluation, TP_hack = multi_root)
}

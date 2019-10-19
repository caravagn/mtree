Suppes_poset = function(x, regions, evaluation = '>=')
{
  samples = x[, regions, drop = FALSE]
  variables = x$cluster
  
  # Marginal p(i)
  marginals = rowSums(samples)/length(regions)
  names(marginals) = variables
  
  # Inequalities between marginals: p(i) >/>= p(j)
  bmarginals = lapply(names(marginals), function(w) {
    if(evaluation == '>=') marginals >= marginals[as.character(w)]
    else marginals > marginals[as.character(w)]
    })
  
  # inequalities TRUE
  bmarginals = lapply(bmarginals, function(w) w[which(w)])
  names(bmarginals) = names(marginals)
  
  # Joint probabilities
  joints = outer(
    1:nrow(samples), 1:nrow(samples),
    FUN = Vectorize( function(i,j) sum(samples[i, ] * samples[j,])/length(regions) ))
  colnames(joints) = rownames(joints) = rownames(samples)
  
  # Inequalities between joints: p(i,j) >/>= p(i)p(j)
  bmarginals = lapply(names(bmarginals), function(b){
    to = b
    bm = bmarginals[[b]]
    bool = sapply(names(bm), function(from)
    {
      if(evaluation == '>=') 
        joints[as.character(to), as.character(from)] >= marginals[as.character(to)] * marginals[as.character(from)]
      else
        joints[as.character(to), as.character(from)] > marginals[as.character(to)] * marginals[as.character(from)]
    })
    names(bool) = names(bm)
    bool
  })
  names(bmarginals) = names(marginals)
  for(b in names(bmarginals)) bmarginals[[b]][b] = FALSE
  
  # inequalities TRUE
  bmarginals = lapply(bmarginals, function(w) w[which(w)])
  bmarginals = bmarginals[sapply(bmarginals, function(w) length(w) > 0)]
  
  return(bmarginals)
}

################### Generate Suppes poset
POSET = bmarginals

# we transform it in the input for the tree sampler
POSET = lapply(POSET, function(w)
{
  n = names(w)
  w = as.numeric(w)/length(w)
  names(w) = n
  w
})

POSET.EDGES = lapply(POSET, names)
POSET.EDGES = lapply(names(POSET.EDGES), function(w){
  expand.grid(from = POSET.EDGES[[w]], to = w, stringsAsFactors = FALSE)
})
POSET.EDGES = Reduce(rbind, POSET.EDGES)



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

colnames(pwise) = c(
  'from',
  'to',
  'p(from)',
  'p(to)',
  'p(from, to)',
  'p(from) * p(to)',
  paste("TP", evaluation),
  paste("PR", evaluation),
  paste0("Suppes", evaluation, ifelse(multi_root, " ~ TP", ''))
)

data('mtree_input')
require(tidyverse)

mtree_input$binary_clusters = mtree_input$CCF_clusters %>%
  mutate(cluster = paste(cluster), nMuts = 4) 

mtree_input$drivers = mtree_input$drivers %>%
  mutate(cluster = paste(cluster), nMuts = 4) 


x = mtree::mtrees(
  binary_clusters = mtree_input$binary_clusters,
  drivers = mtree_input$drivers,
  samples = mtree_input$samples,
  patient = mtree_input$patient,
  sspace.cutoff = mtree_input$sspace.cutoff,
  n.sampling = mtree_input$n.sampling,
  store.max =  mtree_input$store.max
)

x = mtree::mtree(
  binary_clusters = mtree_input$binary_clusters,
  drivers = mtree_input$drivers,
  samples = mtree_input$samples,
  patient = mtree_input$patient,
  M = mtree:::edgesToMatrix(c('4~3', '1~2', '1~4')),
  score = 1
)


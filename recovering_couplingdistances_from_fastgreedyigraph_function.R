
# 16-dic-19 : Set de funciones para recuperar la distancia de acople a partir 
# de community detection using modularity with cluster_fast_greedy of igraph.
# developed in testing_hierarchical_clusteringV12.R

# # # # # # # # # # # # # # # # # GETTING THE NODES OF CLUSTERS # # # # # # # # # # # # #  # # # # # # # # #
# Function that creates a list of nodes of every cluster merge in each iteration
# of the greedy algorithm from cluster_fast_greedy of igraph.
# Inputs:
# 1. merge :  es la lista que viene de cluster_fast_greedy: con dos columnas
#             que indican los nodos que se van fusionando.
# Output: Lista con los nodos de cada cluster de cada iteracion.
#         El numero de clusters va a ser de N-1, siendo N el numero de nodos de la red.
# Output: dataframe merge que contiene V1 / V2 / cluster 
get_each_node_of_each_cluster <- function(merge) {
  Nn <- nrow(merge)
  merge <- as.data.frame(fc$merges)
  max_number_of_nodes = Nn + 1
  # Ahora que le puse los nombres vamos a identificar los nodos reales
  # de los cluster formados
  id <- which(merge$V1 <= max_number_of_nodes)
  merge[id,1] <- -1*merge[id,1] 
  id <- which(merge$V2 <= max_number_of_nodes)
  merge[id,2] <- -1*merge[id,2] 
  #merge$cluster <- seq(from= max_number_of_nodes+1, to= 2*max_number_of_nodes-1)
  merge$cluster <- seq(from= 1, to= max_number_of_nodes-1)
  id <- which(merge$V1 > max_number_of_nodes)
  merge[id,1] <- merge[id,1] - max_number_of_nodes
  id <- which(merge$V2 > max_number_of_nodes)
  merge[id,2] <- merge[id,2] - max_number_of_nodes
  
  # ahora tenemos que crear una lista para cada nombre 
  # de los cluster con los nodos que la componen.
  # esta seria la salida de la funcion
  nodes_of_clusters <- vector(mode = "list", length = Nn-1)
  for ( i in 1:nrow(merge)) {
    nodos <- c(merge$V1[i], merge$V2[i])
    temporal <- nodos
    test <- sum(1*nodos > 0)
    if (test != 0) {
      testeo = 1
      while (testeo > 0 ) {
        idx <- which( temporal > 0 )
        n <- length(idx)
        temporal_original <- temporal
        for (j in 1:n) {
          clu <- temporal_original[idx[j]]
          find_id <- which(merge$cluster == clu)
          temporal <- temporal[- which(temporal == clu) ]
          temporal <- c(temporal, c(merge$V1[find_id], merge$V2[find_id]))
        }
        testeo <- sum(1*temporal > 0) 
      }
      nodos <- temporal
      nodes_of_clusters[[i]] <- nodos
    } else {
      nodes_of_clusters[[i]] <- nodos
    }
  }
  return(list(nodes_of_clusters = nodes_of_clusters, merge = merge))
}
# Ejemplo:
# g <- graph_from_adjacency_matrix(D, weighted=TRUE, mode="undirected", diag=FALSE)
# fc <- cluster_fast_greedy(g)
# is_hierarchical(fc) # TRUE
# plot_dendrogram(fc)
# result <- get_each_node_of_each_cluster(merge = fc$merges)
# nodes_of_clusters <- result$nodes_of_clusters
# merge <- result$merge
# # # # # # # # # # # # # # # # # GETTING THE NODES OF CLUSTERS # # # # # # # # # # # # #  # # # # # # # # #







# # # # # # # # # # # # # # # # # GETTING THE NODES OF CLUSTERS # # # # # # # # # # # # #  # # # # # # # # #
# Con el resultado de la funcion get_heach_node_of_each_cluster, que nos da el merge y la lista
# de nodos de cada cluster, podemos ahora calcular la distancia de acople.
# source("acople_distance_sum_function.R")

# inputs: 
#   merge + nodes_of_clusters from get_heach_node_of_each_cluster
# output: merge con una 4ta columna que es la distancia de acople
report_merge_from_modularity <- function(merge, node_of_clusters) {
  dc <- vector(mode = "numeric", length = nrow(merge))
  for ( i in 1:nrow(merge)) {
    dc[i] <- acople_distance_sum(J, y = as.character(unlist(nodes_of_clusters[i]))  )  
  }
  merge$dc <- dc 
  return(merge)
}
# ejemplo
# merge <- report_merge_from_modularity(merge = result$merge, node_of_clusters = node_of_clusters)
# # # # # # # # # # # # # # # # # GETTING THE NODES OF CLUSTERS # # # # # # # # # # # # #  # # # # # # # # #



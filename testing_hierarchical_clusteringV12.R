# ******** coupling distance computation from modularity optimization  *************
# Intentamos calcular las distancias de acople dc de la 
#solucion de un clustering utilizando modularidad.

# Procedimiento:
# 1. encontrar la solucion jerarquica de comunidades con modularidad
# 2. para cada uno de los clusters que se van formando, computar la distancia de acople



# actual name: testing_hierarchical_clusteringV12.R

#Notas:
# 16-dic-19: creation 




# # # # # # # # # # # # # # # # # # LOADING DATA AND FUNCTIONS # # # # # # # # # # # # # # # # # # # # # # #
# Carga de datos, 
# igual que en clusteragrr_v1.R
rm(list = ls())
#load('new_inferring_parameters_environment270219.RData') # load H and J parameters ising inferred
source("functions_hclust.R") # funciones varias
source("hierarchical_clustering_v2_function.R")
source("hierarchical_clustering_v5_function.R")
source("hierarchical_clustering_greedy_function.R")
source("hierarchical_clustering_probabilistic_greedy_function.R")

source("create_coupling_function.R") # crea matriz J de acoples y D de distancias para simular
source("create_mst_function.R") # crea MST a partir de una matriz de acople.
source("acople_distance_sum_function.R") # suma las distantcias de acople dc
source("find_ady_dc_function.R") # encuentra los nodos adyacentes en un mst dado un nodo
source("get_num_nodes_function.R") # nos dice el numero de nodos que tiene el par de clusters a fusionar
source("get_num_nodes_of_ady_clusters_function.R") # nos dice informacion del numero de nodos que son adyacentes a un cluster
source("get_name_of_the_cluster_function.R") # dado el nombre de un nodo, nos dice a que cluster pertenece
source("get_nodes_of_the_cluster_function.R") # nos da los nodos involucrados en un cluster.
source("simulation_hc_function.R") # Genera una simulacion de hierarchical clustering con algoritmo single linkage y modificado.
source("find_min_distcp_function.R") #permite encontrar las distancias de acople entre un cluster y los restantes
source("pick_a_cluster_function.R") # necesario apra correr algoritmo greedy aprobabilistico

source("find_mst_barrier_function.R")
source("is_integer0_function.R")
source("entropies_functions.R") # nos interesa aqui get_energy function

library(igraph)
library(Matrix)
library(ggplot2)
# # # # # # # # # # # # # # # # # # LOADING DATA AND FUNCTIONS # # # # # # # # # # # # # # # # # # # # # # #


# # # # # # # # # # # # # # # # # GETTING THE NODES OF CLUSTERS # # # # # # # # # # # # #  # # # # # # # # #
# simula 
Jmean = 0
Nn=20
sj=1

set.seed(123)
ot <- create_coupling(Nn=Nn, media=Jmean, sj=sj, lw=-3, up=3)
J <- ot$J
D <- ot$D

# la matriz J debe ser convertida a distancias
D <- sqrt(3 - J)
# tengo que convertir la  matriz J como objeto igraph
g <- graph_from_adjacency_matrix(D, weighted=TRUE, mode="undirected", diag=FALSE)
fc <- cluster_fast_greedy(g)
is_hierarchical(fc) # TRUE
# https://users.dimi.uniud.it/~massimo.franceschet/R/communities.html
plot_dendrogram(fc)


# Function that creates a list of nodes of every cluster merge in each iteration
# of the greedy algorithm from cluster_fast_greedy of igraph.
# Inputs:
# 1. merge :  es la lista que viene de cluster_fast_greedy: con dos columnas
#             que indican los nodos que se van fusionando.
# Output: Lista con los nodos de cada cluster de cada iteracion.
#         El numero de clusters va a ser de N-1, siendo N el numero de nodos de la red.
# Output: dataframe merge que contiene V1 / V2 / cluster 
get_heach_node_of_each_cluster <- function(merge, Nn) {
  merge <- as.data.frame(fc$merges)
  max_number_of_nodes = nrow(merge) + 1
  # Ahora que le puse los nombres vamos a identificar los nodos reales
  # de los cluster formados
  id <- which(merge$V1 <= Nn)
  merge[id,1] <- -1*merge[id,1] 
  id <- which(merge$V2 <= Nn)
  merge[id,2] <- -1*merge[id,2] 
  merge$cluster <- seq(from= max_number_of_nodes+1, to= 2*max_number_of_nodes-1)
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
# result <- get_heach_node_of_each_cluster(merge = fc$merges, Nn=N)
# nodes_of_clusters <- result$nodes_of_clusters
# merge <- result$merge
# # # # # # # # # # # # # # # # # GETTING THE NODES OF CLUSTERS # # # # # # # # # # # # #  # # # # # # # # #




# # # # # # # # # # # # # # # # # GETTING THE NODES OF CLUSTERS # # # # # # # # # # # # #  # # # # # # # # #
# Con el resultado de la funcion get_heach_node_of_each_cluster, que nos da el merge y la lista
# de nodos de cada cluster, podemos ahora calcular la distancia de acople.
source("acople_distance_sum_function.R")

# inputs: 
#   merge + nodes_of_clusters from get_heach_node_of_each_cluster
# output: merge con una 4ta columna que es la distancia de acople
report_merge_from_modularity <- function(merge) {
  dc <- vector(mode = "numeric", length = nrow(merge))
  for ( i in 1:nrow(merge)) {
    dc[i] <- acople_distance_sum(J, y = as.character(unlist(nodes_of_clusters[i]))  )  
  }
  merge$dc <- dc 
  return(merge)
}
# ejemplo
# merge <- report_merge_from_modularity(merge = result$merge)





















































# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 17-nov-19
# Aqui creamos una funcion que para distintas distancias ultrematricas
# corte, determinamos los clusters y los elementos de cada uno. Luego
# a cada cluster sacamos la distancia de acople dc dado la distancia
# de corte h, y luego promediamos.
# Este promedio nos permite comparar con algoritmo SL y greedy con distintos
# niveles de probabilidad.

# funcion que nos permite definir todos los umbrales de corte
# Input: merge: que vienen del clustering jerarquico
#       number_of_intervals = numero de distancias de corte deseadas
# Output: vector con todos los umbrales
get_umbrales <- function(merge, number_of_intervals) {
  # vamos a generar dcut distancias de corte h entre el el minimo y 
  # maximo de distancia ultrametrica
  min_dc <- min(merge$dultr)
  max_dc <- max(merge$dultr)
  
  # b = intervalos deseados - 2 
  b <- number_of_intervals
  dcut <- seq(from=min_dc, to = max_dc, by = (max_dc - min_dc )/b)
  
  # es im portante eliminar el primer y ultima distancia ultrametrica para evitar problemas
  dcut <- dcut[c(-1, -length(dcut) ) ]
  
  return(dcut)
}
# ejemplo:
#dcut <- get_umbrales(merge = merge, number_of_intervals = 10)

# funcion que nos da la media de distancia de acople de los clusters
# encontrados a cierta distancia de corte
# Input: assg que es dataframe con dos columna: nombre de los nodos / solucion de cluster
get_mean_clusters <- function(assg) {
  # aqui calculamos la media de la distancia de acoples para una
  # distancia de corte determinada h.
  uniq_clustrs <- unique(assg$cl)
  inf <- matrix(NA, ncol=2, nrow=length(uniq_clustrs))
  for (cc in 1:length(uniq_clustrs) ) {
    sub <- subset(assg, cl==uniq_clustrs[cc])
    te <- as.numeric(sub$node_names)
    if ( length(te) > 1) {
      coupl_distance <- acople_distance_sum2(J, te)  
    } else {
      coupl_distance <- NA
    }
    
    inf[cc, ] <- c(uniq_clustrs[cc], coupl_distance)
  }
  return(mean(inf[,2], na.rm = TRUE) )
}
# ejemplo:
#lamedia <- get_mean_clusters(assg = assg)


# Funcion que nos permite encontrar la solucion de cluster a una distancia de 
# corte determinada, sin tomar en cuenta el orden de las distancias ultrametricas.
#https://r.789695.n4.nabble.com/Cutting-hierarchical-cluster-tree-at-specific-height-fails-td4693737.html
cutree.h <- function(tree,h) {
  # this line adapted from cutree(...) code
  k <- nrow(tree$merge) + 2L - apply(outer(c(hc$height, Inf), h, ">"), 2, which.max)
  return(cutree(tree,k=k))
}





# Funcion que nos da un dataframe con dos columnas:
# la distancia de corte ultrametrica y la media de la distancia
# de acoples de los cluster a esa distancia
# NOTA: el vector de distancias de cortes dcut es mejor dejarlo 
# afuera, debido a que este vector tiene que ser el mismo para
# todos los demas algortimos de jerarquizacion que utilizamos, aun
# cuando sepamos que el minimo y el maximo de distancia ultrametrica
# siempre va a ser el mismo en todos los algoritmos.
get_vector_mean <- function(J, D, merge, number_of_intervals, dcut) {
  cl = which(colnames(merge)=="dultr")
  hc <- to_dendo(D, merge[,c(1,2,3, cl)], enames=c(1:ncol(D)) )
  
  
  solucion <- matrix(NA, ncol=2, nrow=length(dcut))
  for (n in 1:length(dcut) ) {
    #  testing: 19-nov-19
    #assignments <- cutree(hc, h = dcut[n]) # nos da un vector con el cluster asignado a cada nodo
    assignments <- cutree.h(hc, h = dcut[n]) # nos da un vector con el cluster asignado a cada nodo
    
    node_names <- -1*as.numeric(names(assignments))
    clust <- assignments
    # assg es dataframe con los nodos y su correspondiente asignacion de cluster
    assg <- data.frame(node_names = node_names, cl = clust)
    
    lamedia <- get_mean_clusters(assg = assg)
    
    solucion[n, ] <- c(dcut[n], lamedia )  # este seria el output.
  }
  colnames(solucion) <- c("cut_distance", "cluster_mean_of_dc")
  solucion <- as.data.frame(solucion)
  return(as.numeric(solucion$cluster_mean_of_dc))
}
# ejemplo
dcut <- get_umbrales(merge = merge, number_of_intervals = 10)
get_vector_mean(J=J, D=D, merge=merge, number_of_intervals=10, dcut=dcut)









# evaluacion de clusters.
# https://www.datanovia.com/en/lessons/determining-the-optimal-number-of-clusters-3-must-know-methods/



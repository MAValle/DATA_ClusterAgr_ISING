# **************** Greedy approach to get hierarchical clustering ****************
# Intentamos crear una funcion que va fusionando nodos y clusters 
# minimizando la distancia de acople en cada proceso de fusion.

# Es una aproximacion greedy.

# Procedimiento:
# 1. find min w of E, where E = todos los esges con distancias
#   nodes of e, forma a  cluster Ci 
# 2.1 search for a Ci another node (or cluster) that minimize sum of w
# 2.2 search another edge not in e, with min w
# 2.3 compare 2.1 and 2.2 and merge the one with min sum of w. Put clusters in {include}
# 3 back to 2.1 and repeat until all nodes are in cluster {include} 


# actual name: testing_hierarchical_clusteringV10.R

#Notas:
# 06-noc-19: creation 




# # # # # # # # # # # # # # # # # # LOADING DATA AND FUNCTIONS # # # # # # # # # # # # # # # # # # # # # # #
# Carga de datos, 
# igual que en clusteragrr_v1.R
rm(list = ls())
#load('new_inferring_parameters_environment270219.RData') # load H and J parameters ising inferred
source("functions_hclust.R") # funciones varias
source("hierarchical_clustering_v2_function.R")
source("hierarchical_clustering_v5_function.R")

source("create_coupling_function.R") # crea matriz J de acoples y D de distancias para simular
source("create_mst_function.R") # crea MST a partir de una matriz de acople.
source("acople_distance_sum_function.R") # suma las distantcias de acople dc
source("find_ady_dc_function.R") # encuentra los nodos adyacentes en un mst dado un nodo
source("get_num_nodes_function.R") # nos dice el numero de nodos que tiene el par de clusters a fusionar
source("get_num_nodes_of_ady_clusters_function.R") # nos dice informacion del numero de nodos que son adyacentes a un cluster
source("get_name_of_the_cluster_function.R") # dado el nombre de un nodo, nos dice a que cluster pertenece
source("get_nodes_of_the_cluster_function.R") # nos da los nodos involucrados en un cluster.
source("simulation_hc_function.R") # Genera una simulacion de hierarchical clustering con algoritmo single linkage y modificado.

source("find_mst_barrier_function.R")
source("is_integer0_function.R")
source("entropies_functions.R") # nos interesa aqui get_energy function

library(igraph)
library(Matrix)
library(ggplot2)
# # # # # # # # # # # # # # # # # # LOADING DATA AND FUNCTIONS # # # # # # # # # # # # # # # # # # # # # # #


# simula 
Jmean = 0
Nn=6
sj=1

set.seed(123)
ot <- create_coupling(Nn=Nn, media=Jmean, sj=sj, lw=-3, up=3)
J <- ot$J
D <- ot$D

# inicio - preparacion
N <- ncol(D)
merge <- matrix(NA, ncol=5, nrow=N-1) # col1 y col2 spins merged, col3 numero de cluster, col4 = merge distance
merge <- as.data.frame(merge)
colnames(merge) <- c("node1", "node2", "cluster", "dultr", "dc")
iteraciones <- N - 2
mn <- list() # aqui iremos agregando los nombres de los spins que componen cada clusters 

it=5
# # # # # # # # #  Find the minimum distance among clusters:
rl <- cluster_find_name(D)
id <- rl$id
nombres <- rl$nombres
m <- rl$m
# leaves (original spins) will be negative, merged clusters will be positive starting from 1, 2, ...
nombres_num <- as.numeric(nombres) # nombres de los spins


# # # # # # # # # poniendo los nombres reales de los spins que componen el cluster recien formado.
# detectamos si un nombre de cluster es positivo, en tal caso, extraer los nombres de los spins de ese cluster
te <- nombres_num
hay_cluster_anterior <- which(te > 0)
if( !is.integer0(hay_cluster_anterior) ) {
  el_cluster <- nombres_num[hay_cluster_anterior]
  # rescatar los nombres del cluster en cuestion
  te2 <- c(nombres_num, unlist(mn[el_cluster]) ) 
  te <- te2[-hay_cluster_anterior]
}
mn <- c(mn, list( as.character(te )) ) # y ahora en mn[[2]] cluster 2, se compone de -5 y -6


# calculo de las distancias de acople
d_cpl <- acople_distance_sum2(J, te)  ### > OJo que esta funcion esta en acople_distance_sum_function.R

merge[it,] <- as.numeric( c(nombres_num, it, m, d_cpl) ) #  nuevo: se coloca nodo1 (fusion) / nodo2 (fusion) /nombre del nuevo cluster / distancia a la que se unen


# # # # # # # # # Upgrading the distance matrix
# fila id[1] y columna id[2] se deben borrar
Dold <- D
D <- D[-id, -id]

# pero antes hay que buscar las distancias minimas entre los spins restantes y el cluster
if (length(Dold) <= 9) { # cuando la matriz es de 3X3, no se puede ejecutar directamente la funcion find_min_dist
  # 09-nov-19
  sel <- colno[!(colno %in% as.character(nombres_num) )] # nodo identificado 16-oct-19
  involved_nodes <-  c( mn[[it]], sel )
  tempo <- as.numeric(involved_nodes )
  haypositive <- which(tempo > 0)
  if ( length(haypositive) > 0) {
    involved_nodes <- c( mn[[it]], mn[[ tempo[haypositive] ]]  )
  }
  dit <-  acople_distance_sum2(J, y = involved_nodes) 
  nombres_ <- sel
  # FIN # 09-nov-19
} else {
  #temp <- find_min_dist(D, Dold, nombres)
  temp <- find_min_distcp(J = J, D = D, nombres = as.character(nombres_num) ) # 16-oct-19
  dit <- temp$dit # dit son las minimias distancias de cluster fusionado en la iteracion anterior a todos los demas clusters que quedan
  nombres_ = temp$nombres_p # nombres_ son los nombres de los nodos que estan a minima distancia.
}

#ahora tengo que poner los valores en fila.columna de D
D <- put_dis(D, dit, nombres_, N, it)
print(it)


# Last iteration: when it = N - 1
it <- it + 1
rl <- cluster_find_name(D)
id <- rl$id
nombres <- rl$nombres
m <- rl$m
nombres_num <- as.numeric(nombres)

te <- nombres_num
hay_cluster_anterior <- which(te > 0)
if(!is.integer0(hay_cluster_anterior) ) {
  el_cluster <- nombres_num[hay_cluster_anterior]
  # rescatar los nombres del cluster en cuestion
  te2 <- c(nombres_num, unlist(mn[el_cluster]) ) 
  te <- te2[-hay_cluster_anterior]
}
mn <- c(mn, list( as.character(te )) ) # y ahor

# calculo de las distancias de acople
d_cpl <- acople_distance_sum2(J, te)  ### > OJo que esta funcion esta en acople_distance_sum_function.R

merge[it,] <- as.numeric( c(nombres_num, it, m, d_cpl) ) #  nuevo: se coloca nodo1 (fusion) / nodo2 (fusion) /nombre del nuevo cluster / distancia a la que se unen









# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 10.nov.19
# testeo de la funcion hierarchical_clustering_greedy
source("hierarchical_clustering_greedy_function.R")
source("find_min_distcp_function.R")
# simula 
Jmean = 0
Nn=25
sj=1

#set.seed(123)
ot <- create_coupling(Nn=Nn, media=Jmean, sj=sj, lw=-3, up=3)
J <- ot$J
D <- ot$D
#rm(list=setdiff(ls(), c("J", "D")))
#save.image(file='borra_101119.RData')
#load('borra_101119.RData')

mergeG <-  hierarchical_clustering_greedy(D = D, J = J)
mergeG
hcG <- to_dendo2(mergeG[, c(1,2,3,5) ], enames=c(1:ncol(D)))
plot(hcG)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 













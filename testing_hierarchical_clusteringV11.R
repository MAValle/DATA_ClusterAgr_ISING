# ******** PROBABILISTIC Greedy approach to get hierarchical clustering *************
# Intentamos crear una funcion que utiliza un greedy approach 
# para hierachical clustering igual que en la funcion  hierarchical_clustering_greedy_function.R
# desarrollada en testing_hierarchical_clusteringV11.R.
# Esta aproximacion es un probabilistic greedy en donde se calcula
# la PROBABILIDAD de fusionar dos clusters basado en la distancia
# de acople existente entre los dos clusters.

# Probabilidad de acople entre cluster i e j : P({i,j})
# P({i,j}) ~ -E_ij  
# donde E_ij es la energia de acople de todos los nodos en {i,j}
# Como estamos utilizando metrica de distancia de acople dc, en que
# -E_ij ~ dc_ij, entonces:
# P({i,j}) ~ -dc_ij  para ser mas exactos:
# P({i,j}) ~ 1 - ( dc_ij / sum dc_k,l ) 
# donde sum dc_k,k es la sumatoria de todas las distancias de acople 
# de la mitad superior de la matriz de distancias en la iteracion it.

# Procedimiento:
# 1. find min w of E, where E = todos los edges con distancias
#   nodes of e, forma a  cluster Ci 
# 2.1 search for a Ci another node (or cluster) that minimize sum of w
# 2.2 search another edge not in e, with min w
# 2.3 compare 2.1 and 2.2 and merge the one with min sum of w, 
# with PROBABILITY P({i,j}). Put clusters in {include}
# 3 back to 2.1 and repeat until all nodes are in cluster {include} 


# actual name: testing_hierarchical_clusteringV11.R

#Notas:
# 13-nov-19: creation 




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
source("find_min_distcp_function.R") #permite encontrar las distancias de acople entre un cluster y los restantes

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
merge <- matrix(NA, ncol=6, nrow=N-1) # col1 y col2 spins merged, col3 numero de cluster, col4 = merge distance
merge <- as.data.frame(merge)
colnames(merge) <- c("node1", "node2", "cluster", "prob", "dultr", "dc") # prob es la probabilidad de activacion del cluster seleccionado 
iteraciones <- N - 2
mn <- list() # aqui iremos agregando los nombres de los spins que componen cada clusters 

it=1
# # # # # # # # #  Find probabilies of cluster selection:
rl <- pick_a_cluster(D)
id <- rl$id
nombres <- rl$nombres
m <- rl$m
prob <- rl$p
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

merge[it,] <- as.numeric( c(nombres_num, it, prob, m, d_cpl) ) #  nuevo: se coloca nodo1 (fusion) / nodo2 (fusion) /nombre del nuevo cluster / distancia a la que se unen


# # # # # # # # # Upgrading the distance matrix
# fila id[1] y columna id[2] se deben borrar
Dold <- D
D <- D[-id, -id]

# pero antes hay que buscar las distancias minimas entre los spins restantes y el cluster
if (length(Dold) <= 9) { # cuando la matriz es de 3X3, no se puede ejecutar directamente la funcion find_min_dist
  # 09-nov-19
  colno <- colnames(Dold)
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
  temp <- find_min_distcp(J = J, D = D, nombres = as.character(nombres_num), mn=mn, it=it ) # 16-oct-19
  dit <- temp$dit # dit son las minimias distancias de cluster fusionado en la iteracion anterior a todos los demas clusters que quedan
  nombres_ = temp$nombres_p # nombres_ son los nombres de los nodos que estan a minima distancia.
}

#ahora tengo que poner los valores en fila.columna de D
D <- put_dis(D, dit, nombres_, N, it)
print(it)


# Last iteration: when it = N - 1
it <- it + 1
rl <- pick_a_cluster(D)
id <- rl$id
nombres <- rl$nombres
m <- rl$m
nombres_num <- as.numeric(nombres)
prob <- rl$p

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

merge[it,] <- as.numeric( c(nombres_num, it, prob, m, d_cpl) ) #  nuevo: se coloca nodo1 (fusion) / nodo2 (fusion) /nombre del nuevo cluster / distancia a la que se unen









# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 10.nov.19
# testeo de la funcion hierarchical_clustering_probablistic_greedy
source("hierarchical_clustering_probabilistic_greedy_function.R")
source("find_min_distcp_function.R")
source("pick_a_cluster_function.R")
# simula 
Jmean = 0
Nn=6
sj=1

#set.seed(123)
ot <- create_coupling(Nn=Nn, media=Jmean, sj=sj, lw=-3, up=3)
J <- ot$J
D <- ot$D
#rm(list=setdiff(ls(), c("J", "D")))
#save.image(file='borra_101119.RData')
#load('borra_101119.RData')

mergePG <-  hierarchical_clustering_probabilistic_greedy(D = D, J = J, trh=1)
mergePG
hcPG <- to_dendo2(mergePG[, c(1,2,3,6) ], enames=c(1:ncol(D)))
plot(hcPG)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 













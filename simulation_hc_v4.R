# En este script vamos a simular clusterig jerarquicos de 
# distintos tamanos de nodos de N=  25, 50, 100.
# con matrices de acoples con <J>.

# La intencion es hacer un mapa de calor en que en el eje X va 
# la iteracion de merge (para un determinado numero de nodos)
# y ne el eje Y va <J>.
# Es decir, para un determinado numero de nodos, fabricamos
# una matriz con las iteraciones en un eje y <J> por otro, y la 
# en la superficie graficamos la distancia de acople.

# Podemos hacer el mapa de calor, uno para la distancia de acople
# de merge2, otra de merge5 y otra que sea la diferencia.



# Procedimiento:
# 


# actual name: simulation_hc_V4.R
# creation: 05.nov.19


# Notes:
# 05-nov-19: creation
# 



# # # # # # # # # # # # # # # # # # LOADING DATA AND FUNCTIONS # # # # # # # # # # # # # # # # # # # # # # #
# Carga de datos, 
# igual que en clusteragrr_v1.R
rm(list = ls())
#load('new_inferring_parameters_environment270219.RData') # load H and J parameters ising inferred
source("functions_hclust.R") # funciones varias
source("hierarchical_clustering_v2_function.R")
source("hierarchical_clustering_v5_function.R")
source("hierarchical_clustering_greedy_function.R")

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

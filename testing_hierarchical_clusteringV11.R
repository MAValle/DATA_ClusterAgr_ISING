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
# 16-nov-19: Comparamos para distinto numero de nodos N, las 
#             distancias de acoples que se van formando para distintos
#             valores de trh (nivel de aleatoridad al greedy).
#             Tomamos como referencia el greedy normal 
# 17-nov-19: comenzamos a probar obtener la distancia de acople
#           en distintas iteraciones. El procedimiento es seleccionar
#           distintos umbrales de distancia ultrametrica o de acople dc
#           y calcular la media de la distancia de acople de los clusters
#           encontrados a ese humbral.



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
#   * trh = threshold entre 0 y 1: cuando es 0, es totalmente random, 1 es totalmente deterministico
mergePG <-  hierarchical_clustering_probabilistic_greedy(D = D, J = J, trh=1)
mergePG
hcPG <- to_dendo2(mergePG[, c(1,2,3,6) ], enames=c(1:ncol(D)))
plot(hcPG)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 16-nov-19
# Comparamos para distinto numero de nodos N, las 
# distancias de acoples que se van formando para distintos
# valores de trh (nivel de aleatoridad al greedy).
# Tomamos como referencia el greedy normal 
N <- c(10, 25, 50, 100)
thres <- c(1, 0.5, 0) # 1 greedy deterministico, 0 = greedy random
Jmean = 0
sj = 1
ot <- create_coupling(Nn=N[3], media=Jmean, sj=sj, lw=-3, up=3)
J <- ot$J
D <- ot$D

# Nota:
# mergeG <-  hierarchical_clustering_greedy(D = D, J = J)
# mergePG <-  hierarchical_clustering_probabilistic_greedy(D = D, J = J, trh=1)
# mergeG y mergePG son lo mismo cuendo trh = 1
sims <- length(thres)
df <- matrix(NA, ncol=3, nrow = 0) # threshold utilizado, iteracion, distancia de acople
df <- as.data.frame(df)
for (s in 1:sims) {
  merge <- hierarchical_clustering_probabilistic_greedy(D = D, J = J, trh=thres[s])
  tr <- rep(thres[s], nrow(merge))
  temp <- data.frame(trh = tr, iteration = 1:nrow(merge), dc = merge$dc)
  df <- rbind(df,temp)
}
# vamos a comparar tambien con merge single linkage:
mst_g <- create_mst(J=J, D=D)
merge <- hierarchical_clustering_v2(D=D, J=J, mst_g=mst_g)
tr <- rep("SL", nrow(merge))
temp <- data.frame(trh = tr, iteration = 1:nrow(merge), dc = merge$dc)
df <- rbind(df,temp)
# ahora unimos todos los resultados
df$trh <- as.factor(df$trh)
duplot <- ggplot(df, aes( x = iteration, y = dc, group = trh ) )  + 
  geom_point( aes(shape = trh) )  +
  #geom_line(linetype = "dotted" ) +
  scale_shape_manual(values = c(3, 20, 5, 1, 6, 21)) + 
  #scale_x_continuous(breaks = c(1:nrow(merge) ) ) +
  scale_x_continuous(breaks = seq(1,nrow(merge), 2) ) +
  xlab('Iteration') + ylab('Ultrametric Distance') + theme_light()
duplot

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 















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



# En este script vamos a simular clusterig jerarquicos de 
# distintos tamanos de nodos de N= 20, 30, 50, 100, 250, y 500
# con matrices de acoples con <J>=0.
# Luego haremos una grafica scatterplot en que graficamos en el 
# eje X la distancia de acople del algoritmo normal, y en el eje Y
# la distancia de acople del algoritmo modificado para cada una
# de las iteraciones.

# Procedimiento:



# actual name: simulation_hc_V5.R
# creation: 19.nov.19


# Notes:
# 19-nov-19: creation


  
  
  
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
source("find_min_distcp_function.R") # funcion que encuentra min dist acople en greedy alg.
source("pick_a_cluster_function.R")

source("find_mst_barrier_function.R")
source("is_integer0_function.R")
source("entropies_functions.R") # nos interesa aqui get_energy function
source("multiplot_function.R")

library(igraph)
library(Matrix)
library(ggplot2)
# # # # # # # # # # # # # # # # # # LOADING DATA AND FUNCTIONS # # # # # # # # # # # # # # # # # # # # # # #

# # # # # # # # # # # # # # # # # # # # # FUNCIONES NECESARIAS # # # # # # # # # # # # # # # # # # # # # # #
# FUNCIONES que fueron creadas en testing_hierarchical_clusteringV11.R
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
# dcut <- get_umbrales(merge = merge, number_of_intervals = 10)
# get_vector_mean(J=J, D=D, merge=merge, number_of_intervals=10, dcut=dcut)

# Funcion equivalente a get_vector_mean, pero para calcular la media
# de las distancias de acoples de los clusters segun NUMERO DE CLUSTERS, y no
# distancia de corte.
get_vector_mean2 <- function(J, D, merge, cls) {
  cl = which(colnames(merge)=="dultr")
  hc <- to_dendo(D, merge[,c(1,2,3, cl)], enames=c(1:ncol(D)) )
  
  
  solucion <- matrix(NA, ncol=2, nrow=length(cls))
  for (n in 1:length(cls) ) {
    #  testing: 19-nov-19
    assignments <- cutree(hc, k = cls[n]) # nos da un vector con el cluster asignado a cada nodo
    #assignments <- cutree.h(hc, k = cls[n]) # nos da un vector con el cluster asignado a cada nodo
    
    node_names <- -1*as.numeric(names(assignments))
    clust <- assignments
    # assg es dataframe con los nodos y su correspondiente asignacion de cluster
    assg <- data.frame(node_names = node_names, cl = clust)
    
    lamedia <- get_mean_clusters(assg = assg)
    
    solucion[n, ] <- c(cls[n], lamedia )  # este seria el output.
  }
  colnames(solucion) <- c("cut_distance", "cluster_mean_of_dc")
  solucion <- as.data.frame(solucion)
  return(as.numeric(solucion$cluster_mean_of_dc))
}

# # # # # # # # # # # # # # # # # # # # # FUNCIONES NECESARIAS # # # # # # # # # # # # # # # # # # # # # # #




# # # # # # # # # # # # # # # # # # # # # # # SIMULACION 1 # # # # # # # # # # # # # # # # # # # # # # # #
# PROCEdure
# BASADO EN DISTANCIAS DE CORTE Y SOLO PARA ALGORITMO GREEDY
# 1. inicializamos J, D 
# 2. determinamos numero de cortes deseados: LA SIMULACION SE BASA EN UNA DISTANCIA DE CORTE AL DENDOGRAMA
# 3. se obtiene merge
# 4. se obtienen las medias de las distancias de acoples de lso cluster generados a cada distancia de corte       
# 5. se grafica con boxplot o scatetrplot las medias 

set.seed(123)
N <- c(10, 25, 50, 100)
thres <- c(1, 0.8, 0.5, 0.3, 0) # 1 greedy deterministico, 0 = greedy random
Jmean = 0
sj = 1
ot <- create_coupling(Nn=N[3], media=Jmean, sj=sj, lw=-3, up=3)
J <- ot$J
D <- ot$D

# Para SL: no lo puedo hacer porque la distancia ultrametrica de
# est algoritmo no es la distancia de acople!
# Tendria que hacerlo aparte! pero con distancia MST, no de acople

# Para greedy....
noi = 20 # numero de intervalos deseados para las distancias de corte.
MM <- matrix(, ncol=3)
M <- matrix(NA, ncol=3 , nrow=noi-1)
for (g in 1:length(thres) ) {
  merge <- hierarchical_clustering_probabilistic_greedy(D = D, J = J, trh=thres[g])
  dcut <- get_umbrales(merge = merge, number_of_intervals = noi)
  M[, 1] <- dcut
  M[, 2] <- get_vector_mean(J=J, D=D, merge=merge, number_of_intervals=noi, dcut=dcut)  
  M[, 3] <- rep(thres[g], noi-1)
  MM <- rbind(MM, M)  
}
df <- as.data.frame(MM)
df <- df[complete.cases(df), ]
colnames(df) <- c("Dcut", "Coupling_Distance", "Threshold")

library(ggplot2)
ggplot(df, aes(x=Dcut, y=Coupling_Distance, color=factor(Threshold), shape=factor(Threshold)))+
  geom_point() + 
  #scale_shape_manual(values=c(3, 16, 17)) +
  #scale_color_manual(values=c('#999999','#E69F00', '#56B4E9')) +
  #scale_size_manual(values=c(10,10,12)) +
  theme(legend.position="top") + 
  xlab('Cut Ultrametric distance') + ylab('Mean of Coupling Distance of clusters') + theme_light()

# CONCLUSIONES:
# me doy cuenta que parece mejor hacer la simulacion basada en CLUSTER desdeados Y NO EN DISTANCIAS DE corte.








# # # # # # # # # # # # # # # # # # # # # # # SIMULACION 2 # # # # # # # # # # # # # # # # # # # # # # # #
# 19-nov-19
# PROCEdure
# BASADO EN DISTANCIAS DE CORTE Y SOLO PARA ALGORITMO GREEDY
# 1. inicializamos J, D 
# 2. determinamos numero clusters deseados : LA SIMULACION SE BASA EN NUMERO DE CLUSTERS
# 3. se obtiene merge
# 4. se obtienen las medias de las distancias de acoples de lso cluster generados a cada distancia de corte       
# 5. se grafica con boxplot o scatetrplot las medias 

set.seed(123)
N <- c(10, 25, 50, 100)
thres <- c(1, 0.8, 0.5, 0.3, 0) # 1 greedy deterministico, 0 = greedy random
Jmean = 0
sj = 1
ot <- create_coupling(Nn=N[3], media=Jmean, sj=sj, lw=-3, up=3)
J <- ot$J
D <- ot$D

# determinar numero de clusters para analizar
cls <- 1:N[2]
MM <- matrix(, ncol=3)
M <- matrix(NA, ncol=3 , nrow = length(cls))
for (g in 1:length(thres) ) {
  merge <- hierarchical_clustering_probabilistic_greedy(D = D, J = J, trh=thres[g])
  #dcut <- get_umbrales(merge = merge, number_of_intervals = noi)
  M[, 1] <- cls
  M[, 2] <- get_vector_mean2(J=J, D=D, merge=merge, cls=cls)  
  M[, 3] <- rep(thres[g], length(cls))
  MM <- rbind(MM, M)  
}
df <- as.data.frame(MM)
df <- df[complete.cases(df), ]
colnames(df) <- c("Dcut", "Coupling_Distance", "Threshold")


library(ggplot2)
ggplot(df, aes(x=Dcut, y=Coupling_Distance, color=factor(Threshold), shape=factor(Threshold)))+
  geom_point() + 
  #scale_shape_manual(values=c(3, 16, 17)) +
  #scale_color_manual(values=c('#999999','#E69F00', '#56B4E9')) +
  #scale_size_manual(values=c(10,10,12)) +
  theme(legend.position="top") + 
  xlab('At desired number of clusters') + ylab('Mean of Coupling Distance of clusters') + theme_light()

# # # # # # # # # # # # # # # # # # # # # # # SIMULACION 2 # # # # # # # # # # # # # # # # # # # # # # # #










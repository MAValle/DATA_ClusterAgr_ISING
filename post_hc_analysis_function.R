# 17-dic-19 : Set de funciones para analizar los resultados de 
# varias simulaciones de los algoritmos de hiearchical clustering.



# # # # # # # # # # # # # # # # # # # # # FUNCIONES NECESARIAS # # # # # # # # # # # # # # # # # # # # # # #
# FUNCIONES que fueron creadas en testing_hierarchical_clusteringV11.R
# 17-nov-19
# Aqui creamos una funcion que para distintas distancias ultrametricas
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

# 300120: modificamos la funcion para que nos arreje tambien
# las medias de distancia de acople de cada cluster.
# funcion que nos da la media de distancia de acople de los clusters
# encontrados a cierta distancia de corte
# Input: assg que es dataframe con dos columna: nombre de los nodos / solucion de cluster
# debe venir con nombres de columan node_names y cl (entero)
get_mean_clusters <- function(assg) {
  # aqui calculamos la media de la distancia de acoples para una
  # distancia de corte determinada h.
  uniq_clustrs <- unique(assg$cl)
  # inf tiene numero de cluster, suma de las energias de acople dentro del cluster y numero de elementos de ese cluster
  inf <- matrix(NA, ncol=3, nrow=length(uniq_clustrs))
  for (cc in 1:length(uniq_clustrs) ) {
    sub <- subset(assg, cl==uniq_clustrs[cc])
    te <- as.numeric(sub$node_names)
    if ( length(te) > 1) {
      coupl_distance <- acople_distance_sum2(J, te)  
    } else {
      coupl_distance <- NA
    }
    
    inf[cc, ] <- c(uniq_clustrs[cc], coupl_distance, nrow(sub))
  }
  colnames(inf) <- c("cluster", "sum_coupl_energy", "number_of_elements")
  averages <- as.data.frame(inf) 
  average_d <- mean(inf[,2], na.rm = TRUE)
  return(list(averages = averages, main_average = average_d) )
}
# ejemplo:
#medias <- get_mean_clusters(assg = assg)
#medias$main_average
#medias$averages


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


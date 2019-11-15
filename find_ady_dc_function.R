
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
# creation: 23.oct.19
# viene originalmente de functions_hclust.R 13.sep.19 FUNCION 
#
# # funcion que para un cluster en analisis, nos dice sus nodos adyacentes en el MST y la 
# # distancia de acople dc = f(sumatoria de las energias de acoples entre pares de nodos)
# output: matriz con tres columnas: el cluster bajo analisis, el nodo adyacente y su distancia de acople dc
find_ady_dc <- function(v, n, mst_g, J) {
  registro_interno <- matrix(NA, ncol=3, nrow=0)
  colnames(registro_interno) <- c("cluster_n", "elnodo", "dist_acople") #tengo que ir registrando el cluster (n), el elnodo y distancia_acople entre cluster y elnodo
  # sus nodos adyacentes en el MST que tienen la menor distancia de acople
  ady <- adjacent_vertices(mst_g, as.character(v))
  nodoss <- length(ady) # numero de nodos adyacentes al cluster 
  for (m in 1:nodoss) {
    elnodo <- setdiff(as.character(-1*as.numeric(ady[[m]]) ) , as.character(v)) # encontramos el nodo que esta en ady[[.]] pero no en v
    if (vector.is.empty(elnodo) ) { next } # en caso que elnodo sea vacio
    #elnodo <- -1*as.numeric(elnodo)
    for (x in 1:length(elnodo) ) {
      nodo <- as.numeric(elnodo[x])
      distancia_acople <- acople_distance_sum2(J, c(v, nodo))  ### > OJo que esta funcion esta en acople_distance_sum_function.R
      registro_interno <- rbind(registro_interno, c(n, nodo, distancia_acople))
    }
  }
  return(registro_interno)
}
# ejemplo
# registro_interno <- find_ady_dc(v=lista_nodos_de_clusters[n], n=n, mst_g = mst_g, J = J)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  





# Funcion necesaria para funcion hierarchical_clustering_greedy

# # # # # # # # # # # # # # # # # # NECCESARY FUNCTION # # # # # # # # # # # # # # # # # # # # # # #
# 13-NOV-19
# Creada en testing_gierarchical_clusteringV11.R
# funcion para seleccionar el par de cluster a fusionar en forma
# probabilistica
# https://stats.stackexchange.com/questions/67911/how-to-sample-from-a-discrete-distribution
sampleDist = function(cluster, n=1, probs) { 
  sample(x = cluster, n, replace = T, prob = probs) 
}
# input: matriz de distancia D de la iteracion 
# output: 
#   * $id: vector conla fila y columna de D a seleccionar
#   * $nombres: los nombres de los clusters seleccionados
#   * $m la distancia de acople de los cluster seleccionados
# Nota: para que a mayor distancia la probabilidad sea menor, 
# tenemos que restar 1.
pick_a_cluster <- function(D) {
  D_ <- D
  D_[lower.tri(D_, diag=TRUE)] <- NA
  
  #  Get array indices
  ind <- which( ! is.na(D_) , arr.ind = TRUE ) 
  #  cbind indices to values
  out <- cbind( D[ ! is.na( D_ ) ] , ind )
  out <-  as.data.frame( out ) 
  out$prob <- 1 - out$V1/sum(out$V1)
  
  cluster <- 1:nrow(out)
  idx <- sampleDist(cluster=cluster, n=1, probs=out$prob)
  id <- as.numeric( out[idx, c(2,3)] )
  nombres <- colnames(D)[id]
  m <- D[id[1], id[2]]
  p <- out[idx,4]
  
  return(list(id = id, nombres = nombres, m = m, p = p))
}
# ejemplo
#pick_a_cluster(D = D)
# # # # # # # # # # # # # # # # # # NECCESARY FUNCTION # # # # # # # # # # # # # # # # # # # # # # #



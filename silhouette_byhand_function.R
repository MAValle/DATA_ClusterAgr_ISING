# computing sihouette by hand function
# jan 27, 2020
# silhouette_byhand_function.R

# https://en.wikipedia.org/wiki/Silhouette_(clustering)
# silhouette:
# a(i) = 1/(#Ci - 1) sum d(i,j)  suma de todas las distancias desde i a j, donde
# j son elementos intracluster

# b(i) = min (1/#Ck) sum d(i,j) min de la suma de las distancias desce i a j
# donde j son elementos fuera del cluster donde esta i.


# Procedure
# 1 ver cuantos clusters hay
# Calculo de a(i)
# 2 identificar las posiciones para cluset cl en member
# 3 para cada elemento i , buscar en D la distancia 

# input
# members = vector con los numero de cluster a que corresponde cada obs o cada producto.
# D = distance matrix (debe venor con los nombres de columnas reales)
# output:
# vector con el s(i) de cada elemento.
silhouette_byhand <- function(members, D) {
  num_of_clusters <- length(unique(members))
  number_of_spins <- length(colnames(D))
  clusters <- unique(members)
  sil <-  numeric(number_of_spins)
  A <- numeric(number_of_spins) # valores de a(i)
  B <- numeric(number_of_spins) # valores de b(i)
  b_external <- numeric(number_of_spins) # esto es promedio de acoples entre nodo i, y todo los nodos externos al cluster donde esta i.
  # valores de silhouette
  for (s in 1:number_of_spins) {
    
    cluster_of_spin <- members[s]
    # calculo de a(i)
    id_clusters <- which(members == cluster_of_spin)
    if ( (length (id_clusters) ) == 1 ) {
      sil[s] <- 0
      A[s] <- NA
      next
    }
    dd <- D[s, id_clusters]
    dd <- dd[!is.infinite(dd)]
    a <- mean( dd)
    A[s] <- a
    # calculo de b(i)
    cls <-  setdiff(clusters, cluster_of_spin)
    test <- length(cls)
    b_values <- numeric(test)
    x <- vector(mode="numeric", length=0)
    for(cc in 1:test) {
      id_clusters <- which(members ==cls[cc]) # selecciona los id de los elementos del cluster
      dd <- D[s, id_clusters] # extrae los elementos del cluster
      dd <- dd[!is.infinite(dd)] # saca el termino Inf de la diagonal
      b_values[cc] <- mean( dd ) # calcula la media de los distancias entre pares de elementos del cluster
      x <- c(x,dd)
    }
    b <- min( b_values )
    B[s] <- b
    b_external[s] <- mean(x, na.rm=TRUE) # esto es promedio de acoples entre nodo i, y todo los nodos externos al cluster donde esta i.
    #silhouette computation
    if ( a < b  ) {
      sil[s] <- 1 - (a/b)
    } else if ( a == b ) {
      sil[s] <- 0
    } else {
      sil[s] <- (b/a) - 1 
    }
  }
  output <- cbind(as.numeric(names(members_det)) , a=A, b=B, b_external, sil=sil, as.numeric(members_det))
  colnames(output) <- c("vertex", "a", "b", "b_external",  "sil", "cl")
  return(output)
}



# example
#hc_det <- to_dendo2(merge_greedy_detm[, c(1,2,3,6) ], enames = c(1:ncol(D)))
#members_det <- cutree(hc_det, k = 3)
#colnames(D) <- rownames(D) <- real_names
#silhouette_byhand(members = members_det, D=D)
# da lo mismo que:
# silhouette(as.numeric(members_det), as.dist(D) ) from cluster library




# Feb 01, 2020
# En pag. 273, 274 y 276 redefinimos el indice silhouette para adaptarlo
# a nuestra aplicacion.
# El termino a(i) sige siendo la distancia (energia de acople) promedio del elemento i a 
# otros elementos de su mismo cluster. La idea es que este termino sea 
# el mas pequeno posible y menor que sqrt(3).
# El termino b(i) ahora es el promedio de todas las distancias
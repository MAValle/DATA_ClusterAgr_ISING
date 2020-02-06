# computing my own measure of internal consistency simmilar to the silhouette
# Given a solution with a coupling distance h, we have for each element of the
# system their corresponding cluster.
# Fron this solution we want compute a internal consistency to determine which 
# groups of clusters are more choesive and find whether teh solution has a good
# internal consistency.

# Feb 04, 2020
# my_own_silhouette_function.R

# https://en.wikipedia.org/wiki/Silhouette_(clustering)
# Definiciones
# <a(i)> = 1/(#Ci - 1) sum d(i,j)  suma de todas las distancias desde i a j, donde
# j son elementos intracluster
# <a(i)> = 1/(#Ci - 1) sum d(i,j)  es el promedio de todas las distancias de acople 
# del elemento i al j, donde j son elementos dentro del cluster.
# U*(Ck) es la energia de acople del cluster k que es igual a 1/(#Ck) sum <a(i)>, es decir
# el promedio de todos los promedios de a(i) para todos lo elementos del cluster.
# esto representa el nivel de cohesion de elementos intracluster.

# b(i, Ck) es la energía de acople U* de un cluster k externo a i, cuando agregamos al
# cluster k, el elemento k. Esto se hace para todo cluster k, externo a i.
# min {b(i, Ck)} 
# b(i, Ck) sera un vector de largo M-1, donde M es el numero de clusters.
# <b(i)>  es la media de todos lo minimos de b(i, Ck) para cada elemento i del sistema.

# como resultado, tenemos una dataframe en que cada fila es un nodo o spin 
# con si respectivo U*(Ck) y <b(i, Ck)>. La idea es que U*(Ck) sea lo mas pequeno
# <b(i, Ck)> lo mas grande. Mientras mas elementos tenga un cluster, mas distancia
# de acople tendra, por lo que para comparar U*(Ck) y <b(i, Ck)>, dividimos cada
# termino por el numero de elementos de cada cluster.
# Obviamente para elementos de un mismo cluster, el U*(Ck) sera el mismo. Lo
# que ira cambiando es el <b(i, Ck)>






# Inputs
# members = vector con los numero de cluster a que corresponde cada obs o cada producto o item o spin.
# D = distance matrix (debe venor con los nombres de columnas reales)
# Output:
# dataframe con el nombre del elemento, <a(i)> y <b(i)>.
my_own_silhouette <- function(members, D) {
  
  num_of_clusters <- length(unique(members))
  number_of_spins <- length(colnames(D))
  clusters <- unique(members)
  
  A <- numeric(number_of_spins) # valores de <a(i)>  vector de largo igual al mumero de elementos del sistema
  B <- numeric(number_of_spins) # valores de <b(i)>  vector de largo igual al mumero de elementos del sistema
  
  for (s in 1:number_of_spins) {
    
    # vemos en que cluster esta el spin
    cluster_of_spin <- members[s]
    # calculo de a(i)
    id_clusters <- which(members == cluster_of_spin) # vector de id de los elementos del cluster donde esta s.
    if ( (length (id_clusters) ) == 1 ) {
      sil[s] <- 0
      A[s] <- NA
      next
    }
    sum <- acople_energy_sum(D, y=names(id_clusters))
    A[s] <- sum/((length(id_clusters)*(length(id_clusters)-1))/2)
    
    
    # calculo de b(i)
    cls <-  setdiff(clusters, cluster_of_spin) # vector de clusters en los que NO esta el spin s.
    test <- length(cls)
    b_values <- numeric(test)
    x <- vector(mode="numeric", length=0)
    for(cc in 1:test) {
      id_clusters <- which(members ==cls[cc]) # selecciona los id de los elementos del cluster
      id_clusters2 <- c(id_clusters, s) # agregamos el elemento o spin s, para calcular la energía de acople
      # que habria si es que incluimos al spin s en el cluster.
      names(id_clusters2) <- c(names(id_clusters), colnames(D)[s])
      sum <- acople_energy_sum(D, y= names(id_clusters2) )
      vl <- sum/((length(id_clusters2)*(length(id_clusters2)-1))/2)
      x <- c(x, vl)
    }
    b <- min( x )
    B[s] <- b
    
  }
  
  output <- cbind( as.numeric(names(members_det) ) , a=A, b=B, as.numeric(members_det))
  colnames(output) <- c("vertex", "a", "b", "cl")
  output  <- as.data.frame(output)
  return(output)
  
}
# example
# colnames(D) <- rownames(D) <- real_names
# output <- my_own_silhouette(members = cutree.h(hc_det, h = 20), D = D)
# plot(output$a, output$b)
# abline(a=0, b=1, col="red") 
  


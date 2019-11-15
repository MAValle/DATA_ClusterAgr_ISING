



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# es la misma funcion importada desde functions_hclust.R creada el 23-sep-19.
# Hacemos lo mismo que la funcion hierarchical_clustering_V2, pero a la matriz merge (que es el 
# output), le agregamos dos columnas mas, que son la distancia MST del cluster formado, y la 
# distancia de acople del cluster. 
# cabe aclarar que la distancia mst es un proxy energia de mst, mientras que la distancia de 
# acoples es un proxy de la energia de acoples. La gracia de usar distancias en vez de directamente
# los acoples, es que evitamos tener que sumar numeros negativos y positivos.
# nota: 131019: se incorpora como input el mst_g y matriz de acople J
hierarchical_clustering_v2 <- function(D, J, mst_g) {
  N <- ncol(D)
  merge <- matrix(NA, ncol=6, nrow=N-1) # col1 y col2 spins merged, col3 numero de cluster, col4 = merge distance
  merge <- as.data.frame(merge)
  colnames(merge) <- c("node1", "node2", "cluster", "dultr", "dmst", "dc")
  iteraciones <- N - 2
  mn <- list() # aqui iremos agregando los nombres de los spins que componen cada clusters 
  
  # Inicio de las iteraciones que son N - 2
  for (it in 1:iteraciones) {
    # # # # # # # # #  Find the minimum distance among clusters:
    rl <- cluster_find_name(D)
    id <- rl$id
    nombres <- rl$nombres
    m <- rl$m
    # leaves (original spins) will be negative, merged clusters will be positive starting from 1, 2, ...
    nombres_num <- as.numeric(nombres) # nombres de los spins del cluster recien formado
    
    
    # # # # # # # # # poniendo los nombres reales de los spins que componen el cluster recien formado.
    # detectamos si un nombre de cluster es positivo, en tal caso, extraer los nombres de los spins de ese cluster
    te <- nombres_num
    hay_cluster_anterior <- which(te > 0)
    if(!is.integer0(hay_cluster_anterior) ) {
      el_cluster <- nombres_num[hay_cluster_anterior]
      # rescatar los nombres del cluster en cuestion
      te2 <- c(nombres_num, unlist(mn[el_cluster]) ) 
      te <- te2[-hay_cluster_anterior]
    }
    mn <- c(mn, list( as.character(te )) ) # y ahora en mn[[2]] cluster 2, se compone de -5 y -6
    
    # calculo de la distancia MST
    d_mst <- find_mst_barrier(mst_g, te)
    # calculo de las distancias de acople
    d_cpl <- acople_distance_sum2(J, te)  ### > OJo que esta funcion esta en acople_distance_sum_function.R
    
    merge[it,] <- as.numeric(c(nombres_num, it, m, d_mst, d_cpl)) #  nuevo: se coloca nodo1 (fusion) / nodo2 (fusion) /nombre del nuevo cluster / distancia a la que se unen
    
    
    # # # # # # # # # Upgrading the distance matrix
    # fila id[1] y columna id[2] se deben borrar
    Dold <- D
    D <- D[-id, -id]
    
    # pero antes hay que buscar las distancias minimas entre los spins restantes y el cluster
    if (length(Dold) <= 9) { # cuando la matriz es de 3X3, no se puede ejecutar directamente la funcion find_min_dist
      colno <- colnames(Dold) # nombres de columnas o fila de Dold
      #https://www.youtube.com/watch?v=8hSYEXIoFO8&t=78s
      #sel <- colno[!(colno %in% nombres)] # nodo identificado
      sel <- colno[!(colno %in% as.character(nombres_num) )] # nodo identificado 16-oct-19
      col <- which(colno == sel) # identificamos la columna donde esta el nodo sel.
      dit <- vector()
      dit <- c(dit,min(Dold[, col])) # la distancia minima
      nombres_ <- sel
    } else {
      #temp <- find_min_dist(D, Dold, nombres)
      temp <- find_min_dist(D = D, Dold = Dold, nombres = as.character(nombres_num)) # 16-oct-19
      dit <- temp$dit # dit son las minimias distancias de cluster fusionado en la iteracion anterior a todos los demas clusters que quedan
      nombres_ = temp$nombres_ # nombres_ son los nombres de los nodos que estan a minima distancia.
    }
    
    #ahora tengo que poner los valores en fila.columna de D
    D <- put_dis(D, dit, nombres_, N, it)
    print(it)
  }
  # Last iteration: when it = N - 1
  it <- it + 1
  rl <- cluster_find_name(D)
  id <- rl$id
  nombres <- rl$nombres
  m <- rl$m
  nombres_num <- as.numeric(nombres) # nuevo
  # # # # # # # # # poniendo los nombres reales de los spins que componen el cluster recien formado.
  # detectamos si un nombre de cluster es positivo, en tal caso, extraer los nombres de los spins de ese cluster
  te <- nombres_num
  hay_cluster_anterior <- which(te > 0)
  if(!is.integer0(hay_cluster_anterior) ) {
    el_cluster <- nombres_num[hay_cluster_anterior]
    # rescatar los nombres del cluster en cuestion
    te2 <- c(nombres_num, unlist(mn[el_cluster]) ) 
    te <- te2[-hay_cluster_anterior]
  }
  mn <- c(mn, list( as.character(te )) ) # y ahora en mn[[2]] cluster 2, se compone de -5 y -6
  
  # calculo de la distancia MST
  d_mst <- find_mst_barrier(mst_g, te)
  # calculo de las distancias de acople
  d_cpl <- acople_distance_sum2(J, te)  ### > OJo que esta funcion esta en acople_distance_sum_function.R
  merge[it,] <- as.numeric(c(nombres_num, it, m, d_mst, d_cpl)) #  nuevo: se coloca nodo1 (fusion) / nodo2 (fusion) /nombre del nuevo cluster / distancia a la que se unen
  
  return(merge)
  
}
# ejemplo:
#merge <- hierarchical_clustering_v2(D)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  


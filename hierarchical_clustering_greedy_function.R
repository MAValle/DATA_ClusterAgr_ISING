



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Esta funcion se ha desarrollado en testing_hierarchical_clusteringV10.R
# y lo que hace es simplemente una aproximacion greedy a encontrar 
# clusterig jerarquico basado en greedy approach, minimizando la distancia
# de acople dc en cada fusion de los nodos.
hierarchical_clustering_greedy <- function(D, J) {
  N <- ncol(D)
  Nmerge <- ncol(D) - 1
  merge <- matrix(NA, ncol=5, nrow = Nmerge) # col1 y col2 spins merged, col3 numero de cluster, col4 = merge distance
  merge <- as.data.frame(merge)
  colnames(merge) <- c("node1", "node2", "cluster", "dultr", "dc") # nota: aqui dc y dultr son lo mismo.
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
    
    merge[it,] <- as.numeric( c(nombres_num, it, m, d_cpl) ) #  nuevo: se coloca nodo1 (fusion) / nodo2 (fusion) /nombre del nuevo cluster / distancia a la que se unen
    
    
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
      temp <- find_min_distcp(J = J, D = D, nombres = as.character(nombres_num), mn = mn, it = it) # 10-nov-19
      dit <- temp$dit # dit son las minimias distancias de cluster fusionado en la iteracion anterior a todos los demas clusters que quedan
      nombres_ = temp$nombres_p # nombres_ son los nombres de los nodos que estan a minima distancia.
    }
    
    #ahora tengo que poner los valores en fila.columna de D
    D <- put_dis(D, dit, nombres_, N, it)
    
    #print(m)
    print(paste("Greedy Iteration: ", it))
    
  }
  
  
  # Last iteration: when it = N - 1
  it <- it + 1
  rl <- cluster_find_name(D)
  id <- rl$id
  nombres <- rl$nombres
  m <- rl$m
  nombres_num <- as.numeric(nombres)
  
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
  
  merge[it,] <- as.numeric( c(nombres_num, it, m, d_cpl) ) #  nuevo: se coloca nodo1 (fusion) / nodo2 (fusion) /nombre del nuevo cluster / distancia a la que se unen
  
  
  return(merge)
  
}
# ejemplo:
# mergeG <-  hierarchical_clustering_greedy(D = D, J = J)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  


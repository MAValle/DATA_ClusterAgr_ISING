# funciones para generar el cluster jerarquico usando single linkage
# 22-ago-19

# actual name: functions_hclust.R

#---->cluster_find_name
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 22.ago.19
# pequena funcion que, dado la matriz D, nos encuentre la minima distancia entre dos cluster
# e identifique el nombre de ambos. 
# input: D = matriz de distancias
# output: 
cluster_find_name <- function(D){
  m <- min(D[row(D)!=col(D)]) # nos da el valor minimo 
  id <- which(D == min(D), arr = TRUE) # nos da la locacion del valor minimo fila, columna
  id <- as.vector(id)
  id <- id[c(1,2)] # los otros dos son repetidos
  nombres <- rownames(D)[id] #
  nombres <- nombres[c(1,2)]  # los otros dos son repetidos
  return(list(id=id, nombres=nombres, m=m))
}
# ejemplo
# rl <- cluster_find_name(D)
# id <- rl$id
# nombres <- rl$nombres
# m <- rl$m
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 22.ago.19
# pequena funcion para re-actualizar distancias minimas entre un cluster recien fusionado
# y los demas clusters. 
# input: 
#     D = matriz de distancias actualizada.
#     Dold = matriz de distancias de la iteracion anterior valores en character.
#     nombres : nombres de los nodos fusionados en la iteracion
# output: dit = vector con las minimias distancias de cluster fusionado en la iteracion anterior
#         a todos los demas clusters que quedan. Se utiliza single linkage.
# Nota: esta funcion funciona solamente cuando D es de 4X4 hacia arriba, es decir cuando hay más de 9 o mas elementos en D.
# Nota: el vector "nombres" debe estar cargado en memoria.
find_min_dist <- function(D, Dold, nombres) {
  dit <- vector()
  nombres_ <- colnames(D) # vemos los nombres de las filas o columnas que quedan en D
  n <- length(nombres_)
  for (i in 1:n) {
    nodo <- nombres_[i]
    v <- min(Dold[nombres, nodo])
    #dit <- c(dit,min(Dold[, nodo]))
    dit <- c(dit,v )
  }
  return(list(dit=dit, nombres_ = nombres_))
}
# ejemplo:
# temp <- find_min_dist(D, Dold)
# dit <- temp$dit
# nombres_ = temp$nombres_
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 22.ago.19
# pequena funcion para poner las minimas distancias encontradas entre el cluster
# fusionado en la iteracion anterior y los clusters que quedan, con la funcion find_min_dist:
  # input:
  #      dit <- salida de find_min_dist
  # Nota: La matriz de distancias D debe estar cargada en memoria.
put_dis <- function(D, dit, nombres_, N, it){
    D <- cbind(D, dit)
    #rr <- rep(Inf, ncol(D) )
    rr <- c(dit, Inf)
    D <- rbind(D, rr)
    #colnames(D) <- rownames(D) <- c(nombres_, N+it-1)
    colnames(D) <- rownames(D) <- c(nombres_, it) # nuevo: merged clusters will be positive starting from 1, 2, .
    return(D)
}
# ejemplo:
#D <- put_dis(D, dit, nombres_, N, it)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 





# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 27-ago-19
# Esta es la funcion general. 
# esta funcion esta creada desde testing_hierarchical_clusteringV3.R
# Es necesario todas las funciones ubicadas aqui en este script
## este es el codigo que hace toda la magia.

# Inputs:
# 1. Matriz D de DISTANCIAS de N X N (N=num de spins o elementos) formato matrix. 
#             La diagonal debe venir con inf. 
#             Si la matriz D no viene con nombre de columnas, entonces hay que ponerselas.
# OUTPUTS:
# 1. merge: dataframe con: col1 y col2 spins merged, col3 numero de cluster, col4 = merge distance
# En la primera y segunda columna merge, se ponen en negativo las hojas o nombres originales 
# de los spins que se fusionan, mientras que en positivo los clusters que se forman en ramas
# mas altas del arbol.
hierarchical_clustering <- function(D) {
  N <- ncol(D)
  merge <- matrix(NA, ncol=4, nrow=N-1) # col1 y col2 spins merged, col3 numero de cluster, col4 = merge distance
  merge <- as.data.frame(merge)
  colnames(merge) <- c("node1", "node2", "cluster", "dultr")
  iteraciones <- N - 2
  
  # Inicio de las iteraciones que son N - 2
  for (it in 1:iteraciones) {
    # # # # # # # # #  Find the minimum distance among clusters:
    rl <- cluster_find_name(D)
    id <- rl$id
    nombres <- rl$nombres
    m <- rl$m
    # leaves (original spins) will be negative, merged clusters will be positive starting from 1, 2, ...
    nombres_num <- as.numeric(nombres) # nuevo
    merge[it,] <- as.numeric(c(nombres_num, it, m)) #  nuevo: se coloca nodo1 (fusion) / nodo2 (fusion) /nombre del nuevo cluster / distancia a la que se unen
    
    # # # # # # # # # Upgrading the distance matrix
    # fila id[1] y columna id[2] se deben borrar
    Dold <- D
    D <- D[-id, -id]
    
    # pero antes hay que buscar las distancias minimas entre los spins restantes y el cluster
    if (length(Dold) <= 9) { # cuando la matriz es de 3X3, no se puede ejecutar directamente la funcion find_min_dist
      colno <- colnames(Dold) # nombres de columnas o fila de Dold
      #https://www.youtube.com/watch?v=8hSYEXIoFO8&t=78s
      sel <- colno[!(colno %in% nombres)] # nodo identificado
      col <- which(colno == sel) # identificamos la columna donde esta el nodo sel.
      dit <- vector()
      dit <- c(dit,min(Dold[, col])) # la distancia minima
      nombres_ <- sel
    } else {
      temp <- find_min_dist(D, Dold, nombres)
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
  merge[it,] <- as.numeric(c(nombres_num, it, m)) #  nuevo: se coloca nodo1 (fusion) / nodo2 (fusion) /nombre del nuevo cluster / distancia a la que se unen
  
  return(merge)
  
}

#Example:
# rm(list = ls())
# source("functions_hclust.R")
# D = matrix(c(0, 0.5, 0.4286, 1, 0.25, 0.6250, 0.3750,
#              0.5000, 0, 0.7143, 0.8333, 0.6667, 0.2000, 0.7778,
#              0.4286, 0.7143, 0, 1.0000, 0.4286, 0.6667, 0.3333,
#              1.0000, 0.8333, 1.0000, 0, 1.0000, 0.8000, 0.8571,
#              0.2500, 0.6667, 0.4286, 1.0000, 0, 0.7778, 0.3750,
#              0.6250, 0.2000, 0.6667, 0.8000, 0.7778, 0, 0.7500,
#              0.3750, 0.7778, 0.3333, 0.8571, 0.3750, 0.7500, 0), nrow = 7, ncol = 7)
# colnames(D) <- rownames(D) <- -c(1:ncol(D)) # enumeramos todos las hojas o spines originales con valones negativos.
# diag(D) <- rep(Inf, ncol(D))
# merge <- hierarchical_clustering(D)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 23 sep 19
# Hacemos lo mismo que la funcion hierarchical_clustering, pero a la matriz merge (que es el 
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
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 11, oct, 19
# Hacemos lo mismo que la funcion hierarchical_clustering_V2 en donde se hace single linkage
# pero ahora modificamos el single linkage de tal forma que cuando se intenta unir un 
# cluster de mas de un nodo con otro nodo, se busca si ese nodo es el mejor que minimiza
# la energia de acople (o minimiza la sumatoria de las distancias de acople). 
# Esta modificacion hace que se pierda la propiedad de ultrametricidad, y que no unamos
# cluster a la menor distancia de mst, pero no saseguramos que la energia de acople
# sea lo mas baja posible.
# Cabe senalar que siempre se busca solucion a traves del MST, es decir, se busca
# alternativas de nodos para ser fusionados que siempre esten conectado con el MST.

# Cabe aclarar que la distancia mst es un proxy energia de mst, mientras que la distancia de 
# acoples es un proxy de la energia de acoples. La gracia de usar distancias en vez de directamente
# los acoples, es que evitamos tener que sumar numeros negativos y positivos.
# nota: 131019: se incorpora como input el mst_g y matriz de acople J
hierarchical_clustering_v3 <- function(D, J, mst_g) {
  # # # Inicializacion de la matriz merge que es el output
  N <- ncol(D)
  merge <- matrix(NA, ncol=6, nrow=N-1) # col1 y col2 spins merged, col3 numero de cluster, col4 = merge distance
  merge <- as.data.frame(merge)
  colnames(merge) <- c("node1", "node2", "cluster", "dultr", "dmst", "dc")
  iteraciones <- N - 2
  mn <- list() # aqui iremos agregando los nombres de los spins que componen cada clusters 
  
  # Inicio de las iteraciones que son N - 2
  for (it in 1:iteraciones) {
    rl <- cluster_find_name(D)
    id <- rl$id # 
    nombres <- rl$nombres  # 
    m <- rl$m  #
    # leaves (original spins) will be negative, merged clusters will be positive starting from 1, 2, ...
    nombres_num <- as.numeric(nombres) # nuevo
    
    # segundo, ahora chequeamos si el par de cluster identificados, son 2 nodos, o 
    # hay mas de 2 nodos.
    # # # # este es el cerebro de la funcion
    te <- nombres_num
    hay_cluster_anterior <- which(te > 0)
    if (is.integer0(hay_cluster_anterior)) { # los dos cluster son solo de 1 nodo cada uno
      mn <- c(mn, list( as.character(te )) )
      v <- te
    } else {  # los dos cluster involucran mas de 2 nodos
      # aqui comienza la busqueda de alternativas de fusion que MINIMIZAN dc (distancia de acople)
      num_cluster_anteriores <- length(hay_cluster_anterior) #numero de cluster que tienen mas de 2 nodos
      if (num_cluster_anteriores > 1) { # en caso que se proponga la fusion de dos clusters cada uno con mas de 1 nodos
        # haremos que se fusionen incondicionalmente.
        lista_nodos_de_clusters <- list()
        for (n in 1:num_cluster_anteriores) { # aqui guardamos los nodos de los clusters con mas de 2 nodos
          el_cluster <- te[n]
          #lista_nodos_de_clusters  <- c(lista_nodos_de_clusters, unlist(mn[el_cluster]) )
          lista_nodos_de_clusters <- c(lista_nodos_de_clusters, (mn[el_cluster]) )
        }
        v <- unlist(lista_nodos_de_clusters)
        mn <- c(mn, list( as.character(v )) )
        
      } else { # en caso que se proponga la fusion de dos clusters, en donde uno de ellos tenga solo un nodo
        # comenzaremos a probar si existe otro nodo que minimice la distancia de acople.
        num_cluster_anteriores <- length(hay_cluster_anterior) #numero de cluster que tienen mas de 2 nodos
        lista_nodos_de_clusters <- list()
        for (n in 1:num_cluster_anteriores) { # aqui guardamos los nodos de los clusters con mas de 2 nodos
          #el_cluster <- te[n]   te[which(te > 0)]
          el_cluster <- te[which(te > 0)]
          #lista_nodos_de_clusters  <- c(lista_nodos_de_clusters, unlist(mn[el_cluster]) )
          lista_nodos_de_clusters <- c(lista_nodos_de_clusters, (mn[el_cluster]) )
        }
        registro <- matrix(NA, ncol=3, nrow=0)
        colnames(registro) <- c("cluster_n", "elnodo", "dist_acople") #tengo que ir registrando el cluster (n), el elnodo y distancia_acople entre cluster y elnodo
        for (n in 1:num_cluster_anteriores) { # ahora buscamos para cada cluster con mas de dos nodos, 
          # sus nodos adyacentes en el MST que tienen la mayor distancia de acople
          #registro_interno <- find_ady_dc(v=unlist(lista_nodos_de_clusters[n]), n=te[n])
          registro_interno <- find_ady_dc(v=unlist(lista_nodos_de_clusters[n]), n=el_cluster, mst_g = mst_g, J = J)
          registro <- rbind(registro, registro_interno)
        }
        # veamos cual seria la proxima fusion de acuerdo a la minima distancia de acople
        id2 <- which.min(registro[,3])
        proxima_fusion <- registro[id2, c(1,2)]
        
        #ahora tendriamos que ver si el nodo que esta en registro[id, 2] pertenece ya a otro cluster mas grande
        #nombres_num <- final_fusion(nombre_columnas_D=colnames(D), proxima_fusion=proxima_fusion, elregistro=registro, id=id2) # esto ya no es necesario
        nombres_num <- as.numeric(proxima_fusion)
        #proxima_fusion_numeric <- as.numeric(nombres_num)
        #nombres <- nombres_num
        rm(registro)
        
        
        # # # # ESTO YA NO ES NECESARIO  # # # # 
        # ahora tenemos que formar el vector v que son todos los vectores originales que componen los clusters
        # rescatar los nombres del cluster en cuestion
        #clusters_con_mas_de_un_nodo <- proxima_fusion_numeric[proxima_fusion_numeric > 0]
        #for (n in 1:length(clusters_con_mas_de_un_nodo)) {
        #  temp <- c(nombres_num,  unlist(   mn[ clusters_con_mas_de_un_nodo[n]   ])  )
        #}
        # # # # ESTO YA NO ES NECESARIO  # # # # 
        
        #temp <- c(nombres_num, unlist(mn[el_cluster]) ) 
        #v <- temp[-hay_cluster_anterior]
        v <-  c( unlist(lista_nodos_de_clusters[n]), nombres_num[which(nombres_num < 0)] )
        
        
        # actualizacion del id, necesario para borrar las filas y columnas en D
        #id <- c( which(colnames(D) == proxima_fusion_numeric[1]) ,   which(colnames(D) == proxima_fusion_numeric[2])   )
        id <- c( which(colnames(D) == nombres_num[1]) ,   which(colnames(D) == nombres_num[2])   )
        
        # actualizacion de la distancia ultrametrica m
        m <- D[which(colnames(D)==nombres_num[1]), which(colnames(D)==nombres_num[2])]
        
        # actualizacion de mn
        mn <- c(mn, list( v ) )
      }
      
    }
    # # # # este es el cerebro de la funcion
    
    # calculo de la distancia MST
    d_mst <- find_mst_barrier(mst_g, v)
    # calculo de las distancias de acople
    d_cpl <- acople_distance_sum2(J, v)  ### > OJo que esta funcion esta en acople_distance_sum_function.R
    
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
    
    #ahora tengo que poner los valores en fila.columna de D.
    D <- put_dis(D, dit, nombres_, N, it)
    print(it)
  }
  
  # # # # LAST ITERATION # # # # 
  # Last iteration: when it = N - 1
  it <- it + 1
  rl <- cluster_find_name(D)
  id <- rl$id
  nombres <- rl$nombres
  m <- rl$m
  nombres_num <- as.numeric(nombres) # nuevo
  clusters_con_mas_de_un_nodo <- nombres_num[nombres_num > 0]
  if (length(clusters_con_mas_de_un_nodo ) > 1) { # En caso que se vayan a fusionar dos clusters que a su vez contienen varios nodos cada uno.
    temp <- c()
    for (n in 1:length(clusters_con_mas_de_un_nodo)) {
      temp <- c(temp,  unlist(   mn[ clusters_con_mas_de_un_nodo[n] ] )  )
    }
    v <- temp
  } else {
    temp <- c()
    for (n in 1:length(clusters_con_mas_de_un_nodo)) {
      temp <- c(nombres_num[nombres_num < 1],  unlist(   mn[ clusters_con_mas_de_un_nodo ] )  )
    }
    #temp <- c(nombres_num, unlist(mn[el_cluster]) ) 
    v <- temp
  }
  # calculo de la distancia MST
  d_mst <- find_mst_barrier(mst_g, v)
  # calculo de las distancias de acople
  d_cpl <- acople_distance_sum2(J, v)  ### > OJo que esta funcion esta en acople_distance_sum_function.R
  merge[it,] <- as.numeric(c(nombres_num, it, m, d_mst, d_cpl)) 
  # # # # LAST ITERATION # # # # 
  
  # retorno de los resultados matriz merge
  return(merge)
  
}
# ejemplo:
#merge <- hierarchical_clustering_v3(D)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 17-oct-19
# Esta version es la misma que la hierarchical_clustering_v3, pero en este caso, cuando hay dos 
# cluster que se van a fusionar en donde los dos tienen más de 1 nodo, 
# se buscar para cada cluster, otro nodo u otro cluster que logre la menor 
# distancia de acople. En la version hierarchical_clustering_v3, cuando teníamos este caso, los dos 
# cluster se unían sin buscar minima distancia de acople.
# nota: 131019: se incorpora como input el mst_g y matriz de acople J
hierarchical_clustering_v4 <- function(D, J, mst_g) {
  # # # Inicializacion de la matriz merge que es el output
  N <- ncol(D)
  merge <- matrix(NA, ncol=6, nrow=N-1) # col1 y col2 spins merged, col3 numero de cluster, col4 = merge distance
  merge <- as.data.frame(merge)
  colnames(merge) <- c("node1", "node2", "cluster", "dultr", "dmst", "dc")
  iteraciones <- N - 2
  mn <- list() # aqui iremos agregando los nombres de los spins que componen cada clusters 
  
  for (it in 1:iteraciones) {
    #comenzamos a iterar
    rl <- cluster_find_name(D)
    id <- rl$id # 
    nombres <- rl$nombres  # 
    m <- rl$m  #
    # leaves (original spins) will be negative, merged clusters will be positive starting from 1, 2, ...
    nombres_num <- as.numeric(nombres) # nuevo
    
    # segundo, ahora chequeamos si el par de cluster identificados, son 2 nodos, o 
    # hay mas de 2 nodos.
    # # # # este es el cerebro de la funcion
    te <- nombres_num
    hay_cluster_anterior <- which(te > 0)
    if (is.integer0(hay_cluster_anterior)) { # los dos cluster son solo de 1 nodo cada uno
      mn <- c(mn, list( as.character(te )) )
      v <- te
    } else {  # los dos cluster involucran mas de 2 nodos
      # aqui comienza la busqueda de alternativas de fusion que MINIMIZAN dc (distancia de acople)
      num_cluster_anteriores <- length(hay_cluster_anterior) #numero de cluster que tienen mas de 2 nodos
      if (num_cluster_anteriores > 1) { # en caso que se proponga la fusion de dos clusters cada uno con mas de 1 nodos
        # haremos que se fusionen incondicionalmente.
        # # # # # # # Aqui viene lo nuevo 16-oct-19  # # # # # # #  # # # # # # #  # # # # # # #  # # # # # # # 
        # Step 1: identificar nodos de cada cluster
        nodos_cluster1 <- mn[[te[1]]]
        nodos_cluster2 <- mn[[te[2]]]
        
        # step 2: determinar nodos adyacentes de cada cluster
        nodos_ady_cluster1 <- find_ady_dc(v=nodos_cluster1, n=te[1], mst_g = mst_g, J = J)
        nodos_ady_cluster2 <- find_ady_dc(v=nodos_cluster2, n=te[2], mst_g = mst_g, J = J)
        
        # Step 3: extraer nodos adyacentes de cluster 1 y 2, son considerar los nodos
        # del propio cluster 1 y 2.
        id_drop <- which(nodos_ady_cluster1[,2] %in% as.numeric(nodos_cluster2))
        nodos_ady_cluster1 <- nodos_ady_cluster1[-id_drop,]
        id_drop <- which(nodos_ady_cluster2[,2] %in% as.numeric(nodos_cluster1))
        nodos_ady_cluster2 <- nodos_ady_cluster2[-id_drop,]
        nodos_ady <- rbind(nodos_ady_cluster1, nodos_ady_cluster2)
        
        # Step 4: deterinar el proximo nodo a fusionar eligiendo aquel que tenma menor distancia de acople
        id_sel <- which.min(nodos_ady[,3])
        proxima_fusion_tentativa <- as.numeric(nodos_ady[id_sel, c(1,2)])
        
        # Step 5: verificar si nodo seleccionado para fusion pertenece a otro cluster
        nodo_a_fusionar <- as.character(proxima_fusion_tentativa[2])
        for (n in 1:length(mn) ) {
          if ( nodo_a_fusionar %in% mn[[n]] ) {
            proxima_fusion <- c(proxima_fusion_tentativa[1], n)
            break
          } else { proxima_fusion  <- proxima_fusion_tentativa }
        }
        
        # Step 6: determinar todos los nodos de los clusters fusionados
        v <- c( mn[[ proxima_fusion[1] ]], proxima_fusion[2])
        mn <- c(mn, list( as.character(v )) )
        nombres_num <- as.numeric(proxima_fusion)
        
        # step 7: 
        # actualizacion del id, necesario para borrar las filas y columnas en D
        id <- c( which(colnames(D) == nombres_num[1]) ,   which(colnames(D) == nombres_num[2])   )
        
        # actualizacion de la distancia ultrametrica m
        m <- D[which(colnames(D)==nombres_num[1]), which(colnames(D)==nombres_num[2])]
        
        # # # # # ## # ## # ## # ## # ## # ## # #lo antiguo
        # lista_nodos_de_clusters <- list()
        # for (n in 1:num_cluster_anteriores) { # aqui guardamos los nodos de los clusters con mas de 2 nodos
        #   el_cluster <- te[n]
        #   #lista_nodos_de_clusters  <- c(lista_nodos_de_clusters, unlist(mn[el_cluster]) )
        #   lista_nodos_de_clusters <- c(lista_nodos_de_clusters, (mn[el_cluster]) )
        # }
        # v <- unlist(lista_nodos_de_clusters)
        # mn <- c(mn, list( as.character(v )) )
        # # # FIN FIN FIN FIN FIN FIN FIN FIN FIN FIN FIN FIN 
      } else { # en caso que se proponga la fusion de dos clusters, en donde uno de ellos tenga solo un nodo
        # comenzaremos a probar si existe otro nodo que minimice la distancia de acople.
        num_cluster_anteriores <- length(hay_cluster_anterior) #numero de cluster que tienen mas de 2 nodos
        lista_nodos_de_clusters <- list()
        for (n in 1:num_cluster_anteriores) { # aqui guardamos los nodos de los clusters con mas de 2 nodos
          #el_cluster <- te[n]   te[which(te > 0)]
          el_cluster <- te[which(te > 0)]
          #lista_nodos_de_clusters  <- c(lista_nodos_de_clusters, unlist(mn[el_cluster]) )
          lista_nodos_de_clusters <- c(lista_nodos_de_clusters, (mn[el_cluster]) )
        }
        registro <- matrix(NA, ncol=3, nrow=0)
        colnames(registro) <- c("cluster_n", "elnodo", "dist_acople") #tengo que ir registrando el cluster (n), el elnodo y distancia_acople entre cluster y elnodo
        for (n in 1:num_cluster_anteriores) { # ahora buscamos para cada cluster con mas de dos nodos, 
          # sus nodos adyacentes en el MST que tienen la mayor distancia de acople
          #registro_interno <- find_ady_dc(v=unlist(lista_nodos_de_clusters[n]), n=te[n])
          registro_interno <- find_ady_dc(v=unlist(lista_nodos_de_clusters[n]), n=el_cluster, mst_g = mst_g, J = J)
          registro <- rbind(registro, registro_interno)
        }
        # veamos cual seria la proxima fusion de acuerdo a la minima distancia de acople
        id2 <- which.min(registro[,3])
        proxima_fusion <- registro[id2, c(1,2)]
        
        #ahora tendriamos que ver si el nodo que esta en registro[id, 2] pertenece ya a otro cluster mas grande
        #nombres_num <- final_fusion(nombre_columnas_D=colnames(D), proxima_fusion=proxima_fusion, elregistro=registro, id=id2) # esto ya no es necesario
        nombres_num <- as.numeric(proxima_fusion)
        #proxima_fusion_numeric <- as.numeric(nombres_num)
        #nombres <- nombres_num
        rm(registro)
        v <- c( unlist(lista_nodos_de_clusters[n]), nombres_num[which(nombres_num < 0)] )
        
        # actualizacion del id, necesario para borrar las filas y columnas en D
        #id <- c( which(colnames(D) == proxima_fusion_numeric[1]) ,   which(colnames(D) == proxima_fusion_numeric[2])   )
        id <- c( which(colnames(D) == nombres_num[1]) ,   which(colnames(D) == nombres_num[2])   )
        
        # actualizacion de la distancia ultrametrica m
        m <- D[which(colnames(D)==nombres_num[1]), which(colnames(D)==nombres_num[2])]
        
        # actualizacion de mn
        mn <- c(mn, list( v ) )
      }
    }
    # # # # este es el cerebro de la funcion
    
    # calculo de la distancia MST
    d_mst <- find_mst_barrier(mst_g, v)
    # calculo de las distancias de acople
    d_cpl <- acople_distance_sum2(J, v)  ### > OJo que esta funcion esta en acople_distance_sum_function.R
    
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
    
    #ahora tengo que poner los valores en fila.columna de D.
    D <- put_dis(D, dit, nombres_, N, it)
    print(it)
  }
  
  # # # # LAST ITERATION # # # # 
  # Last iteration: when it = N - 1
  it <- it + 1
  rl <- cluster_find_name(D)
  id <- rl$id
  nombres <- rl$nombres
  m <- rl$m
  nombres_num <- as.numeric(nombres) # nuevo
  clusters_con_mas_de_un_nodo <- nombres_num[nombres_num > 0]
  if (length(clusters_con_mas_de_un_nodo ) > 1) { # En caso que se vayan a fusionar dos clusters que a su vez contienen varios nodos cada uno.
    temp <- c()
    for (n in 1:length(clusters_con_mas_de_un_nodo)) {
      temp <- c(temp,  unlist(   mn[ clusters_con_mas_de_un_nodo[n] ] )  )
    }
    v <- temp
  } else {
    temp <- c()
    for (n in 1:length(clusters_con_mas_de_un_nodo)) {
      temp <- c(nombres_num[nombres_num < 1],  unlist(   mn[ clusters_con_mas_de_un_nodo ] )  )
    }
    #temp <- c(nombres_num, unlist(mn[el_cluster]) ) 
    v <- temp
  }
  # calculo de la distancia MST
  d_mst <- find_mst_barrier(mst_g, v)
  # calculo de las distancias de acople
  d_cpl <- acople_distance_sum2(J, v)  ### > OJo que esta funcion esta en acople_distance_sum_function.R
  merge[it,] <- as.numeric(c(nombres_num, it, m, d_mst, d_cpl)) 
  
  # retorno de los resultados matriz merge
  return(merge)
}
# ejemplo:
#merge <- hierarchical_clustering_v4(D, J=J, mst_g)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 




# # # # # # # # # # # # # # # # # # NECESARY FUNCTIONS # # # # # # # # # # # # # # # # # # # #
# Funciones necesarioa desarrolladas en testing_hierarchica_clusteringV7b.R para ejecutar
# la funcion hierarchical_clustering_v3
vector.is.empty <- function(x) return(length(x) ==0 )
# # funcion que  verifica si el nodo que esta en registro[id, 2] pertenece ya a otro cluster mas grande
# output: vector "nombre" con los nombres definitivos de los cluster que van a fusionar
final_fusion <- function(nombre_columnas_D, proxima_fusion, elregistro, id) {
  esta_elnodo <- elregistro[id, 2] %in% nombre_columnas_D
  if (!esta_elnodo) { #si no esta el nodo, hay que buscarlo
    cluster_del_nodo <- which(registro[id, 2] %in% mn)
    fusion_sera_entre <- c(proxima_fusion[1], cluster_del_nodo)
  } else {
    fusion_sera_entre <- c(proxima_fusion[1], proxima_fusion[2])
    nombres <- as.character(fusion_sera_entre)
  }
  return(nombres)
}
# ejemplo
# nombre <- final_fusion(nombre_columnas_D=colnames(D), proxima_fusion=proxima_fusion, elregistro=registro)
# # funcion que para un cluster en analisis, nos dice sus nodos adyacentes en el MST y la 
# # distancia de acople dc = f(sumatoria de las energias de acoples entre pares de nodos)
# output: matriz con tres columnas: el cluster bajo analisis, el nodo adyacente y su distancia de acople dc
# find_ady_dc <- function(v, n) {
#   registro_interno <- matrix(NA, ncol=3, nrow=0)
#   colnames(registro_interno) <- c("cluster_n", "elnodo", "dist_acople") #tengo que ir registrando el cluster (n), el elnodo y distancia_acople entre cluster y elnodo
#   # sus nodos adyacentes en el MST que tienen la menor distancia de acople
#   ady <- adjacent_vertices(mst_g, as.character(v))
#   nodoss <- length(ady) # numero de nodos adyacentes al cluster 
#   for (m in 1:nodoss) {
#     elnodo <- setdiff(as.character(-1*as.numeric(ady[[m]]) ) , as.character(v)) # encontramos el nodo que esta en ady[[.]] pero no en v
#     if (vector.is.empty(elnodo) ) { next } # en caso que elnodo sea vacio
#     #elnodo <- -1*as.numeric(elnodo)
#     for (x in 1:length(elnodo) ) {
#       nodo <- as.numeric(elnodo[x])
#       distancia_acople <- acople_distance_sum2(J, c(v, nodo))  ### > OJo que esta funcion esta en acople_distance_sum_function.R
#       registro_interno <- rbind(registro_interno, c(n, nodo, distancia_acople))
#     }
#   }
#   return(registro_interno)
# }
# ejemplo
# registro_interno <- find_ady_dc(v=lista_nodos_de_clusters[n], n=n)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
# 13.sep.19 FUNCION 
# Funcion que calcula la suma de distancias de acoples
# Por ejemplo, para sumar acoples a, b y c, se hace: f(a) + f(b) + f(c)
# donde f(*) es funcion que comvierte el acoples en una medidda de distancia.
# por ejemplo f(x) = sqrt(3 - x). 
# input: J = matriz de acoples J
# input: y = vector con todos los nombres de los spins 
# output: suma de f(a) + f(b) + f(c) + .... Si son N spins la suma sobre n(n-1)/2 acoples.
# acople_distance_sum <- function(J, y) {
#   sumdist <- 0
#   combinaciones <- combn(y,2)
#   lk <- ncol(combinaciones)
#   for (i in 1:lk) {
#     idx <- combinaciones[, i]
#     cupl <- J[idx[1], idx[2]]
#     dcupl <- sqrt(3 - cupl)
#     sumdist <- sumdist + dcupl
#   }
#   return(sumdist)
# }
# ejemplo
# si y = c("-4" "-2" "-1")
# sumdist <- acople_distance_sum(J, y)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
# FUNCION para detectar si el largo de un vector es nulo
is.integer0 <- function(x){
  is.integer(x) && length(x) == 0L
}
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
# # # # # # # # # # # # # # # # # # NECESARY FUNCTIONS # # # # # # # # # # # # # # # # # # # #











# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 27-ago-19
# Funcion que transforma la informacion del output de la funcion hierarchical_clustering
# en un objeto de dendograma para poder graficarlo.
# INPUTS: 
# 1. Matriz D de DISTANCIAS de N X N (N=num de spins o elementos) formato matrix. 
#             La diagonal debe venir con inf. 
#             Si la matriz D no viene con nombre de columnas, entonces hay que ponerselas.
#             Es importante que los nonbres de las columnas sean numeros negativos.
# 2. merge: dataframe que es el output de hierarchical_clustering
# 3. enames: vector con los combres que quiero ponerle a cada elemento. El largo debe ser igual
#           al numero de columnas o numero de filas de D
# OUTPUT
# hc: lista con 4 objetos con toda la información  del dendograma que necesita R para plotear
to_dendo <- function(D, merge, enames) {
  # 1
  # para genera rl hc$order, necesitamos identificar de merge[,c(1,2)] todos los
  # nodos negativos y ordenarlos en un vector en su orden de aparicion.
  temp <- (merge[,c(1,2)])
  temp <- as.numeric(rbind(temp$node1, temp$node2))
  temp <- temp[temp < 0]
  temp <- -1*temp
  # # # # # # # # # # 
  
  hc <- list()  # initialize empty object
  # define merging pattern: 
  #    negative numbers are leaves, 
  #    positive are merged clusters (defined by row number in $merge)
  hc$merge <- as.matrix(merge[,c(1,2)])
  hc$height <- as.vector(merge$dultr)    # define merge heights
  hc$order <- temp   
  #hc$labels <- LETTERS[1:7]    # labels of leaves -----#poner nombres originales
  hc$labels <- enames
  class(hc) <- "hclust"        # make it an hclust object
  dd <- as.dendrogram(hc)
  
  temp2 <- order.dendrogram(dd)  # esto es iportante porque nos permite
  # colocar en orden los merges para poder plotear.
  
  # Luego repetimos el proceso nuevamente con temp2.
  a <- list()  # initialize empty object
  # define merging pattern: 
  #    negative numbers are leaves, 
  #    positive are merged clusters (defined by row number in $merge)
  a$merge <- as.matrix(merge[,c(1,2)])
  a$height <- as.vector(merge$dultr)    # define merge heights
  a$order <- temp2   
  #hc$labels <- LETTERS[1:7]    # labels of leaves -----#poner nombres originales
  #a$labels <- c(1:ncol(D))
  a$labels <- enames
  class(a) <- "hclust"        # 
  
  return(a)
}
# EJEMPLO
# hc <- to_dendo(D, merge, enames=c(1:nolc(D)) )
# plot(hc)
# # si fuese necesario para dexdend package:
# hc <- as.dendogram(hc)


# 05-nov-19
# funcion creada para hacer el dendograma, version mejorada y que sirve 
# tambien para hacer el dendograma sobre la base de cualquier distancia
# es importante que en merge, la distancia que se utilziara para el dendograma
# por lo general la ultrametrica, venga en la ultima columna de merge
to_dendo2 <- function(merge, enames) {
  col = ncol(merge)
  # 1
  # para genera rl hc$order, necesitamos identificar de merge[,c(1,2)] todos los
  # nodos negativos y ordenarlos en un vector en su orden de aparicion.
  temp <- (merge[,c(1,2)])
  temp <- as.numeric(rbind(temp$node1, temp$node2))
  temp <- temp[temp < 0]
  temp <- -1*temp
  # # # # # # # # # # 
  
  hc <- list()  # initialize empty object
  # define merging pattern: 
  #    negative numbers are leaves, 
  #    positive are merged clusters (defined by row number in $merge)
  hc$merge <- as.matrix(merge[,c(1,2)])
  hc$height <- as.vector(merge[, col])    # define merge heights
  hc$order <- temp   
  #hc$labels <- LETTERS[1:7]    # labels of leaves -----#poner nombres originales
  hc$labels <- enames
  class(hc) <- "hclust"        # make it an hclust object
  dd <- as.dendrogram(hc)
  
  temp2 <- order.dendrogram(dd)  # esto es iportante porque nos permite
  # colocar en orden los merges para poder plotear.
  
  # Luego repetimos el proceso nuevamente con temp2.
  a <- list()  # initialize empty object
  # define merging pattern: 
  #    negative numbers are leaves, 
  #    positive are merged clusters (defined by row number in $merge)
  a$merge <- as.matrix(merge[,c(1,2)])
  a$height <- as.vector(merge[, col])   # define merge heights
  a$order <- temp2   
  #hc$labels <- LETTERS[1:7]    # labels of leaves -----#poner nombres originales
  #a$labels <- c(1:ncol(D))
  a$labels <- enames
  class(a) <- "hclust"        # 
  
  return(a)
}
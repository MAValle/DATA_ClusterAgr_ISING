# hacemos algo muy similar que en testing_hierarchical_clusteringV5.R, utilizando
# la funcion hierarchical_clustering_v2, pero en forma modificada de acuerdo 
# a lo expresado en pag. 243, 244 y 250 de los apuntes.

# La idea consiste en unir todos los clusters  clusters de un solo spin o nodo, basado
# en al distancia mst (dmst) con single linkage (igual que siempre), pero a partir
# de la busqueda de un tercer nodo o cluster para unir, se busca el que minimice
# la energia de acople, siempre y cuando haya una conexion a traves del MST.
# Esto nos forzara a unir clusters que no son necesariamente de minima distancia 
# de mst, pero si de minima energia de acople, pero al menos que estan unidos
# via el MST y en consecuencia de baja distancia MST.

# Que diferencia hay entre este algoritmo y el single_linkage de 
# testing_hierarchical_clusteringV5.R?
# En escencia es lo mismo, solo que en la primera parte cuando encontramos el par de 
# clusters que la distancia minima (el par de cluster (i,j) con menor distancia), 
# es decir encontramos d = min d[(i), (j)] entre el cluster i y j, comenzamos
# a buscar para cluster i y j, otros nodos adyacentes a ellos que tienen
# la menor distancia de acople (ir, Ec energia de acople). Puede que como 
# resultado de esta busqueda, el cluster i y j ya no se unan, y ahora se
# una el i con clusrer k, donde k es adyacente a cluster i, pero no a j.
# Luego, el proceso sigue de manera normal.


# actual name: testing_hierarchical_clusteringV7.R
# 30.sep.19


# # # # # # # # # # # # # # # # # # # # DATA LOADING # # # # # # # # # # # # # # # # # # # # #

# Carga de datos, 
# igual que en clusteragrr_v1.R
rm(list = ls())
load('new_inferring_parameters_environment270219.RData') # load H and J parameters ising inferred
source("functions_hclust.R")
#source("network_functions_v1.R") # aqui nos interesa la funcion find_mst_barrier
source("find_mst_barrier_function.R")
source("acople_distance_sum_function.R")
source("is_integer0_function.R")
source("entropies_functions.R") # nos interesa aqui get_energy function
library(igraph)

# # # # # # # # # # # # # # # # # # NECESARY FUNCTIONS # # # # # # # # # # # # # # # # # # # #
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
# # funcion que para un cluster en analisis, nos dice sus nodos adyacentes en el MST y el 
# # que produce la mayor distancia de acople (que en realidad es la menor energi de acople)
# output: matriz con tres columnas: el cluster bajo analisis, el nodo adyacente y su distancia de acople dc
find_ady_dc <- function(v, n) {
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
# registro_interno <- find_ady_dc(v=lista_nodos_de_clusters[n], n=n)
# # # # # # # # # # # # # # # # # # NECESARY FUNCTIONS # # # # # # # # # # # # # # # # # # # #


# # # # # # # # # # # # # # # # # # # # DATA LOADING # # # # # # # # # # # # # # # # # # # # #
# de los 25 nodos de J, seleccionemos los primeros 6
#set.seed(123)
# seleccionemos 5 nodos cualquiera
nd <- c(4,8,12,17,19,22) # 01oct19
nd <- c(2,4,6,8,12,15,17,19, 22) # 09oct19
J <- J[nd, nd]
# # # # # # # # # # # # # # # # # # # # DATA LOADING # # # # # # # # # # # # # # # # # # # # #



# # # # # # # # # # # # # # # # # # # # # FIND MST # # # # # # # # # # # # # # # # # # # # # #
# Conforming the COUPLING Network
colnames(J) <- colnames(wb[, nd])
rownames(J) <- colnames(wb[, nd])

# 12-sep-19
# Vamos a transformar los nombres de los nodos originales de una vez en nombres para el dendograma
colnames(J) <- rownames(J) <- -c(1:6) # ojo -c(1:Nn) tiene el mismo orden que colnames(dis)

#http://www.shizukalab.com/toolkits/sna/plotting-networks-pt-2
net <- graph_from_adjacency_matrix(as.matrix(J), mode="upper", weighted=TRUE, diag=FALSE)
#E(net)$coupling <- E(net)$weight
# Conforming the MST Network
dis <- sqrt(-J + 3) 
# Convertir la matriz de distancia en objeto igraph
g <- graph_from_adjacency_matrix(dis, mode = "undirected", weighted = TRUE, diag = FALSE, add.colnames = NULL, add.rownames = NA)
E(g)$coupling <- E(net)$weight # asignamos las energias de acoples de net a E(g)$coupling
mst_g <- minimum.spanning.tree(g, algorithm="prim")
edg <- as_edgelist(mst_g, names = TRUE)
# # # # # # # # # # # # # # # # # # # # # FIND MST # # # # # # # # # # # # # # # # # # # # # #


# # # # # # # # # # # # # # # # # # #  PREPARATION # # # # # # # # # # # # # # # # # # # # # #
D <- dis
Nn <- ncol(D)
diag(D) <- rep(Inf, Nn)
colnames(D) <- rownames(D) <- -c(1:Nn) # ojo -c(1:Nn) tiene el mismo orden que colnames(dis)
# equivalencia entre nombre de nodos reales y nombres en el dendograma:
names_spins <- data.frame(real = colnames(J), dend = colnames(D))
# por ejemplo: encontrar las filas en names_spins que calzan con los valores en nombres.
# ii = which(names_spins$dend %in% nombres)
# names_spins[ii, "real"] # aqui obtenemos los nombres reales de los spins.
# # # # # # # # # # # # # # # # # # #  PREPARATION # # # # # # # # # # # # # # # # # # # # # #


# # # # # # # # # # # # # # # # # # # # INICIO # # # # # # # # # # # # # # # # # # # # # # # #

# VAMOS A EJECUTAR AL FUNCION hierarchical_clustering por PEDAZOS.
# # # PARTE I
N <- ncol(D)
merge <- matrix(NA, ncol=6, nrow=N-1) # col1 y col2 spins merged, col3 numero de cluster, col4 = merge distance
merge <- as.data.frame(merge)
colnames(merge) <- c("node1", "node2", "cluster", "dultr", "dmst", "dc")
iteraciones <- N - 2
mn <- list() # aqui iremos agregando los nombres de los spins que componen cada clusters 

# Inicio de las iteraciones que son N - 2
# # # # # # # # # # # # # # # # # # # #  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>iteracion 1
# Primero encontramos el par de cluster con minima distancia
it = 1 # iteracion iteracion iteracion iteracion iteracion iteracion iteracion 
rl <- cluster_find_name(D)
id <- rl$id # (4 2)
nombres <- rl$nombres  #  "-4" "-2"
m <- rl$m  # 1.261858
# leaves (original spins) will be negative, merged clusters will be positive starting from 1, 2, ...
nombres_num <- as.numeric(nombres) # nuevo

# segundo, ahora chequeamos si el par de cluster identificados, son 2 nodos, o 
# hay mas de 2 nodos.
te <- nombres_num
hay_cluster_anterior <- which(te > 0)
if (is.integer0(hay_cluster_anterior)) { # los dos cluster son solo dos nodos
  mn <- c(mn, list( as.character(te )) )
  v <- te
} else {  # los dos cluster involucran mas de 2 nodos
  # aqui comienza la busqueda de alternativas de fusion que MINIMIZAN dc (distancia de acople)
  num_cluster_anteriores <- length(hay_cluster_anterior) #numero de cluster que tienen mas de 2 nodos
  lista_nodos_de_clusters <- list()
  for (n in 1:num_cluster_anteriores) { # aqui guardamos los nodos de los clusters con mas de 2 nodos
    el_cluster <- te[n]
    #lista_nodos_de_clusters  <- c(lista_nodos_de_clusters, unlist(mn[el_cluster]) )
    lista_nodos_de_clusters <- c(lista_nodos_de_clusters, (mn[el_cluster]) )
  }
  registro <- matrix(NA, ncol=3, nrow=0)
  colnames(registro) <- c("cluster_n", "elnodo", "dist_acople") #tengo que ir registrando el cluster (n), el elnodo y distancia_acople entre cluster y elnodo
  for (n in 1:num_cluster_anteriores) { # ahora buscamos para cada cluster con mas de dos nodos, 
    # sus nodos adyacentes en el MST que tienen la menor distancia de acople
    registro_interno <- find_ady_dc(v=unlist(lista_nodos_de_clusters[n]), n=te[n])
    registro <- cbind(registro_interno)
  }
  # veamos cual seria la proxima fusion de acuerdo a la minima distancia de acople
  id2 <- which.max(registro[,3])
  proxima_fusion <- registro[id2, c(1,2)]
  
  
  #ahora tendriamos que ver si el nodo que esta en registro[id, 2] pertenece ya a otro cluster mas grande
  nombres_num <- final_fusion(nombre_columnas_D=colnames(D), proxima_fusion=proxima_fusion, elregistro=registro, id=id2)
  proxima_fusion_numeric <- as.numeric(nombres_num)
  nombres <- nombres_num
  
  # ahora tenemos que formar el vector v que son todos los vectores originales que componen los clusters
  # rescatar los nombres del cluster en cuestion
  clusters_con_mas_de_un_nodo <- proxima_fusion_numeric[proxima_fusion_numeric > 0]
  for (n in 1:length(clusters_con_mas_de_un_nodo)) {
    temp <- c(nombres_num,  unlist(   mn[ clusters_con_mas_de_un_nodo[n]   ])  )
  }
  #temp <- c(nombres_num, unlist(mn[el_cluster]) ) 
  v <- temp[-hay_cluster_anterior]
  
  # actualizacion del id, necesario apra borrar las filas y columnas en D
  id <- c( which(colnames(D) == proxima_fusion_numeric[1]) ,   which(colnames(D) == proxima_fusion_numeric[2])   )
  
  # actualizacion de mn
  mn <- c(mn, list( v ) )
  
}



# calculo de la distancia MST
d_mst <- find_mst_barrier(mst_g, nombres_num)
# calculo de las distancias de acople
d_cpl <- acople_distance_sum2(J, nombres_num)  ### > OJo que esta funcion esta en acople_distance_sum_function.R

merge[it,] <- as.numeric(c(nombres_num, it, m, d_mst, d_cpl)) #  nuevo: se coloca nodo1 (fusion) / nodo2 (fusion) /nombre del nuevo cluster / distancia a la que se unen


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


#ahora tengo que poner los valores en fila.columna de D.
D <- put_dis(D, dit, nombres_, N, it)
print(it)








# # # # # # # # # # # # # # # # # # # #  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>iteracion 2
it = it+1 # iteracion iteracion iteracion iteracion iteracion iteracion iteracion 
# ahora viene el weveo de calcular gamma, las Ec y las Emst
rl <- cluster_find_name(D)
id <- rl$id
nombres <- rl$nombres
m <- rl$m
# leaves (original spins) will be negative, merged clusters will be positive starting from 1, 2, ...
nombres_num <- as.numeric(nombres) # nuevo


# segundo, ahora chequeamos si el par de cluster identificados, son 2 nodos, o 
# hay mas de 2 nodos.
te <- nombres_num
hay_cluster_anterior <- which(te > 0)
if (is.integer0(hay_cluster_anterior)) { # los dos cluster son solo dos nodos
  mn <- c(mn, list( as.character(te )) )
  v <- te
} else {  # los dos cluster involucran mas de 2 nodos
  # aqui comienza la busqueda de alternativas de fusion que MINIMIZAN dc (distancia de acople)
  num_cluster_anteriores <- length(hay_cluster_anterior) #numero de cluster que tienen mas de 2 nodos
  lista_nodos_de_clusters <- list()
  for (n in 1:num_cluster_anteriores) { # aqui guardamos los nodos de los clusters con mas de 2 nodos
    el_cluster <- te[n]
    #lista_nodos_de_clusters  <- c(lista_nodos_de_clusters, unlist(mn[el_cluster]) )
    lista_nodos_de_clusters <- c(lista_nodos_de_clusters, (mn[el_cluster]) )
  }
  registro <- matrix(NA, ncol=3, nrow=0)
  colnames(registro) <- c("cluster_n", "elnodo", "dist_acople") #tengo que ir registrando el cluster (n), el elnodo y distancia_acople entre cluster y elnodo
  for (n in 1:num_cluster_anteriores) { # ahora buscamos para cada cluster con mas de dos nodos, 
    # sus nodos adyacentes en el MST que tienen la menor distancia de acople
    registro_interno <- find_ady_dc(v=unlist(lista_nodos_de_clusters[n]), n=te[n])
    registro <- cbind(registro_interno)
  }
  # veamos cual seria la proxima fusion de acuerdo a la minima distancia de acople
  id2 <- which.max(registro[,3])
  proxima_fusion <- registro[id2, c(1,2)]
  
  
  #ahora tendriamos que ver si el nodo que esta en registro[id, 2] pertenece ya a otro cluster mas grande
  nombres_num <- final_fusion(nombre_columnas_D=colnames(D), proxima_fusion=proxima_fusion, elregistro=registro, id=id2)
  proxima_fusion_numeric <- as.numeric(nombres_num)
  nombres <- nombres_num
  
  # ahora tenemos que formar el vector v que son todos los vectores originales que componen los clusters
  # rescatar los nombres del cluster en cuestion
  clusters_con_mas_de_un_nodo <- proxima_fusion_numeric[proxima_fusion_numeric > 0]
  for (n in 1:length(clusters_con_mas_de_un_nodo)) {
    temp <- c(nombres_num,  unlist(   mn[ clusters_con_mas_de_un_nodo[n]   ])  )
  }
  #temp <- c(nombres_num, unlist(mn[el_cluster]) ) 
  v <- temp[-hay_cluster_anterior]
  
  # actualizacion del id, necesario apra borrar las filas y columnas en D
  id <- c( which(colnames(D) == proxima_fusion_numeric[1]) ,   which(colnames(D) == proxima_fusion_numeric[2])   )
  
  # actualizacion de mn
  mn <- c(mn, list( v ) )
  
}



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


#ahora tengo que poner los valores en fila.columna de D.
D <- put_dis(D, dit, nombres_, N, it)
print(it)





# # # # # # # # # # # # # # # # # # # #  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>iteracion 3
it = 3 # iteracion iteracion iteracion iteracion iteracion iteracion iteracion 
# ahora viene el weveo de calcular gamma, las Ec y las Emst
rl <- cluster_find_name(D)
id <- rl$id
nombres <- rl$nombres
m <- rl$m
# leaves (original spins) will be negative, merged clusters will be positive starting from 1, 2, ...
nombres_num <- as.numeric(nombres) # nuevo


# segundo, ahora chequeamos si el par de cluster identificados, son 2 nodos, o 
# hay mas de 2 nodos.
te <- nombres_num
hay_cluster_anterior <- which(te > 0)
if (is.integer0(hay_cluster_anterior)) { # los dos cluster son solo dos nodos
  mn <- c(mn, list( as.character(te )) )
  v <- te
} else {  # los dos cluster involucran mas de 2 nodos
  # aqui comienza la busqueda de alternativas de fusion que MINIMIZAN dc (distancia de acople)
  num_cluster_anteriores <- length(hay_cluster_anterior) #numero de cluster que tienen mas de 2 nodos
  lista_nodos_de_clusters <- list()
  for (n in 1:num_cluster_anteriores) { # aqui guardamos los nodos de los clusters con mas de 2 nodos
    el_cluster <- te[n]
    #lista_nodos_de_clusters  <- c(lista_nodos_de_clusters, unlist(mn[el_cluster]) )
    lista_nodos_de_clusters <- c(lista_nodos_de_clusters, (mn[el_cluster]) )
  }
  registro <- matrix(NA, ncol=3, nrow=0)
  colnames(registro) <- c("cluster_n", "elnodo", "dist_acople") #tengo que ir registrando el cluster (n), el elnodo y distancia_acople entre cluster y elnodo
  for (n in 1:num_cluster_anteriores) { # ahora buscamos para cada cluster con mas de dos nodos, 
    # sus nodos adyacentes en el MST que tienen la menor distancia de acople
    registro_interno <- find_ady_dc(v=unlist(lista_nodos_de_clusters[n]), n=te[n])
    registro <- cbind(registro_interno)
  }
  # veamos cual seria la proxima fusion de acuerdo a la minima distancia de acople
  id2 <- which.max(registro[,3])
  proxima_fusion <- registro[id2, c(1,2)]
  
  
  #ahora tendriamos que ver si el nodo que esta en registro[id, 2] pertenece ya a otro cluster mas grande
  nombres_num <- final_fusion(nombre_columnas_D=colnames(D), proxima_fusion=proxima_fusion, elregistro=registro, id=id2)
  proxima_fusion_numeric <- as.numeric(nombres_num)
  nombres <- nombres_num
  
  # ahora tenemos que formar el vector v que son todos los vectores originales que componen los clusters
  # rescatar los nombres del cluster en cuestion
  clusters_con_mas_de_un_nodo <- proxima_fusion_numeric[proxima_fusion_numeric > 0]
  for (n in 1:length(clusters_con_mas_de_un_nodo)) {
    temp <- c(nombres_num,  unlist(   mn[ clusters_con_mas_de_un_nodo[n]   ])  )
  }
  #temp <- c(nombres_num, unlist(mn[el_cluster]) ) 
  v <- temp[-hay_cluster_anterior]
  
  # actualizacion del id, necesario apra borrar las filas y columnas en D
  id <- c( which(colnames(D) == proxima_fusion_numeric[1]) ,   which(colnames(D) == proxima_fusion_numeric[2])   )
  
  # actualizacion de mn
  mn <- c(mn, list( v ) )
  
}


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


#ahora tengo que poner los valores en fila.columna de D.
D <- put_dis(D, dit, nombres_, N, it)
print(it)




















# # # # # # # # # # # # # # # # # # # #  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>iteracion 4
it = 4
# ahora viene el weveo de calcular gamma, las Ec y las Emst
rl <- cluster_find_name(D)
id <- rl$id
nombres <- rl$nombres
m <- rl$m
# leaves (original spins) will be negative, merged clusters will be positive starting from 1, 2, ...
nombres_num <- as.numeric(nombres) # nuevo


# segundo, ahora chequeamos si el par de cluster identificados, son 2 nodos, o 
# hay mas de 2 nodos.
te <- nombres_num
hay_cluster_anterior <- which(te > 0)
if (is.integer0(hay_cluster_anterior)) { # los dos cluster son solo dos nodos
  mn <- c(mn, list( as.character(te )) )
  v <- te
} else {  # los dos cluster involucran mas de 2 nodos
  # aqui comienza la busqueda de alternativas de fusion que MINIMIZAN dc (distancia de acople)
  num_cluster_anteriores <- length(hay_cluster_anterior) #numero de cluster que tienen mas de 2 nodos
  lista_nodos_de_clusters <- list()
  for (n in 1:num_cluster_anteriores) { # aqui guardamos los nodos de los clusters con mas de 2 nodos
    el_cluster <- te[n]
    #lista_nodos_de_clusters  <- c(lista_nodos_de_clusters, unlist(mn[el_cluster]) )
    lista_nodos_de_clusters <- c(lista_nodos_de_clusters, (mn[el_cluster]) )
  }
  registro <- matrix(NA, ncol=3, nrow=0)
  colnames(registro) <- c("cluster_n", "elnodo", "dist_acople") #tengo que ir registrando el cluster (n), el elnodo y distancia_acople entre cluster y elnodo
  for (n in 1:num_cluster_anteriores) { # ahora buscamos para cada cluster con mas de dos nodos, 
    # sus nodos adyacentes en el MST que tienen la menor distancia de acople
    registro_interno <- find_ady_dc(v=unlist(lista_nodos_de_clusters[n]), n=te[n])
    registro <- cbind(registro_interno)
  }
  # veamos cual seria la proxima fusion de acuerdo a la minima distancia de acople
  id2 <- which.max(registro[,3])
  proxima_fusion <- registro[id2, c(1,2)]
  
  
  #ahora tendriamos que ver si el nodo que esta en registro[id, 2] pertenece ya a otro cluster mas grande
  nombres_num <- final_fusion(nombre_columnas_D=colnames(D), proxima_fusion=proxima_fusion, elregistro=registro, id=id2)
  proxima_fusion_numeric <- as.numeric(nombres_num)
  nombres <- nombres_num
  
  # ahora tenemos que formar el vector v que son todos los vectores originales que componen los clusters
  # rescatar los nombres del cluster en cuestion
  clusters_con_mas_de_un_nodo <- proxima_fusion_numeric[proxima_fusion_numeric > 0]
  for (n in 1:length(clusters_con_mas_de_un_nodo)) {
    temp <- c(nombres_num,  unlist(   mn[ clusters_con_mas_de_un_nodo[n]   ])  )
  }
  #temp <- c(nombres_num, unlist(mn[el_cluster]) ) 
  v <- temp[-hay_cluster_anterior]
  
  # actualizacion del id, necesario apra borrar las filas y columnas en D
  id <- c( which(colnames(D) == proxima_fusion_numeric[1]) ,   which(colnames(D) == proxima_fusion_numeric[2])   )
  
  # actualizacion de mn
  mn <- c(mn, list( v ) )
  
}


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


#ahora tengo que poner los valores en fila.columna de D.
D <- put_dis(D, dit, nombres_, N, it)
print(it)



# # # # # # # # # # # # # # # # # # # #  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>iteracion 5
# Last iteration: when it = N - 1
it <- it + 1
rl <- cluster_find_name(D)
id <- rl$id
nombres <- rl$nombres
m <- rl$m
nombres_num <- as.numeric(nombres) # nuevo
clusters_con_mas_de_un_nodo <- nombres_num[nombres_num > 0]
for (n in 1:length(clusters_con_mas_de_un_nodo)) {
  temp <- c(nombres_num,  unlist(   mn[ clusters_con_mas_de_un_nodo[n]   ])  )
}
#temp <- c(nombres_num, unlist(mn[el_cluster]) ) 
v <- unique(temp[-hay_cluster_anterior])
# calculo de la distancia MST
d_mst <- find_mst_barrier(mst_g, v)
# calculo de las distancias de acople
d_cpl <- acople_distance_sum2(J, v)  ### > OJo que esta funcion esta en acople_distance_sum_function.R
merge[it,] <- as.numeric(c(nombres_num, it, m, d_mst, d_cpl)) 




# # # # # # # # # # # # # # # # # # # #  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>PLOTEO
Dinicial <- sqrt(-J + 3) 
hc <- to_dendo(D = Dinicial, merge = merge[,c(1:4)], enames = colnames(dis) )
plot(hc)







# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
# 13.sep.19 FUNCION 
# Funcion que calcula la suma de distancias de acoples
# Por ejemplo, para sumar acoples a, b y c, se hace: f(a) + f(b) + f(c)
# donde f(*) es funcion que comvierte el acoples en una medidda de distancia.
# por ejemplo f(x) = sqrt(3 - x). 
# input: J = matriz de acoples J
# input: y = vector con todos los nombres de los spins 
# output: suma de f(a) + f(b) + f(c) + .... Si son N spins la suma sobre n(n-1)/2 acoples.
acople_distance_sum <- function(J, y) {
  sumdist <- 0
  combinaciones <- combn(y,2)
  lk <- ncol(combinaciones)
  for (i in 1:lk) {
    idx <- combinaciones[, i]
    cupl <- J[idx[1], idx[2]]
    dcupl <- sqrt(3 - cupl)
    sumdist <- sumdist + dcupl
  }
  return(sumdist)
}
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

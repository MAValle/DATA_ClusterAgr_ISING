# Esta version es la misma que la V7b, pero en este caso, cuando hay dos 
# cluster que se van a fusionar en donde los dos tienen más de 1 nodo, 
# se buscar para cad acluster, otro nodo u otro cluster que logre la menor 
# distancia de acople. En la version V7b, cuando teníamos este caso, los dos 
# cluster se unían sin buscar minima distancia de acople.



# actual name: testing_hierarchical_clusteringV8.R
# 16.oct.19


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
# # funcion que para un cluster en analisis, nos dice sus nodos adyacentes en el MST y la 
# # distancia de acople dc = f(sumatoria de las energias de acoples entre pares de nodos)
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


# # # # # # # # # # # # # # # # # # # # CREATING J SUBSET # # # # # # # # # # # # # # # # # # # # #
# de los 25 nodos de J, seleccionemos los primeros 6
#set.seed(123)
# seleccionemos 5 nodos cualquiera
nd <- c(4,8,12,17,19,22) # 01oct19
nd <- c(2,4,6,8,12,15,17,19, 22) # 09oct19
nd <- c(1,4,11,13,20,25) # 10oct19
nd <- c( 1, 3, 7, 11, 13, 14, 15, 17, 25) # 11oct19a
nd <- c(1, 2, 6, 7, 8, 10, 11, 13, 14, 19, 20, 21) # 11oct19b
nd <- 1:25 # 11oct19c
# creando matriz J a mano oct, 13
J <- matrix(c(Inf, 2.8, 2.5, -10, -30,   2.8, Inf, -50, -25, -28,   2.5, -50, Inf, 1.5, -26,    -10, -25, 1.5, Inf, 1.8,   -30, -28, -26, 1.8, Inf), nrow = 5, ncol = 5)
# luego ejecutamos directamente la L107
J <- J[nd, nd]
# # # # # # # # # # # # # # # # # # # # CREATING J SUBSET # # # # # # # # # # # # # # # # # # # # #

# # # # # # # # # # # # # # # # # # # # CREATING A J MATRIX # # # # # # # # # # # # # # # # # # # # #
# creating artifitial matrix of couplings
Nn = 25 # number of spins
library(Matrix)
set.seed(134)
x <- Matrix(rnorm(Nn*Nn),Nn)  # los J vienen con media 0 y desv estandar 1
J <- forceSymmetric(x) # couplings
D <- sqrt(3 - as.matrix(J)) # distance
colnames(D) <- rownames(D) <- -c(1:Nn) # enumeramos todos las hojas o spines originales con valones negativos.
diag(D) <- rep(Inf, Nn)
colnames(J) <- rownames(J) <- -c(1:ncol(J)) 
# # # # # # # # # # # # # # # # # # # # CREATING A J MATRIX # # # # # # # # # # # # # # # # # # # # #




# # # # # # # # # # # # # # # # # # # # # FIND MST # # # # # # # # # # # # # # # # # # # # # #
# Conforming the COUPLING Network
# 12-sep-19
# Vamos a transformar los nombres de los nodos originales de una vez en nombres para el dendograma
colnames(J) <- rownames(J) <- -c(1:ncol(J)) # ojo -c(1:Nn) tiene el mismo orden que colnames(dis)

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
#names_spins <- data.frame(real = colnames(J), dend = colnames(D))
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
it = 0

#comenzamos a iterar
it = it + 1
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
    nodos_ady_cluster1 <- find_ady_dc(v=nodos_cluster1, n=te[1])
    nodos_ady_cluster2 <- find_ady_dc(v=nodos_cluster2, n=te[2])
    
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
      registro_interno <- find_ady_dc(v=unlist(lista_nodos_de_clusters[n]), n=el_cluster)
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






















# esto se hara una vez que tengamos la funcion hierarchical_clustering_v4
# the magic
merge2 <- hierarchical_clustering_v2(D)
hc <- to_dendo(dis, merge2[,c(1:4)], enames=c(1:ncol(D)) )
plot(hc)

merge2 <- hierarchical_clustering_v3(D)
hc <- to_dendo(D, merge, enames=c(1:ncol(D)) )
plot(hc)

merge4 <- hierarchical_clustering_v4(D)
hc <- to_dendo(D, merge4[,c(1:4)], enames=c(1:ncol(D)) )
plot(hc)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


















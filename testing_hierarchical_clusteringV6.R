# Ahora que ya tenemos en funtions_hlcust.R un clustering jerarquico basado en 
# single linkage funcionando, lo que comenzamos a intentar aqui es un 
# clustering de la misma forma, pero en que ahora la matriz de distancia se 
# vaya actualizando con distancias del MST multiplicadas por el factor gamma.
# Var pag. 242, 246.

# El factor gamma es la relacion que hay entre la energia de acoples del cluster
# y la energia de MST en recorrer lso nodos del cluster por el MST.
# Ec = sumatoria de todos los pares de acoples entre nodos del cluster que se 
# pretende formar.
# Emst es la energia de recorrer los nodos del cluster que se pretende formar
# a traves del MST.
# gamma = d(c)/d(mst)
# donde d(c) es la energia de acople transformada en distancia (dik + dik + ...)
# y d(mst) es la energia de mst transformada en distancia.
# es decir ahora la distancia ultrametrica d<' = d< X gamma.



# actual name: testing_hierarchical_clusteringV6.R
# 09.sep.19


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
# # # # # # # # # # # # # # # # # # # # DATA LOADING # # # # # # # # # # # # # # # # # # # # #
# de los 25 nodos de J, seleccionemos los primeros 6
#set.seed(123)
J <- J[c(1:6),(1:6)]
# # # # # # # # # # # # # # # # # # # # DATA LOADING # # # # # # # # # # # # # # # # # # # # #



# # # # # # # # # # # # # # # # # # # # # FIND MST # # # # # # # # # # # # # # # # # # # # # #
# Conforming the COUPLING Network
colnames(J) <- colnames(wb[, 1:6])
rownames(J) <- colnames(wb[, 1:6])

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
merge <- matrix(NA, ncol=4, nrow=N-1) # col1 y col2 spins merged, col3 numero de cluster, col4 = merge distance
merge <- as.data.frame(merge)
colnames(merge) <- c("node1", "node2", "cluster", "dultr")
iteraciones <- N - 2
mn <- list() # aqui iremos agregando los nombres de los spins que componen cada clusters 

# Inicio de las iteraciones que son N - 2
# # # # # # # # # # # # # # # # # # # #  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>iteracion 1
it = 1 # iteracion iteracion iteracion iteracion iteracion iteracion iteracion 
rl <- cluster_find_name(D)
id <- rl$id # (4 2)
nombres <- rl$nombres  #  "-4" "-2"
m <- rl$m  # 1.261858
# leaves (original spins) will be negative, merged clusters will be positive starting from 1, 2, ...
nombres_num <- as.numeric(nombres) # nuevo
merge[it,] <- as.numeric(c(nombres_num, it, m)) #  nuevo: se coloca nodo1 (fusion) / nodo2 (fusion) /nombre del nuevo cluster / distancia a la que se unen

# # # # # # # # # Upgrading the distance matrix
# fila id[1] y columna id[2] se deben borrar
Dold <- D
D <- D[-id, -id]

# # # # # # # # # poniendo los nombres reales de los spins que componen el cluster recien formado.
#ii = which(names_spins$dend %in% nombres) # indices de los nombres de lso spins del nuevo cluster
te <- c(merge[it,1], merge[it,2]  ) # nombres de los spins del nuevo cluster
# detectamos si un nombre de cluster es positivo, en tal caso, extraer los nombres de los spins de ese cluster
hay_cluster_anterior <- which(te > 0)
if(!is.integer0(hay_cluster_anterior) ) {
  el_cluster <- te[hay_cluster_anterior]
  # rescatar los nombres del cluster en cuestion
  te2 <- c(te, mn[[el_cluster]])
  te <- te2[-hay_cluster_anterior]
}
mn <- c(mn, list( as.character(te )) ) # y ahora en mn[[2]] cluster 2, se compone de -5 y -6


# AQUI EN DONDE TENGO QUE HACER LAS MODIFICACIONES, function find_min_dist
# la funcion find_min_dist(D, Dold, nombres) es para re-actualizar la matriz D, utilizando 
# single linkage method entre el cluster recien funsionado y los demas clusters.
#temp <- find_min_dist(D, Dold, nombres)   #$dit 1.696616 1.273414 1.753150 1.741162   $nombres_ [1] "-1" "-3" "-5" "-6"
#dit <- temp$dit # dit son las minimias distancias de cluster fusionado en la iteracion anterior a todos los demas clusters que quedan
#nombres_ = temp$nombres_ # nombres_ son los nombres de los nodos que estan a minima distancia.
# abriendo la funcion find_min_dist
dit <- vector()
nombres_ <- colnames(D) # 
n <- length(nombres_) 
for (i in 1:n) {
  nodo <- nombres_[i]
  v <- min(Dold[nombres, nodo])  # ojo que en "nombres" estan los nombres de los clusters que se fusionaron
  #dit <- c(dit,min(Dold[, nodo]))
  # vector con nombres de los spins involucrados
  y <- c(nombres, nodo)  # los spins o nodos del cluster recien formado con los nodos de los otros clusters
  # como "y" puede contener cluster de antes, hay que detectar los nombres de los spins de ese cluster:
  hay_cluster_anterior <- which(as.numeric(y) >0)
  if (!is.integer0(hay_cluster_anterior)  ) { # en caso que hayan clusters anteriores
    el_cluster <- as.numeric(y[hay_cluster_anterior])
    # rescatar los nombres del cluster en cuestion
    y2 <- c(y, unlist(mn[el_cluster]) ) 
    y <- y2[-hay_cluster_anterior]
  }
  # calculo de la distancia MST
  d_mst <- find_mst_barrier(mst_g, y)
  # calculo de las distancias de acople
  d_cpl <- acople_distance_sum(J, y)
  # calculo de gamma
  gamma <-  d_cpl/d_mst
  #
  vprima <- v*gamma
  dit <- c(dit,vprima )
}
# el resultado es:                            1.876042 1.932618 1.875510 2.755208
# el resultado original con find_min_dist es: 1.681610 1.273414 1.473493 1.473493

#ahora tengo que poner los valores en fila.columna de D.
D <- put_dis(D, dit, nombres_, N, it)
print(it)

# # # # # # # # # # # # # # # # # # # #  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>iteracion 2
it = 2 # iteracion iteracion iteracion iteracion iteracion iteracion iteracion 
# ahora viene el weveo de calcular gamma, las Ec y las Emst
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

# # # # # # # # # poniendo los nombres reales de los spins que componen el cluster recien formado.
te <- c(merge[it,1], merge[it,2]  ) # nombres de los spins del nuevo cluster
# detectamos si un nombre de cluster es positivo, en tal caso, extraer los nombres de los spins de ese cluster
hay_cluster_anterior <- which(te > 0)
if(!is.integer0(hay_cluster_anterior) ) {
  el_cluster <- te[hay_cluster_anterior]
  # rescatar los nombres del cluster en cuestion
  te2 <- c(te, mn[[el_cluster]])
  te <- te2[-hay_cluster_anterior]
}
mn <- c(mn, list( as.character(te )) ) # y ahora en mn[[2]] cluster 2, se compone de -5 y -6


# AQUI EN DONDE TENGO QUE HACER LAS MODIFICACIONES, function find_min_dist
# la funcion find_min_dist(D, Dold, nombres) es para re-actualizar la matriz D, utilizando 
# single linkage method entre el cluster recien funsionado y los demas clusters.
#temp <- find_min_dist(D, Dold, nombres)   #$dit 1.696616 1.273414 1.753150 1.741162   $nombres_ [1] "-1" "-3" "-5" "-6"
#dit <- temp$dit # dit son las minimias distancias de cluster fusionado en la iteracion anterior a todos los demas clusters que quedan
#nombres_ = temp$nombres_ # nombres_ son los nombres de los nodos que estan a minima distancia.
# abriendo la funcion find_min_dist
dit <- vector()
nombres_ <- colnames(D) # 
n <- length(nombres_) # 
for (i in 1:n) {
  nodo <- nombres_[i]
  v <- min(Dold[nombres, nodo])  # ojo que en "nombres" estan los nombres de los clusters que se fusionaron
  #dit <- c(dit,min(Dold[, nodo]))
  # calculo de la distancia MST
  y <- c(nombres, nodo)  # los spins o nodos del cluster recien formado con los nodos de los otros clusters
  d_mst <- find_mst_barrier(mst_g, y)
  # calculo de las distancias de acople
  # como y puede contener cluster de antes, hay que detectar los nombres de los spins de ese cluster:
  hay_cluster_anterior <- which(as.numeric(y) >0)
  if (!is.integer0(hay_cluster_anterior)  ) { # en caso que hayan clusters anteriores
    el_cluster <- as.numeric(y[hay_cluster_anterior])
    # rescatar los nombres del cluster en cuestion
    y2 <- c(y, mn[[el_cluster]])
    y <- y2[-hay_cluster_anterior]
  }
  d_cpl <- acople_distance_sum(J, y)
  # calculo de gamma
  gamma <-  d_cpl/d_mst
  #
  vprima <- v*gamma
  dit <- c(dit,vprima )
}
# el resultado es:                            1.443342  1.966347 12.413364
# el resultado original con find_min_dist es: 1.767269 1.762505 1.875510

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
merge[it,] <- as.numeric(c(nombres_num, it, m)) #  nuevo: se coloca nodo1 (fusion) / nodo2 (fusion) /nombre del nuevo cluster / distancia a la que se unen

# # # # # # # # # Upgrading the distance matrix
# fila id[1] y columna id[2] se deben borrar
Dold <- D
D <- D[-id, -id]

# # # # # # # # # poniendo los nombres reales de los spins que componen el cluster recien formado.
#ii = which(names_spins$dend %in% nombres) # indices de los nombres de lso spins del nuevo cluster
te <- c(merge[it,1], merge[it,2]  ) # nombres de los spins del nuevo cluster
# detectamos si un nombre de cluster es positivo, en tal caso, extraer los nombres de los spins de ese cluster
hay_cluster_anterior <- which(te > 0)
if(!is.integer0(hay_cluster_anterior) ) {
  el_cluster <- te[hay_cluster_anterior]
  # rescatar los nombres del cluster en cuestion
  te2 <- c(te, mn[[el_cluster]])
  te <- te2[-hay_cluster_anterior]
}
mn <- c(mn, list( as.character(te )) ) # y ahora en mn[[2]] cluster 2, se compone de -5 y -6


# AQUI EN DONDE TENGO QUE HACER LAS MODIFICACIONES, function find_min_dist
# la funcion find_min_dist(D, Dold, nombres) es para re-actualizar la matriz D, utilizando 
# single linkage method entre el cluster recien funsionado y los demas clusters.
#temp <- find_min_dist(D, Dold, nombres)   #$dit 1.696616 1.273414 1.753150 1.741162   $nombres_ [1] "-1" "-3" "-5" "-6"
#dit <- temp$dit # dit son las minimias distancias de cluster fusionado en la iteracion anterior a todos los demas clusters que quedan
#nombres_ = temp$nombres_ # nombres_ son los nombres de los nodos que estan a minima distancia.
# abriendo la funcion find_min_dist
dit <- vector()
nombres_ <- colnames(D) # 
n <- length(nombres_) 
for (i in 1:n) {
  nodo <- nombres_[i]
  v <- min(Dold[nombres, nodo])  # ojo que en "nombres" estan los nombres de los clusters que se fusionaron
  #dit <- c(dit,min(Dold[, nodo]))
  # vector con nombres de los spins involucrados
  y <- c(nombres, nodo)  # los spins o nodos del cluster recien formado con los nodos de los otros clusters
  # como "y" puede contener cluster de antes, hay que detectar los nombres de los spins de ese cluster:
  hay_cluster_anterior <- which(as.numeric(y) >0)
  if (!is.integer0(hay_cluster_anterior)  ) { # en caso que hayan clusters anteriores
    el_cluster <- as.numeric(y[hay_cluster_anterior])
    # rescatar los nombres del cluster en cuestion
    y2 <- c(y, unlist(mn[el_cluster]) ) 
    y <- y2[-hay_cluster_anterior]
  }
  # calculo de la distancia MST
  d_mst <- find_mst_barrier(mst_g, y)
  # calculo de las distancias de acople
  d_cpl <- acople_distance_sum(J, y)
  # calculo de gamma
  gamma <-  d_cpl/d_mst
  #
  vprima <- v*gamma
  dit <- c(dit,vprima )
}
# el resultado es:                            2.794855 4.220567
# el resultado original con find_min_dist es: 1.681610 1.876042

#ahora tengo que poner los valores en fila.columna de D.
D <- put_dis(D, dit, nombres_, N, it)
print(it)


# # # # # # # # # # # # # # # # # # # #  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>iteracion 4
it = 4
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

# # # # # # # # # poniendo los nombres reales de los spins que componen el cluster recien formado.
#ii = which(names_spins$dend %in% nombres) # indices de los nombres de lso spins del nuevo cluster
te <- c(merge[it,1], merge[it,2]  ) # nombres de los spins del nuevo cluster
# detectamos si un nombre de cluster es positivo, en tal caso, extraer los nombres de los spins de ese cluster
hay_cluster_anterior <- which(te > 0)
if(!is.integer0(hay_cluster_anterior) ) {
  el_cluster <- te[hay_cluster_anterior]
  # rescatar los nombres del cluster en cuestion
  te2 <- c(te, mn[[el_cluster]])
  te <- te2[-hay_cluster_anterior]
}
mn <- c(mn, list( as.character(te )) ) # y ahora en mn[[2]] cluster 2, se compone de -5 y -6

if (length(Dold) <= 9) { # cuando la matriz es de 3X3, no se puede ejecutar directamente la funcion find_min_dist
  colno <- colnames(Dold) # nombres de columnas o fila de Dold
  #https://www.youtube.com/watch?v=8hSYEXIoFO8&t=78s
  sel <- colno[!(colno %in% nombres)] # nodo identificado
  col <- which(colno == sel) # identificamos la columna donde esta el nodo sel.
  v <- min(Dold[nombres, nodo]) 
  y <- mn[[it]]
  # calculo de la distancia MST
  d_mst <- find_mst_barrier(mst_g, y)
  # calculo de las distancias de acople
  d_cpl <- acople_distance_sum(J, y)
  # calculo de gamma
  gamma <-  d_cpl/d_mst
  #
  vprima <- v*gamma
  dit <- vprima 
  nombres_ <- sel
} else {
  temp <- find_min_dist_con_acoples(D, Dold, nombres) #esta es la que hay que crear.
  dit <- temp$dit # dit son las minimias distancias de cluster fusionado en la iteracion anterior a todos los demas clusters que quedan
  nombres_ = temp$nombres_ # nombres_ son los nombres de los nodos que estan a minima distancia.
}

#ahora tengo que poner los valores en fila.columna de D
D <- put_dis(D, dit, nombres_, N, it)


# # # # # # # # # # # # # # # # # # # #  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>iteracion 5
# Last iteration: when it = N - 1
it <- it + 1
rl <- cluster_find_name(D)
id <- rl$id
nombres <- rl$nombres
m <- rl$m
nombres_num <- as.numeric(nombres) # nuevo
merge[it,] <- as.numeric(c(nombres_num, it, m)) 


# # # # # # # # # # # # # # # # # # # #  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>PLOTEO
Dinicial <- sqrt(-J + 3) 
hc <- to_dendo(D = Dinicial, merge = merge, enames = colnames(dis) )
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

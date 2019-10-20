# testing the algorithm for hierarchical agglomeration clustering using simple linkage.

# En esta version 5, probamos las funciones de functions_hclust.R para hacer
# HAC sobre la matriz de distancias del MST de los acoples para el paper y verificamos
# que el resultado del dendograma sea equivalente al del MST. Es decir, deseamos
# verificar lo indicado en pag. 237 (ver tambien pag. 241 de los apuntes.)

# El esquema de operacion es similar a lo hecho en testing_hierarchical_clusteringV4.R.


# actual name: testing_hierarchical_clusteringV5.R
# 28.ago.19


# # # # # # # # # # # # # # # # # # # # DATA LOADING # # # # # # # # # # # # # # # # # # # # #

# Carga de datos, 
# igual que en clusteragrr_v1.R
rm(list = ls())
load('new_inferring_parameters_environment270219.RData') # load H and J parameters ising inferred
source("functions_hclust.R")
source("find_mst_barrier_function.R")
source("acople_distance_sum_function.R")
source("is_integer0_function.R")
library(igraph)
# # # # # # # # # # # # # # # # # # # # DATA LOADING # # # # # # # # # # # # # # # # # # # # #

# # # # # # # # # # # # # # # # # # # # DATA LOADING # # # # # # # # # # # # # # # # # # # # #
# de los 25 nodos de J, seleccionemos los primeros 6
#set.seed(123)
# seleccionemos 5 nodos cualquiera
nd <- c(4,8,12,17,19,22) # 01oct19
nd <- c(2,4,6,8,12,15,17,19, 22) # 09oct19
nd <- c(1,4,11,13,20,25) # 10oct19
nd <- c( 1, 3, 7, 11, 13, 14, 15, 17, 25) # 11oct19a
nd <- c(1, 2, 6, 7, 8, 10, 11, 13, 14, 19, 20, 21) # 11oct19b
nd <- 1:25 # 11oct19c
J <- J[nd, nd]
# # # # # # # # # # # # # # # # # # # # DATA LOADING # # # # # # # # # # # # # # # # # # # # #


# # # # # # # # # # # # # # # # # # # # # FIND MST # # # # # # # # # # # # # # # # # # # # # #
# Conforming the COUPLING Network
#colnames(J) <- colnames(wb)
#rownames(J) <- colnames(wb)
# Vamos a transformar los nombres de los nodos originales de una vez en nombres para el dendograma
colnames(J) <- rownames(J) <- -c(1:ncol(J)) # ojo -c(1:Nn) tiene el mismo orden que colnames(dis)

#rownames(h) <- c(1:dim(J)[1])
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



D <- dis
Nn <- ncol(D)
diag(D) <- rep(Inf, Nn)
colnames(D) <- rownames(D) <- -c(1:Nn) # ojo -c(1:Nn) tiene el mismo orden que colnames(dis)


# the magic
merge <- hierarchical_clustering(D)
hc <- to_dendo(D, merge, enames = colnames(dis) )
plot(hc)
plot(merge$cluster, merge$dultr, type = "S", pch = 19, 
     col = "red", xlab = "Cluster Number", ylab = "Ultrametric Distance")





# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 22 sept:
# a merge, necesito agregar dos columnas: la energia Emst asociada a cluster, y 
# la energia de acoples Ec asociada a cada cluster.
# the magic  --- viene de la L52
merge <- hierarchical_clustering_v2(D)
hc <- to_dendo(D, merge[,c(1:4)], enames = colnames(dis) )
plot(hc)
plot(merge$cluster, merge$dultr, type = "S", pch = 19, 
     col = "red", xlab = "Cluster Number", ylab = "Ultrametric Distance")
plot(merge$cluster, merge$dmst, type = "S", pch = 19, 
     col = "red", xlab = "Cluster Number", ylab = "dmst")
plot(merge$cluster, merge$dc, type = "S", pch = 19, 
     col = "red", xlab = "Cluster Number", ylab = "dc")
plot(merge$dmst, merge$dc, pch=19, xlab="dmst", ylab="dc")
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 



















##### SOLUCIONADO!!!! 
# vemos que tenemos un problema con hc$order.........

# vienen de testing_hierarchical_clusteringV3.R
# # # # # # # # # # # # # # 
# para genera rl hc$order, necesitamos identificar de merge[,c(1,2)] todos los
# nodos negativos y ordenarlos en un vector en su orden de aparicion.
temp <- (merge[,c(1,2)])
temp <- as.numeric(rbind(temp$node1, temp$node2))
temp <- temp[temp < 0]
temp <- -1*temp
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

dd <- as.dendrogram(hc)
temp2 <- order.dendrogram(dd)


a <- list()  # initialize empty object
# define merging pattern: 
#    negative numbers are leaves, 
#    positive are merged clusters (defined by row number in $merge)
a$merge <- as.matrix(merge[,c(1,2)])
a$height <- as.vector(merge$dultr)    # define merge heights
a$order <- temp2   
#hc$labels <- LETTERS[1:7]    # labels of leaves -----#poner nombres originales
a$labels <- c(1:ncol(D))
class(a) <- "hclust"        # make it an hclust object
plot(a)


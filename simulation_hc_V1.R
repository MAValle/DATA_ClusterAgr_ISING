# En este script vamos a utilizar la funcion hierarchical_clustering_v2
# y funcion hierarchical_clustering_v4 que estan en functions_hclust.R
# para siguiente simulacion:

# Generar una matriz J con J distribuido N(0,1) truncado en +-3 con dimension NXN
# calcular las distancias siguiendo la regla de sqrt(3-J)
# calcular el MST.
# calcular single linkage normal V2, y single linkage modificado V4
# Comparar la distribucion de las dmst y dc para cada caso
# comparar los energy path de cada caso.


# actual name: simulation_hc_V1.R
# 16.oct.19


# Notes:
# 19-oct-19: creation
#23-oct-19: creamos funcion create_coupling y create_mst 
# 27-oct-19: este script ya esta obsoleto porque la utilizacion de la funcion hierarchical_clustering_v4
#             ya esta obsoleta.




# # # # # # # # # # # # # # # # # # # # # LOADING DATA # # # # # # # # # # # # # # # # # # # # # # # # #
# Carga de datos, 
# igual que en clusteragrr_v1.R
rm(list = ls())
#load('new_inferring_parameters_environment270219.RData') # load H and J parameters ising inferred
source("functions_hclust.R")
#source("network_functions_v1.R") # aqui nos interesa la funcion find_mst_barrier
source("find_mst_barrier_function.R")
source("acople_distance_sum_function.R")
source("is_integer0_function.R")
source("entropies_functions.R") # nos interesa aqui get_energy function
library(igraph)
library(Matrix)
# # # # # # # # # # # # # # # # # # # # # LOADING DATA # # # # # # # # # # # # # # # # # # # # # # # # #



# # # # # # # # # # # # # # # # # # # # CREATING A J MATRIX # # # # # # # # # # # # # # # # # # # # # # #
# Funcion que crea matriz de acople y de distancia con media mu, desv estandar de acople sj y con Nn nodos
create_copupling <- function(Nn, media, sj) {
  library(Matrix)
  x <- Matrix(rnorm(Nn*Nn, mean=media, sd=sj),Nn)  # los J vienen con media 0 y desv estandar 1
  J <- forceSymmetric(x) # couplings
  D <- sqrt(3 - as.matrix(J)) # distance
  colnames(D) <- rownames(D) <- -c(1:Nn) # enumeramos todos las hojas o spines originales con valones negativos.
  diag(D) <- rep(Inf, Nn)
  colnames(J) <- rownames(J) <- -c(1:ncol(J)) 
  J <- as.matrix(J)
  return(list(J=J, D=D))
}
# ejemplo
# ot <- create_copupling(Nn=6, media=0, sj=1)
# J <- ot$J
# D <- ot$D

# creating artifitial matrix of couplings
Nn = 6 # number of spins
#library(Matrix)
set.seed(134)
ot <- create_copupling(Nn=6, media=0, sj=1)
J <- ot$J
D <- ot$D

# x <- Matrix(rnorm(Nn*Nn, mean=0, sd=1),Nn)  # los J vienen con media 0 y desv estandar 1
# J <- forceSymmetric(x) # couplings
# D <- sqrt(3 - as.matrix(J)) # distance
# colnames(D) <- rownames(D) <- -c(1:Nn) # enumeramos todos las hojas o spines originales con valones negativos.
# diag(D) <- rep(Inf, Nn)
# colnames(J) <- rownames(J) <- -c(1:ncol(J)) 
# # # # # # # # # # # # # # # # # # # # CREATING A J MATRIX # # # # # # # # # # # # # # # # # # # # #



# # # # # # # # # # # # # # # # # # # # # FIND MST # # # # # # # # # # # # # # # # # # # # # #
# Funcion que crea como salida la red de acoples y el mst de la red de acoples.
# imput: matriz de distancia D y de acoples J
create_mst <- function(J, D) {
  library(igraph)
  net <- graph_from_adjacency_matrix(as.matrix(J), mode="upper", weighted=TRUE, diag=FALSE)
  dis <- D
  # Convertir la matriz de distancia en objeto igraph
  g <- graph_from_adjacency_matrix(dis, mode = "undirected", weighted = TRUE, diag = FALSE, add.colnames = NULL, add.rownames = NA)
  E(g)$coupling <- E(net)$weight # asignamos las energias de acoples de net a E(g)$coupling
  mst_g <- minimum.spanning.tree(g, algorithm="prim")
  return(mst_g)
}
# ejemplo:
mst_g <- create_mst(J=J, D=D)

# Conforming the COUPLING Network
# 12-sep-19
# Vamos a transformar los nombres de los nodos originales de una vez en nombres para el dendograma
#colnames(J) <- rownames(J) <- -c(1:ncol(J)) # ojo -c(1:Nn) tiene el mismo orden que colnames(dis)

# #http://www.shizukalab.com/toolkits/sna/plotting-networks-pt-2
# net <- graph_from_adjacency_matrix(as.matrix(J), mode="upper", weighted=TRUE, diag=FALSE)
# #E(net)$coupling <- E(net)$weight
# # Conforming the MST Network
# dis <- D
# # Convertir la matriz de distancia en objeto igraph
# g <- graph_from_adjacency_matrix(dis, mode = "undirected", weighted = TRUE, diag = FALSE, add.colnames = NULL, add.rownames = NA)
# E(g)$coupling <- E(net)$weight # asignamos las energias de acoples de net a E(g)$coupling
# mst_g <- minimum.spanning.tree(g, algorithm="prim")
# edg <- as_edgelist(mst_g, names = TRUE)
# # # # # # # # # # # # # # # # # # # # # FIND MST # # # # # # # # # # # # # # # # # # # # # #

# # # # # # # # # # # # # # # # # # # # # PLOT MST # # # # # # # # # # # # # # # # # # # # # #
# https://mariliagaiarsa.weebly.com/uploads/3/8/6/2/38628397/igraphtutorialeng.html
# https://kateto.net/netscix2016.html
par(mfrow=c(2,1), oma=c(0,0,0,0), mar=c(1,1,1,1))#
l <- layout_in_circle(g)
plot(g, layout = l, main="Coupling net",
     vertex.color="yellow",
     vertex.color="lightgrey",
     #vertex.label=NA,
     vertex.size=15,
     edge.curved=0.3,
     edge.width=1/log(E(g)$weight))

plot(mst_g, layout = l, main="Coupling MSTnet",
     vertex.color="yellow",
     vertex.color="lightgrey",
     #vertex.label=NA,
     vertex.size=15,
     edge.curved=0.3,
     edge.width=1/log(E(mst_g)$weight))
# # # # # # # # # # # # # # # # # # # # # PLOT MST # # # # # # # # # # # # # # # # # # # # # #




# # # # # # # # # # # # # # # # # # #  PREPARATION # # # # # # # # # # # # # # # # # # # # # #
#D <- dis
Nn <- ncol(D)
#diag(D) <- rep(Inf, Nn)
#colnames(D) <- rownames(D) <- -c(1:Nn) # ojo -c(1:Nn) tiene el mismo orden que colnames(dis)
# equivalencia entre nombre de nodos reales y nombres en el dendograma:
#names_spins <- data.frame(real = colnames(J), dend = colnames(D))
# por ejemplo: encontrar las filas en names_spins que calzan con los valores en nombres.
# ii = which(names_spins$dend %in% nombres)
# names_spins[ii, "real"] # aqui obtenemos los nombres reales de los spins.
# # # # # # # # # # # # # # # # # # #  PREPARATION # # # # # # # # # # # # # # # # # # # # # #


# # # # # # # # # # # # # # # # # # # # # # # HC # # # # # # # # # # # # # # # # # # # # # # #
merge2 <- hierarchical_clustering_v2(D, J, mst_g)
hc2 <- to_dendo(D, merge2[,c(1:4)], enames=c(1:ncol(D)) )
plot(hc2)
library("ape")
plot(as.phylo(hc2), type = "unrooted", cex = 0.6, no.margin = TRUE)
# ultrametric distancias: http://www.bioinfo.org.cn/lectures/index-68.html
# https://stats.stackexchange.com/questions/103767/extract-ultrametric-distances-from-hclust-or-dendrogram
# https://en.wikipedia.org/wiki/Cophenetic_correlation
cophon2 <- cophenetic(hc2)

merge5 <- hierarchical_clustering_v5(D, J, mst_g)
hc5 <- to_dendo(D, merge5[,c(1:4)], enames=c(1:ncol(D)) )
plot(hc5)
plot(as.phylo(hc5), type = "unrooted", cex = 0.6, no.margin = TRUE)
cophon5 <- cophenetic(hc5)
# # # # # # # # # # # # # # # # # # # # # # # HC # # # # # # # # # # # # # # # # # # # # # # #

# # # # # # # # # # # # # # # # # # # # COPHONETICS DISIIMILARITIES # # # # #  # # # # # # # # # # # # # # #
# Calculamos la correlacion entre disimilirades cofoneticas (distancias ultrametricas entre nodos)
# las distancias MST originales.
# para single linkage normal
# https://stats.stackexchange.com/questions/33066/on-cophenetic-correlation-for-dendrogram-clustering?rq=1
temp <- (as.matrix(cophon2))
plot( D[lower.tri(D)] , temp[lower.tri(temp)])
cor(D[lower.tri(D)] , temp[lower.tri(temp)]) #0.702476

# para single linkage modificado cophenetics correlation
temp <- (as.matrix(cophon4))
plot( D[lower.tri(D)] , temp[lower.tri(temp)])
cor(D[lower.tri(D)] , temp[lower.tri(temp)]) # 0.6179028

# matrix plots
#https://stackoverflow.com/questions/3081066/what-techniques-exists-in-r-to-visualize-a-distance-matrix
#http://ichthyology.usm.edu/courses/multivariate/coldiss.R
library(lattice)
trellis.par.set(standard.theme(color = FALSE))
levelplot(D, xlab = "Node", ylab = "Node")
temp <- as.matrix(cophon2) 
diag(temp) <- Inf
levelplot(temp, xlab = "-Node", ylab = "-Node")
# otra alternativa:
require(gclus)
spe.color = dmat.color(1-D, byrank=FALSE, cm.colors(4))
spe.o = order.single(1-D)
speo.color = spe.color[spe.o,spe.o]
plotcolors(spe.color, rlabels=attributes(D)$dimnames[[1]], 
           main="Dissimilarity Matrix")
plotcolors(speo.color, rlabels=attributes(D)$dimnames[[1]], 
           main="Ordered Dissimilarity Matrix")
# 
# https://towardsdatascience.com/hierarchical-clustering-on-categorical-data-in-r-a27e578f2995
# # # # # # # # # # # # # # # # # # # # COPHONETICS DISIIMILARITIES # # # # #  # # # # # # # # # # # # # # #





# # # # # # # # # # # # # # # # # # # # SIMPLE GRAPHS # # # # #  # # # # # # # # # # # # # # #
# graficas simples recuros:
# https://www.datamentor.io/r-programming/subplot/
# https://stackoverflow.com/questions/2564258/plot-two-graphs-in-same-plot-in-r
# https://www.statmethods.net/advgraphs/parameters.html
# http://r-statistics.co/Top50-Ggplot2-Visualizations-MasterList-R-Code.html

# comparacion de la distancia ultrametrica
plot(merge2$cluster, merge2$dultr, type="l", col="black", pch=21, lwd=2)
lines(merge5$cluster, merge5$dultr , col="red", pch=22, lwd=2)
title(main = "d_ultr")
# comparacion de la distancia MST
plot(merge2$cluster, merge2$dmst, type="l", col="black", pch=21, lwd=2)
lines(merge4$cluster, merge4$dmst , col="red", pch=22, lwd=2)
title(main = "d_mst")
# comparacion de la distancia de acople
plot(merge2$cluster, merge2$dc, type="l", col="black", pch=21, lwd=2)
lines(merge4$cluster, merge4$dc , col="red", pch=22, lwd=2)
title(main = "d_c")
# https://stats.stackexchange.com/questions/103767/extract-ultrametric-distances-from-hclust-or-dendrogram
# # # # # # # # # # # # # # # # # # # # SIMPLE GRAPHS # # # # #  # # # # # # # # # # # # # # #





# # # # # # # # # # # # # # # # # # # # # MODULARITY # # # # #  # # # # # # # # # # # # # # # #


# modularity.
# https://stackoverflow.com/questions/25287032/modularity-calculation-for-weighted-graphs-in-igraph
# https://en.wikipedia.org/wiki/Louvain_modularity
# https://igraph.org/r/doc/modularity.igraph.html
# https://igraph.org/r/doc/cluster_fast_greedy.html
# https://www.geeksforgeeks.org/introduction-hill-climbing-artificial-intelligence/
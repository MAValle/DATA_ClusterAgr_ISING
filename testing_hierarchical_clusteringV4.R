# testing the algortith for hierarchical agglomeration clustering using simple linkage.

# En esta version 4, hacemos un test utilizando una matriz de distancia mas grande
# y con las funciones generales en functions_hclust.R desarrolladas en el 27.ago.19
# basados en testing_hierarchical_clusteringV3.R

# actual name: testing_hierarchical_clusteringV4.R
# 27.ago.19

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

rm(list = ls())
source("functions_hclust.R")

# creating artifitial matrix of couplings
Nn = 25 # number of spins
library(Matrix)
set.seed(134)
x <- Matrix(rnorm(Nn*Nn),Nn)  # los J vienen con media 0 y desv estandar 1
J <- forceSymmetric(x) # couplings
D <- sqrt(3 - as.matrix(J)) # distance
colnames(D) <- rownames(D) <- -c(1:Nn) # enumeramos todos las hojas o spines originales con valones negativos.
diag(D) <- rep(Inf, Nn)


# the magic
merge <- hierarchical_clustering(D)
hc <- to_dendo(D, merge, enames=c(1:ncol(D)) )
plot(hc)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 





















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


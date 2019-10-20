# testing the algortith for hierarchical agglomeration clustering using simple linkage.

# En esta version 1, hacemos un test utilizando solo N=4 nodos.


# actual name: testing_hierarchical_clustering.R
# 22.ago.19



rm(list = ls())
source("functions_hclust.R")

N = 4 # number of spins
merge <- matrix(NA, ncol=4, nrow=N-1) # col1 y col2 spins merged, col3 numero de cluster, col4 = merge distance
merge <- as.data.frame(merge)
colnames(merge) <- c("node1", "node2", "cluster", "dultr")

library(Matrix)
set.seed(13)
x <- Matrix(rnorm(N*N),N)  # los J vienen con media 0 y desv estandar 1
J <- forceSymmetric(x) # couplings
D <- sqrt(3-J) # distance
diag(D) <- rep(0, N)
colnames(D) <- rownames(D) <- c(1:N)
D[lower.tri(D, diag = TRUE)] <- Inf # nos quedamos solo con la parte superior de la matgriz

# # # primera iteracion
it <- 1
rl <- cluster_find_name(D)
id <- rl$id
nombres <- rl$nombres
m <- rl$m
merge[it,] <- as.numeric(c(nombres, N+it, m))
it <- it + 1




# ahora tengo fusionar los spins indicados en id en la matriz D
# nombre del nuevo cluster = N+1

# fila id[1] y columna id[2] se deben borrar
Dold <- D
D <- D[-id, -id]

# ahora hay que agregar fila y columna 
# pero antes hay que buscar las distancias minimas entre los spins restantes y el cluster
if (length(Dold) <= 9) { # cuando la matriz es de 3X3, no se puede ejecutar directamente la funcion find_min_dist
  colno <- colnames(Dold) # nombres de columnas o fila de Dold
  #https://www.youtube.com/watch?v=8hSYEXIoFO8&t=78s
  sel <- colno[!(colno %in% nombres)] # nodo identificado
  col <- which(colno == sel) # identificamos la columna donde esta el nodo sel.
  dit <- vector()
  dit <- c(dit,min(Dold[, col])) # la distancia minima
} else {
  temp <- find_min_dist(D, Dold)
  dit <- temp$dit # dit son las minimias distancias de cluster fusionado en la iteracion anterior a todos los demas clusters que quedan
  nombres_ = temp$nombres_ # nombres_ son los nombres de los nodos que estan a minima distancia.
}


#aqui
#ahora tengo que poner los valores en fila.columna de D
D <- put_dis(D, dit, nombres_, N, it)


# # # segunda iteracion
rl <- cluster_find_name(D)
id <- rl$id
m <- rl$m
nombres <- rl$nombres
merge[it,] <- as.numeric(c(nombres, N+it, m))
it <- it + 1


Dold <- D
D <- D[-id, -id]  # como es el ultimo nodo, aqui la matriz desaparece 

# como lo hago aqui para identificar el nodo 4?
# solucion:
if (length(Dold) <= 9) { # cuando la matriz es de 3X3, no se puede ejecutar directamente la funcion find_min_dist
  colno <- colnames(Dold) # nombres de columnas o fila de Dold
  #https://www.youtube.com/watch?v=8hSYEXIoFO8&t=78s
  nombres_ <- colno[!(colno %in% nombres)] # nodo identificado
  col <- which(colno == nombres_) # identificamos la columna donde esta el nodo sel.
  dit <- vector()
  dit <- c(dit,min(Dold[, col])) # la distancia minima
} else {
  temp <- find_min_dist(D, Dold)
  dit <- temp$dit # dit son las minimias distancias de cluster fusionado en la iteracion anterior a todos los demas clusters que quedan
  nombres_ = temp$nombres_ # nombres_ son los nombres de los nodos que estan a minima distancia.
}


#aqui
D <- put_dis(D, dit, nombres_, N, it)





# # # tercera y ultima iteracion (son 4 nodos)
rl <- cluster_find_name(D)
id <- rl$id
nombres <- rl$nombres
m <- rl$m

# sp <- (c( as.numeric(rownames(D)[id[1]]), as.numeric(colnames(D)[id[2]])))
# merge[3,] <- as.numeric(c(sp, N+it, m))
merge[it,] <- as.numeric(c(nombres, N+it, m))
it <- it + 1

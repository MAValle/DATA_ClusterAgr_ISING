# este es un test intentando generar el algoritmo de aglomeracion con single linkage 
# con una prueba simple.

# Lo que hacemos aqui es  utilizar el concepto de clustering aggregation basado en aglomeración. 
# El procedimiento es así:
#   
#   \1. Generamos matriz E = \alpha*E_{c} + (1-\alpha)*E_{mst} para cada par de nodos.
#   
#   \2. Search for a pair of nodes i,j que produzca la menor E y lo agrupamos.
#   
#   \3. Se repite el proceso de agregacion aglomerativo hasta que se terminen todos los nodos.
#   
#   Si se hace el merge o aglomeracion con sigle linkage, obtendria el MST tree (con alpha=0).
#   
# 

# Recursos de agglomerative clustering:
# Paper de ARISTIDES GIONIS "Clustering Aggregation"
# https://en.wikipedia.org/wiki/Particle_aggregation   

# Recursos de -clusterign aggregation-
# https://www.youtube.com/watch?v=XJ3194AmH40

# recursos de plot hclust
# https://atrebas.github.io/post/2019-06-08-lightweight-dendrograms/
# https://www.gastonsanchez.com/visually-enforced/how-to/2012/10/03/Dendrograms/

# Nombre del acutual archivo: clusteragrr_v2.R
# ubicacion: dropbox->Research->Project_cluster_aggregation->DATA_ClusterAgr_ISING
# creation date: 16-aug-19
# companion file: 

#Notas:
# 16-aug-19: creation
# 21.ago-19: he terminado el ejemplo con N=4 nodos.

rm(list = ls())

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
m <- min(D[row(D)!=col(D)]) # nos da el valor minimo 
id <- which(D == min(D), arr = TRUE) # nos da la locacion del valor minimo fila, columna
id <- as.vector(id)
nombres <- rownames(D)[id] # nos da los nombres del nodo fila y nodo columna de D que se fusionan
#sp <- (c( as.numeric(rownames(D)[id[1]]), as.numeric(colnames(D)[id[2]])))
#merge[1,] <- as.numeric(c(id, N+it, m))
merge[1,] <- as.numeric(c(nombres, N+it, m))
it <- it + 1




# ahora tengo fusionar los spins indicados en id en la matriz D
# nombre del nuevo cluster = N+1

# fila id[1] y columna id[2] se deben borrar
Dold <- D
D <- D[-id, -id]

# ahora hay que agregar fila y columna 
# pero antes hay que buscar las distancias minimas entre los spins restantes y el cluster
dit <- vector() # dit son las minimias distancias de cluster fusionado en la iteracion anterior a todos los demas clusters que quedan
nombres_ <- colnames(D) # vemos los nombres de las filas o columnas que quedan en D
n <- length(nombres_)
for (i in 1:n) {
  nodo <- nombres_[i]
  dit <- c(dit,min(Dold[, nodo]))
}



#ahora tengo que poner los valores en fila.columna de D
D <- cbind(D, dit)
D <- rbind(D, rep(Inf, it+1))
colnames(D) <- rownames(D) <- c(nombres_, N+it-1)


# # # segunda iteracion
m <- min(D[row(D)!=col(D)]) # nos da el valor minimo 
id <- which(D == min(D), arr = TRUE) # nos da la locacion del valor minimo
id <- as.vector(id)
nombres <- rownames(D)[id] # nos da los nombres del nodo fila y nodo columna de D que se fusionan
# sp <- (c( as.numeric(rownames(D)[id[1]]), as.numeric(colnames(D)[id[2]])))
# merge[2,] <- as.numeric(c(sp, N+it, m))
merge[2,] <- as.numeric(c(nombres, N+it, m))
it <- it + 1


Dold <- D
D <- D[-id, -id]  # como es el ultimo nodo, aqui la matriz desaparece 

# como lo hago aqui para identificar el nodo 4?
# solucion:
if (length(Dold) <= 9) {
  colno <- colnames(Dold) # nombres de columnas o fila de Dold
  #https://www.youtube.com/watch?v=8hSYEXIoFO8&t=78s
  sel <- colno[!(colno %in% nombres)] # nodo identificado
  col <- which(colno == sel) # identificamos la columna donde esta el nodo sel.
}
dit <- vector()
# nodo <- 4
# col = 2
dit <- c(dit,min(Dold[, col])) # la distancia minima


D <- cbind(D, dit)
D <- rbind(D, rep(Inf, ncol(D) ) )
colnames(D) <- rownames(D) <- c(sel,  N+it-1)



# # # tercera y ultima iteracion (son 4 nodos)
m <- min(D[row(D)!=col(D)]) # nos da el valor minimo 
id <- which(D == min(D), arr = TRUE) # nos da la locacion del valor minimo
id <- as.vector(id)
nombres <- rownames(D)[id] # nos da los nombres del nodo fila y nodo columna de D que se fusionan

# sp <- (c( as.numeric(rownames(D)[id[1]]), as.numeric(colnames(D)[id[2]])))
# merge[3,] <- as.numeric(c(sp, N+it, m))
merge[3,] <- as.numeric(c(nombres, N+it, m))
it <- it + 1


# testing the algortith for hierarchical agglomeration clustering using simple linkage.

# En esta version 3, hacemos un test utilizando un ejemplo de matriz de distancias
# que se encuentra en:
# http://84.89.132.1/~michael/stanford/maeb7.pdf


# actual name: testing_hierarchical_clusteringV3.R
# 23.ago.19



rm(list = ls())
source("functions_hclust.R")
N = 7 # number of spins
merge <- matrix(NA, ncol=4, nrow=N-1) # col1 y col2 spins merged, col3 numero de cluster, col4 = merge distance
merge <- as.data.frame(merge)
colnames(merge) <- c("node1", "node2", "cluster", "dultr")
D = matrix(c(0, 0.5, 0.4286, 1, 0.25, 0.6250, 0.3750,
             0.5000, 0, 0.7143, 0.8333, 0.6667, 0.2000, 0.7778,
             0.4286, 0.7143, 0, 1.0000, 0.4286, 0.6667, 0.3333,
             1.0000, 0.8333, 1.0000, 0, 1.0000, 0.8000, 0.8571,
             0.2500, 0.6667, 0.4286, 1.0000, 0, 0.7778, 0.3750,
             0.6250, 0.2000, 0.6667, 0.8000, 0.7778, 0, 0.7500,
             0.3750, 0.7778, 0.3333, 0.8571, 0.3750, 0.7500, 0), nrow = N, ncol = N)
colnames(D) <- rownames(D) <- -c(1:N) # enumeramos todos las hojas o spines originales con valones negativos.
#D[lower.tri(D, diag = TRUE)] <- Inf # nos quedamos solo con la parte superior de la matgriz
diag(D) <- rep(Inf, N)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
## este es el codigo que hace toda la magia.
iteraciones <- N - 2
for (it in 1:iteraciones) {
  # # # # # # # # #  Find the minimum distance among clusters:
  rl <- cluster_find_name(D)
  id <- rl$id
  nombres <- rl$nombres
  m <- rl$m
  # leaves (original spins) will be negative, merged clusters will be positive starting from 1, 2, ...
  nombres_num <- as.numeric(nombres) # nuevo
  #r <- which(nombres_num <= N) # nuevo
  #nombres_num[r] <- -1*nombres_num[r] # nuevo
  # originalmerge[it,] <- as.numeric(c(nombres, N+it, m)) # se coloca nodo1 (fusion) / nodo2 (fusion) /nombre del nuevo cluster / distancia a la que se unen
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
    temp <- find_min_dist(D, Dold)
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
#merge[it,] <- as.numeric(c(nombres, N+it, m))
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# como creamos el dendograma a partir de merge
# https://stackoverflow.com/questions/2310913/how-do-i-manually-create-a-dendrogram-or-hclust-object-in-r
hc <- list()  # initialize empty object
# define merging pattern: 
#    negative numbers are leaves, 
#    positive are merged clusters (defined by row number in $merge)
hc$merge <- as.matrix(merge[,c(1,2)])
hc$height <- as.vector(merge$dultr)    # define merge heights
hc$order <- temp   
#hc$labels <- LETTERS[1:7]    # labels of leaves -----#poner nombres originales
hc$labels <- 1:7
class(hc) <- "hclust"        # make it an hclust object
plot(hc)


# # # # # # # # # # # # # # 
# para genera rl hc$order, necesitamos identificar de merge[,c(1,2)] todos los
# nodos negativos y ordenarlos en un vector en su orden de aparicion.
temp <- (merge[,c(1,2)])
temp <- as.numeric(rbind(temp$node1, temp$node2))
temp <- temp[temp < 0]
temp <- -1*temp
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 








# solo un ejemplo...................
a <- list()  # initialize empty object
# define merging pattern: 
#    negative numbers are leaves, 
#    positive are merged clusters (defined by row number in $merge)
a$merge <- matrix(c(-1, -2,
                    -3, -4,
                    1,  2), nc=2, byrow=TRUE ) 
a$height <- c(1, 1.5, 3)    # define merge heights
a$order <- c(4,3,2,1)             # order of leaves(trivial if hand-entered)
a$labels <- 1:4    # labels of leaves
class(a) <- "hclust"        # make it an hclust object
plot(a)
#convert to a dendrogram object if needed
ad <- as.dendrogram


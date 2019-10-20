# testing the algortith for hierarchical agglomeration clustering using simple linkage.

# En esta version 2, hacemos un test utilizando un ejemplo de matriz de distancias
# que se encuentra en:
# https://home.deib.polimi.it/matteucc/Clustering/tutorial_html/hierarchical.html


# actual name: testing_hierarchical_clusteringV2.R
# 22.ago.19



rm(list = ls())
source("functions_hclust.R")
N = 6 # number of spins
merge <- matrix(NA, ncol=4, nrow=N-1) # col1 y col2 spins merged, col3 numero de cluster, col4 = merge distance
merge <- as.data.frame(merge)
colnames(merge) <- c("node1", "node2", "cluster", "dultr")
D = matrix(c(0, 662, 877, 255, 412, 996,
             662, 0, 295, 468, 268, 400,
             877, 295, 0, 754, 564, 138, 
             255, 468, 754, 0, 219, 869,
             412, 268, 564, 219, 0, 669,
             996, 400, 138, 869, 669, 0), nrow = N, ncol = N)
colnames(D) <- rownames(D) <- c(1:N)
#D[lower.tri(D, diag = TRUE)] <- Inf # nos quedamos solo con la parte superior de la matgriz
diag(D) <- rep(Inf, N)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
iteraciones <- N - 2
for (it in 1:iteraciones) {
  # # # # # # # # #  Find the minimum distance among clusters:
  rl <- cluster_find_name(D)
  id <- rl$id
  nombres <- rl$nombres
  m <- rl$m
  merge[it,] <- as.numeric(c(nombres, N+it, m))
  #it <- it + 1
  
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
  D <- put_dis(D, dit, nombres_, N, it+1)
  print(it)
}
# Last iteration: when it = N - 1
it <- it + 1
rl <- cluster_find_name(D)
id <- rl$id
nombres <- rl$nombres
m <- rl$m
merge[it,] <- as.numeric(c(nombres, N+it, m))
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


# como creamos el dendograma a partir de merge
# https://stackoverflow.com/questions/2310913/how-do-i-manually-create-a-dendrogram-or-hclust-object-in-r
# ejemplo

a <- list()  # initialize empty object
# define merging pattern: 
#    negative numbers are leaves, 
#    positive are merged clusters (defined by row number in $merge)
a$merge <- matrix(c(-1, -2,
                    -3, -4,
                    1,  2), nc=2, byrow=TRUE ) 
a$height <- c(1, 1.5, 3)    # define merge heights
a$order <- 1:4              # order of leaves(trivial if hand-entered)
a$labels <- LETTERS[1:4]    # labels of leaves
class(a) <- "hclust"        # make it an hclust object
plot(a)
#convert to a dendrogram object if needed
ad <- as.dendrogram(a)
# Lo que hacemos aqui es  utilizar el concepto de clustering aggregation basado en aglomeración. El procedimiento es así:
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

# Nombre del acutual archivo: clusteragrr_v1.R
# ubicacion: dropbox->Research->Project_cluster_aggregation->DATA_ClusterAgr_ISING
# creation date: 21-jun-19
# companion file: 

#Notas:
# 21-jun-19: creation
# 08-ago-19: para generar la energia de mst (distancia mst) Emst, estamos utilizando 
#     la funcion find_mst_barrier(mst_g, items) en L129, ya que esta produce 
#     coincidencia entre el E(mst_g)$weight y el Emst, que es lo que corresponde.




# # # # # # # # # # # # # # # # # # # # CARGA de datos # # # # # # # # # # # # # # # # # # # # 
# 040319
rm(list = ls())
library(dplyr)
library(ggplot2)
library(latex2exp)
library(gdata)
library(scales)  # makes pretty labels on the x-axis
library(colorspace)
library(purrr)
library(igraph)
library(svMisc) # progress bar
library(purrr)
source("entropies_functions.R")
source("network_functions_v1.R")
source("find_mst_barrier_function.R")
source("find_starting_vertex_function.R")

load('new_inferring_parameters_environment270219.RData') # load H and J parameters ising inferred
# # # # # # # # # # # # # # # # # # # # CARGA de datos # # # # # # # # # # # # # # # # # # # # 



# # # # # # # # # # # # # # # # # # # # FUNCIONES # #  # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # FUNCIONES # #  # # # # # # # # # # # # # # # # # # # # 




# # # # # # # # # # # # # # # # # # # # # FIND MST # # # # # # # # # # # # # # # # # # # # # #
# Conforming the COUPLING Network
colnames(J) <- colnames(wb)
rownames(J) <- colnames(wb)
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


# Starting -----------------------21, jun, 2019
# Create the matrix of energies E = \alpha*E_{c} + (1-\alpha)*E_{mst} para cada par de nodos
nodes <- colnames(wb)
matrix_energy <- matrix(NA, ncol=length(nodes), nrow=length(nodes))
colnames(matrix_energy) <-  rownames(matrix_energy) <- nodes
# creation of dataframe de trabajo, matriz a dataframe
temp <- combn(nodes,2)
temp <- t(temp)
df <- as.data.frame(temp)
df <- cbind(df, rep(NA, nrow(df)))
colnames(df) <- c("node1", "node2", "E")
# creation of dataframe de trabajo, matriz a dataframe
### dataframe de los indexs row, col de matrix_energy. como n=25 nodos, entonces la matriz superior tiene n(n+1)/2 
# elementos sin ingluir la diagonal.
indexes <- as.data.frame(which(lower.tri(matrix_energy, diag = FALSE), arr.ind=T))
# vamos a trabajar con la diagonal inferior.
# in indexes tenemos los ids de df
cluster = paste(df$node2, df$node1, sep=",")
#cluster <- paste0("{", cluster,"}")
df <- cbind(indexes, cluster, df) # df tiene N(N+1)/2 filas
N <- length(nodes)
rm(indexes)
# con df vamos a trabajar en vez de matrix_energy, pero en cualquier momento podemos convertir df a una matriz.

# ----------------------- inicio del llenado de la matrix df (proxy de distancia en caso de clustering aglomerativo)
alpha = 0.0 # minimizar E=alphaX*Ec + (1-alpha)*Emst
# resultado: llenado de df
for (i in 1:nrow(df) ) {
  node1 <- matrix(0, ncol=N, nrow=1)
  node2 <- matrix(0, ncol=N, nrow=1)
  colnames(node1) <- colnames(node2) <- nodes
  x <- which(colnames(node1) == df[i,"node2"])
  node1[x] <- 1
  x <- which(colnames(node2) == df[i,"node1"])
  node2[x] <- 1
  v <- node1 + node2
  # subrutina para calcular energia de acople
  Ec <- get_energy(v, mode="c") # funcion de entropies_functions.R / Matriz J de acoples debe estar cargada.
  # subrutina para calcular energia de mst
  t1 <- which(v==1)
  items <- nodes[t1]
  items <- as.factor(items)
  items <- as.numeric(as.character(items))
  # # usando solo d_mst_distances function
  #r <- d_mst_distances(mst_g, vec = items)  ###########OJO cambiarla por find_mst_barrier en find_mst_barrier_function.R
  #dmst <- r[[1]] #distancia
  #Emst <-  -1*r[[2]] #energia mst (multiplicamos por -1 porque si bien los acoples son positivos, la energia es negativa)
  # # usando solo d_mst_distances function
  Emst <- find_mst_barrier(mst_g, items)
  df[i, "E"] <- alpha*Ec + (1-alpha)*Emst
}  # OJO como hclust single-link maximiza similitud, tenemos que multiplicar la energia E por -1.
# ahora identificamos el minimo y extraemos informacion
id <- which(df$E == min(df$E)) #170, o sea 72 con 221 se unen, ahora este {72,221} es un cluster y hau quue calcular la enrgia
x <- df[id, "cluster"]
x <- as.numeric(strsplit(x,split=",",fixed=TRUE)[[1]]) # extraemos los nodos involucrados
# de df tengo que quitar row 170 que es el cluster {72,221} y todas aquellas filas en que sale el 72 y el 221


# de este cluster a todos los demas nodos
# ----------------------- Actualizacion de matrix_energy (no se si esto es necesario, ya que trabajaremos con df)
for (i in 1:nrow(df)) {
  matrix_energy[df[i, "row"], df[i, "col"]] <- 1*df[i, "E"]
  # 
}
matrix_energy[upper.tri(matrix_energy)] = t(matrix_energy)[upper.tri(matrix_energy)]
# ahora esta seria mi matriz de distancia.
var_HC <- hclust(as.dist(matrix_energy), method="single") 
var_HC_c <- hclust(as.dist(matrix_energy), method="complete") 
var_HC_av <- hclust(as.dist(matrix_energy), method="average") 
# en var_HC$merge (cuando alpha=0.5) vemos que en la primera iteracion se unio el vector nodes[9] con el node[15]
# que corresponde justamente row 170 de df es decir, cluster {72,221}. Luego en la segunda iteracion
# se une a este cluster el nodo nodes[7], es decir, {72,221,100}, etc.....
plot(var_HC)
coef(var_HC)
plot(var_HC_c)
plot(var_HC_av)
cutree(var_HC, h = 1.5)
library(dendextend)
library(ape)
plot(as.phylo(var_HC), cex = 0.9)
# En var_HC$height  estan las energias de cada cluster fusionado secuencialmente.
plot(1:length(var_HC$height), var_HC$height, xlab="merge steps", ylab="Energy Potential - U", pch=19)
# esto nos da un elemento distintivo unico del sistema, lo cual puedo comparar con la agregacion desde un punto de partida.

plot(1:length(var_HC$height), log(var_HC$height), xlab="merge steps", ylab="Energy Potential - U", pch=19)
hist(var_HC$height,15)
temp <- as.data.frame(var_HC$merge)
temp$height <- var_HC$height
temp # nodos que se van aglomerando en col1 y col2, y col3 es la distancia ultrametrica a la cual se unen.

# ----------------------- Recalculo de df y energy_matrix con el nuevo cluster fusionado
# https://www.youtube.com/watch?v=aXsaFNVzzfI <---- Lance_williams algorithm (LW)
# el LW es la base sobre la cual funciona el clusterign aglomerativo. En nuestro caso lo podemos
# utilizar de manera modificada puesto que el metodo de fusion no es ninguno de los conocidos (ward, single, average, etc)
# sino uno modificado para m¡nuestrois propositos porque utilizamos energia y no distancia.
# En nuestro caso, esto se reduce a utilizar el single link, que es el min{D_ki, D_kj}.
Si lo quisiera programa todo yo, sin hclust, para la segunda iteracion tendria que hacer:
En df tengo que quitar row 170 que es el cluster {72,221} y todas aquellas filas en que sale el 72 y el 221
creo otra lista o matriz de cluster Cl donde coloco cluster {72,221}, E=energia (Cl es la lista de cluster que se van generando)
df <- df + fila con cluster {72,221, X} donde X son todos los demas nodos.

# hierarchical Agglomerative clustering HAC
# https://nlp.stanford.edu/IR-book/html/htmledition/hierarchical-agglomerative-clustering-1.html
# 22-jun-19 
x <- "60,234"
as.numeric(strsplit(x,split=",",fixed=TRUE)[[1]])
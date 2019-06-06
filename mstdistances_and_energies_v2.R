# Hacemos lo mismo que esta en mstdistances_and_energies.R, utilizando sus mismos codigos
# pero con datos de acoples inferidos con 25 productos: (ver new_inferring_parameters.R)
# Las inferencias se encuentran en: new_inferring_parameters_environment270219.RData

# Dado que son 25 productos, el numero de estados es de app 33.5 millones, lo cual 
# hace que generar la matriz de active_vertex sea impractico. Hay que hacer otro enfoque:

# 1. se comienza por un vertice (producto) de posea alto strength en el MST
# (equivale a un producto de alta rotacion)
# 2. A este vertice le encontramos su energia de acople Ec y distancia MST 
# (creamos matriz tableau con columnas: vertex, Ec , Emst, Et=Ec+Emst)
# 3. inicio de loop: buscar otro vertex adyacente a v en el mst, que minimice DEc (ver paper congreso)
# 4. agregar ese vertex a tableau
# 5 volver a paso 3 hasta desired number of size of the cluster
# cluster seria el conjunto de todos los elementos de la columna vertgex en tableau.


# Nombre del acutual archivo: mstdistances_and_energies_v2.R
# ubicacion: dropbox->Research->PAPER MST-Fondecyt_2->data
# creation date: 05-mar-19
# companion file: 

#Notas:
# 05-mar-19: creation
# 26-may-19: agregamos al algoritmo de cluster aggregation, la posibilidad de ir agregando
#         no solo por alpha*Ec + (1-alpha)*Emst, sino tambien por solo minimizando Ec


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
source("entropies_functions.R")
source("network_functions_v1.R")

load('new_inferring_parameters_environment270219.RData')
# # # # # # # # # # # # # # # # # # # # CARGA de datos # # # # # # # # # # # # # # # # # # # # 


# # # # # # # # # # # # # # # # # # # # # FIND MST # # # # # # # # # # # # # # # # # # # # # #
# Conforming the COUPLING Network
colnames(J) <- colnames(wb)
rownames(J) <- colnames(wb)
#rownames(h) <- c(1:dim(J)[1])
#http://www.shizukalab.com/toolkits/sna/plotting-networks-pt-2
net <- graph_from_adjacency_matrix(as.matrix(J), mode="upper", weighted=TRUE, diag=FALSE)
#signos_de_J <- as.matrix(sign(upperTriangle(as.matrix(J), diag = FALSE, byrow=TRUE))) #obtenemos los signos de los acoples
#E(net)$signos <- signos_de_J
E(net)$coupling <- E(net)$weight
#E(net)$color[E(net)$signos > 0] <- 'orange' #acople positivo.
#E(net)$color[E(net)$signos < 0] <- 'blue' #acople negativo
#V(net)$magn <- abs(t(H))
#V(net)$signs <- sign(H)
#V(net)$color <- 'lightgrey'
#V(net)$color[V(net)$signs > 0] <- 'green'

# Conforming the MST Network
dis <- sqrt(-J + 3) 
# Convertir la matriz de distancia en objeto igraph
g <- graph_from_adjacency_matrix(dis, mode = "undirected", weighted = TRUE, diag = FALSE, add.colnames = NULL, add.rownames = NA)
#E(g)$signos <- signos_de_J
E(g)$coupling <- E(net)$weight
#E(g)$color[E(net)$signos > 0] <- 'orange' #acople positivo.
#E(g)$color[E(net)$signos < 0] <- 'blue' #acople negativo
#V(g)$magn <- abs(t(H))
#V(g)$signs <- sign(H)
#V(g)$color <- 'lightgrey'
#V(g)$color[V(g)$signs > 0] <- 'green'
mst_g <- minimum.spanning.tree(g, algorithm="prim")
edg <- as_edgelist(mst_g, names = TRUE)
# # # # # # # # # # # # # # # # # # # # # FIND MST # # # # # # # # # # # # # # # # # # # # # #

#convertir de distancia mst a energia J
dist_to_j <- function(vect) {
  value <- 3 - vect^2
  return(value)
}

# # FUNCTION
# # # # # # # # # # # #  # # # # # # FIND strongest vertex # # # # # # # # ## #  # # # # # # #
# Input: igraph object the mst
# output: the name of teh starting vertex in numeric format
find_starting_vertex <- function(mst_g) {
  # strengh de los nodos en el MST
  #st <- strength(mst_g, weights = E(mst_g)$coupling)
  st <- strength(mst_g, weights = E(mst_g)$weight)
  id <- which(st == max(st))
  v <- as.numeric(as_ids(V(mst_g)[id]))
  v <- c(v)
  return(v)
}
# example:
#v <- find_starting_vertex(mst_g)
# Result: node 80
# # # # # # # # # # # #  # # # # # # FIND strongest vertex # # # # # # # # ## #  # # # # # # #
# # FUNCTION
# # # # # # # # # # # # # # # # # # # # # FIND Emst # # # # # # # # # # #  # # # # # # # # # #
# Inspirado en la idea de la L164 en adelante en mstdistances_and_energies.R
# input: objeto mst igraph 
# input: vector de nombres de spins activos en el MST
# output: E_mst : energia o distancia de MST entre los spins activos en el MST
find_mst_barrier <- function(mst_g, v) {
  temp <- combn(v,2)
  calce <- vector(mode="numeric", length=0)
  for (j in 1:ncol(temp) ) {
    # Primero vemos si los dos nodos son adyacentes:
    node1 <- which(as.character(temp[1,j]) == V(mst_g)$name)
    node2 <- which(as.character(temp[2,j]) == V(mst_g)$name)
    
    #are_adjacent(mst_g, as.character(temp[1,j]), as.character(temp[2,j]) )
    #ed <- get.edge.ids(mst_g, vp=c(as.character(temp[1,j]), as.character(temp[2,j]) ) ) 
    
    if ( ( are_adjacent(mst_g, node1, node2 ) ) ) {
      ed <- get.edge.ids(mst_g, vp=c(node1, node2 ) )
      calce <- c(calce, ed)
    } else {
      # identifica todos los edges para ir de nodo temp[1,j] a nodo temp[2,j]
      #sp <- shortest_paths(mst_g, from=as.character(temp[1,j]), to=as.character(temp[2,j]), output="epath")
      sp <- shortest_paths(mst_g, from=node1, to=node2, output="epath")
      
      #ed <- sp$epath[[1]]
      ed <- unlist(sp$epath) #09-dic-18
      #ed <- all_simple_paths(mst_g, from = as.character(temp[1,j]), to = as.character(temp[2,j]) )
      #ed <- ed[[1]]
      
      # identifica los id de los edges en ed
      #calce <- c(calce, match(ed, E(mst_g), nomatch = NA_integer_, incomparables = NULL) )
      calce <- c(calce, ed ) # 09-dic-18
    }
    
    distance <- sum(E(mst_g)$weight[unique(calce)]) # esta es la suma de energias de acoples en el MST.
    #distance <- sum( dist_to_j(E(mst_g)$weight[unique(calce)]) ) # esta es la suma de energias de acoples en el MST.
    #E_mst <- sum(E(mst_g)$coupling[unique(calce)]) 
  }
  return(distance)
}
# Ejemplo
#E_mst <- find_mst_barrier(mst_g, v)
# # # # # # # # # # # # # # # # # # # # # FIND Emst # # # # # # # # # # #  # # # # # # # # # #





# Starting
v <- find_starting_vertex(mst_g) # basado en min alpha*Ec + (1-alpha)*Emst
number_of_vertexs <- 20
alpha <- 1 #nivel de importancia a la energia de acople Ec 
tableau <- matrix(NA, ncol=4, nrow=number_of_vertexs ) # vertex name, Ec, Emst, Etotal
# generar el vector de estado:
vertex_positions <- as.numeric(rownames(H))
state <-  rep(0, length(H))
id <- which(vertex_positions == v)
state[id] <- 1
#energia_total <- get_energy(as.matrix(state), mode="t")
#energia_field <- get_energy(as.matrix(state), mode="f")
energia_coupl <- get_energy(t(as.matrix(state)), mode="c")
Ec <- 0 # la distancia MST de un unico vertex es nula.
tableau[1, ] <- c(v[length(v)], energia_coupl, 0, alpha*energia_coupl + (1-alpha)*0)
# End Starting 

# # # # ALGORITMO de CLUSTER AGGREGATION
# Find the next vertexs
# buscar estado (80, u) todas las N-(L-1) combinaciones siendo L el numero de spins activos, 
# que (ya no: produzcan la minima energia) que tengan la minima distancia MST.
# Algoritmo:
# step 1: get the starting vertex and put that vertex into vector v
# step 2: find all the posible companions to v and for each companion compute Ec and Emst
# step 3: select the companion that minimize alpha*Ec + (1-alpha)*Emst
# step 4: put the selected vertex into v
# step 5: back to step 2 until all nodes are aggregated to v.
its <- number_of_vertexs - 1
for (vertex in 1:its) {
  cat(sprintf("Node: %s", vertex), "\n")
  id <- which(vertex_positions %in% v)
  a_combinar <- vertex_positions[-id]
  n <- length(a_combinar)
  temp <- matrix(NA, ncol=3, nrow=n) # primera columna Ec, segunda es Emst, tercera: suma de ambas
  for (i in 1:n) {
    new_vertex <- c(v, a_combinar[i])
    # generar el vector de estado:
    state <-  rep(0, length(H))
    id <- which(vertex_positions %in% new_vertex)
    state[id] <- 1 
    energia_coupl <- get_energy(t(as.matrix(state)), mode="c")
    #temp[i,] <- energia_coupl
    E_mst <- find_mst_barrier(mst_g, new_vertex)
    temp[i,] <- c(energia_coupl, E_mst, alpha*energia_coupl+(1-alpha)*E_mst)
  }
  id_ <- which(temp[,3] == min(temp[,3])) # encontramos id de temp donde hay minima energia de acople y de barreras (mst)
  E_mst <- temp[id_,2]
  v <- c(v, a_combinar[id_])
  
  #luego le encontramos a todos estos este estados el Ec
  state <-  rep(0, length(H))
  id <- which(vertex_positions %in% v)
  state[id] <- 1
  energia_coupl <- get_energy(t(as.matrix(state)), mode="c")
  tableau[vertex+1, ] <- c(v[length(v)], energia_coupl, E_mst, energia_coupl + E_mst)
}
# Find the next vertexs

# Analysis
tableau <- as.data.frame(tableau)
colnames(tableau) <- c("item", "Ec", "Emst", "Etotal")
#tableau$real_Emst <- dist_to_j(tableau$Emst)
#tableau$name <- row.names(tableau)
#myarrow=arrow(angle = 15, ends = "both", type = "closed")
myarrow=arrow(angle = 15, length=unit(0.3, "cm"))
#myarrow=arrow()
ggplot(subset(tableau), aes(-Ec, Emst)) + geom_point() + geom_path() +
  geom_segment(aes(xend=c(tail(-Ec, n=-1), NA), yend=c(tail(Emst, n=-1), NA)), arrow=myarrow, size=0.5, color="black") +
  geom_text(aes(label=item),hjust=1, vjust=-1.0, size=3) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     panel.background = element_blank(), axis.line = element_line(colour = "black")) 
# Analysis





# 07-mar-19
# Analysis but with alpha = 0.5 and alpha=0.7
tableau1 <- tableau
tableau2 <- tableau

tableau <- rbind(tableau1, tableau2)
tableau$alpha <- rep(0.5, nrow(tableau))
tableau[c(21:40), "alpha"] <- 0.7

ggplot(tableau, aes(-Ec, Emst, colour=factor(alpha))) + geom_point() + geom_path() +
  scale_color_manual(values=c("black", "blue", "#E69F00")) + 
  #geom_segment(aes(xend=c(tail(-Ec, n=-1), NA), yend=c(tail(Emst, n=-1), NA)), arrow=myarrow, size=0.5) +
  geom_text(aes(label=item),hjust=1, vjust=-1.0, size=2) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  guides(fill=FALSE, color=FALSE)





# 09-mar-19
# Analysis but with different starting points
tableau1 <- tableau # vinitial=80
tableau2 <- tableau # vinitial=72
tableau3 <- tableau # vinitial=214
tableau <- rbind(tableau1, tableau2, tableau3)
tableau$case <- c(rep("A", nrow(tableau)/3), rep("B", nrow(tableau)/3), rep("C", nrow(tableau)/3))

ggplot(tableau, aes(-Ec, Emst, colour=factor(case))) + geom_point() + geom_path() +
  scale_color_manual(values=c("black", "darkgoldenrod4", "darkorange3")) + 
  #geom_segment(aes(xend=c(tail(-Ec, n=-1), NA), yend=c(tail(Emst, n=-1), NA)), arrow=myarrow, size=0.5) +
  geom_text(aes(label=item),hjust=1, vjust=-1.0, size=2) +
  scale_y_continuous(breaks=seq(0,40,5)) + 
  scale_x_continuous(breaks=seq(0,40,5)) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  guides(fill=FALSE, color=FALSE)
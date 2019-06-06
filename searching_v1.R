# Lo que hacemos aqui es ir buscando las energias de acople Ec y las barreras de
# energia Emst, de tal forma de minimizar alphaX*Ec + (1-alpha)*Emst.
# 

# Para detalles ver el paper de conferencia a ICANN, y pag. 225-228 de los apuntes.



# Nombre del acutual archivo: searching<-v1.R
# ubicacion: dropbox->Research->Project_cluster_aggregation->DATA_ClusterAgr_ISING
# creation date: 06-mar-19
# companion file: 

#Notas:
# 06-jun-19: creation


# Algoritmo:
# 1. Se selecciona un vertice o nodo (o producto) de partida S = {nodo inicial}
# nota: el nodo inicial puede ser aquel con el mayor strength en el MST.
# 2. subrutina energy_coupling_search:
#         2.1 Para cada uno de los N-1 nodos restantes not in S, calcular la energia de acople
#         y tener el vector de energias de acoples para cada agregacion de producto. (temp_Ec)
# 3. subrutina de energy_mst_search:
#         3.1 Para cada uno de los N-1 nodos restantos not in S, calcular la energia de barreras
#         o distancia a lo largo del MST que una S con el nodo not in S. Obtener el vector
#         de energias mst para cada uno de los nodos. (temp_Emst)
# 4. sumar alpha*temp_Ec + (1-alpha)*temp_Ec para cada nodo y seleccionar nodo x que corresponda 
#         al minimo de alphaX*Ec + (1-alpha)*Emst.
# 5. S = S Union {nodo x}
# 6. repeat to 2 until S = {esten todos los nodos o hasta un numero determinado de nodos}




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



# # # # # # # # # # # # # # # # # # # # FUNCIONES # #  # # # # # # # # # # # # # # # # # # # # 
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

# # # # # # # # # # # # # # # # # # # # FUNCIONES # #  # # # # # # # # # # # # # # # # # # # # 




# # # # # # # # # # # # # # # # # # # # # FIND MST # # # # # # # # # # # # # # # # # # # # # #
# Conforming the COUPLING Network
colnames(J) <- colnames(wb)
rownames(J) <- colnames(wb)
#rownames(h) <- c(1:dim(J)[1])
#http://www.shizukalab.com/toolkits/sna/plotting-networks-pt-2
net <- graph_from_adjacency_matrix(as.matrix(J), mode="upper", weighted=TRUE, diag=FALSE)
E(net)$coupling <- E(net)$weight
# Conforming the MST Network
dis <- sqrt(-J + 3) 
# Convertir la matriz de distancia en objeto igraph
g <- graph_from_adjacency_matrix(dis, mode = "undirected", weighted = TRUE, diag = FALSE, add.colnames = NULL, add.rownames = NA)
E(g)$coupling <- E(net)$weight
mst_g <- minimum.spanning.tree(g, algorithm="prim")
edg <- as_edgelist(mst_g, names = TRUE)
# # # # # # # # # # # # # # # # # # # # # FIND MST # # # # # # # # # # # # # # # # # # # # # #


# Starting
v <- find_starting_vertex(mst_g) # basado en min alpha*Ec + (1-alpha)*Emst  v=80
###### subroutine energy_coupling_search








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

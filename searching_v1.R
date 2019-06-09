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
# 09 jun: he creado la funcion energy_coupling_and_mst_search que hace de una vez la 
#       subrutina energy_coupling_search y energy_mst_search y el resultado es una matriz
#       con todas las energias de acoples y de mst para la agregacion de cada uno de los 
#       nodos disponibles al vector de cluster v. 

# Procedure before the algorithm
# load de wb = matriz de canastas 179,610 X 25
# load J matriz de acoples 256 X 25
# create net = red completa de acoples con 25 nodos y 24X24/2 = 300 edges
# create mst_g = red mst de la red de acoples en que dis = sqrt (-J + 3)
# create nodo de partida.

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

load('new_inferring_parameters_environment270219.RData') # load H and J parameters ising inferred
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

# # # # # # # # # # # # # # # # # # energy_coupling_search# # # # # # # # ## # # # # # # # # #
# Esta funcion lo que hace es buscar el producto que incrementa la minima energia al cluster
# dado o contenido en v.

###### FUNCTION - subroutine coupling_and_mst_search
# Created on Jun 7, 2019.

# Lo que hace esta funcion es que para un vecttor de estado v (o cluster), se calculan la energia
# de acoples y energia mst de la agregacion de un producto mas al vector n. El resultado es una
# matriz con esas energias para cada uno de los restantes productos o nodos para postprocesamiento.
# input: vector de cluster v
#        tes = es un vector con los id de columnas de v que se utiliza despues. Siempre es el mismo
#        option = "all" si queremos la matriz de resultados con todos las combinaciones o "one" si queremos solo el minimo.
# output: item o nodo que produce la menor energia y la energia de acople correspondiente.
# Nota: la matriz de acoples J debe estar cargada en memoria.
energy_coupling_and_mst_search <- function(v, tes, option="one") {
  node_names <- colnames(J)
  ocupados <- which(v==1)    #ocupados = c(3,4)
  a_probar <- tes[! tes %in% ocupados] # dejamos los cols id de v, que estan vacios con ceros.
  a_probar_names <- node_names[!node_names %in% node_names[ocupados]] # nombre de los nodos que vamos a probar
  almacen <- matrix(NA, ncol=5, nrow=length(a_probar)) # aqui vamos dejando las energias de acople y el producto correspondiente
  almacen[,2] <- a_probar
  almacen[,1] <- as.numeric(a_probar_names)
  colnames(almacen) <- c("vertex_name", "vertex", "coupling_energy", "mst_distance", "mst_energy")
  for (i in 1:length(a_probar)) {
    v_temporal <- v # esto despues es el input
    id <- a_probar[i]
    v_temporal[id] <- 1
    #print(v_temporal)
    # subrutina para calcular la energia de acople
    almacen[i,3] <- get_energy(v_temporal, mode="c")
    # subrutina para calcular las distancias mst y energias mst
    id_temporal <- which(v_temporal==1)
    items <- as.factor(names[id_temporal])
    items <- as.numeric(as.character(items))
    r <- d_mst_distances(mst_g, vec = items)  
    almacen[i,4] <- r[[1]]
    almacen[i,5] <- r[[2]]
  }
  if (option == "one") {
    mini <- which(almacen[,3] == min(almacen[,3])) # retorna el minimo de energia de acople
    mini2 <- which(almacen[,5] == min(almacen[,5])) # retorna el minimo de energia mst 
    return( list( min_of_coupling =  almacen[mini,], min_of_mstenergy = almacen[mini2,]) )
  } else {
    return(almacen)
  }
}
# Ejemplo
#min_energy_coupling_search(v = v, tes = tes, option = "one")
#min_energy_coupling_search(v = v, tes = tes, option = "all")
# En el caso de "all" el resultado es la energia de acople de vector de entrada v con cardinalidad
# N, con los N-1 nodos posibles para combinar con v.
# # # # # # # # # # # # # # # # # # energy_coupling_search# # # # # # # # ## # # # # # # # # #
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


# Starting ------------------------07, jun, 2019
# create the vector-clustering
v <- matrix(0, ncol=ncol(wb), nrow=1)
colnames(v) <- colnames(wb)

# Create de initial node to start
item <- find_starting_vertex(mst_g) # basado en min alpha*Ec + (1-alpha)*Emst  v=80
i <- which(colnames(v) == item)
v[i] <- 1
Ec <- get_energy(v, mode="c") # funcion de entropies_functions.R / Matriz J de acoples debe estar cargada.
tes <- seq(1:ncol(v)) # un vector con los id de columnas de v que se utiliza despues.

# subrutina energy_coupling_search y subrutina y subrutina de energy_mst_search en la funcion energy_coupling_and_mst_search
r_energies <- energy_coupling_and_mst_search(v = v, tes = tes, option = "all")


# Step 4: sumar alpha*temp_Ec + (1-alpha)*temp_Ec para cada nodo y seleccionar nodo x que corresponda 
#         al minimo de alphaX*Ec + (1-alpha)*Emst.
alpha = 0.5
ponderation <- alpha*r_energies[,4] + (1-alpha)*r_energies[,5]
r_energies <-  cbind(r_energies, ponderation)
min_id <- which(r_energies[,6] == min(r_energies[,6]))
resultado <- r_energies[min_id, ]

# Step 5. S = S Union {nodo x}
v[resultado[2]] <- 1





























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

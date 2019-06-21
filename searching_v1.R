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
# 10 jun: repetimos el algoritmo para un punto de partida, y luego para cada uno de los
#       restantes nodos.

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
library(svMisc) # progress bar
source("entropies_functions.R")
source("network_functions_v1.R")
source("find_mst_barrier_function.R")
source("find_starting_vertex_function.R")

load('new_inferring_parameters_environment270219.RData') # load H and J parameters ising inferred
# # # # # # # # # # # # # # # # # # # # CARGA de datos # # # # # # # # # # # # # # # # # # # # 



# # # # # # # # # # # # # # # # # # # # FUNCIONES # #  # # # # # # # # # # # # # # # # # # # # 
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
  # tenemos que chequear si el vectro v de cluster ya esta lleno de puros unos
  # en ese caso hay que parar la funcion.
  if (length(ocupados) == length(v) ) {
    stop("The vector of clusters is complete. There are no more nodes to add")
    #e_acople <- get_energy(v, mode="c")
    #id_temporal <- which(v==1)
    #items <- as.factor(colnames(v)[id_temporal])
    #items <- as.numeric(as.character(items))
    #r <- d_mst_distances(mst_g, vec = items)
    #d_mst <- r[[1]]
    #e_mst <- r[[2]]
  } 
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
    items <- as.factor(colnames(v_temporal)[id_temporal])
    items <- as.numeric(as.character(items))
    r <- d_mst_distances(mst_g, vec = items)  ###########OJO cambiarla por find_mst_barrier en find_mst_barrier_function.R
    almacen[i,4] <- r[[1]] #distancia
    almacen[i,5] <- -1*r[[2]] #energia mst (multiplicamos por -1 porque si bien los acoples son positivos, la energia es negativa)
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
#E(net)$coupling <- E(net)$weight
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



# Starting the loop ------------------------10, jun, 2019
# Debemos cargar los datos, funciones, librerias y encontrar el MST.
# create the vector-clustering
v <- matrix(0, ncol=ncol(wb), nrow=1)
colnames(v) <- colnames(wb)

# Create de initial node to start
item <- find_starting_vertex(mst_g) # basado en min alpha*Ec + (1-alpha)*Emst  v=80
i <- which(colnames(v) == item)
v[i] <- 1
Ec <- get_energy(v, mode="c") # funcion de entropies_functions.R / Matriz J de acoples debe estar cargada.
tes <- seq(1:ncol(v)) # un vector con los id de columnas de v que se utiliza despues.
alpha = 0.5 # minimizar alphaX*Ec + (1-alpha)*Emst

# Comienza el loop recorriendo cada uno de los N-1 nodos
nodes <- ncol(v) - 1
final_result <- matrix(NA, ncol=3, nrow=nodes+1)
colnames(final_result) <- c("node", "Ec", "Emst")
final_result[1,] <- c(item, 0, 0) # we put here on the first row, the starting vertex.
for (node in 1:nodes) {
  # https://medium.com/human-in-a-machine-world/displaying-progress-in-a-loop-or-function-r-664796782c24
  progress(node, progress.bar = TRUE)
  Sys.sleep(0.01)
  if (node == nodes) cat("Done!\n")
  #print(node)
  # subrutina energy_coupling_search y subrutina y subrutina de energy_mst_search en la funcion energy_coupling_and_mst_search
  r_energies <- energy_coupling_and_mst_search(v = v, tes = tes, option = "all")
  #ponderation <- alpha*r_energies[,3] + (1-alpha)*r_energies[,5] # en col3 energia de acople, en col5 energia MST
  ponderation <- alpha*r_energies[,3] + (1-alpha)*r_energies[,4] # en col3 energia de acople, en col4 distancia MST
  r_energies <-  cbind(r_energies, ponderation)
  min_id <- which(r_energies[,6] == min(r_energies[,6]))
  resultado <- r_energies[min_id, ]
  v[resultado[2]] <- 1
  final_result[node+1,] <- c(resultado[1], resultado[3], resultado[5])
}
final_result


# como se interpreta final_result?
# Cada fila representa el nodo que se va agregando al cluster con su 
# respectiva Ec y Emst



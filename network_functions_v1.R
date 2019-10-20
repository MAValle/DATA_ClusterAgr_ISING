
# funciones necesarias para post-procesar signed networks, es decir,
# redes que edges positivos y negativos, como por ejemplo, la que 
# se encuentra de la maquina de boltzmann (boltzmann_machine_realdata_V1.R).
# Todas estas funciones requieren como input una signed netqork.

# Nombre del acutual archivo: network_functions_v1.R
# ubicacion: dropbox->Research->PAPER MST-Fondecyt_2->data
# date: 11.may.18


# Modificaciones: 

#Notas:
# 11.may.18:  Agregamos funciones:
#             get_signed_net, check_zero_degree, retrieve_attribute
#             que fueron creadas en BM_toy_asonam.R
# 07-jun-18: agregamos la funcion sampling_elist que viene de physicaA_v1.R
#             que nos permite generar sampleos de un edge_list maestro. Tambien
#             agregamos la funcion from_edges_to_widebasketV2 que viene 
#             de physicaA_v1.R que nos permite obtener la matriz de transaccion
#             de 0 y 1nos (wide_basket)
# 09-dic-18: agregamos la funcion d_mst_distances que calcula 
#             Calculamos D_mst de acuerdo a lo indicado en pg. 145 de mis apuntes como:
#           D_mst = sumatoria  (unique (edges(d12), edges(d13), edges(d23), ..... ) )
#           esto se utiliza en MST_for_physicaA.R.

library(igraph)

# get_signed_net -----
# Funcion
# Input: net: red net completa con acoples positivos y negativos
# Input: type = "p", o "n"
# Output: objeto igraph de red original con solo los pesos positivos o negativos
get_signed_net <- function(net, type="p") {
  if (type =="n") {
    red <- delete_edges(net, which(E(net)$weight>0)) # red es la red de acoples con solo pesos negativos
  } else {
    red <- delete_edges(net, which(E(net)$weight<0)) # red es la red de acoples con solo pesos positivos
  }
  # Ahora descomponemos la red en sus partes (pueden quedas subgrafos o nodos aislados)
  dg <- decompose.graph(red) # returns a list of three graphs
  n <- length(dg) # n es el numero de subgrafos
  vec <- check_zero_degree(dg) # chequeo de cual subgrafo de dg consisten en solo 1 nodo
  dg <- dg[vec]
  return(dg)
}
# Ejemplo
# p_net <- get_signed_net(net, "p")
# n_net <- get_signed_net(net, "n")


# NOTA: 05-ene-18: esta funcion no funciona, pero si la dejo como parte de un loop (sin funcion), ahi funciona.
# ver L.154 mst_distances_and_energies.R
# d_mst_distances -----
# Creada y testeada el 09-dic-18 en MST_for_physicaA.R
# Funcion para calular el D_mst segun pag.145
# Input: mst_g: el MST
#         Es importante que mst_g venga con nombres de nodos V(mst_g)$name y con los acoples E(mst_g)$coupling
# Input: vec:  es un dataframe tipo "factor" que identifica el estado de subcategorias de productos activas
d_mst_distances <- function(mst_g, vec) {
  # como vec es una dataframe de 1 fila y con valores factor, tenemos que pasarlo a numero:
  indx <- sapply(vec, is.factor)
  vec[indx] <- lapply(vec[indx], function(x) as.numeric(as.character(x)))
  vec <- as.numeric(vec) # lo pasamos a vector numerico
  temp <- combn(vec,2)
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
    
    D_mst <- sum(E(mst_g)$weight[unique(calce)])
    E_mst <- sum(E(mst_g)$coupling[unique(calce)]) # esta es la suma de energias de acoples en el MST.
  }
  return(list(mst_distance=D_mst, mst_energy=E_mst))
}
#ejemplo
# d_mst_distances(mst_g, vec)
#d_mst_distances(mst_g, active_vertex_sample)


# 13-sep-19
# find_mst_barrier -----
# # # # # # # # # # # # # # # # # # # # # FIND Emst # # # # # # # # # # #  # # # # # # # # # #
# Inspirado en la idea de la L164 en adelante en mstdistances_and_energies.R
# esta funcion fue tomada de mst_distances_and_energies_v2.R
# Es una funcion que nos entrega la sumatoria de todas las distancias en el MST dado un set
#     de spins activos.
# Reemplaza la funcion d_mst_distances que no funciona bien (ver mas arriba)
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







# check_zero_degree -----
# Funcion para borrar subgrafos con nodos con grado 0.
# input: dg es una lista son subgrafos
# output: vector de largo igual al numero de subgrafos de dg que indica cual hay que eliminar
check_zero_degree <- function(dg) {
  n <- length(dg) # n es el numero de subgrafos
  salida <- !logical(length=n)
  for (i in seq_along(1:n) ) {
    subg <- dg[[i]]
    elist <- get.edgelist(dg[[i]])
    dimen <- dim(elist)
    if (dimen[1] == 0) {
      salida[i] <- FALSE
    }
  }
  return(salida)
}
# Ejemplo:
# vec <- check_zero_degree(dg)


# retrieve_attribute -----
# Funcion para rescatar los weigths de los edges y las magnitudes de los vertex
# de la union de varios subgrafos
# Input: la red con varios subgrafos unidas: totalnet
# Output: una lista con vector de weigths y magnitudes 
retrieve_attribute <- function(totalnet) {
  # Retrieving the weights of the edges
  atributos_edges <- edge_attr_names(totalnet)
  atributos_edges_de_interes <- atributos_edges[grepl("weight_", atributos_edges)]
  n <- length(atributos_edges_de_interes)
  num_edges <- length(edge_attr(totalnet, atributos_edges_de_interes[1]))
  temporal <- matrix(, ncol=num_edges, nrow=n)
  for (i in seq_along(1:n)) {
    temporal[i,] <- edge_attr(totalnet, atributos_edges_de_interes[i])
  }
  salida_edges <- colSums( temporal, na.rm=TRUE)
  # Retrieving the magnitudes of the vertexs
  atributos_vertex <- vertex_attr_names(totalnet)
  atributos_vertex_de_interes <- atributos_vertex[grepl("magn_", atributos_vertex)]
  n <- length(atributos_vertex_de_interes)
  num_vertex <- length(vertex_attr(totalnet, atributos_vertex_de_interes[1]))
  temporal <- matrix(, ncol=num_vertex, nrow=n)
  for (i in seq_along(1:n)) {
    temporal[i,] <- vertex_attr(totalnet, atributos_vertex_de_interes[i])
  }
  salida_vertex <- colSums( temporal, na.rm=TRUE)
  # Retrieving the colors of the vertexs
  atributos_vertex_de_interes <- atributos_vertex[grepl("color_", atributos_vertex)]
  n <- length(atributos_vertex_de_interes)
  num_vertex <- length(vertex_attr(totalnet, atributos_vertex_de_interes[1]))
  temporal <- matrix(, ncol=num_vertex, nrow=n)
  for (i in seq_along(1:n)) {
    temporal[i,] <- vertex_attr(totalnet, atributos_vertex_de_interes[i])
  }
  bindeo <- pmin(temporal[1,], temporal[2,], na.rm = TRUE)
  if (n > 2) {
    N <- n - 2
    for (i in seq_along(1:N)){
      cont <- i + 2
      bindeo <- pmin(bindeo, temporal[cont, ], na.rm=TRUE)
    }
  }
  salida_vertex_color <- bindeo
  return(list(salida_edges, salida_vertex, salida_vertex_color))
}
# Ejemplo
# att <- retrieve_attribute(p_net_union)
# E(p_net_union)$weight <- att[[1]]
# V(p_net_union)$magn <- att[[2]]
# V(p_net_union)$color <- att[[3]]


# get_signed_net_ready -----
# Funcion que nos entrega la red positiva o negativa, con los atributos
# weight y magn listos para ser graficados y analizados topologicamente
# Input: signed network con edge attributes weight y vertex atributes magn
# Output: positive or negative network with the weight and magn attributes
get_signed_net_ready <- function(net, type="p") {
  signed_net <- get_signed_net(net, type=type)
  if (length(signed_net) == 1 ) {
    signed_net <- signed_net[[1]]
  } else {
    signed_net <- graph.union(signed_net, byname=T)
    att <- retrieve_attribute(signed_net)
    E(signed_net)$weight <- att[[1]]
    V(signed_net)$magn <- att[[2]]
    V(signed_net)$color <- att[[3]]
  }
  return(signed_net)
}
# Example:
# p_net <- get_signed_net_ready(net, type="p")




# sampling_elist -----
# Esta funcion few creada en physicaA_v1.R
# Function sampling_list: genera una lista con varios edges_list con sampleos 
# a un subconjunto de productos que vienen de una edge list maestro.
# inputs
# elist_basket : el edge list maestro de toda la base de datos transaccional
#   vector items = vecttor de productos a analizar (usualmente los 21 de mayor frecuencia de compra)
#   n = tamano de la muestra Mi que deseamos
#   m = numero de muestras que queremos (el maximo sera mCn )
# outputs
#   una lista con las las muestras
sampling_elist <- function(elist_basket, items, n=10, m=15) {
  library(dplyr)
  library(progress)
  # 1. primero filtrar elist_basket para que solo contenga transacciones en items
  # 2. luego seleccionar aleatoriamente de items N productos 
  # 3. luego volver a filtrar para que aparezcon solo los N productos. guardar ese elist.
  # 4. volver a 2 cuantas veces sea necesario (poner semilla para replicar resultados)
  result = vector("list", m) 
  #pb <- progress_bar$new(format = "  Progress [:bar] :percent eta: :eta", total = m, clear = FALSE, width= 60) 
  for (muestra in seq_along(1:m)) {
    # progress bar
    #pb$tick()
    #Sys.sleep(1/m)
    # 1. https://blog.exploratory.io/filter-data-with-dplyr-76cf5f1a258e
    elist_basket_ <- elist_basket %>%
      filter(P1 %in% items) %>%
      filter(P2 %in% items)
    # Seleccionar de items aleatoriamente n productos de items
    #set.seed(123)
    id <- sample(1:length(items), rep = FALSE, n)
    sample_items <- items[id] # algo importante a guardar
    # 3. Volvemos a filtrar
    elist_basket_sample <- elist_basket_ %>%  #algo importante a guardar.
      filter(P1 %in% sample_items) %>%
      filter(P2 %in% sample_items)
    result[[muestra]] <- elist_basket_sample
    #result <- c(result, elist_basket_sample)
  }
  return(result)
  #https://stackoverflow.com/questions/2436688/append-an-object-to-a-list-in-r-in-amortized-constant-time-o1
}
# ejemplo:
#set.seed(123)
#elists <- sampling_elist(elist_basket, items, n=12, m=80)
# head(elists)
# head(elists[[2]])




# from_edges_to_widebasketV2 -----
# Esta funcion few creada en physicaA_v1.R
#  Encontrar el wide_basket (matriz de transacciones) y correlaciones para cada Mi 
# Funcion para transformar el edge_list en un wide_basket 
# Tengo de input una edge listo del tipo:
# P1    P2    newid
# 28    196   1
# 28    25    1
# 196   25    1
# y necesito convertirlo a una transaction data matrix del tipo:
# newid   28  25  196   33    50    ....
# 1       1   1   1     0     0   ..
# 2 .....
# Step 1:
# subset newid del edge list
# Step 2:
# para cada newid del edge list, unir vector de P1 y P2  = F
# y luego Encontrar los unique de F  = productos_union
# Step 3:
# poner valor 1 a todos las columnas donde corresponda J
# Inputs: 
#     elists: lista de los edge list que vienen de la funcion sampling_elist
#     items:  numero total de productos a aser considerados.
# Outputs: wide_basket (matriz de transacciones de o y 1nos)
# version de la funcion from_edges_to_widebasket con datatable
from_edges_to_widebasketV2 <- function(items, edge_list) {
  library(data.table)
  losids <- unique(edge_list$newid)
  filas <- length(losids)
  # El numero de columnas de wide_basket esta dado por el vector de entrada "items"
  columnas <- length(items)
  wide_basket <- matrix(0, ncol=columnas+1, nrow=filas)
  wide_basket <- as.data.frame(wide_basket)
  colnames(wide_basket) <- c("newid", items)
  #rownames(wide_basket) <- c(1:filas)
  wide_basket <- setDT(wide_basket) #convertimos la matrix a datatable
  edge_list <- setDT(edge_list) #convertimos el dataframe a datatable
  setkey(edge_list, newid) # sorts the datatable by P1 and P2 in ascending order
  for (i in seq_along(1:filas)) {
    # Step 1:
    # subset newid del edge list
    idx <- losids[i] #esto se hace para cada newid
    #sub_edges <- subset(edge_list, newid==idx)
    sub_edges <- edge_list[newid == idx]
    
    # Step 2:
    # para cada newid del edge list, unir vector de P1 y P2  = F
    #productos_union <- c(sub_edges$P1, sub_edges$P2) 
    productos_union <- unique(c(sub_edges$P1, sub_edges$P2) )
    
    # Step 3:
    # poner valor 1 a todos las columnas donde corresponda 
    #rownames(wide_basket)[i] <- idx
    #w_cols <- which(colnames(wide_basket) %in% productos_union )
    #wide_basket[i, w_cols] <- 1
    productos_union <- as.character(productos_union)
    #wide_basket[i, c(productos_union) := list(rep(1, length(productos_union))) ]
    wide_basket[i,1] <- idx
    wide_basket[i, c(productos_union) := 1]
  }
  return(wide_basket)
}
# Ejemplo
#system.time( { wide_basket <- from_edges_to_widebasketV2(items=items, edge_list = elists[[50]]) } )# demora app 6.3 minutos!

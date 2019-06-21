# # # # # # # # # # # #  # # # # # # FIND strongest vertex # # # # # # # # ## #  # # # # # # #
# # FUNCTION
# June, 11, 2019
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


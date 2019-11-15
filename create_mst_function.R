

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
# 23.sep.19 FUNCION creada en simulation_hc_v1.R
# # # # # # # # # # # # # # # # # # # # # FIND MST # # # # # # # # # # # # # # # # # # # # # #
# Funcion que crea como salida la red de acoples y el mst de la red de acoples.
# imput: matriz de distancia D y de acoples J
create_mst <- function(J, D) {
  library(igraph)
  net <- graph_from_adjacency_matrix(as.matrix(J), mode="upper", weighted=TRUE, diag=FALSE)
  dis <- D
  # Convertir la matriz de distancia en objeto igraph
  g <- graph_from_adjacency_matrix(dis, mode = "undirected", weighted = TRUE, diag = FALSE, add.colnames = NULL, add.rownames = NA)
  E(g)$coupling <- E(net)$weight # asignamos las energias de acoples de net a E(g)$coupling
  mst_g <- minimum.spanning.tree(g, algorithm="prim")
  return(mst_g)
}
# ejemplo:
# mst_g <- create_mst(J=J, D=D)


# # # # # # # # # # # # # # # # # # # # # FIND MST # # # # # # # # # # # # # # # # # # # # # #
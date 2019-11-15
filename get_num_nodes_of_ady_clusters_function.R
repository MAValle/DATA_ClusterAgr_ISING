

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
# 27.oct.19 FUNCION 
# funcion que nos arroja el numero de nodos que tiene los
# clusters adyacentes a un par de clusters:
# identificamos los nodos adyacentes en el MST de los clusters A y B y su numero de nodos
# creamos matriz temporal que tiene: cluster / node_ady / cluster_ady / num_nodes_clust_ady
# donde cluster_ady es el nombre del cluster en donde se encuentra node_ady, y num_nodes_clust_ady
# es el numero de nodos de cluster_ady .

# input: vector con los nombres de los clusters a los cuales
# queremos encontrar el numero de clusters adyacentes a ellos.
# input: mst_g
# input: lista con los nodos que integran cada clusters
# output: matriz con: cluster / node_ady / cluster_ady / num_nodes_clust_ady
# donde cluster_ady es el nombre del cluster en donde se encuentra node_ady, y num_nodes_clust_total
# num_nodes_clust_total = numero de nodos del cluster + numero de nodos del cluster adyacente a cluster.
# es el numero de nodos de cluster_ady .


get_num_nodes_of_ady_clusters <- function(name_cls, lista, mst_g) {
  register <- matrix(NA, ncol = 4, nrow = 0)
  for (i in 1:2) {
    cl <- name_cls[i]
    number_of_nodes_cluster <- get_num_nodes(name_cls = cl, lista = lista)
    if (cl > 0) {
      # identificar nodos adyacentes a cl (considerando cl > 0):
      nodos_de_cl <-  lista[[ cl ]] 
      ady <- adjacent_vertices(mst_g, nodos_de_cl )
      ady <- -1*as.numeric(  unlist(ady)   )
      # botamos los nodos de ady que son parte de cl
      #id_drop <-  which(as.character(ady)  %in% nodes_de_cl )   
      #ady <- ady[-id_drop]
      num_nodos_ady <- length(ady) # numero de nodos adyacentes al cluster 
      for (j in 1:num_nodos_ady) {
        namen <- get_name_of_the_cluster( nady = ady[j], lista = lista )
        numen <- get_num_nodes(name_cls = namen, lista = lista) + number_of_nodes_cluster
        filadata <- c(cl, ady[j], namen, numen)
        register <- rbind(register, filadata)
      }
    }
    if (cl < 0) {
      # identificar nodos adyacentes a cl (considerando cl < 0)
      ady <- adjacent_vertices(mst_g, as.character(cl) )
      ady <- -1*as.numeric(  unlist(ady)   )
      num_nodos_ady <- length(ady) # numero de nodos adyacentes al cluster 
      for (j in 1:num_nodos_ady) {
        namen <- get_name_of_the_cluster( nady = ady[j], lista = lista )
        numen <- get_num_nodes(name_cls = namen, lista = lista) + number_of_nodes_cluster
        filadata <- c(cl, ady[j], namen, numen)
        register <- rbind(register, filadata)
      }
    }
  }
  colnames(register) <- c("cluster", "node_ady", "cluster_ady", "num_nodes_clust_total")
  return(register)
}

# Ejemplo:
# numeros <- get_num_nodes(name_cls = nombres_num, lista=mn)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  

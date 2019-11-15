

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
# 28.oct.19 FUNCION 
# funcion que nos arroja todos los nodos de un determinado cluster

# input: name_cluster = numero positivo o megativo que indica el cluster
# input: lista con los nodos de cada cluster agregado
# output: numero que indica todos los nodos involucrados en el cluster (formato numerico)

get_nodes_of_the_cluster <- function(name_cluster, lista) {
  largo <- length(lista)
  thenodes <- c()
  if (name_cluster < 0) {
    thenodes <- name_cluster
    return(thenodes);
  }
  
  thenodes <- as.numeric( lista[[name_cluster]] )
  
  
  return(thenodes)
}

# Ejemplo:
# numeros <- get_nodes_of_the_cluster(name_cluster = 1, lista=mn)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  

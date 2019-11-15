

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
# 27.oct.19 FUNCION 
# funcion que nos arroja el numero de nodos que tienen dos
# clusters:
# input: name_cls = vector con el nombre de los "clusters".  Debe ser positivo.
# input: lista con los nodos de cada cluster agregado
# output: vector con el numero de nodos de cada clusters

get_num_nodes <- function(name_cls, lista) {
  # if (is.na (name_cls) ) { # no 
  #   num_cls <- 0
  # } 
  # if ( name_cls < 0) {
  #   num_cls <- 0
  # }
  largo <- length(name_cls)

  num_cls <- vector(mode="numeric", length=largo)
  for (i in 1:largo) {
    cl <- name_cls[i]
    if (cl > 0) { # el cluster tiene mas de 1 nodo
      num_cls[i] <- length( lista[[ cl ]] )
    } else {
      num_cls[i] <- 1
    }
  }
  
  return(num_cls)
}

# Ejemplo:
# numeros <- get_num_nodes(name_cls = nombres_num, lista=mn)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
# 27.oct.19 FUNCION 
# funcion que dado el nombre de un nodo, nos dice a que cluster pertenece
# Es importante que el nombre del cluster debe ser positivo
get_name_of_the_cluster <- function(nady, lista) {
  if (length(lista) == 0 ) {
    thename <- nady 
    return(thename);
  }
    save_r <- c()
    for (n in 1:length(lista) ) {
      if ( nady %in% lista[[n]] ) {
        save_r <- c(save_r, n)
      } 
    }
    if ( length(save_r) == 0 ) { # el cluster adyacente no pertenece a un cluster.
      thename <- nady 
    } else {
      thename <- max(save_r)
    }

  return(thename)
}

# Example
# namen <- get_name_of_the_cluster( nady = 3, lista = mn )
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  

# Version antigua.
# get_name_of_the_cluster <- function(nady, lista) {
#   if ( nady > 0) {
#     thename <- nady 
#   }
#   
#   if ( nady < 0 ) {
#     save_r <- c()
#     for (n in 1:length(lista) ) {
#       if ( nady %in% lista[[n]] ) {
#         save_r <- c(save_r, n)
#       } 
#     }
#     if ( length(save_r) == 0 ) { # el cluster adyacente no pertenece a un cluster.
#       thename <- nady 
#     } else {
#       thename <- max(save_r)
#     }
#   }
#   return(thename)
# }

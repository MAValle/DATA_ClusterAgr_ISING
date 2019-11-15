

# Funcion necesaria para funcion hierarchical_clustering_greedy
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 06-nov-19
# Esta funcion es analoga a find_min_dist, pero se acondiciona para trabajar con las
# DISTANCIAS DE ACOPLES, en vez de las distancias MST.
# pequena funcion para re-actualizar distancias minimas entre un cluster recien fusionado
# y los demas clusters. 
# input: 
#     J = matriz de acoples
#     D = matriz de distancias actualizada.
#     Dold = matriz de distancias de la iteracion anterior valores en character.
#     nombres : nombres de los nodos fusionados en la iteracion
#     mn = lista de los nodos de cada cluster que se han ido generando
#     it = numero de iteracion actual
# output: dit = vector con las minimias distancias de cluster fusionado en la iteracion anterior
#         a todos los demas clusters que quedan. Se utiliza single linkage.
# Nota: esta funcion funciona solamente cuando D es de 4X4 hacia arriba, es decir cuando hay más de 9 o mas elementos en D.
# Nota: el vector "nombres" debe estar cargado en memoria.
find_min_distcp <- function(J, D, nombres, mn, it) {
  dit <- vector()
  nombres_p <- colnames(D) # vemos los nombres de las filas o columnas que quedan en D
  n <- length(nombres_p)
  for (i in 1:n) {
    nodo <- nombres_p[i]
    nods <- c(nodo, mn[[it]] )
    # tenemos que chequear si en nods hay algun cluster con mas de un nodo:
    tempo <- as.numeric(nods)
    haypositive <- which(tempo > 0)
    if ( length(haypositive) > 0) {
      masnodos <- vector()
      for (i in 1:length(haypositive) ) {
        masnodos <- c( masnodos, mn[[ tempo[ haypositive[i] ] ]] )
      }
      nods <- c(masnodos, mn[[it]])
    } 
    # fin chequeo # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    v <-  acople_distance_sum2(J, y = nods) 
    dit <- c(dit,v )
  }
  return(list(dit=dit, nombres_p = nombres_p))
}
# ejemplo:
# temp <- find_min_dist(J, D, Dold,  nombres = as.character(nombres_num), mn = mn)
# dit <- temp$dit
# nombres_ = temp$nombres_
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
# 13.sep.19 FUNCION 
# Funcion que calcula la suma de distancias de acoples
# Por ejemplo, para sumar acoples a, b y c, se hace: f(a) + f(b) + f(c)
# donde f(*) es funcion que comvierte el acoples en una medidda de distancia.
# por ejemplo f(x) = sqrt(3 - x). 
# input: J = matriz de acoples J
# input: y = vector con todos los nombres de los spins 
# output: suma de f(a) + f(b) + f(c) + .... Si son N spins la suma sobre n(n-1)/2 acoples.
acople_distance_sum <- function(J, y) {
  sumdist <- 0
  combinaciones <- combn(y,2)
  lk <- ncol(combinaciones)
  for (i in 1:lk) {
    idx <- combinaciones[, i]
    cupl <- J[idx[1], idx[2]]
    dcupl <- sqrt(3 - cupl)
    sumdist <- sumdist + dcupl
  }
  return(sumdist)
}
# ejemplo
# si y = c("-4" "-2" "-1")
# sumdist <- acople_distance_sum(J, y)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
#sep, 23, 2019
# input: matriz de acople J, y y = vector con el nombre de los nodos en formato numerico
# lo mismo que en acople_distance_sum , dado y, el vector de nombres con los spins
# utilizamos el acceso a la matriz J, pro los nombres.
acople_distance_sum2 <- function(J, y) {
  sumdist <- 0
  combinaciones <- combn(y,2)
  lk <- ncol(combinaciones)
  for (i in 1:lk) {
    idx <- combinaciones[, i]
    cupl <- J[as.character(idx[1]), as.character(idx[2])]
    dcupl <- sqrt(3 - cupl)
    sumdist <- sumdist + dcupl
  }
  return(sumdist)
}
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  

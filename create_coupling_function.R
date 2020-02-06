

# 24-oct-19 : modificamos funcion apra que simula J truncados entre -inf e inf
# 18-dic-19: modificamso funcion para incorporar generar numero aleatorios de una distribucion uniforme

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
# 23.sep.19 FUNCION creada en simulation_hc_v1.R
# # # # # # # # # # # # # # # # # # # # CREATING A J MATRIX # # # # # # # # # # # # # # # # # # # # # # #
# Funcion que crea matriz de acople y de distancia con media mu, desv estandar de acople sj y con Nn nodos
# inputs: 
# kind = tipo de distribucion : puede ser "normal" o "uniform"
# lw = limite inferior para J
# up = limite superior para J
# Nn = numero de nodos de la matriz o numero de nodos de la red

create_coupling <- function(Nn, media, sj, lw, up, kind="normal") {
  library(Matrix)
  library(msm)
  if (kind == "normal") {
    #x <- Matrix(rnorm(Nn*Nn, mean=media, sd=sj),Nn)  # los J vienen con media 0 y desv estandar 1
    x <- Matrix(rtnorm(Nn*Nn, mean=media, sd=sj, lower=lw, upper=up), Nn)
  }
  if (kind == "uniform") {
    x <- Matrix(runif(Nn*Nn, min=lw, max=up ), Nn)
  }
  J <- forceSymmetric(x) # couplings
  D <- sqrt(3 - as.matrix(J)) # distance
  colnames(D) <- rownames(D) <- -c(1:Nn) # enumeramos todos las hojas o spines originales con valones negativos.
  diag(D) <- rep(Inf, Nn)
  colnames(J) <- rownames(J) <- -c(1:ncol(J)) 
  J <- as.matrix(J)
  return(list(J=J, D=D))
}
# ejemplo
# ot <- create_copupling(Nn=6, media=0, sj=1, lw=-3, up=3)
# J <- ot$J
# D <- ot$D

# # # # # # # # # # # # # # # # # # # # CREATING A J MATRIX # # # # # # # # # # # # # # # # # # # # #


# funciones necesarias para ejecutar inverse ISING problem.
# mediante una maquina de Boltzmann without hidden units.
# en compania de annealing_ising_toy2.R

# Nombre del acutual archivo: ising_functions_v3.R
# ubicacion: dropbox->Research->PAPER MST-Fondecyt_2->data
# date: 05-dic-17

# Modificaciones: 

#Notas:
# 11-dic-17:  Agregamos funcion kullback, que computa el kullback_leibler 
#             divergence en su forma tradicional.
# 13-dic-17: Agregamos funcion contrastive_divergence3 que calcula el J(e+1)
#            de acuerdo al paper de Acklof de 1985.
# 20-dic-17: Creamos funcion perf que nos grafica Cij real en eje X contra Cij estimado 
#           en el eje Y. 
# 22-dic-17: creamso funcion kullback3, que hace lo mismo que kullback, pero ahora
#           tenemos la opcion de considerar un porcentaje de los estados (ej.0.8) 
#           para poder apurar el culculo de KL especialmente cuando hay mas de 15 nodos.
# 25-dic-17: Hemos modificado las funciones goodness_activation-energy, sweep y metropolis-hasting
#            para incorporar las magnetizaciones.
# 28-dic-17: funcion que genera un grafico en ggplot entre las correlaciones Cij (eje X) y
#            los acoples Jij.
# 20-abr-18: agregamos funcion para calcular los two-body connections Cij = <sisj> - <si><sj>
#            se llama pairwise_cij. Tambien agregamos la funcion rmse que extrae el RMSE de
#            una regresion lineal.
# 20-abr-18 Creamos una funcion para traspasar la matriz de transacciones al edge-list.
#           La funcion se llama 
# 19-may-18: Ponemos la funcion bmachine que lleva a cabo el boltzmann learning
#           utilizando varias otras funciones (inicital_conditions, metropolis-sampling
#           pairproducts, contrastive_divergence2 y kullback3) que estan en este script
#           y que fue tambien utilizado en BM_toy_asonam.R
# 23-jun-18: creamos frecuencyV2 que es lo mismo que la funcion frecuency, pero mas rapida.


# # # # # # # # # # # # # # # # # CONDICIONES INCIALES # # # # # # # # # # # # # # # # # # # 
#04-ago-17
# Alternativas de condicionales iniciales de los acoples y/o magnetizaciones - tipe:
# a1 = todos igual a cero.
# a2 = todos igual al vector de medias para h0 y a la mtriz de correlaciones para J
# a3 = las magnetizaciones iniciales h0 son igual a log(pi/(1-pi)) donde pi es la proporcion
#      de spines 1 en la variable i. Los acoples J0 igual a la matriz de correlaciones
# a4 = las magnetizaciones iniciales son iguales a log(pi((1-pi))) donde pi es la prob. que 
#      spin este activada (Si=1). Los acoples son iguales a 0. Se necesita como input la 
#      la matriz de datos 
# a5 = la matriz de acople son valores aleatorios de una gaussiana con media igual a 0 y 
#      desv. estandar igual a desv. y las magnetizaciones todas igual a cero.

#INPUTS:
#  * randomon = obtiene un valor aleatorio de un valor de una gaussiana con media igual al valor 
#  y desv. estandar igual a 0.05 (on, off)
#  * si y sij = vector sigma_i_data y matriz sigma_ij_data de magnetizaciones y correlaciones de la data
#  * df = data original con spines 1 y -1.
#  * N = numero de nodos o variables
# * desv = desviacion estandar en caso que tipe = a5

#OUTPUTS:
#  * ho y J0 : vector y matriz inicial de magnetizacion y acoples.
initial_conditions <- function(df=NULL, N, tipe, si, sij, desv) {
  require(Matrix)
  if (tipe == 'a1') {
    h0 <- t(as.matrix(matrix(0L, nrow=1, ncol=N)))
    J0 <- matrix(0L, nrow=N, ncol=N) 
  } else if (tipe == 'a2') {
    h0 <- si
    J0 <- sij
  } else if (tipe == 'a3') {
    df = 0.5*df+0.5
    p = colSums(df)/nrow(df)
    h0 <- log(p/(1-p))
    J0 <- matrix(0L, nrow=N, ncol=N) 
  } else if (tipe == 'a4'){
    h0 <-matrix(0L, nrow=1, ncol=N)
    for (cl in c(1:N)) {
      s <- sum(df[,cl] == 1)
      p <- s/nrow(df)
      h0[1, cl] <- log(p/(1-p))
    }
    h0 <- t(as.matrix(h0))
    J0 <- matrix(0L, nrow=N, ncol=N) 
  } else if (tipe == 'a5') {
    h0 <- t(as.matrix(matrix(0L, nrow=1, ncol=N)))
    J0 <- matrix(runif(N*N, -1, 1),N )
    ind <- lower.tri(J0) 
    J0[ind] <- t(J0)[ind] 
    diag(J0) <- NA
  }
  return(list(h0, J0))
}
# # # # # # # # # # # # # # # # # CONDICIONES INCIALES # # # # # # # # # # # # # # # # # # # 





# # # # # # # # # # # ## # # # TRANSITION PROBABILITY # ## # # # # ## # # # # ## # # # # ##  # # # # ## 
# Funcion para calcular la probabilidad de que un spin Si cambie de estado P(Si = 1) mientras
# que P(Si = 0) = 1 - P(Si = 1)
# Date: 05-12-17
# # # # # # # # # # # ## # # # TRANSITION PROBABILITY # ## # # # # ## # # # # ## # # # # ##  # # # # ## 
#Inputs: 
#  vector de entrada = spins
#  J matriz de acoples = J
#  h vector de magnetizaciones = H
#  spin <- numero de neurona que se desea calcular la probabilidad
#Output: activation energy  

# 11-DIC-17 ACLARACION:
# Si Delta E = E(xi=1) - E(xi=0). Este Delta E es la energia de activacion.
# Si Delta E = E(xi=0) - E(xi=1). Este Delta E es la energia de desactivacion.
activation_energy <- function(H, acoples, spins, idspin) {
  # computo de energia de activacion:
  rows <- length(spins)
  id <- combn(1:rows,2)
  #N_par <- nrow(id) #numero de parametros
  N_par = dim(id)[1] #25-12-17
  #N_units <- nrow(acoples)
  N_units <- dim(acoples)[1] #25-12-17
  E <- H[idspin] #valor de la energia de activacion
  losj <- seq_along(1:N_units)
  losj <- losj[-(idspin)]
  for (j in losj) {
    E <- E + acoples[idspin,j]*spins[j]
  }
  return(E)
}
# # # # # # # # # # # ## # # # TRANSITION PROBABILITY # ## # # # # ## # # # # ## # # # # ##  # # # # ## 
#input: 
#  energy = activation or desactivation energy
activation_probability <- function(energy, Temp) {
  p <- 1/( 1 + exp(-energy/Temp) )
  return(p)
}
# # # # # # # # # # # ## # # # TRANSITION PROBABILITY # ## # # # # ## # # # # ## # # # # ##  # # # # ## 
# Ejemplo
#activation_probability(activation_energy(H=c(0,0,0,0,0,0,0), J=I, spins=c(0,1,1,0,1,0,0), idspin=2), T=NULL)
# # # # # # # # # # # ## # # # TRANSITION PROBABILITY # ## # # # # ## # # # # ## # # # # ##  # # # # ## 

# # # # # # # # # # # ## # # # METROPOLIS-HASTING # ## # # # # ## # # # # ## # # # # ##  # # # # ## 
#INPUT: 
#   ACOPLES: matriz de acoples J del momento
#   initial_spin (pueden ser todos cero, todos 1, o aleatorios)
#   kmax = maximo numero de iteraciones sweep de MH
#   To = temperatura inicial. Si To = NULL, entonces T = 1. 
#   ks = numero de iteraciones adicionales en temperatura T = 1.
#OUTPUTS:
  # matriz de spines
  # vector de energias de cada sweep

# sweep function
#input: vector de estados del momentos (spines) y matriz de acoples J
sweep <- function(magnetizaciones, acoples, vector, Temp) {
  N_nodos <- length(vector)
  v <- vector
  for (i in seq_along(1:N_nodos)) {
    logico <- as.logical(v)
    #energy_activation <- activation_energy(H=rep(0, N_nodos), acoples=acoples, spins=v, idspin=i)
    energy_activation <- activation_energy(H=magnetizaciones, acoples=acoples, spins=v, idspin=i) #25-12-17
    Pon <- activation_probability(energy_activation, Temp=Temp)
    Poff <- 1 - Pon
    u <- runif(1)
    if (logico[i]) {
      if (Poff > u) { 
        logico[i] <- !logico[i] 
      } 
    } else {
      if ( Pon > u) { 
        logico[i] <- !logico[i] } 
    }
    v <- 1*logico
    #print(v)
  }
  return(v)
}
metropolis.sampling <- function(magnetizaciones, acoples, initial_spin, kmax, ks, To=NULL, alpha) {
  output <- matrix(NA, ncol=length(initial_spin), nrow=kmax+ks)
  energies <- matrix(NA, ncol=2, nrow=kmax+ks)
  k <- 1
  vec <- initial_spin
  if ( is.null(To) ) { To = 1 ; alpha = 1 } 
  Temp <- To
  #sampleo annealing
  while (k <= kmax) {
    spines <- sweep(magnetizaciones=magnetizaciones, acoples = acoples, vector = vec, Temp = Temp)
    output[k, ] <- spines
    energies[k, ] <- c(k, goodness(spins = spines, magnetizaciones=magnetizaciones, acoples = acoples))
    vec <- spines
    Temp <- tempfunc(Temp, alpha = alpha)
    k <- k + 1
  }
  #sampleo con temperatura T = 1.
  if ( ks > 0 ) {
    while (k <= kmax + ks) {
      spines <- sweep(magnetizaciones=magnetizaciones, acoples = acoples, vector = vec, Temp=1)
      output[k, ] <- spines
      energies[k, ] <- c(k, goodness(spins = spines, magnetizaciones=magnetizaciones, acoples = acoples))
      vec <- spines
      #Temp <- tempfunc(Temp, alpha=0)
      k <- k + 1
    }
  }

  return(list(output, energies))
}
#ejemplo
#mh <- metropolis.sampling(acoples = I, initial_spin = rep(1,7), kmax=20000, ks=0, To=NULL, alpha=NULL)
# si quiero las ultimas Kcut muestras: muestras <- ss[[1]]; muestras <- muestras[-(1:Kcut), , drop=FALSE]
#ss[[2]] #energias al final de cdada sweep
# # # # # # # # # # # ## # # # METROPOLIS-HASTING # ## # # # # ## # # # # ## # # # # ##  # # # # ## 






# # # # # # # # # # # ## # # # ENERGY EVALUATION# ## # # # # ## # # # # ## # # # # ##  # # # # ## 
# FUNCION DE ENERGIA (O COSTO) A SER MINIMIZADA
# INPUT: vector de 1XN de spines (N=numero de nodos de la red)
#spins <- vector(mode = "numeric", length = node_numbers)
# OUTPUT: valor de la energia de la red para un spin dado y matiz de acople I dado.
#goodness <- function(spins, acoples) {
#  value = 0
  #acoples <- I[upper.tri(I)]
#  rows <- length(spins)
#  id <- which(upper.tri(matrix(, rows, rows)) == TRUE, arr.ind=T)
#  n <- nrow(id)
#  for (i in seq_along(1:n)) {
#    value <- value + acoples[id[i,1], id[i,2]]*spins[id[i,1]]*spins[id[i,2]]
    #print(value)
#  }
#  return(value)
#}
goodness <- function(spins, magnetizaciones, acoples) {
  simples <- t(magnetizaciones)%*%(as.matrix(spins)) #25-12-17
  value = 0
  rows <- length(spins)
  colu <- rev(abs(sequence(seq.int(rows - 1)) - rows) + 1)
  fila <- rep.int(seq.int(rows - 1), rev(seq.int(rows - 1)))
  id <- cbind(fila, colu)
  n <- dim(id)[1] #15-12-17
  for (i in seq_along(1:n)) {
    value <- value + acoples[id[i,1], id[i,2]]*spins[id[i,1]]*spins[id[i,2]]
    #print(value)
  }
  value <-  value + simples #25-12-17
  return(value)
}
# # # # # # # # # # # ## # # # ENERGY EVALUATION# ## # # # # ## # # # # ## # # # # ##  # # # # ## 






# # # # # # # # # # # ## # # # LEARNING PARAMETER # ## # # # # ## # # # # ## # # # # ##  # # # # ## 
# FUNCION PARA GENERAR TPARAMETRO DE APRENDIZAJE EN LA MAQUINA DE BOLTZMANN
# INPUT: valor de alpha usualmente entre 0.8 y 0.9999
#        temperatura_antes
#        ciclo de iteracion actual
# OUTPUT: una prametro para el ciclo m
learning_par <- function(ciclo, decay, LearningPar) {
  LearningPar <- LearningPar * (1/(1 + decay*ciclo))
  return(LearningPar)
}
# ejemplo
#LearningRate <- learning_par(ciclo = 3, decay = 0.0008, LearningPar = 0.001)
# # # # # # # # # # # ## # # # LEARNING PARAMETER # ## # # # # ## # # # # ## # # # # ##  # # # # ## 





# # # # # # # # # # # ## # # # SPINES ALEATORIOS # ## # # # # ## # # # # ## # # # # ##  # # # # ## 
# FUNCION PARA GENERAR VECTOR DE SPINES ALEATORIOS                                               #
#Input: N que es el numero de nodos de la red                                                    #
#OUTPUT: vector de 1XN de spines aleatorios 0 o 1.                                               #
random_spin <- function(N) {
  espines <- floor(runif(N,0,2))
  return(espines)
}
# # # # # # # # # # # ## # # # SPINES ALEATORIOS # ## # # # # ## # # # # ## # # # # ##  # # # # ## 


# # # # # # # # # # # ## # # # FUNCION PARA GENERAR TEMPERATURA # # # # # # # # # # # ## # # # # #
# INPUT: valor de alpha udualmente entre 0.8 y 0.9999                                            #
#        temperatura_antes                                                                       #
#        ciclo de iteracion actual                                                               #
# OUTPUT: una temperatura para el ciclo k.                                                       #
tempfunc <- function(temperatura_antes, alpha) {
  temp <- temperatura_antes*alpha
  return(temp)
}
# # # # # # # # # # # ## # # # FUNCION PARA GENERAR TEMPERATURA # # # # # # # # # # # ## # # # # #


























# # # # # # # # # # # # # # # # # CONTRASTIVE DIVERGENCE # # # # # # # # # # # # # # # # # # # 
#OUTPUTS:  h matriz de Nx1 con los parametros de magnetizacion
#          J matriz NxN con los acoples
# INPUTS: sigma_i_data  medias de los spines de la data empirica
#         sigma_ij_data correlaciones entre i y j de la data empirica
#         nu parametro de aprendizaje 
#         h, J parametros de magnetizacion y acoples en la epoca t-1
contrastive_divergence <- function(sigma_i_data, sigma_ij_data, sigma_i_sample, sigma_ij_sample, nu, h, J) {
  require(gdata)
  h_old <- h
  J_old <- J
  #contrastive divergence
  h <- h_old + nu*(sigma_i_data - sigma_i_sample)
  J <- J_old + nu*(sigma_ij_data - sigma_ij_sample)
  #trackeo de convergencia de parametros h i J
  diff_hpar <- h - h_old
  diff_Jpar <- J - J_old
  diff_Jpar <- as.matrix(upperTriangle(diff_Jpar, diag = FALSE))
  #mean square of differences
  msd_hpar <- (t(diff_hpar)%*%diff_hpar)/length(diff_hpar) #  suma de las diferencias al cuadrado de hi
  msd_Jpar <- (t(diff_Jpar)%*%diff_Jpar)/length(diff_Jpar) # suma de las diferencias al cuadrado de los Jij
  cat(sprintf("hi square sum diffs %.7f -- Jij square sum diffs: %.7f \n", msd_hpar, msd_Jpar))
  
  
  #trackeo de convergencia de los sigma
  diff_h <- sigma_i_data - sigma_i_sample
  diff_J <- sigma_ij_data - sigma_ij_sample
  diff_J <- as.matrix(upperTriangle(diff_J, diag = FALSE))
  #mean square of differences
  msd_h <- (t(diff_h)%*%diff_h)/length(diff_h) #  suma de las diferencias al cuadrado de hi
  msd_J <- (t(diff_J)%*%diff_J)/length(diff_J) # suma de las diferencias al cuadrado de los Jij
  cat(sprintf("sigmas_i square sum diffs %.7f -- sigmas_ij square sum diffs: %.7f \n", msd_h, msd_J))
  
  
  return(list(h, J, msd_h, msd_J, msd_hpar, msd_Jpar))
}
#El mismo que antes, pero el error en las correlaciones entre el modelo y cada iteracion se calcula distinto.
contrastive_divergence2 <- function(sigma_i_data, sigma_ij_data, sigma_i_sample, sigma_ij_sample, nu, h, J) {
  #sigma_ij_data son los <sisj>+ de la fase posotiva, que en nuestra caso se calcula solo 1 vez a partir de la data de entranmiento
  #sigma_ij_sample son los <sisj>- de la fase negativa que viene de la funcion annealing.
  require(gdata)
  h_old <- h
  J_old <- J
  #contrastive divergence
  h <- h_old + nu*(sigma_i_data - sigma_i_sample)
  J <- J_old + nu*(sigma_ij_data - sigma_ij_sample)
  #J <- as-matrix(J) #puesto el 031217
  
  #trakeo de <sisj>+ - <sisj>-
  dif <- sigma_ij_data - sigma_ij_sample
  N <- ncol(J)
  N <- N*(N+1)/2 - N
  trk <- matrix(NA, nrow=1, ncol=N)  #trakeo de <sisj>+ - <sisj>-
  idx <- combn(1:ncol(J),2)
  for (i in seq_along(1:N)) {
    trk[ 1, i ] <- dif[ idx[1,i], idx[2,i] ]
  }
  trk <- as.vector(trk)
  
  #trackeo de convergencia de parametros h i J
  diff_hpar <- h - h_old
  diff_Jpar <- J - J_old
  diff_Jpar <- as.matrix(upperTriangle(diff_Jpar, diag = FALSE))
  #mean square of differences
  msd_hpar <- (t(diff_hpar)%*%diff_hpar)/length(diff_hpar) #  suma de las diferencias al cuadrado de hi
  msd_Jpar <- (t(diff_Jpar)%*%diff_Jpar)/length(diff_Jpar) # suma de las diferencias al cuadrado de los Jij
  #cat(sprintf("hi square sum diffs %.7f -- Jij square sum diffs: %.7f \n", msd_hpar, msd_Jpar))
  
  #trackeo de convergencia de los sigma como the absolute correlation difference metric (paper de Broderick: Faster solutions of the inverse pairwise ising problem)
  sigmas_i_diff <- abs(sigma_i_data - sigma_i_sample)
  sigmas_ij_diff <- abs(sigma_ij_data - sigma_ij_sample)
  sigmas_ij_diff <- as.matrix(upperTriangle(sigmas_ij_diff, diag = FALSE))
  #mean 
  mean_sigmas_i_diff <- (mean(sigmas_i_diff))    # log(mean(abs(sigma_i_data - sigma_i_sample))
  mean_sigmas_ij_diff <- (mean(sigmas_ij_diff))  # log(mean(abs(sigma_ij_data - sigma_ij_sample))
  #cat(sprintf("sigmas_i abs means diffs %.7f -- sigmas_ij abs means  diffs: %.7f \n", mean_sigmas_i_diff, mean_sigmas_ij_diff))
  
  return(list(h, J, mean_sigmas_i_diff, mean_sigmas_ij_diff, msd_hpar, msd_Jpar, trk))
}
#contrastive_divergence3 aplica el siguiente criterio para calcular el nuevo J:
#si <sisj>+ > <sisj>- entonces Jij = Jij + nu,
#si <sisj>+ < <sisj>- entonces Jij = Jij - nu,
contrastive_divergence3 <- function(sigma_i_data, sigma_ij_data, sigma_i_sample, sigma_ij_sample, nu, h, J) {
  #sigma_ij_data son los <sisj>+ de la fase posotiva, que en nuestra caso se calcula solo 1 vez a partir de la data de entranmiento
  #sigma_ij_sample son los <sisj>- de la fase negativa que viene de la funcion annealing.
  require(gdata)
  h_old <- h
  J_old <- J
  #contrastive divergence (segun paper Ackley 1985)
  tempJ <- ((sigma_ij_data > sigma_ij_sample)*1)*2-1
  J <- J_old + nu*tempJ
  tempH <- ((sigma_i_data > sigma_i_sample)*1)*2-1
  h <- h_old + nu*tempH
  #J <- as-matrix(J) #puesto el 031217
  
  #trakeo de <sisj>+ - <sisj>-
  dif <- sigma_ij_data - sigma_ij_sample
  N <- ncol(J)
  N <- N*(N+1)/2 - N
  trk <- matrix(NA, nrow=1, ncol=N)  #trakeo de <sisj>+ - <sisj>-
  idx <- combn(1:ncol(J),2)
  for (i in seq_along(1:N)) {
    trk[ 1, i ] <- dif[ idx[1,i], idx[2,i] ]
  }
  trk <- as.vector(trk)
  
  #trackeo de convergencia de parametros h i J
  diff_hpar <- h - h_old
  diff_Jpar <- J - J_old
  diff_Jpar <- as.matrix(upperTriangle(diff_Jpar, diag = FALSE))
  #mean square of differences
  msd_hpar <- (t(diff_hpar)%*%diff_hpar)/length(diff_hpar) #  suma de las diferencias al cuadrado de hi
  msd_Jpar <- (t(diff_Jpar)%*%diff_Jpar)/length(diff_Jpar) # suma de las diferencias al cuadrado de los Jij
  #cat(sprintf("hi square sum diffs %.7f -- Jij square sum diffs: %.7f \n", msd_hpar, msd_Jpar))
  
  #trackeo de convergencia de los sigma como the absolute correlation difference metric (paper de Broderick: Faster solutions of the inverse pairwise ising problem)
  sigmas_i_diff <- abs(sigma_i_data - sigma_i_sample)
  sigmas_ij_diff <- abs(sigma_ij_data - sigma_ij_sample)
  sigmas_ij_diff <- as.matrix(upperTriangle(sigmas_ij_diff, diag = FALSE))
  #mean 
  mean_sigmas_i_diff <- log(mean(sigmas_i_diff))    # log(mean(abs(sigma_i_data - sigma_i_sample))
  mean_sigmas_ij_diff <- log(mean(sigmas_ij_diff))  # log(mean(abs(sigma_ij_data - sigma_ij_sample))
  #cat(sprintf("sigmas_i abs means diffs %.7f -- sigmas_ij abs means  diffs: %.7f \n", mean_sigmas_i_diff, mean_sigmas_ij_diff))
  
  return(list(h, J, mean_sigmas_i_diff, mean_sigmas_ij_diff, msd_hpar, msd_Jpar, trk))
}
# # # # # # # # # # # # # # # # # CONTRASTIVE DIVERGENCE # # # # # # # # # # # # # # # # # # # 



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# FUNCION PARA CALCULAR LOS <sigma_{i}*sigma{j}>                                                        #
# CREAtion date: 10-ago-17                                                                              #
#Inputs:                                                                                                #
# M = matriz de spines                                                                                  #
# Output: matriz NXN simetrica con los productos #Calculo de <sigma_{i}*sigma{j}> con valores NA        #
# en la diagonal                                                                                        #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
pairproducts <- function(M) {
  require(Matrix)
  N <- ncol(M)
  P <- matrix(NA, ncol=N, nrow=N)
  ids <- combn(c(1:N),2) #15-12-17
  #rep <- ncol(ids)
  rep <- dim(ids)[2]
  for (i in seq_along(1:rep)) {
    fila <- ids[1,i]
    colu <- ids[2,i]
    P[fila, colu] <- mean(M[,fila]*M[,colu])
  }
  P <- forceSymmetric(P,uplo="U")
  return(P)
}               
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# FUNCION PARA determinar la frecuencia de los estados de una muestra                   #
# Fecha creacion: 03-dic-17                                                             #
# Fecha modificacion:                                                                   #
# Input: matriz muestra (posiblemente datos simulados de metropolis hasting o datos de compra)                                           #                     
# Output: matriz de todos los estados encontrados en la muestra                         #
# Outputs:  vector con la fecuencia de cada uno de los estados encontrados              #
# Nota: su hay NAs en las salidas es porque No todos los estados estaban presentes      #
# en la muestra                                                                         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
frecuency <- function(muestra) {
  if (!is.matrix(muestra)) {
    muestra <- as.matrix(muestra)
  }
  #N = nrow(muestra) #muestra es una matriz con lso estados, simulados de un metropolis hasting
  N = dim(muestra)[1] #25-12-17
  #n = ncol(muestra)
  n = dim(muestra)[2] #15-12-17
  n_states = 2^n
  frqq <- matrix(NA, nrow=n_states, ncol=n+1) #colocamos todos los posibles estados y la frecuencia
  
  it = 1
  matriz_states <- muestra
  #while ( !(is.null(nrow(matriz_states) ) ) | (nrow(matriz_states) != 0)  ) { 
  while ( nrow(matriz_states) != 0) { 
    selected_state <- matriz_states[1,]
    row.is.a.match <- apply(matriz_states, 1, identical, selected_state)
    total.matches <- sum(row.is.a.match)
    frqq[it, ] <- c(selected_state, total.matches) #gradamos estado y su frecuencia
    match.idx <- which(row.is.a.match)
    #ahora borramos todos los selected_state de matriz_state
    matriz_states <- matriz_states[-match.idx, ]
    if ( !(is.matrix(matriz_states)) ) {
      matriz_states <- t(as.matrix(matriz_states))
    }
    #print(nrow(matriz_states))
    it <- it + 1
  }
  frqq  <- frqq[order(-frqq[,n+1]), ]
  states <- frqq[,c(1:n)] #matriz con los estados contados
  counts <- frqq[, n+1] #matriz con el conteo de cada estado
  
  return(list(states, counts))
}
#ejemplo:
#ff <- frecuency(muestra = mh[[1]]) #muestra generada del metropolis hasting
#ff[[1]] #matriz con los estados contados
#ff[[2]] #vector con el conteo de cada estado
# frecuencyV2 corre marginalmente un poco mas rapido que frecuency.
frecuencyV2 <- function(muestra) {
  if (!is.matrix(muestra)) {
    muestra <- as.matrix(muestra)
  }
  # chequeo si muestra viene con columna newid
  if ("newid" %in% colnames(muestra))
  {
    borrar <- which( colnames(muestra)== "newid" )
    muestra <- muestra[, -c(borrar)]
  }
  
  #N = nrow(muestra) #muestra es una matriz con lso estados, simulados de un metropolis hasting
  N = dim(muestra)[1] #25-12-17 numero de filas
  #n = ncol(muestra)
  n = dim(muestra)[2] #15-12-17 numero de variables
  
  # identificamos los estados de muestra
  unique_wb <- as.matrix(unique(muestra[,c(1:ncol(muestra))])) 
  #n_states = 2^n
  n_states <- nrow(unique_wb)
  frqq <- matrix(NA, nrow=n_states, ncol=n+1) #colocamos todos los posibles estados y la frecuencia
  
  # comenzamos:
  matriz_states <- muestra
  it = 1
  while ( nrow(matriz_states) != 0) { 
    selected_state <- matriz_states[1,]
    row.is.a.match <- apply(matriz_states, 1, identical, selected_state)
    total.matches <- sum(row.is.a.match)
    frqq[it, ] <- c(selected_state, total.matches) #gradamos estado y su frecuencia
    match.idx <- which(row.is.a.match)
    #ahora borramos todos los selected_state de matriz_state
    matriz_states <- matriz_states[-match.idx, ]
    if ( !(is.matrix(matriz_states)) ) {
      matriz_states <- t(as.matrix(matriz_states))
    }
    #print(nrow(matriz_states))
    it <- it + 1
  }
  frqq  <- frqq[order(-frqq[,n+1]), ]
  states <- frqq[,c(1:n)] #matriz con los estados contados
  counts <- frqq[, n+1] #matriz con el conteo de cada estado
  
  return(list(states, counts))
}
#Ejemplo
#ff <- frecuencyV2(muestra = mh[[1]]) 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# FUNCION PARA KL                                                                             #
# KL = sum P * sum log (P/Q)    P es al pdf original, Q es la aproximacion                    #
# Teniendo la pdf o frecuencia de base , ahora debemos calcular la frecuencia de los          #
# estados de la muestra (que queremso comparar, posiblemente la que                           #
# viene de BM) en el orden de los estados contados anteriormente)                             #
#                                                                                             #
# Fecha creacion: 03-dic-17                                                                   #
# Fecha modificacion:                                                                         #
# Inputs:                                                                                     #
#       sample es la muestra a la cual se sacamos la pdf y luego comparamos con base_states   #
#       base_states son los estados de la muestra original  ff[[1]]                           #
#       base_counts es la frecuencia de los estados de la muestra original ff[[2]]            #
# Output: Kullback--Leibler divergence y matriz con estados/conteo muestra original/          #
# conteo muestra a ser comparada con la original                                              #
# NOtas:                                                                                      #
# Dado que empiricamente vemos que la distribucion P y Q suelen tener conteo igual a ceros,   #
# no e sposible calcular KL porque esta asume que todo P y Q son distintos de cero. Entonces  #
# tenemos dos opciones: 1) eliminar todos los estados de P y Q en los cuales el conteo es cero#
## 2) calcular Jensen???Shannon divergence. Nosotros optamos por esta ultima.                 #
# Ver: https://en.wikipedia.org/wiki/Jensen%E2%80%93Shannon_divergence                        #
# https://stats.stackexchange.com/questions/97938/calculate-the-kullback-leibler-divergence-in-practice
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
KL <- function(base_states, base_counts, sample) {
  require(entropy)
  base_states <- na.omit(base_states)
  base_cpunts <- na.omit(base_counts)
  #nodos <- ncol(base_states)
  nodos <- dim(base_states)[2] #25-12-17
  #filas <- nrow(base_states)
  filas <- dim(base_states)[1] #25-12-17
  # matriz con: estados de muestra original, frequency de pdf original, frequency de la pdf de sample
  M <- matrix(NA, ncol=nodos + 2, nrow=filas)
  #tenemos que obtener la frecuencia de sample
  ff <- frecuency(muestra = sample)
  conteo_estados_sample <- na.omit(ff[[2]])
  estados_encontrados_sample <- na.omit(ff[[1]])
  
  for (i in seq_along(1:filas)) {
    vec <- base_states[i,]
    M[i,c(1:nodos)] <- vec
    M[i,nodos+1] <- base_counts[i]
    id <- which(apply(estados_encontrados_sample, 1, function(x) identical(x[1:nodos], vec )))
    if (length(id) == 0) {
      M[i,nodos+2] <- 0
    } else {
      M[i,nodos+2] <- conteo_estados_sample[id]
    }
  }
  #calculo de Jensen Shannon divergence:
  p <- M[,nodos+1]/sum(M[,nodos+1])
  q <- M[,nodos+2]/sum(M[,nodos+2])
  m <- (p + q)
  
  kl_pm <- KL.empirical(p, m)
  kl_qm <- KL.empirical(q, m)
  JSD <- 0.5*kl_pm + 0.5*kl_qm
  return(list(JSD, M))
} 
#Ejemplo
#mh <-  metropolis_hasting(10000, couplings = I, perf = 'off')
#ff <- frecuency(muestra = mh[[1]])
#J <- as.matrix(J) #para que sea mas rapido
#Q <- metropolis_hasting(10000, couplings = I, perf = 'off')
#rr <- KL(base_states = ff[[1]], base_counts = ff[[2]], sample = Q[[1]])
#rr[[1]] #Jensen-Shannon divergence (mientras mas chico mejor: mide how many bit is expected to lose approximating P with Q)
#rr[[2]] #matriz con los estados el el conteo de l amuestra original y la comparada.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 





# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# FUNCION PARA determinar KUllback-Leibler tradicional.                                 #
# #https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence                    #
# Fecha creacion: 12-dic-17                                                             #
# Inputs:                                                                               #
# P distribucion real empirica                                                          #
# Q distribucion aproximada (ej de un metropolis hasting)                               #
# Si P es NULL entonces imputs adicionales: P_states y P_frequency                      #
# Output: KUllback-Leibler tradicional index                                            #
#  kullback2 es lo mismo pero se asume que P=NULL                                       #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
kullback <- function(P=NULL, states=NULL, frequencies=NULL, Q) {
  if ( !is.null(P) ) {
    P_info <- frecuency(muestra = P )
    P_states <- na.omit(P_info[[1]])
    P_frequency <- na.omit(P_info[[2]])
  } else {
    P_states <- states
    P_frequency <- frequencies
  }
  Q_info <- frecuency(muestra = Q )
  Q_states <- na.omit(Q_info[[1]])
  Q_frequency <- na.omit(Q_info[[2]])

  # Politica:Todos los estados en Q_states que NO estan en P_states hay que borrarlos.
  #filas <- nrow(P_states )
  filas <- dim(P_states )[1] #25-12-17
  #nodos <- ncol(P_states )
  nodos <- dim(P_states )[2] #25-12-17
  # matriz con: estados de muestra original, frequency de pdf original (P_frequency), frequency de la pdf de sample (Q_frequency)
  M <- matrix(NA, ncol=nodos + 2, nrow=filas)
  for (i in seq_along(1:filas)) {
    vec <- P_states[i,]
    M[i,c(1:nodos)] <- vec
    M[i,nodos+1] <- P_frequency[i]
    id <- which(apply(Q_states, 1, function(x) identical(x[1:nodos], vec )))
    if (length(id) == 0) {
      M[i,nodos+2] <- 0
    } else {
      M[i,nodos+2] <- Q_frequency[id]
    }
  }
  require(entropy)
  # The Kullback-Leibler divergence from Q to P is often denoted DKL(P||Q).
  KL <- KL.empirical(M[,nodos+2], M[,nodos+1])
  
  return(KL)
}
kullback2 <- function(states=NULL, frequencies=NULL, Q) {
 
  P_states <- states
  P_frequency <- frequencies

  Q_info <- frecuency(muestra = Q )
  Q_states <- na.omit(Q_info[[1]])
  Q_frequency <- na.omit(Q_info[[2]])
  
  # Politica:Todos los estados en Q_states que NO estan en P_states hay que borrarlos.
  #filas <- nrow(P_states )
  filas <- dim(P_states )[1] #25-12-17
  #nodos <- ncol(P_states )
  nodos <- dim(P_states )[2] #25-12-17
  # matriz con: estados de muestra original, frequency de pdf original (P_frequency), frequency de la pdf de sample (Q_frequency)
  M <- matrix(NA, ncol=nodos + 2, nrow=filas)
  for (i in seq_along(1:filas)) {
    vec <- P_states[i,]
    M[i,c(1:nodos)] <- vec
    M[i,nodos+1] <- P_frequency[i]
    id <- which(apply(Q_states, 1, function(x) identical(x[1:nodos], vec )))
    if (length(id) == 0) {
      M[i,nodos+2] <- 0
    } else {
      M[i,nodos+2] <- Q_frequency[id]
    }
  }
  require(entropy)
  # The Kullback???Leibler divergence from Q to P is often denoted DKL(P||Q).
  KL <- KL.empirical(M[,nodos+2], M[,nodos+1])
  
  return(KL)
}
#Ejemplo: kl <- kullback(P=muestra, Q=mh[[1]])
# ejemplo: kullback(P=NULL, states=ff[[1]], frequencies=ff[[2]], Q=mh[[1]] )
# 22-dc-17
# Vemos que el calulo de KL cuando el numero de nodos es de 20, es extremamemnte largo 
# y por tanto impractico. Por lo tanto creamos kullback3. kullback3, calcula kl pero 
# para un porcentaje acumulado preestablecido (ej. 0.8) de la distribucion de P.
# Inputs : * P (la muestra real) -- opcional
#          * states de ff[1]] que son los estados ya cortados al umbral acumulado preestablecido de P
#          * frequencies  de ff[[2]] que es la frecuencia  de los estados de states de P.
#          * Q muestra en proceso.
kullback3 <- function(P=NULL, states=NULL, frequencies=NULL, Q, umbral=0.8) {
  if ( !is.null(P) ) {
    P_info <- frecuency(muestra = P )
    P_states <- na.omit(P_info[[1]])
    P_frequency <- na.omit(P_info[[2]])
  } else {
    P_states <- na.omit(states)
    P_frequency <- na.omit(frequencies)
  }
  Q_info <- frecuency(muestra = Q )
  Q_states <- na.omit(Q_info[[1]])
  Q_frequency <- na.omit(Q_info[[2]])
  
  #Cortamos Q_states y Q_frequency al umbral preestablecido
  pacum <- cumsum(Q_frequency/nrow(Q))
  corte <- which.min(abs(pacum-umbral))
  corte <- corte + 1
  Q_states <-  Q_states[-c(corte:nrow(Q_states)), ]
  Q_frequency <- Q_frequency[-c(corte:length(Q_frequency))]
  
  # Politica:Todos los estados en Q_states que NO estan en P_states hay que borrarlos.
  #filas <- nrow(P_states )
  filas <- dim(P_states )[1] #25-12-17
  #nodos <- ncol(P_states )
  nodos <- dim(P_states )[2] #25-12-17
  # matriz con: estados de muestra original, frequency de pdf original (P_frequency), frequency de la pdf de sample (Q_frequency)
  M <- matrix(NA, ncol=nodos + 2, nrow=filas)
  for (i in seq_along(1:filas)) {
    vec <- P_states[i,]
    M[i,c(1:nodos)] <- vec
    M[i,nodos+1] <- P_frequency[i]
    id <- which(apply(Q_states, 1, function(x) identical(x[1:nodos], vec )))
    if (length(id) == 0) {
      M[i,nodos+2] <- 0
    } else {
      M[i,nodos+2] <- Q_frequency[id]
    }
  }
  require(entropy)
  # The Kullback-Leibler divergence from Q to P is often denoted DKL(P||Q). P es la original, Q es la aproximada
  KL <- KL.empirical(M[,nodos+2], M[,nodos+1])
  
  return(KL)
} 
#ejemplo:
#kl <- kullback3(states=states, frequencies=fr, Q=Q, umbral=1 )
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 



# # # # # # # # # # # ## # # # BOLTZMANN-LEARNING # ## # # # # ## # # # # ## # # # # ##  # # # # ## 
# Creation date: 19-may-18                                                                        #
# Creamos la funcion bmachine (utilizada en BM_toy_asonam.R) para llevar a cabo                   #
# el Boltzmann-learning datos transaccionales.                                                    #
# Bolstzamnn machine sacado de boltzmann_machine_realdata_V1.R                                    #
# INPUTS:                                                                                         #
# condInitial: puede ser -a2- para la funcion initial_condition                                   #
# muestra: base de datos con columnas (productos) y filas (compras) en formato 0 y 1 matrix       #
# Kmax: Numero maximo de iteraciones de la maquina de Boltzmann                                   #
# parametros aprendizaje                                                                          #
# LearningRate: parametro de aprendizaje 0.9 el original funcionando                              #
# decay: decaimiento del parametro de aprendizaje (en 0.02)                                       #
# N_sa: Numero de muestras para la fase negativa (numero de veces que repetimoes el annealing)    #
# umbral: umbral para el calculo del KL                                                           #
# every: cada cuantas iteraciones medir el KL.                                                    #
# metrop_its: numero iteraciones de metropolis hasting                                            #
# OUTPUTS:                                                                                        #
# just_for_track: trackeo de el ajuste <sisj>+ - <sisj>-)                                         #
# diver : divergencia KL                                                                          #
# cd: diferencias <sisj>+ - <sisj>-                                                               #
# inferences[Kmak, ] : parametros inferidos y los sigmas_{ij} y Hi                                #
# # # # # # # # # # # ## # # # BOLTZMANN-LEARNING # ## # # # # ## # # # # ## # # # # ##  # # # # ## 
bmachine <- function(condInitial, muestra, Kmax, metrop_its, LearningRate, decay, N_sa, umbral, every) {
  library(progress)
  #require(ggplot2)
  #library(latex2exp)
  #library(ggthemes)
  # ADQUISICION DE ALGUNOS PARAMETROS PREVIOS
  ff <- frecuency(muestra = muestra ) #calcula frecuencias de los estados de la muestra 
  sigma_i_sample <-  as.matrix(colMeans(muestra)) #15-abr-18
  sigmas_sample <- as.matrix(pairproducts(muestra))
  ic <- initial_conditions(df=NULL, N=ncol(muestra), tipe=condInitial, si=sigma_i_sample, sij=sigmas_sample, desv=0.5) 
  H0 <- ic[[1]] 
  J0 <- ic[[2]]
  numero_mediciones <- Kmax/every
  # generacion de vector senal para medir kl
  if (every == 1) {
    senal <- rep(TRUE, Kmax)
  } else {
    inic <- c(TRUE, rep(FALSE, every-2), TRUE)
    v <- c(rep(FALSE, every-1), TRUE)
    vv <- rep(v, numero_mediciones-1)
    senal <- c(inic, vv)
    rm(inic, v, vv)
  }
  
  # PPRESETEO DE VARIABLES Y MATRICES
  #comenzamos:
  pb <- progress_bar$new(format = "  Progress [:bar] :percent eta: :eta", total = Kmax, clear = FALSE, width= 60) 
  # Annealing network
  initial_spins=rep(1,ncol(sigmas_sample))
  J <- as.matrix(J0)
  H <- as.matrix(H0)
  epoch <- 1
  just_for_track <- matrix(NA, ncol=3, nrow=Kmax) #col1:epoch, col2:log(mean(abs(<sisj>+ - <sisj>-)), col3:suma de las diferencias al cuadrado de los Jij
  numbers_of_nodes = ncol(muestra)
  numbers_of_coupl = (numbers_of_nodes*(numbers_of_nodes+1)/2) - numbers_of_nodes
  #matriz inferences: epoch de la BA / Jij / hi. El orden es J12, J13, J23, J14, J24, J34, J15, J25, J35, etc..
  inferences <- matrix(NA, nrow=Kmax, ncol=numbers_of_nodes+numbers_of_coupl+1) #guardando los parametros h y J de cada epoch de la BM.
  #matriz de las diferencias <sisj>+ - <sisj>- para cada epoch
  cd <-  matrix(NA, nrow=Kmax, ncol=numbers_of_coupl+1)
  #matriz de los KUllback-Leibler Divergence
  diver <- matrix(NA, nrow=Kmax, ncol=2) # numero de epoch / KL
  
  # MAQUINA
  while (epoch <= Kmax) {
    # progress bar
    pb$tick()
    Sys.sleep(1 / Kmax)
    # progress bar fin 
    # Iniciamos la fase negativa: network is run free by using metropolis hasting
    sampling <- metropolis.sampling(magnetizaciones=H, acoples = J, initial_spin = initial_spins, kmax=metrop_its, ks=0, To=NULL, alpha=NULL)
    # en caso de que quisiera operar la maquina sin las magnetizaciones.
    #sampling <- metropolis.sampling(magnetizaciones=as.matrix(rep(0, ncol(muestra))), acoples = J, initial_spin = initial_spins, kmax=20000, ks=0, To=NULL, alpha=NULL)
    sigmas <- as.matrix(pairproducts( sampling[[1]] ))
    rescue <- contrastive_divergence2(sigma_i_data = sigma_i_sample, 
                                      sigma_ij_data = sigmas_sample, 
                                      sigma_i_sample = as.matrix( colMeans( sampling [[1]]) ), 
                                      sigma_ij_sample = sigmas,
                                      nu=learning_par(ciclo = epoch, decay = decay, LearningPar = LearningRate),
                                      h=H,
                                      J=J)
    J <- as.matrix(rescue[[2]])
    H <- rescue[[1]]
    put_all <- c(epoch, rescue[[4]],  rescue[[6]]) #epoch / log(mean(abs(<sisj>+ - <sisj>-)) / msd_Jpar
    just_for_track[epoch,] <- put_all
    #Colocamos los acoples Jij y hi en matriz inferences
    inferences[epoch, ] <- c(epoch, upperTriangle(J, diag = FALSE), rescue[[1]])
    #colocamos las diferencias <sisj>+ - <sisj>- de cada epoch
    cd[epoch, ] <- c(epoch, rescue[[7]])
    #calculamos Jensen-Shannon divergence y la colocamos en matrix diver
    # vamos a ir calculando kl solo cuando sea necesario y no en cada iteracion. Dato de entrada -every- cada cuanto medimos
    if (senal[epoch] == TRUE) {
      kl <- kullback3(states=ff[[1]], frequencies=ff[[2]], Q=sampling[[1]], umbral=umbral ) #22-dic-17
      diver[epoch, ] <- c(epoch, log2(kl) )
    }
    epoch <- epoch + 1
  }
  # OUTPUTS
  return(list(just_for_track, diver, cd, H, J, inferences, sampling[[1]] ))
}
# EJEMPLO:
#system.time(
#  machine_results <- bmachine(condInitial = "a5", muestra = wbasket, Kmax = 100, LearningRate = 0.8, decay = 0.02, N_sa = 10, umbral = .7, every = 10)
#)
# Recuperacion de la informacion
#just_for_track <- machine_results[[1]]
#diver <- machine_results[[2]]
#cd <- machine_results[[3]]
#H <- machine_results[[4]]
#J <- machine_results[[5]]
#inferences <- machine_results[[6]]
#last_mh_sampling <- machine_results[[7]] # equivale a sampling[[1]] del ultimo metropolis hasting efectuado.
# # # # # # # # # # # ## # # # BOLTZMANN-LEARNING # ## # # # # ## # # # # ## # # # # ##  # # # # ## 





# # # # # # # # # # # ## # # # PLOT OF Cij# ## # # # # ## # # # # ## # # # # ##  # # # # ## 
# Modification date: 20-12-17
# Calculo de los Cij = <SiSj> - <Si><Sj> y plot 
# INPUT = muestra1 USUALMENTE la muestra real
#         muestra2, USUALMENTE la muestra simulada con los acoples J estimados con BM
# OUTPUT: * matriz con: fila, columna, Cij de muestar1, Cij de muestra2
#         * plot de Cij real en eje X, y de Cij estimado en el eje Y.
# # # # # # # # # # # ## # # # PLOT OF Cij# ## # # # # ## # # # # ## # # # # ##  # # # # ## 
perf <- function(muestra1, muestra2) {
  require(gdata)
  #para muestra1
  nodos <- ncol(muestra1)
  C <- matrix(NA, ncol=nodos, nrow=nodos)
  R <- matrix(NA, nrow=nodos*(nodos-1)/2, ncol=4)
  Si <- colMeans(muestra1)
  SiSj <- pairproducts(muestra1)
  id <- combn(1:nodos,2)  #fila 1 es fila de la matriz, fila 2 es la columna de la matriz
  n <- ncol(id)
  for (i in seq_along(1:n)) {
    C[ id[1,i], id[2,i] ] <- SiSj[ id[1,i], id[2,i] ] - Si[ id[1,i] ]*Si[ id[2,i] ]
  }
  R[, c(1,2)] <- t(id)
  R[, 3] <- upperTriangle(C, diag = FALSE, byrow=TRUE)
  #para muestra2
  nodos <- ncol(muestra2)
  Si <- colMeans(muestra2)
  SiSj <- pairproducts(muestra2)
  id <- combn(1:nodos,2)  #fila 1 es fila de la matriz, fila 2 es la columna de la matriz
  n <- ncol(id)
  for (i in seq_along(1:n)) {
    C[ id[1,i], id[2,i] ] <- SiSj[ id[1,i], id[2,i] ] - Si[ id[1,i] ]*Si[ id[2,i] ]
  }
  R[, 4] <- upperTriangle(C, diag = FALSE, byrow=TRUE)
  #Generacion de la grafica de Cij real en eje X contra Cij estimado en el eje Y
  library(ggplot2)
  library(latex2exp)
  d <- as.data.frame(R) #V3 es el Cij original, V4 es el Cij estimado
  #https://drsimonj.svbtle.com/pretty-scatter-plots-with-ggplot2
  d$pc <- predict(prcomp(~V3+V4, d))[,1]
  pl <- ggplot(d, aes(V3, V4, color = pc)) +
    geom_point(shape = 16, size = 2, show.legend = FALSE, alpha = .999) +
    theme_minimal() + 
    geom_abline(slope=1, intercept=0, alpha=0.7) + 
    scale_color_gradient(low = "#f0650e", high = "#0091ff") +
    labs(x=TeX(' $C_{ij}$  '), y=TeX( '$ \\hat{C_{ij} } $' ) ) +
    ggtitle('Connected two-body averages Cij' )
  return(list(R, pl))
}
#Ejemplo
#rr <- perf(muestra, sampling[[1]])
#rr[[1]] #da la matriz con los datos
#rr[[2]] #da el plot
# # # # # # # # # # # ## # # # PLOT OF Cij# ## # # # # ## # # # # ## # # # # ##  # # # # ## 

# # # # # # # # # # # ## # # # PLOT OF Cij contra los Jij# ## # # # # ## # # # # ## # # # # 
# Modification date: 27-12-17
# Plot en ggplot2 de Cij (correlaciones en el eje X) contra Jij (los acoples en el eje Y)
# INPUT = Matriz Cij
#         Matrix Jij
# OUTPUT: * matriz con: fila, columna, Cij , Jij
#         * plot de Cij real en eje X, y de Jij  en el eje Y.
# # # # # # # # # # # ## # # # PLOT OF Cij# ## # # # # ## # # # # ## # # # # ##  # # # # ## 
plotCJ <- function(correlations, couplings) {
  S <- correlations
  J <- couplings
  require(gdata)
  nodos <- ncol(J)
  R <- matrix(NA, nrow=nodos*(nodos-1)/2, ncol=4)
  id <- combn(1:nodos,2)  #fila 1 es fila de la matriz, fila 2 es la columna de la matriz
  n <- ncol(id)
  
  R[, c(1,2)] <- t(id)
  R[, 3] <- upperTriangle(S, diag = FALSE, byrow=TRUE)
  R[, 4] <- upperTriangle(J, diag = FALSE, byrow=TRUE)
  
  #Generacion de la grafica de Cij real en eje X contra Cij estimado en el eje Y
  library(ggplot2)
  library(latex2exp)
  d <- as.data.frame(R) #V3 es el Cij , V4 es el Jij estimado
  #https://drsimonj.svbtle.com/pretty-scatter-plots-with-ggplot2
  d$pc <- predict(prcomp(~V3+V4, d))[,1]
  pl <- ggplot(d, aes(V3, V4, color = pc)) +
    geom_point(shape = 16, size = 2, show.legend = FALSE, alpha = .999) +
    theme_minimal() + 
    geom_abline(slope=1, intercept=0, alpha=0.7) + 
    scale_color_gradient(low = "#f0650e", high = "#0091ff") +
    labs(x=TeX(' $C_{ij}$  '), y=TeX( '$ \\hat{J_{ij} } $' ) ) +
    ggtitle('Connected two-body averages Cij vs Estimated couplings' )
  return(list(R, pl))
}
#ejemplo:
#hh <- plotCJ(correlations=as.matrix(pairproducts(muestra)), couplings = J)
#hh[[1]] #data frame
#hh[[2]]
# # # # # # # # # # # ## # # # PLOT OF Cij# ## # # # # ## # # # # ## # # # # ##  # # # # ## 


# # # # # # # # # # # ## # # # COMPUTATION OF Cij # ## # # # # ## # # # # ## # # # # 
# CREADO el 18-abr-18 
# Calculo de los Cij = <SiSj> - <Si><Sj> y plot 
# INPUT = muestra1 USUALMENTE la muestra real
#         muestra2, USUALMENTE la muestra simulada con los acoples J estimados con BM
# OUTPUT: * matriz con: fila, columna, Cij de muestar1, Cij de muestra2
# # # # # # # # # # # ## # # # PLOT OF Cij# ## # # # # ## # # # # ## # # # # ##  # # # # ## 
pairwise_cij <- function(muestra1, muestra2) {
  require(gdata)
  #para muestra1
  nodos <- ncol(muestra1)
  C <- matrix(NA, ncol=nodos, nrow=nodos)
  R <- matrix(NA, nrow=nodos*(nodos-1)/2, ncol=4)
  Si <- colMeans(muestra1)
  SiSj <- pairproducts(muestra1)
  id <- combn(1:nodos,2)  #fila 1 es fila de la matriz, fila 2 es la columna de la matriz
  n <- ncol(id)
  for (i in seq_along(1:n)) {
    C[ id[1,i], id[2,i] ] <- SiSj[ id[1,i], id[2,i] ] - Si[ id[1,i] ]*Si[ id[2,i] ]
  }
  R[, c(1,2)] <- t(id)
  R[, 3] <- upperTriangle(C, diag = FALSE, byrow=TRUE)
  #para muestra2
  nodos <- ncol(muestra2)
  Si <- colMeans(muestra2)
  SiSj <- pairproducts(muestra2)
  id <- combn(1:nodos,2)  #fila 1 es fila de la matriz, fila 2 es la columna de la matriz
  n <- ncol(id)
  for (i in seq_along(1:n)) {
    C[ id[1,i], id[2,i] ] <- SiSj[ id[1,i], id[2,i] ] - Si[ id[1,i] ]*Si[ id[2,i] ]
  }
  R[, 4] <- upperTriangle(C, diag = FALSE, byrow=TRUE)
  d <- as.data.frame(R) #V3 es el Cij original, V4 es el Cij estimado
  return(d)
}
# ejemplo:
# d <- pairwise_cij(wbasket, last_mh_sampling)
# # # # # # # # # # # ## # # # COMPUTATION OF Cij # ## # # # # ## # # # # ## # # # # 




# # # # # # # # # # # # # # # # EDGE-LIST EXTRACTION # # # # # # # # # # # # # # # # 
#funcion para crear el edge list a partir de la data transactional.
# INPUT: m: matrix of N products and M rows. 
#output: edge vertex1 -  vertex 2
transaction_to_edges <- function(m) {
  #check if the matrix without colnames
  if (is.null(colnames(m))) {
    colnames(m) <- seq(1:ncol(m))
    rownames(m)  <- seq(1:nrow(m))
  } 
  #start
  edges <- NULL
  i = 1
  descartes <- as.matrix(rowSums(m))
  id <- which(descartes[,1] == 1 | descartes[,1] == 0 )
  if ( length(id) != 0) {
    m <- m[-c(id), ] #aqui eliminamos las canastas con solo 1 producto.
  }
  filas = nrow(m)
  while (i <= filas) {
    #fla <- as.matrix(m[i,])
    fla <- m[i, , drop=FALSE]
    # Me quedo con solo los productos en la canasta:
    fla <- fla[1, colSums(fla) > 0, drop=FALSE]
    
    # identifico las columnas con cero:
    #id <- as.numeric(which(fla[,1] == 0))
    #if ( length(id) != 0) {
    #  fla <-  as.matrix(fla[-c(id),])
    #}
    productos <- colnames(fla)
    ccb <- (combn(productos,2)) 
    # formando los miniedges
    miniedges <- matrix(NaN, nrow=ncol(ccb), ncol=2) #guardamos  vertex1 y vertex2
    for (j in seq_along(1:ncol(ccb))) {
      miniedges[j,1] <- as.character(ccb[1,j])
      miniedges[j,2] <- as.character(ccb[2,j])
    }
    edges <- rbind(edges,miniedges)   
    i = i + 1 
  }
  return(edges)
}
# ejemplo:
# edges <- transaction_to_edges(last_mh_sampling)
# # # # # # # # # # # # # # # # EDGE-LIST EXTRACTION # # # # # # # # # # # # # # # # 




# # # # # # # # # # # # # # # # RMSE EXTRACTION # # # # # # # # # # # # # # # # 
# Funcion para extraer el RMSE de los plots a continuacion.
# INPUTS: el modelo lineal: model (ejm model <- lm (...))
# Output: el rmse
#https://stackoverflow.com/questions/43123462/how-to-obtain-rmse-out-of-lm-result
rmse <- function(model) {
  # Mean squared error:
  RSS <- c(crossprod(model$residuals))
  MSE <- RSS / length(model$residuals)
  #Root MSE:
  rMSE <- sqrt(MSE)
  return(rMSE)
}
# # # # # # # # # # # # # # # # RMSE EXTRACTION # # # # # # # # # # # # # # # # 







# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
# ejemplo
#multiplot(p1, p2, p3, p4, cols=2)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


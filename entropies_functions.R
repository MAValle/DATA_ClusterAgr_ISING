
# funciones necesarias para computar entropias
# Estas funciones fueron desarrolladas en one_and_two_body_analysisV2.R

# Nombre del acutual archivo: entropies_functions.R
# ubicacion: dropbox->Research->PAPER MST-Fondecyt_2->data
# date: 12-08-18

# Modificaciones: 

#Notas:
# 12-08-18: creation
# 03-sep-18: creamos funcion get_energy que evalua la energia de un vector de estados.
# 08-sep-18: agregamos funcion custom_rbind que permite agregar a una dataframe X1, las
#           filas de otra dataframe X2, en los numero de columnas correspondientes a X1.
# 18-dic-18: agregamos la funcion get_couplings_from_state nos reporta los acoples 
#           involucrados en un estado.






# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 26-jun-18
# Comentario:_
# La funcion entropy_S1 esta mal calculada debido a que toma la sumatoria de las entropias
# de cada variable por separado. Lo correcto es tomar la probabilidad de ocurrencia de un 
# estado cualquiera como la multiplicacion de las probabilidades de ocurrencia o no ocurrencia
# de cada variable por separado. Esto es lo que hacemos en entropy_S1v2
# Funcion para calcular la entropia S1 
# Inputs:
#       wbasket: matriz transaction 
# Output:
#       S1: entropia
entropy_S1v2 <- function(wbasket) {
  #step 0: preparar el input: wbasket
  borrar <- as.numeric(which(colMeans(wbasket) == 0)) #hay que quitar las columnas que tienen puros ceros.
  wbasket <- wbasket[, -c(1, borrar)] # hay que quitar la columna con el id
  rm(borrar)
  #step 1: calcular la P(si=1) para todo i., luego calcular P(si=0).
  probs1 <- as.matrix(colMeans(wbasket)) #matriz de 12X1
  probs0 <- 1 - probs1
  #step 2: encontrar y almacenar todos los estados la muestra
  # identificar y dejar los estados unicos.
  unique_wb <- as.matrix(unique(wbasket[,c(1:ncol(wbasket))]))
  if (nrow(unique_wb) == nrow(wbasket) ) {
    # https://stackoverflow.com/questions/36699272/why-is-message-a-better-choice-than-print-in-r-for-writing-a-package
    message("Number of states is equal to number of samples!!")
  }
  # starting
  num_of_states <- nrow(unique_wb)
  suma_acumulada = 0
  state_probs <- matrix(NA, ncol=1, nrow = num_of_states)
  
  #step 3: para cada estado, calcular su probabilidad = P(s1)*P(s2)*....*P(sn)
  for (i in c(1:num_of_states)) { # vamos recorriendo todos los estados de la muestra
    estado <- unique_wb[i, , drop=FALSE] #matriz de 1X12
    id1 <- which(estado == 1) #sacamos el column id de estado en que valores son unos
    id0 <- which(estado == 0) #sacamos el column id de estado en que valores son ceros.
    valor1 <- prod(probs1[id1,])
    valor0 <- prod(probs0[id0,])
    state_probs[i,1] <- valor1*valor0
  }
  state_probs <- state_probs/sum(state_probs) #matriz de 1526 X 1
  #step 4: calcular la entropia = entropia + entropia del estado anterior
  logFd <- log2(state_probs)
  S1 <- -1*t((state_probs))%*%(as.matrix(logFd))
  
  # step 5: reportar entropia.
  return(S1)
}
# Ejemplo
#entropy_S1v2(all_widebaskets[[1]])
# Utilizando modelo independiente de poisson (viene de correlations_support.R)
entropy_S1poisson <- function(wbasket) {
  #step 0: preparar el input: wbasket
  borrar <- as.numeric(which(colMeans(wbasket) == 0)) #hay que quitar las columnas que tienen puros ceros.
  wbasket <- wbasket[, -c(1, borrar)] # hay que quitar la columna con el id
  rm(borrar)
  #step 1: calcular lambda
  spikes_sums <- rowSums(wbasket)
  N = ncol(wbasket)
  mean_spikes <- spikes_sums/N
  lambda <- sum(mean_spikes)/nrow(wbasket)
  ##step 2: encontrar y almacenar todos los estados la muestra
  # identificar y dejar los estados unicos.
  unique_wb <- as.matrix(unique(wbasket[,c(1:ncol(wbasket))]))
  if (nrow(unique_wb) == nrow(wbasket) ) {
    # https://stackoverflow.com/questions/36699272/why-is-message-a-better-choice-than-print-in-r-for-writing-a-package
    message("Number of states is equal to number of samples!!")
  }
  # step3: para cada elemento de unique_wb, encontrar su probabilidad poisson
  num_of_states <- nrow(unique_wb)
  state_probs <- matrix(NA, ncol=1, nrow = num_of_states)
  for (i in seq_along(1:num_of_states)) { # vamos recorriendo todos los estados de la muestra
    estado <- unique_wb[i, , drop=FALSE] #matriz de 1X12
    k <- sum(estado)
    state_probs[i,1] <-dpois(k, lambda)
  }
  state_probs <- state_probs/sum(state_probs) #matriz de 1526 X 1
  #step 4: calcular la entropia = entropia + entropia del estado anterior
  logFd <- log2(state_probs)
  S1 <- -1*t((state_probs))%*%(as.matrix(logFd))
  
  # step 5: reportar entropia.
  return(S1)
}
# ejemplo
#entropy_S1poisson(all_wbaskets[[3]])
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 03-sep-18
# Comentario:_
# La funcion get_energy calcula la energia total, de acople y de campo para un estado cuaqluiera
# Inputs:
# * vector y de spines en formato matrix
# H &  J vectro de fields y de acoples respectivamente
# modo = "f" o "c" o "t" para sacar energia de fields, energia de coupligs o energia total del estado
get_energy <- function(y, mode="f") {
  library(gdata)
  # step 1: get the state
  #vect <- as.matrix(wbasket[i,])
  vect <- y
  #if ( is.matrix(vect) == FALSE ) { #si NO es una matriz
  #  vect <- as.matrix(vect)
  #}
  # step 2: compute sum hi*si y guardarlo 
  Efield <- -1*vect%*%H
  #Efield <- H %*% vect
  
  # step 3: compute sum Jij*si*sj y guardarlo
  # https://stackoverflow.com/questions/41744915/convert-upper-triangular-part-of-a-matrix-to-3-column-long-format
  ind <- which(upper.tri(J, diag = FALSE), arr.ind = TRUE)
  dfw <- as.data.frame(cbind(ind, J[ind]))
  dfw$E <- vect[dfw$row]*vect[dfw$col]*dfw$V3
  Ecoupl <- -1*sum(dfw$E)
  
  # step 4: sumar lo del step 2 y step 3 y guardarlo.
  Etotal <- Efield + 0.5*Ecoupl
  
  if (mode == "f") {
    return(Efield)
  }
  if (mode == "c") {
    return(Ecoupl)
  } 
  if (mode == "t") {
    return(Etotal)
  }
}
# ejemplo:
# Efield <- apply(df[,c(1:ncol(df))], 1, get_energy, mode="f")
# Ecoupl <- apply(df[,c(1:ncol(df))], 1, get_energy, mode="c")
# Etotal <- apply(df[,c(1:ncol(df))], 1, get_energy, mode="t")
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 






# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 18-dic-18
# Comentario:_
# La funcion get_couplings_from_state nos reporta los acoples involucrados en un estado.
# Esta funcion nos sirve despues para saber si los estados son fustrados o no 
# multiplicando los acoples y si el signo es negativo, es estado es fustrado.
# Inputs:
# * vector y de spines en formato matrix
# * mode = "acoples" si quiero que retorne los coples, mode="sign" si quiero que retorne
# si es fustrado o no. -1 significa fustrado, 1 significa no fustrado.
# H &  J vectro de fields y de acoples respectivamente
# Salida
# * un vector con los acoples involucrados.
get_couplings_from_state <- function(y, mode="acoples") {
  library(gdata)
  vect <- y
  # step 1: compute the couplings involved
  # https://stackoverflow.com/questions/41744915/convert-upper-triangular-part-of-a-matrix-to-3-column-long-format
  ind <- which(upper.tri(J, diag = FALSE), arr.ind = TRUE)
  dfw <- as.data.frame(cbind(ind, J[ind]))
  temp <- vect[dfw$row]*vect[dfw$col]*dfw$V3
  # solo nos quedamos con los acoples de los estados removiendo todo lo que es cero.
  temp <- temp[temp != 0]
  if (mode == "acoples") {
    return(temp)
  } 
  if (mode == "sign") {
    multiplicacion <- prod(temp)
    return(sign(multiplicacion))
  }
  
}
# ejemplo:
#cpl <- get_couplings_from_state(y=as.matrix(df[100,]), mode="acoples")
#cpl <- get_couplings_from_state(y=as.matrix(df[100,]), mode="sign")
#cpl <- get_couplings_from_state(y=c(0, 0, 1, 0 , 0, 1, 1, 0, 0, 0, 0, 0))
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 






# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 15-jun-18
# Funcion para calcular la entropia S2 utilizando la PDF de boltzmann (SIN sampleo metropolis hasting)
# Inputs:
#       wbasket: matriz transaction 
#       J y H que son los acoples y magnetizaciones
# Output:
#       S2: entropia
entropy_S2 <- function(wbasket, J, H) {
  borrar <- as.numeric(which(colMeans(wbasket) == 0))
  wbasket <- wbasket[, -c(1, borrar)]
  rm(borrar)
  # identificar y dejar los estados unicos.
  unique_wb <- as.matrix(unique(wbasket[,c(1:ncol(wbasket))]))
  if (nrow(unique_wb) == nrow(wbasket) ) {
    # https://stackoverflow.com/questions/36699272/why-is-message-a-better-choice-than-print-in-r-for-writing-a-package
    message("Number of states is equal to number of samples!!")
  }
  #transformamos a spines {-1,1}
  #unique_wb <- 2*unique_wb - 1 
  # starting
  num_of_states <- nrow(unique_wb)
  suma_acumulada = 0
  state_probs <- matrix(NA, ncol=1, nrow = num_of_states)
  #idx <- combn(x = colnames(J), m = 2)
  idx <- combn(x = nrow(H), m = 2)
  largo <- ncol(idx)
  for (i in c(1:num_of_states)) { # vamos recorriendo todos los estados de la muestra
    estado <- t(as.matrix(unique_wb[i, ]))
    #sum_magnetizaciones <- -1*sum(estado%*%H)
    sum_magnetizaciones <- estado %*% H
    sum_acoples = 0
    for (j in c(1:largo)) { # calculamos el sum_{i<>j} Jij*Si*Sj
      fila <- idx[1,j]
      colu <- idx[2,j]
      sum_acoples <- sum_acoples + J[fila,colu]*estado[fila]*estado[colu]
    }
    #sum_acoples <- -0.5*sum_acoples
    #sum_acoples <- sum_acoples
    state_probs[i,1] <- exp(sum_magnetizaciones + sum_acoples )
  }
  Z <- sum(state_probs) # funcion de particion
  # P({S})
  state_probs <- state_probs/Z
  S2 <- -t(state_probs) %*% log2(state_probs)
  return(S2)
}
# Ejemplo
#entpy <- entropy_S2(wbasket = all_widebaskets[[1]], J = a[[2]], H =  a[[1]] )
# Esta version toma las 2^n combinaciones
entropy_S2v2 <- function(J, H) {
  num_of_nodes <- nrow(H)
  num_of_states <- 2^nrow(H)
  Estados <- expand.grid(replicate(num_of_nodes, 0:1, simplify = FALSE))
  #transformamos a spines {-1,1}
  #Estados <- as.matrix(2*Estados - 1)
  # starting
  suma_acumulada = 0
  state_probs <- matrix(NA, ncol=1, nrow = num_of_states)
  #idx <- combn(x = colnames(J), m = 2)
  idx <- combn(x = nrow(H), m = 2)
  largo <- ncol(idx)
  for (i in c(1:num_of_states)) { # vamos recorriendo todos los estados de la muestra
    estado <- as.matrix(Estados[i, ])
    #sum_magnetizaciones <- -1*sum(estado%*%H)
    sum_magnetizaciones <- estado %*% H
    sum_acoples = 0
    for (j in c(1:largo)) { # calculamos el sum_{i<>j} Jij*Si*Sj
      fila <- idx[1,j]
      colu <- idx[2,j]
      sum_acoples <- sum_acoples + J[fila,colu]*estado[fila]*estado[colu]
    }
    #sum_acoples <- -0.5*sum_acoples
    #sum_acoples <- sum_acoples
    state_probs[i,1] <- exp(sum_magnetizaciones + sum_acoples )
  }
  Z <- sum(state_probs) #funcion de particion
  # P({S})
  state_probs <- state_probs/Z
  S2 <- -t(state_probs)%*% log2(state_probs)
  return(S2)
}
# Calcula entropia de S2 a partir de un sampleo de metropolis hastion
entropy_S2_mh <- function(wbasket, J, H) {
  borrar <- as.numeric(which(colMeans(wbasket) == 0))
  wbasket <- wbasket[, -c(1, borrar)]
  rm(borrar)
  filas <- nrow(wbasket)
  rm(wbasket)
  
  initial_spins=rep(1,ncol(J))
  sampling <- metropolis.sampling(magnetizaciones=H, acoples = J, initial_spin = initial_spins, kmax=filas, ks=0, To=NULL, alpha=NULL)
  wbasket <- sampling[[1]]
  
  S2 <- entropy_SN(wbasket, drop_zeros = "off", drop_newid = "off" )
  return(S2)
}
# ejemplo:
#entropy_S2_mh(wbasket = all_wbaskets[[12]], J = couplings[[muestra]]$J, H = couplings[[muestra]]$H)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 






# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Funcion para calcular la entropia SN del model entero
# Esta funcion utiliza la funcion frecuency que esta en ising_functions_v3.R
# Inputs:
#       wbasket: matriz transaction 
#       drop_zeros = "on" si deseo quitar columnas de wbasket que tengan putos ceros.  
#       drop_newide = "on" si widebasket viene con columna newid, hay que quitarla
# Output:
#       SN: entropia
entropy_SN <- function(wbasket, drop_zeros = "on", drop_newid = "on") {
  if (drop_newid == "on") {
    wbasket <- wbasket[, c(-1)]
  }
  if (drop_zeros == "on") {
    borrar <- as.numeric(which(colMeans(wbasket) == 0)) #hay que quitar las columnas que tienen puros ceros.
    wbasket <- wbasket[, -c(borrar)] # hay que quitar la columna con el newid
    rm(borrar)
  }
  ff <- frecuency(muestra = wbasket )
  #ff[[1]] #matriz con los estados contados
  #ff[[2]] #vector con el conteo de cada estado
  N <- nrow(wbasket) 
  sum(ff[[2]], na.rm=TRUE)  #42077 muestras que es igual al numero de filas de wbasket.
  Fd <- ff[[2]]/N
  #quitamos los NA y ceros de Fd. Esto ocurre porque frecuency genera 2^n estados. Algunos sobran.
  idx <- (Fd == 0 | is.na(Fd))
  Fd <- Fd[!idx]
  logFd <- log2(Fd)
  SN <- -1*t(as.matrix(Fd))%*%(as.matrix(logFd)) 
  return(SN)
}
#Ejemplo
#entpy <- entropy_SN(wide_baskets[[1]], drop_zeros = "on, drop_newid = "on)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 





# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#02-jun-18
#Ahora vamos a calcular la divergencia KL(PN||P1) que equivale a IN segun el paper 
# "market structure explained by pairwise interactions" de T. Bury y calcular el I2/IN.
# PN es la pdf de wbasket real, P1 es la pdf del modelo independiente.
# Input: wbasket
# Output: IN que equivale al KL(PN||P1) = S1 - SN
# Lo que hace esta funcion es calcular KL(PN||P1) a partir de wbasket.
multi_information <- function(wbasket) {
  if (!is.matrix(wbasket)) {
    wbasket <- as.matrix(wbasket)
  }
  # chequeo si muestra viene con columna newid
  if ("newid" %in% colnames(wbasket))
  {
    borrar <- which( colnames(wbasket)== "newid" )
    muestra <- wbasket[, -c(borrar)]
  }
  borrar <- as.numeric(which(colMeans(wbasket) == 0)) #hay que quitar las columnas que tienen puros ceros.
  wbasket <- wbasket[, -c(1,borrar)] # hay que quitar la columna con el newid
  
  # Step 1: buscar todos los estados de PN
  # identificar y dejar los estados unicos.
  # Step 2: encontrar el conteo  de todos los estados de PN
  ff <- frecuencyV2(muestra = wbasket ) #frecuencyV2 es 1 segundo mas rapido que frecuency pq no recorre los 2^n estados
  P <- ff[[2]]/sum(ff[[2]])
  unique_wb <- ff[[1]]
  num_of_states <- nrow(unique_wb)
  #sum(ff[[2]], na.rm=TRUE)  #42077 muestras que es igual al numero de filas de wbasket. ###Importante
  # Step 3: Encontrar las probabilidades de P1 (modelo independiente)
  probs1 <- as.matrix(colMeans(wbasket)) #matriz de 12X1
  probs0 <- 1 - probs1
  state_probs <- matrix(NA, ncol=1, nrow = num_of_states) 
  #
  for (i in c(1:num_of_states)) { # vamos recorriendo todos los estados de la muestra
    estado <- unique_wb[i, , drop=FALSE] #matriz de 1X12
    id1 <- which(estado == 1) #sacamos el column id de estado en que valores son unos
    id0 <- which(estado == 0) #sacamos el column id de estado en que valores son ceros.
    valor1 <- prod(probs1[id1,])
    valor0 <- prod(probs0[id0,])
    state_probs[i,1] <- valor1*valor0
  }
  Q <- state_probs/sum(state_probs) #matriz de 1526 X 1   
  # Step 3: calculo de KL.
  m <- matrix(NA, ncol=3, nrow=nrow(Q))
  m[,1] <- as.matrix(P)
  m[,2] <- as.matrix(Q)
  m[,3] <- as.matrix(log2(m[,1]/m[,2]))
  colnames(m) <- c("P", "Q", "cuo")
  library(entropy)
  KL <- KL.empirical(m[,1], m[,2], unit="log2") #da 0.5 para all_widebaskets[[1]]
  return(as.numeric(KL))
}
# Ejemplo:
#multi_information(all_widebaskets[[1]])
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 





# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 08-sep-18
# Comentario:_
#Funcion para hacer rbind de acuerdo a los colnames de los dataframes
#https://stackoverflow.com/questions/16962576/how-can-i-rbind-vectors-matching-their-column-names
# inputs: x1, x2 : x1 es dataframe grande al cual deseamos anexar la data frame x2 en los cumeros
# de columna correspondiente.
# output : es la dataframe x1 con las filas de x2 agregadas en las columnas correspondientes de x1.
custom_rbind <- function(x1,x2){
  c1 = setdiff(colnames(x1),colnames(x2))
  c2 = setdiff(colnames(x2),colnames(x1))
  for(c in c2){##Adding missing columns from 2 in 1
    x1[[c]]=NA
  }
  for(c in c1){##Similiarly ading missing from 1 in 2
    x2[[c]]=NA
  }
  x2 = x2[colnames(x1)]
  rbind(x1,x2)
}
# ejemplo:
#M <- custom_rbind(M, df_f)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# ******** Basado en los resultados de simulation_hc_V6.R y 
# simulation_hc_V7.R, vemos que al agregar mas randomness al 
# algoritmo greedy, la media de distancia de acople d_c va 
# aumentando en relación al greedy deterministico. Asimismo, la 
# varianza tambien aumenta, lo cual indica que a veces, el greedy 
# probabilistico hace muy bien el trabajo buscando una distancia de
# acople menor al greedy deterministico, pero a veces tambien lo hace 
# muy mal. Por lo tanto, es de suponer, que el greedy probabilistico 
# necesita una ayuda. Lo que hacemos en este script es dar un punto de 
# partida basado en el greedy deterministico, para mejor aun mas el 
# desempeno del greedy probabilistico, disminuyendo su variabilidad 
# en los resultado y tratando de estar siempre mejor que el deterministico.

# Para implementar lo anterior, crearemos una nueva funcion que se llama
# hierarchical_clustering_probabilistic_greedy_improved.

# Procedimiento:
# 1. Se comienza por encontrar una solucion con greedy deterministico.
#     se guardan los resultados de la distancia de acople lograda 
#     por greedy deterministico en cada iteracion. Esto sera de benchmark
#     para el greedy probabilistico.
#     merge_greedy_detm <- hierarchical_clustering_probabilistic_greedy(D = D, J = J, trh=1 ) 
# en cada iteracion:
# 2.  en la Iteracion X, partimos con la solucion de greedy deterministico.
#     X puede ser igual a 1, y luego seguimos con greedy probabilistico.
# 3. Podemos generar N paths (simulaciones) y elegir aquella que supere
#   al path de benchmark deterministico por mayor numero de iteraciones
#   o aquella que se mantenga con la menor pendiente hasta cierto numero
#   de iteraciones.



# actual name: testing_hierarchical_clusteringV13.R

#Notas:
# 26-dic-19: creation 
# 30-dic-19: 



# # # # # # # # # # # # # # # # # # LOADING DATA AND FUNCTIONS # # # # # # # # # # # # # # # # # # # # # # #
# Carga de datos, 
# igual que en clusteragrr_v1.R
rm(list = ls())
#load('new_inferring_parameters_environment270219.RData') # load H and J parameters ising inferred
source("functions_hclust.R") # funciones varias
source("hierarchical_clustering_v2_function.R")
source("hierarchical_clustering_v5_function.R")
source("hierarchical_clustering_greedy_function.R")
source("hierarchical_clustering_probabilistic_greedy_function.R")

source("create_coupling_function.R") # crea matriz J de acoples y D de distancias para simular
source("create_mst_function.R") # crea MST a partir de una matriz de acople.
source("acople_distance_sum_function.R") # suma las distantcias de acople dc
source("find_ady_dc_function.R") # encuentra los nodos adyacentes en un mst dado un nodo
source("get_num_nodes_function.R") # nos dice el numero de nodos que tiene el par de clusters a fusionar
source("get_num_nodes_of_ady_clusters_function.R") # nos dice informacion del numero de nodos que son adyacentes a un cluster
source("get_name_of_the_cluster_function.R") # dado el nombre de un nodo, nos dice a que cluster pertenece
source("get_nodes_of_the_cluster_function.R") # nos da los nodos involucrados en un cluster.
source("simulation_hc_function.R") # Genera una simulacion de hierarchical clustering con algoritmo single linkage y modificado.
source("find_min_distcp_function.R") #permite encontrar las distancias de acople entre un cluster y los restantes
source("pick_a_cluster_function.R") # necesario apra correr algoritmo greedy aprobabilistico

source("find_mst_barrier_function.R")
source("is_integer0_function.R")
source("entropies_functions.R") # nos interesa aqui get_energy function

library(igraph)
library(Matrix)
library(ggplot2)
# # # # # # # # # # # # # # # # # # LOADING DATA AND FUNCTIONS # # # # # # # # # # # # # # # # # # # # # # #


# # # # # # # # # # # # # # # # # GETTING THE NODES OF CLUSTERS # # # # # # # # # # # # #  # # # # # # # # #
# simula 
Jmean = 0
Nn=20
sj=1

# set.seed(123)
ot <- create_coupling(Nn=Nn, media=Jmean, sj=sj, lw=-3, up=3)
J <- ot$J
D <- ot$D
# # # # # # # # # # # # # # # # # GETTING THE NODES OF CLUSTERS # # # # # # # # # # # # #  # # # # # # # # #




# dic 27, 2019
# modification to function in hierarchical_clustering_probabilistic_greedy_function.R
# modification:
# 1. In the first iteration, the algorithm runs with trh=1, ie, deterministic
#     so, the first merge is at minimum coupling energy.
# 2. After first iteration, the algorithm follows as normal.
# #

#   * trh = threshold entre 0 y 1: cuando es 0, es totalmente random, 1 es totalmente deterministico
hc_probabilistic_greedy_withelp <- function(D, J, trh=0) {
  # inicio - preparacion
  N <- ncol(D)
  merge <- matrix(NA, ncol=6, nrow=N-1) # col1 y col2 spins merged, col3 numero de cluster, col4 = merge distance
  merge <- as.data.frame(merge)
  colnames(merge) <- c("node1", "node2", "cluster", "prob", "dultr", "dc") # prob es la probabilidad de activacion del cluster seleccionado 
  iteraciones <- N - 2
  mn <- list() # aqui iremos agregando los nombres de los spins que componen cada clusters 
  
  # Inicio de las iteraciones que son N - 2
  for (it in 1:iteraciones) {
    if ( ( runif(1) > trh )  && ( it > 1 ) ) {
      # # # # # # # # #  Find probabilies of cluster selection:
      rl <- pick_a_cluster(D)
      id <- rl$id
      nombres <- rl$nombres
      m <- rl$m
      prob <- rl$p
      # leaves (original spins) will be negative, merged clusters will be positive starting from 1, 2, ...
      nombres_num <- as.numeric(nombres) # nombres de los spins
    } else {
      # # # # # # # # #  Find the minimum distance among clusters:
      rl <- cluster_find_name(D)
      id <- rl$id
      nombres <- rl$nombres
      m <- rl$m
      # leaves (original spins) will be negative, merged clusters will be positive starting from 1, 2, ...
      nombres_num <- as.numeric(nombres) # nombres de los spins
      prob <- 1
    }
    
    # # # # # # # # # poniendo los nombres reales de los spins que componen el cluster recien formado.
    # detectamos si un nombre de cluster es positivo, en tal caso, extraer los nombres de los spins de ese cluster
    te <- nombres_num
    hay_cluster_anterior <- which(te > 0)
    if( !is.integer0(hay_cluster_anterior) ) {
      el_cluster <- nombres_num[hay_cluster_anterior]
      # rescatar los nombres del cluster en cuestion
      te2 <- c(nombres_num, unlist(mn[el_cluster]) ) 
      te <- te2[-hay_cluster_anterior]
    }
    mn <- c(mn, list( as.character(te )) ) # y ahora en mn[[2]] cluster 2, se compone de -5 y -6
    
    
    # calculo de las distancias de acople
    d_cpl <- acople_distance_sum2(J, te)  ### > OJo que esta funcion esta en acople_distance_sum_function.R
    
    merge[it,] <- as.numeric( c(nombres_num, it, prob, m, d_cpl) ) #  nuevo: se coloca nodo1 (fusion) / nodo2 (fusion) /nombre del nuevo cluster / distancia a la que se unen
    
    
    # # # # # # # # # Upgrading the distance matrix
    # fila id[1] y columna id[2] se deben borrar
    Dold <- D
    D <- D[-id, -id]
    
    # pero antes hay que buscar las distancias minimas entre los spins restantes y el cluster
    if (length(Dold) <= 9) { # cuando la matriz es de 3X3, no se puede ejecutar directamente la funcion find_min_dist
      # 09-nov-19
      colno <- colnames(Dold)
      sel <- colno[!(colno %in% as.character(nombres_num) )] # nodo identificado 16-oct-19
      involved_nodes <-  c( mn[[it]], sel )
      tempo <- as.numeric(involved_nodes )
      haypositive <- which(tempo > 0)
      if ( length(haypositive) > 0) {
        involved_nodes <- c( mn[[it]], mn[[ tempo[haypositive] ]]  )
      }
      dit <-  acople_distance_sum2(J, y = involved_nodes) 
      nombres_ <- sel
      # FIN # 09-nov-19
    } else {
      #temp <- find_min_dist(D, Dold, nombres)
      temp <- find_min_distcp(J = J, D = D, nombres = as.character(nombres_num), mn=mn, it=it ) # 16-oct-19
      dit <- temp$dit # dit son las minimias distancias de cluster fusionado en la iteracion anterior a todos los demas clusters que quedan
      nombres_ = temp$nombres_p # nombres_ son los nombres de los nodos que estan a minima distancia.
    }
    
    #ahora tengo que poner los valores en fila.columna de D
    D <- put_dis(D, dit, nombres_, N, it)
    
    
    #print(m)
    print(paste("Probabilistic Greedy Iteration: ", it))
    
  }
  
  # Last iteration: when it = N - 1
  it <- it + 1
  rl <- cluster_find_name(D)
  id <- rl$id
  nombres <- rl$nombres
  m <- rl$m
  nombres_num <- as.numeric(nombres)
  prob <- 1
  
  te <- nombres_num
  hay_cluster_anterior <- which(te > 0)
  if(!is.integer0(hay_cluster_anterior) ) {
    el_cluster <- nombres_num[hay_cluster_anterior]
    # rescatar los nombres del cluster en cuestion
    te2 <- c(nombres_num, unlist(mn[el_cluster]) ) 
    te <- te2[-hay_cluster_anterior]
  }
  mn <- c(mn, list( as.character(te )) ) # y ahor
  
  # calculo de las distancias de acople
  d_cpl <- acople_distance_sum2(J, te)  ### > OJo que esta funcion esta en acople_distance_sum_function.R
  
  merge[it,] <- as.numeric( c(nombres_num, it, prob, m, d_cpl) ) #  nuevo: se coloca nodo1 (fusion) / nodo2 (fusion) /nombre del nuevo cluster / distancia a la que se unen
  
  
  return(merge)
  
}
# ejemplo:


# deterministico 
invisible(capture.output( merge_greedy_detm <- hc_probabilistic_greedy_withelp(D = D, J = J, trh= 1 ) ) )
invisible(capture.output( merge_greedy_detm <- hierarchical_clustering_probabilistic_greedy(D = D, J = J, trh=1 ) ) )
# probabilistico 100% con ayuda en la primera iteracion
invisible(capture.output( merge_greedy_help <- hc_probabilistic_greedy_withelp(D = D, J = J, trh = 0 ) ) )
# probabilistico 100% sin ayuda
invisible(capture.output( merge_greedy_prob <- hierarchical_clustering_probabilistic_greedy(D = D, J = J, trh = 0 ) ) )

# transformacion a objeto dendograma para ver el dendograma
hc_ghelp <- to_dendo2(merge_greedy_help[, c(1,2,3,6) ], enames = c(1:ncol(D)))
hc_gprob <- to_dendo2(merge_greedy_prob [, c(1,2,3,6) ], enames = c(1:ncol(D)))
plot(hc_ghelp)
plot(hc_gprob)

# comparacion de desempeno tomando siempre como benchmarck el deterministico
plot(merge_greedy_detm[,3], merge_greedy_detm[,6], type="l", col="red")
par(new=T)
plot(merge_greedy_prob[,3], merge_greedy_prob[,6], type="l", col="blue")
plot(merge_greedy_detm[,3], merge_greedy_prob[,6]-merge_greedy_detm[,6], type="l", col="black")
abline(h = 0, lty=2, col="red")
# la idea es que la suma sea negativa o lo mas chica. 
sum(merge_greedy_prob[,6]-merge_greedy_detm[,6])
# la idea es que sea grande este numero.
100*sum(merge_greedy_prob[,6]-merge_greedy_detm[,6] <= 0)/nrow(merge_greedy_prob)

plot(merge_greedy_detm[,3], merge_greedy_detm[,6], type="l", col="red")
par(new=T)
plot(merge_greedy_help[,3], merge_greedy_help[,6], type="l", col="green")
plot(merge_greedy_detm[,3], merge_greedy_help[,6]-merge_greedy_detm[,6], type="l", col="black")
abline(h = 0, lty=2, col="red")
sum(merge_greedy_help[,6]-merge_greedy_detm[,6])
100*sum(merge_greedy_help[,6]-merge_greedy_detm[,6] <= 0)/nrow(merge_greedy_prob)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


# dic 30, 2019
# 1. we generate a matrix J for N 25 nodes
# 2. we get the path (coupling distance vs iterations) solution for greedy deterministic
# 3. we get I path-simulations using greedy probabilistic with help
# 4. we get I path-simulations using greedy probabilistic 100% probabilistic
# 5. we get I path-simulations using greedy probabilistic 50% probabilistic
# 6. from 3, 4 and 5 we compute the means of coupling distance for each iteration
# 7. graph 6
# 8. we build a function that computes sum Delta  and (num de veces que  Delta dc es negativo sea maximo)
#     and the we choose the best one
# 9. graph 8.

# por hacer:
# una simulacion para N = 25, para un J dado, y repetir muchas
# veces el greedy probabilistico 100%, 50% y con 100% con ayuda 
# para obtener varios paths. Cada uno de esos paths compararlos con 
# el greedy deterministico, y elegir uno que:
# * sum Delta dc sea minimo
# * num de veces que  Delta dc es negativo sea maximo
# * aquel que tenga Delta dc negativo mas grande.

# 2. we get the path (coupling distance vs iterations) solution for greedy deterministic
invisible(capture.output( merge_greedy_detm <- hierarchical_clustering_probabilistic_greedy(D = D, J = J, trh=1 ) ) )

I = 500
ware <- list()
ware[["greedy_prob"]] = list()
ware[["greedy_help"]] = list()
ware[["d1"]] = list()  # delta distancia de acople entre probabilistico y deterministico
ware[["d2"]] = list()  # delta distancia de acople entre probabilistico con help y deterministico

for (i in 1:I) {
  # creating solutions 
  # greedy probabilistico 100%
  invisible(capture.output( merge_greedy_prob <- hierarchical_clustering_probabilistic_greedy(D = D, J = J, trh = 0 ) ) )
  # probabilistic greedy with help
  invisible(capture.output( merge_greedy_help <- hc_probabilistic_greedy_withelp(D = D, J = J, trh = 0 ) ) )
  
  # savings
  ware[["greedy_detm"]][[i]] <- merge_greedy_detm[, c(1,2,6)]
  ware[["greedy_prob"]][[i]] <- merge_greedy_prob[, c(1,2,6)]
  ware[["greedy_help"]][[i]] <- merge_greedy_help[, c(1,2,6)]
  #
  
  # computing delta distance coupling
  ware[["d1"]][[i]] <- (merge_greedy_prob[,6]-merge_greedy_detm[,6])
  ware[["d2"]][[i]] <- (merge_greedy_help[,6]-merge_greedy_detm[,6])
  
  cat("\n Iteration number :", i) 
  
}


# getting the means and creating the final dataframe to plot
# https://stackoverflow.com/questions/31839996/mean-of-each-row-of-data-frames-column-nested-in-list
Gprob <- do.call(cbind, lapply(ware[["greedy_prob"]], "[", "dc"))
dc_mean_Gprob <- rowMeans(Gprob)
#dc_sdev_Gsemi <- apply(Gsemi,1,sd)
Ghelp <- do.call(cbind, lapply(ware[["greedy_help"]], "[", "dc"))
dc_mean_Ghelp <- rowMeans(Ghelp)
#dc_sdev_Gprob <- apply(Gprob,1,sd)
# delta distancia de acople entre probabilistico y deterministico
d1 <- unlist(lapply(ware[["d1"]], sum)) 
nd1 <- unlist(lapply(ware[["d1"]], function(x)  100*sum(x <= 0)/(Nn-1)))

# delta distancia de acople entre probabilistico con help y deterministico
d2 <- unlist(lapply(ware[["d2"]], sum))
nd2 <- unlist(lapply(ware[["d2"]], function(x)  100*sum(x <= 0)/(Nn-1)))




# Preparing dataframe to ggplot
iteration <- rep(seq(from = 1, to = Nn-1),3)
algo <- c(rep("Gdet", Nn-1), rep("Gprob", Nn-1), rep("Ghelp", Nn-1))
medias <- c(merge_greedy_detm$dc, dc_mean_Gprob, dc_mean_Ghelp)
#desvs <- c(rep(NA, N-1), dc_sdev_Gsemi, dc_sdev_Gprob)
df <- data.frame(iteration = iteration, Algorithm = algo, 
                 mean = medias)

# doing the ggplot
f <- ggplot(df, aes(x = iteration, y = log(mean), colour = Algorithm)) + 
  geom_point(shape=20, size=1.5) +
  #geom_smooth(show.legend = FALSE, fill = "gray", span = 0.3) +
  #geom_errorbar(aes(ymin=min, ymax=max ), width=.2, position=position_dodge(.01)) +
  #geom_errorbar(aes(ymin=log(mean-sd/sqrt(3)), ymax=log(mean+sd/sqrt(3))), width=.2, position=position_dodge(.01)) +
  labs(x ="Iteration", y = "mean of log(U*)") +
  theme_linedraw() + theme(legend.position="bottom") +
  #scale_x_continuous(breaks = seq(from=1, to=N-1, by=2) ) +
  scale_y_continuous(breaks = c(-1, 0, 1, 2, 3,4,5,6,7,8)) +
  #scale_y_continuous(breaks = seq(from=1, to=8) ) +
  theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "white"))
f

# graph of delta coupling distances
plot(1:I, d1, type="l", col="red") # delta distancia de acople entre probabilistico y deterministico
par(new=T)
plot(1:I, d2, type="l", col="blue") # delta distancia de acople entre probabilistico con help y deterministico
par(mfrow=c(2,1))
hist(d1)
abline(v=mean(d1),col="red")
hist(d2)
abline(v=mean(d2),col="red")

# graph of delta coupling distances <= 0
plot(1:I, nd1, type="l", col="red") # delta distancia de acople entre probabilistico y deterministico
par(new=T)
plot(1:I, nd2, type="l", col="blue") # delta distancia de acople entre probabilistico con help y deterministico
par(mfrow=c(2,1))
hist(nd1)
abline(v=mean(nd1),col="red")
hist(nd2)
abline(v=mean(nd2),col="red")


which()





# evaluacion de clusters.
# https://www.datanovia.com/en/lessons/determining-the-optimal-number-of-clusters-3-must-know-methods/



# En este script vamos a utilizar la funcion hierarchical_clustering_v2
# y funcion hierarchical_clustering_v4 que estan en functions_hclust.R
# para siguiente simulacion:

# Generar una matriz J con J distribuido N(0,1) truncado en +-3 con dimension NXN
# calcular las distancias siguiendo la regla de sqrt(3-J)
# calcular el MST.
# calcular single linkage normal V2, y single linkage modificado V4
# Comparar la distribucion de las dmst y dc para cada caso
# comparar los energy path de cada caso.

# Lo anterior se repite muchas veces para distinto numero de N, mu y sj
# donde N es el numero de nodos, mu es la media de J, y sj es la desv estandar de J.

# actual name: simulation_hc_V2.R
# 22.oct.19


# Notes:
# 22-oct-19: creation
# 23-oct-19: incorporamos a la simulacion las funciones create_copupling_function.R
#             y create_mst_function.R
# 29-oct-19: se modifica script para considerar version V9 de la funcion hierarchical_clustering.
# 04-nov-19: vemos que tomar el promedio de las distancias de acople de todas las iteraciones
#         de merge2 y merge5, no nos dice mucho, porque vemos que las medias de estas 
#         distancias de acoples son muy similares. Los datos quedaron guardados en 
#         results_031119_from_simulation_hv_V3.RData. 
#         Dado estos resultados, sera mejor calcular como ejemplo, una matriz J aleatoria
#         y graficar la distancia de acople de merge2 y merge5. Y hacer esto para 
#         N = N <- c(10, 20, 30, 50, 100, 200).
# 12-nov-19: en el analisis de L339, incorporamos las funciones para clustering jerarquico
#           con greedy algorithm.



# # # # # # # # # # # # # # # # # # # # # LOADING DATA # # # # # # # # # # # # # # # # # # # # # # # # #
# Carga de datos, 
# igual que en clusteragrr_v1.R
rm(list = ls())
#load('new_inferring_parameters_environment270219.RData') # load H and J parameters ising inferred
source("functions_hclust.R") # funciones varias
source("hierarchical_clustering_v2_function.R")
source("hierarchical_clustering_v5_function.R")
source("hierarchical_clustering_greedy_function.R")

source("create_coupling_function.R") # crea matriz J de acoples y D de distancias para simular
source("create_mst_function.R") # crea MST a partir de una matriz de acople.
source("acople_distance_sum_function.R") # suma las distantcias de acople dc
source("find_ady_dc_function.R") # encuentra los nodos adyacentes en un mst dado un nodo
source("get_num_nodes_function.R") # nos dice el numero de nodos que tiene el par de clusters a fusionar
source("get_num_nodes_of_ady_clusters_function.R") # nos dice informacion del numero de nodos que son adyacentes a un cluster
source("get_name_of_the_cluster_function.R") # dado el nombre de un nodo, nos dice a que cluster pertenece
source("get_nodes_of_the_cluster_function.R") # nod da los nodos involucrados en un cluster.
source("simulation_hc_function.R") # Genera una simulacion de hierarchical clustering con algoritmo single linkage y modificado.
source("find_min_distcp_function.R") # funcion que encuentra min dist acople en greedy alg.


source("find_mst_barrier_function.R")
source("is_integer0_function.R")
source("entropies_functions.R") # nos interesa aqui get_energy function
source("multiplot_function.R")

library(igraph)
library(Matrix)
library(ggplot2)

# Funcion:
# Input: df, que es merge en formato data.frame
fit_mean <- function(df) {
  hasta <- nrow(df)-1 
  apply(df[c(2:hasta), 4:6], 2, mean)
}  # calcula la media de las ultimas tres columnas y de la fila 2 a la 4 de merge
# # # # # # # # # # # # # # # # # # # # # LOADING DATA # # # # # # # # # # # # # # # # # # # # # # # # #




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# SIMULATION TO GET A LOT OF MERGE MATRICES FROM A RANDOM J COUPLING MATRIX
rep = 1000 # numero de repeticiones: la menos 1000
Nn = 10 # numero de nodos de la red de acoples
listm2 <- vector(mode = "list", length = rep) # guardamos los merges de V2
listm5 <- vector(mode = "list", length = rep) # guardamos los merges de V5
for (i in 1:rep) {
  
  obj <- simulation_hc_v1(N=Nn, mu=Jmean, sj=1)
  m2 <- obj$merge2
  m5 <- obj$merge5
  
  print(paste("ending iteration", i))
  #analysis <- rbind(analysis, difs)
  listm2[[i]] <- m2
  listm5[[i]] <- m5
  
}
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


# OBSERVE THE COUPLING (and all) DISTANCES BETWEEN TWO ALGORITHMS
# aqui calculamos el incremento de distancias de acople para cada algoritmo
library(purrr)
R <- map(listm2, fit_mean)
M2 <-  matrix(NA, ncol=3, nrow=rep)
for (i in 1:rep) {
  M2[i, ] <- R[[i]]
}
colnames(M2) <- c("mean_dultr", "mean_dmst", "mean_dcoupl")
M2 <- as.data.frame(M2)

R <- map(listm5, fit_mean)
M5 <-  matrix(NA, ncol=3, nrow=rep)
for (i in 1:rep) {
  M5[i, ] <- R[[i]]
}
colnames(M5) <- c("mean_dultr", "mean_dmst", "mean_dcoupl")
M5 <- as.data.frame(M5)

M <- rbind(M2, M5)
M$alg <- c(rep("Normal", rep), rep("Modified", rep))






library(ggplot2)
#violin plot
dp <- ggplot(M, aes(x=alg, y=mean_dcoupl, fill=alg)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="Plot of Coupling Distance",x="Algorithm", y = "Dc")
dp + theme_classic()

# boxplot
# https://owi.usgs.gov/blog/boxplots/
library(cowplot)
legend_plot <- ggplot_box_legend()
bx <- ggplot(data = M, aes(x = alg, y = mean_dcoupl)) +
  geom_boxplot() 
  #boxplot_framework(upper_limit = 70) + 
  xlab(label = "Algorithm") 
  #ylab(label = cl_name) +
  #labs(title = cl_site)

plot_grid(bx, legend_plot, nrow = 1, rel_widths = c(.6,.4))

# scatterplot entre d_ultr y d_c 
# https://stackoverflow.com/questions/11838278/plot-with-conditional-colors-based-on-values-in-r
# https://www.biostars.org/p/296296/
ggplot(M) + geom_point(aes(x = mean_dultr, y = mean_dcoupl, colour = alg )) +
  scale_colour_manual(name = 'Algorithm', values = setNames(c('red','lightgrey'), c("Normal", "Modified"))) +
  xlab('<dultr>') + ylab('<dcoupl>')


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 30.10.19
# vamos a fabricar una funcion que nos de la media y varianza de las distancias
# tomando los "rep" merges para algoritmo V2 y V5 
# Inputs:
#   listm2 y listm5 que son las listas con todas las matrices merges 
# NOTA: la funcion fit_mean debe estar cargada.
# rep: numero de merges que hay en cada lista
compute_mv_distances <- function(lista2, lista5, rep ) {
  library(purrr)
  library(dplyr)
  
  R <- map(lista2, fit_mean)
  M2 <-  matrix(NA, ncol=3, nrow=rep)
  for (i in 1:rep) {
    M2[i, ] <- R[[i]]
  }
  colnames(M2) <- c("mean_dultr", "mean_dmst", "mean_dcoupl")
  M2 <- as.data.frame(M2)
  
  R <- map(lista5, fit_mean)
  M5 <-  matrix(NA, ncol=3, nrow=rep)
  for (i in 1:rep) {
    M5[i, ] <- R[[i]]
  }
  colnames(M5) <- c("mean_dultr", "mean_dmst", "mean_dcoupl")
  M5 <- as.data.frame(M5)
  
  M <- rbind(M2, M5)
  M$alg <- c(rep("Normal", rep), rep("Modified", rep))
  
  medias <- M %>% 
    group_by(alg) %>% 
    summarise(avgdu = mean(mean_dultr), avgdm = mean(mean_dmst), avgco = mean(mean_dcoupl))
  
  return(medias)
}
medias <- compute_mv_distances(lista2 = listm2, lista5 = listm5, rep = rep)
medias$N <- rep(Nn, 2)
medias$Jmean <- rep(Jmean, 2)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 














# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# SIMULATION TO OBSERVE DIFFERENCES BETWEEN dMST AND COUPLING DISTANCES.

# Generacion de las diferencias en distancias entre el signle linkage y single linkage modificado
N=6
NN <- N - 1
difs <- matrix(NA, ncol=4, nrow=N-1)
colnames(difs) <- c("merge", "Dc", "Dmst", "Ddultr")
# vamos calculando Ddc =  dc(single linkage) - dc(single linkage modificado)  x$dc - y$dc
# vamos calculando Dmst =  dmst(single linkage) - dmst(single linkage modificado)  x$dmst - y$dmst
# vamos calculando Dultr =  dultr(single linkage) - dultr(single linkage modificado)  x$dultr - y$dultr
difs[, ] <-  cbind(c(1:NN), x$dultr - y$dultr, x$dmst - y$dmst, x$dc - y$dc)


# Aqui calculamos la diferencia de las distancias enter los dos algoritmos
# numero de simulaciones
rep = 10 # numero de repeticiones: la menos 1000
Nn = 6 # numero de nodos de la red de acoples
NN <- Nn - 1
# aqui guardamos una lista que contiene todas las diferencias entre distancias de mst y 
# de acople para cada simulacion en cada iteracion del proceso de acople.
empty_list <- vector(mode = "list", length = rep)
for (i in 1:rep) {
  
  obj <- simulation_hc_v1(N=Nn, mu=0, sj=1)
  m2 <- obj$merge2
  m5 <- obj$merge5
  
  # Generacion de las diferencias en distancias entre el signle linkage y single linkage modificado
  difs <- matrix(NA, ncol=4, nrow=NN)
  colnames(difs) <- c("merge", "Dc", "Dmst", "Ddultr")
  # m5$dultr - m5$dultr > 0 por lo general
  # m5$dmst - m2$dmst > 0 por lo general
  # m5$dc - m2$dc < 0 por lo general
  difs[, ] <-  cbind(c(1:NN), m5$dultr - m5$dultr, m5$dmst - m2$dmst, m5$dc - m2$dc)
  #analysis <- rbind(analysis, difs)
  empty_list[[i]] <- difs
  #i=i+1
}











# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 04-nov-19
# # para un determinado numero de nodos N, simulamos varias veces (rep=500)
# una matriz de acople y calculamos la media de dc, dmst desde merge
# para hacer la pdf de esas medias y compararlas entre el algoritmo normal
# y modificado.

# procedimiento:
# SImular matriz de acople J con <J>=0 para N numero de nodos, y calcular la matriz de distancia
# Obtenemos marge2 y merge5
# calculamos <dc> y <dmst> de merge 2 y merge5
# repetimos el proceso muchas veces (rep veces), y tendremos rep medias de dc y dmst para un N.
# Luego para todos los N, hacemos una grafica de pdf de esas medias.

# 
# definition of number of nodes for each simulation:
N <- c(10, 25, 50, 100)
# Parametros
Jmean = 0 # media de la matriz de acoples
rep = 50
contenedor <- data.frame()
for (i in 1:length(N) ) {
  vect_c1 <- vector(mode="numeric", length=rep)
  vect_c2 <- vector(mode="numeric", length=rep)
  for (j in 1:rep) {
    obj <- simulation_hc_v1(N=N[i], mu=0, sj=1)
    m2 <- obj$merge2
    m5 <- obj$merge5
    
    # L134
    medias2 <- fit_mean(m2)
    medias5 <- fit_mean(m5)
    
    vect_c1[j] <- medias2[3]
    vect_c2[j] <- medias5[3]
    

  }
  alg <- rep(c("Normal", "Modified"), rep)
  nodes <- rep(N[i], 2*rep)
  df <-  data.frame(distances = c(vect_c1, vect_c2), alg = alg, nodes = nodes)
  
  contenedor <- rbind(contenedor, df)
  

}
rm(list=setdiff(ls(), c("contenedor")))
save.image(file='results_041119_from_simulation_hv_V2.RData')


df <- subset(contenedor, nodes == 100)
library(plyr)
mu <- ddply(df, "alg", summarise, grp.mean=mean(distances))
head(mu)

p <- ggplot(df, aes(x=distances, color=alg)) +
  geom_density() +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=alg),
             linetype="dashed") +
  labs(title="Weight density curve", x="Weight(kg)", y = "Density")
p

#violin plot
dp <- ggplot(df, aes(x=alg, y=distances, fill=alg)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="Plot of Coupling Distance", x="Algorithm", y = "Dc")
dp + theme_classic()
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 









# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 05-nov-19
# Simulamos para distinto numero de nodos, una matriz J cualquiera y obtenemos el clustering
# jerarquico normal y modificado (merge2 y merge5 respectivamentr)
# Luego grficamos la distancia de acople para cada iteracion de cada algortimo
# para chequear que la distancia de acople con el algoritmo modificado es siempre menor

# funcion para ticks marks en ggplot2
interleave <- function(x,y){
  lx <- length(x)
  ly <- length(y)
  n <- max(lx,ly)
  as.vector(rbind(rep(x, length.out=n), rep(y, length.out=n)))
}


# procedimiento:
# Para cada valor de numero de nodos N do:
#     Calcular J y las distancia
#     carcular merge 2 y merge 5
#     Plotear dc como scatterplot con linea de distinto color.


# definition of number of nodes for each simulation:
N <- c(10, 25, 50, 100)
# Parametros
Jmean = 0 # media de la matriz de acoples
listm2 <- vector(mode = "list", length = length(N)) # guardamos los merges de V2
listm5 <- vector(mode = "list", length = length(N)) # guardamos los merges de V5
listmG <- vector(mode = "list", length = length(N)) # guardamos los merges de V5
for (i in 1:length(N)) {
  obj <- simulation_hc_v1(Nn=N[i], mu=Jmean, sj=1)
  m2 <- obj$merge2
  m5 <- obj$merge5
  mG <- obj$mergeG
  listm2[[i]] <- m2
  listm5[[i]] <- m5
  listmG[[i]] <- mG
}
rm(list=setdiff(ls(), c("listm2", "listm5", "listmG")))
save.image(file='results_121119_from_simulation_hv_V2.RData')
load('results_121119_from_simulation_hv_V2.RData')


i=4
variable <- "dc"
merge2 <- listm2[[i]]
merge5 <- listm5[[i]]
mergeG <- listmG[[i]]
merge2$alg <- rep("SL", nrow(merge2))
merge5$alg <- rep("ModSL", nrow(merge5))
mergeG$alg <- rep("Greedy", nrow(mergeG))
merge <- rbind(subset(merge2, select = c("alg", "cluster",  variable)), 
               subset(merge5, select = c("alg", "cluster",  variable)),
               subset(mergeG, select = c("alg", "cluster",  variable)))

# https://ggplot2.tidyverse.org/reference/ggtheme.html
# https://stackoverflow.com/questions/46803260/assigning-40-shapes-or-more-in-scale-shape-manual/46814238
my_breaks <- seq(1,nrow(merge2),by=1)
my_labs <- interleave(seq(1,nrow(merge2) ,by=2), "")
dcplot <- ggplot(merge, aes( x = cluster, y = dc, group = alg ) )  + 
            geom_point( aes(shape = alg) )  +
            #geom_line(linetype = "dotted" ) +
            scale_shape_manual(values = c(3, 20, 5)) + 
            #scale_x_continuous(breaks = seq(1:nrow(merge2), by=1), labels=my_labs ) +
            scale_x_continuous(breaks = c(1:nrow(merge2) ) ) +
            xlab('Iteration') + ylab('Coupling Distance') + theme_light()
dcplot

variable <- "dultr"
merge2 <- listm2[[i]]
merge5 <- listm5[[i]]
mergeG <- listmG[[i]]
merge2$alg <- rep("SL", nrow(merge2))
merge5$alg <- rep("ModSL", nrow(merge5))
mergeG$alg <- rep("Greedy", nrow(mergeG))
merge <- rbind(subset(merge2, select = c("alg", "cluster",  variable)), 
               subset(merge5, select = c("alg", "cluster",  variable)),
               subset(mergeG, select = c("alg", "cluster",  variable)))

duplot <- ggplot(merge, aes( x = cluster, y = dultr, group = alg ) )  + 
  geom_point( aes(shape = alg) )  +
  #geom_line(linetype = "dotted" ) +
  scale_shape_manual(values = c(3, 20, 5)) + 
  scale_x_continuous(breaks = c(1:nrow(merge2) ) ) +
  xlab('Iteration') + ylab('Ultrametric Distance') + theme_light()
multiplot(dcplot, duplot, cols=2)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 05-nov-19
# Simulacion de  merge3 y merge5 simuladas para ver los dendogramas basados en 
# la distancia de acople, y no del MST.

#from simulacion_hc_V1.R:
Jmean = 0 
obj <- simulation_hc_v1(Nn=30, mu=Jmean, sj=1)
merge2 <- obj$merge2
merge5 <- obj$merge5
D <- obj$D
hc2 <- obj$hc2
hc5 <- obj$hc5

hc5b <- to_dendo2(merge5[, c(1,2,3,6) ], enames=c(1:ncol(D)))
par(mfrow=c(1,2))
plot(hc2)
plot(hc5b)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 



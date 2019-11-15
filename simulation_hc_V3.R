# En este script vamos a simular clusterig jerarquicos de 
# distintos tamanos de nodos de N= 20, 30, 50, 100, 250, y 500
# con matrices de acoples con <J>=0.
# Luego haremos una grafica scatterplot en que graficamos en el 
# eje X la distancia de acople del algoritmo normal, y en el eje Y
# la distancia de acople del algoritmo modificado para cada una
# de las iteraciones.

# Procedimiento:
# SImular matriz de acople J con <J>=0 para N numero de nodos, y calcular la matriz de distancia
# Obtenemos marge2 y merge5
# graficamos el scatterplot 


# actual name: simulation_hc_V3.R
# creation: 01.nov.19


# Notes:
# 01-nov-19: creation
# 12-nov-19: he modificado el analisis para incluir el algoritmo greedy


# # # # # # # # # # # # # # # # # # LOADING DATA AND FUNCTIONS # # # # # # # # # # # # # # # # # # # # # # #
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
source("get_nodes_of_the_cluster_function.R") # nos da los nodos involucrados en un cluster.
source("simulation_hc_function.R") # Genera una simulacion de hierarchical clustering con algoritmo single linkage y modificado.
source("find_min_distcp_function.R") # funcion que encuentra min dist acople en greedy alg.

source("find_mst_barrier_function.R")
source("is_integer0_function.R")
source("entropies_functions.R") # nos interesa aqui get_energy function
source("multiplot_function.R")

library(igraph)
library(Matrix)
library(ggplot2)
# # # # # # # # # # # # # # # # # # LOADING DATA AND FUNCTIONS # # # # # # # # # # # # # # # # # # # # # # #

# definition of number of nodes for each simulation:
N <- c(10, 25, 50, 100)
# Parametros
Jmean = 0 # media de la matriz de acoples
listas <- vector(mode = "list", length = length(N)) # guardamos los merges de V2
for (i in 1:length(N)) {
  obj <- simulation_hc_v1(Nn = N[i], mu = Jmean, sj = 1)
  m2 <- obj$merge2
  m5 <- obj$merge5
  mG <- obj$mergeG
  
  x <- m2$dc # distancia de acople obtenido en cada fusion del algoritmo single linkage.
  y5 <- m5$dc # distancia de acople obtenido en cada fusion del algoritmo single linkage modificado.
  yG <- mG$dc # distancia de acople obtenido en cada fusion del algoritmo single linkage modificado.
  
  nodes <- rep((length(x)+1), 2*(N[i]-1) )
  alg <- rep(c("ModSL", "Greedy"), (N[i]-1) )
  df <- data.frame(x = rep(x, 2), dc = c(y5, yG), alg = alg, nodes = nodes, it=1:length(x))
  listas[[i]] <- df
}

data <- listas[[1]]
con10 <- ggplot(data, aes(x = x, y = dc, group = alg) ) +
  #geom_point(shape=16, alpha=0.8) +    # 1/4 opacity 
  geom_point( aes(shape = alg) )  +
  #geom_line(linetype = "dotted" ) +
  scale_shape_manual(values = c(3, 20)) +
  geom_abline(intercept = 0, slope=1, colour="red",size=0.5, linetype="dotted") + 
  geom_text(aes(label = it), size = 3, hjust=-0.5, vjust=0.3) + 
  xlab('Coupling Distance SL') + ylab('Coupling Distance of Variations') + theme_light()

data <- listas[[2]]
con25 <- ggplot(data, aes(x = x, y = dc, group = alg) ) +
  #geom_point(shape=16, alpha=0.8) +    # 1/4 opacity 
  geom_point( aes(shape = alg) )  +
  #geom_line(linetype = "dotted" ) +
  scale_shape_manual(values = c(3, 20)) +
  geom_abline(intercept = 0, slope=1, colour="red",size=0.5, linetype="dotted") + 
  geom_text(aes(label = it), size = 3, hjust=-0.5, vjust=0.3) + 
  xlab('Coupling Distance SL') + ylab('Coupling Distance of Variations') + theme_light()


data <- listas[[3]]
con50 <- ggplot(data, aes(x = x, y = dc, group = alg) ) +
  #geom_point(shape=16, alpha=0.8) +    # 1/4 opacity 
  geom_point( aes(shape = alg) )  +
  #geom_line(linetype = "dotted" ) +
  scale_shape_manual(values = c(3, 20)) +
  geom_abline(intercept = 0, slope=1, colour="red",size=0.5, linetype="dotted") + 
  geom_text(aes(label = it), size = 3, hjust=-0.5, vjust=0.3) + 
  xlab('Coupling Distance SL') + ylab('Coupling Distance of Variations') + theme_light()


data <- listas[[4]]
con100 <- ggplot(data, aes(x = x, y = dc, group = alg) ) +
  #geom_point(shape=16, alpha=0.8) +    # 1/4 opacity 
  geom_point( aes(shape = alg) )  +
  #geom_line(linetype = "dotted" ) +
  scale_shape_manual(values = c(3, 20)) +
  geom_abline(intercept = 0, slope=1, colour="red",size=0.5, linetype="dotted") + 
  geom_text(aes(label = it), size = 3, hjust=-0.5, vjust=0.3) + 
  xlab('Coupling Distance SL') + ylab('Coupling Distance of Variations') + theme_light()
multiplot(con10,  con50, con25, con100, cols=2)
# 12-nov-19
rm(list=setdiff(ls(), c("listas")))
save.image(file='results_121119_from_simulation_hc_V3.RData')




# Esto ya esta desactualizado!!!!!!!!!!!!!!!!!!!  se debe actualizar!!!!!!
# https://stackoverflow.com/questions/31993704/storing-ggplot-objects-in-a-list-from-within-loop-in-r
# funcion para generar graficas de una vez en cada elemento de la lista.
plot_data_from_lists = function (data, title="") {
  ggplot(data, aes(x = x, y = y) ) +
    geom_point(shape=16, alpha=0.8) +    # 1/4 opacity 
    geom_abline(intercept = 0, slope=1, colour="red",size=0.5, linetype="dotted") + 
    geom_text(aes(label = it), size = 3, hjust=-0.5, vjust=0.3) + 
    xlab('Coupling Distance SL') + ylab('Coupling Distance modified SL') +
    ggtitle(title)
}

myplots <- lapply(listas, plot_data_from_lists, title = "")

rm(list=setdiff(ls(), c("listas", "myplots")))
save.image(file='results_031119_from_simulation_hv_V3.RData')

load(file='results_031119_from_simulation_hv_V3.RData')

source("multiplot_function.R")
library(ggplot2)
multiplot(myplots[[1]], myplots[[2]], myplots[[3]], myplots[[4]], myplots[[6]], myplots[[6]], cols=2)

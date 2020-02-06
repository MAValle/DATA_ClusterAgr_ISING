# En este script simulamos varias matrices de acople
# N(0,1) y con N=20 para graficar las distancias de acople
# en cada iteraciones en cada uno de los algoritmos de prueba.

# Prcedure:
# 1 create a coupling matrix of size N
# 2 get the merge dataframe with coupling distances of each iteration for:
#     * single linkage
#     * modularity method
#     * deterministic greedy
#     * probabilistic greedy with gamma = 0.5
#     * probabilistic greedy with gamma = 0 (100% probabilistic)
# 3 save coupling distances
# 4 back to 1  I=1000 times
# 5 get the means of coupling distances for each iteration for each algorthm
# 6 ggplot of means of each iteration for each algorithm with confidence.



# actual name: simulation_hc_V6.R
# creation: 16.dic.19


# Notes:
# 16-dic-19: creation


  
  
  
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
source("find_min_distcp_function.R") # funcion que encuentra min dist acople en greedy alg.
source("pick_a_cluster_function.R")

# funciones para recuperar distancias de acople para los clusters
# que vienen de fast_greedy para calcular clustering con modularidad.
source("recovering_couplingdistances_from_fastgreedyigraph_function.R")

# Funciones para hacer post analysis con las simulaciones de HC.
source("post_hc_analysis_function.R")

source("find_mst_barrier_function.R")
source("is_integer0_function.R")
source("entropies_functions.R") # nos interesa aqui get_energy function
source("multiplot_function.R")

library(igraph)
library(Matrix)
library(ggplot2)
# # # # # # # # # # # # # # # # # # LOADING DATA AND FUNCTIONS # # # # # # # # # # # # # # # # # # # # # # #


# # # # # # # # # # # # # # # # # # # # # # # SIMULACION 1 # # # # # # # # # # # # # # # # # # # # # # # #
# dic 17, 2019
# Using J distributed with N(0,1)
set.seed(123)
#N <- c(10, 25, 50, 100)
N <- 50
thres <- c(1, 0.5, 0) # 1 greedy deterministico, 0 = greedy random
Jmean = 1.5
sj = 1

# Begining the simulations
I = 10
ware <- list()
ware[["singlelink"]] = list()
ware[["modularity"]] = list()
ware[["greedy_detm"]] = list()
ware[["greedy_semi"]] = list()
ware[["greedy_prob"]] = list()
for (i in 1:I){
  # CREATING COUPLINGS
  ot <- create_coupling(Nn = N, media = Jmean, sj=sj, lw=-3, up=3, kind="normal") # para J normal
  #ot <- create_coupling(Nn = N, media = Jmean, sj=sj, lw=-3, up=3, kind="uniform") # para J uniforme
  J <- ot$J
  D <- ot$D
  #J
  # # # # # # SINGLE LINKAGE
  # Conforming the COUPLING Network
  mst_g <- create_mst(J=J, D=D)
  # dic, 18, 2019
  invisible(capture.output( merge_singlelink <- hierarchical_clustering_v2(D=D, J=J, mst_g=mst_g) ) )
  # original
  #merge_singlelink <- hierarchical_clustering_v2(D=D, J=J, mst_g=mst_g)
  #hc2 <- to_dendo(D, merge_singlelink[,c(1:4)], enames=c(1:ncol(D)) )
  
  # # # # # # MODULARITY
  g <- graph_from_adjacency_matrix(D, weighted=TRUE, mode="undirected", diag=FALSE)
  fc <- cluster_fast_greedy(g)
  # is_hierarchical(fc) # TRUE
  # plot_dendrogram(fc)
  result <- get_each_node_of_each_cluster(merge = fc$merges)
  nodes_of_clusters <- result$nodes_of_clusters
  merge <- result$merge
  merge_modularity <- report_merge_from_modularity(merge = result$merge, node_of_clusters = node_of_clusters)
  
  # # # # # # GREEDY
  # dic, 18, 2019
  invisible(capture.output( merge_greedy_detm <- hierarchical_clustering_probabilistic_greedy(D = D, J = J, trh= thres[1] ) ) )
  invisible(capture.output( merge_greedy_semi <- hierarchical_clustering_probabilistic_greedy(D = D, J = J, trh= thres[2] ) ) )
  invisible(capture.output( merge_greedy_prob <- hierarchical_clustering_probabilistic_greedy(D = D, J = J, trh= thres[3] ) ) )
  # original
  # merge_greedy_detm <- hierarchical_clustering_probabilistic_greedy(D = D, J = J, trh= thres[1] )
  # merge_greedy_semi <- hierarchical_clustering_probabilistic_greedy(D = D, J = J, trh= thres[2] )
  # merge_greedy_prob <- hierarchical_clustering_probabilistic_greedy(D = D, J = J, trh= thres[3] )
  
  # # # # # # SAVING merges IN LISTS
  ware[["singlelink"]][[i]] = merge_singlelink[, c(1,2,6)]
  ware[["modularity"]][[i]] = merge_modularity[, c(1,2,4)]
  ware[["greedy_detm"]][[i]] = merge_greedy_detm[, c(1,2,6)]
  ware[["greedy_semi"]][[i]] = merge_greedy_semi[, c(1,2,6)]
  ware[["greedy_prob"]][[i]] = merge_greedy_prob[, c(1,2,6)]
  
  #
  cat("\n Iteration number :", i) 
  #flush.console()
}

# getting the means and creating the final dataframe to plot
# https://stackoverflow.com/questions/31839996/mean-of-each-row-of-data-frames-column-nested-in-list
SL <- do.call(cbind, lapply(ware[["singlelink"]], "[", "dc"))
dc_mean_SL <- rowMeans(SL)
dc_sdev_SL <- apply(SL,1,sd)
MOD <- do.call(cbind, lapply(ware[["modularity"]], "[", "dc"))
dc_mean_MOD <- rowMeans(MOD)
dc_sdev_MOD <- apply(MOD,1,sd)
Gdet <- do.call(cbind, lapply(ware[["greedy_detm"]], "[", "dc"))
dc_mean_Gdet <- rowMeans(Gdet)
dc_sdev_Gdet <- apply(Gdet,1,sd)
Gsemi <- do.call(cbind, lapply(ware[["greedy_semi"]], "[", "dc"))
dc_mean_Gsemi <- rowMeans(Gsemi)
dc_sdev_Gsemi <- apply(Gsemi,1,sd)
Gprob <- do.call(cbind, lapply(ware[["greedy_prob"]], "[", "dc"))
dc_mean_Gprob <- rowMeans(Gprob)
dc_sdev_Gprob <- apply(Gprob,1,sd)
# 
# iteration <- rep(seq(from = 1, to = N-1),5)
# algo <- c(rep("SL", N-1), rep("MOD", N-1), rep("Gdet", N-1), 
#           rep("Gsemi", N-1), rep("Gprob", N-1))
# medias <- c(dc_mean_SL, dc_mean_MOD, dc_mean_Gdet, dc_mean_Gsemi, dc_mean_Gprob)
# desvs <- c(dc_sdev_SL, dc_sdev_MOD, dc_sdev_Gdet, dc_sdev_Gsemi, dc_sdev_Gprob) 

iteration <- rep(seq(from = 1, to = N-1),3)
algo <- c(rep("SL", N-1), rep("MOD", N-1), rep("Gdet", N-1))
medias <- c(dc_mean_SL, dc_mean_MOD, dc_mean_Gdet)
desvs <- c(dc_sdev_SL, dc_sdev_MOD, dc_sdev_Gdet) 

df <- data.frame(iteration = iteration, Algorithm = algo, 
                      mean = medias, sd = desvs)

# doing the ggplot
f <- ggplot(df, aes(x = iteration, y = log(mean), colour = Algorithm)) + 
  geom_point(shape=20, size=1.5) +
  geom_smooth(show.legend = FALSE, fill = "gray", span = 0.3) +
  #geom_errorbar(aes(ymin=log(mean-sd/sqrt(I)), ymax=log(mean+sd/sqrt(I))), width=.2, position=position_dodge(.01)) +
  #geom_errorbar(aes(ymin=log(mean-sd/sqrt(1)), ymax=log(mean+sd/sqrt(1))), width=.2, position=position_dodge(.01)) +
  labs(x ="Iteration", y = "mean of log(U*)") +
  theme_linedraw() + theme(legend.position="bottom") +
  scale_x_continuous(breaks = seq(from=1, to=N-1, by=2) ) +
  scale_y_continuous(breaks = c(-1, 0, 1, 2, 3,4,5,6,7,8)) +
  #scale_y_continuous(breaks = seq(from=1, to=8) ) +
  theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "white"))
cairo_ps("figure.eps")
f + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
ggsave("figure2eps_10samplings.eps")
save(f5, file = "figure5_from_simulation_hc_V6.Rdata")
# # # # # # # # # # # # # # # # # # # # # # # SIMULACION 1 # # # # # # # # # # # # # # # # # # # # # # # #




# # # # # # # # # # # # # # # # # # # # # # # SIMULACION 2 # # # # # # # # # # # # # # # # # # # # # # # #
# dic 18, 2019
# # Using J distributed with U(-2.5, 2.5)


# # # # # # # # # # # # # # # # # # # # # # # SIMULACION 2 # # # # # # # # # # # # # # # # # # # # # # # #

# add confidence interval to splines from quantile regression
# https://stackoverflow.com/questions/47464691/add-confidence-interval-to-splines-from-quantile-regression
# https://stackoverflow.com/questions/34742555/adding-uncertainty-bands-to-a-smooth-spline-in-a-scatterplot




# https://stackoverflow.com/questions/5142842/export-a-graph-to-eps-file-with-r


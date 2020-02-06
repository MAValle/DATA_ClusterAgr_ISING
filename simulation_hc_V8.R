# En este script analizamos el ejemplo real de la bases de 
# datos transaccional que he utilizado en los papers de ICANN 
# con 25 nodos, para calcular los dendogramas con MST y con greedy. 
# Adicionalmente, vemos las simulaciones con greedy probabilistico.

# este script va de la mano con figure_netcoupling_for_IEEE.R utilizado
# para hacer la red de acoples y MST con los datos reales de 25 nodos.



# actual name: simulation_hc_V8.R
# creation: 10.ene.20


# Notes:
# 10-ene-20: creation
# # jan 18, 2020:  script para simular varios greedy probabilisticos con gamma=0, 0.5 y 0.75 y calculamos el cophonetic correlations



# # # # # # # # # # # # # # # # # # LOADING DATA AND FUNCTIONS # # # # # # # # # # # # # # # # # # # # # # #
# Carga de datos, 
# igual que en clusteragrr_v1.R
rm(list = ls())
# LOADING DATA BASE REAL DATA
load('new_inferring_parameters_environment270219.RData')
rm(list=setdiff(ls(), c("J", "wb" ) ))

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
source("cophonetic_cor_function.R") # para calcular las correlaciones cofoneticas

source("entropies_functions.R")
source("network_functions_v1.R")

library(igraph)
library(Matrix)
library(ggplot2)

# Conforming the DISTANCES Network
real_names <- colnames(wb)
D <- sqrt(3 - as.matrix(J)) # distance
# Por ahora, tenemos que enumerar los productos en forma consecutiva
# despues, para el dendograma, le ponemos los nombres reales
colnames(D) <- rownames(D) <- -c(1:length(colnames(wb)) ) # enumeramos todos las hojas o spines originales con valones negativos.
colnames(J) <- rownames(J) <- -c(1:length(colnames(wb)) ) # enumeramos todos las hojas o spines originales con valones negativos.
diag(D) <- rep(Inf, length(colnames(J)))
# # # # # # # # # # # # # # # # # # LOADING DATA AND FUNCTIONS # # # # # # # # # # # # # # # # # # # # # # #




# # # # # # # # # # # # # # # # # # # # # # # SIMULACION # # # # # # # # # # # # # # # # # # # # # # # # #
# jan 10 , 2020
# tomando como ejemplo L203 de simulation_hc_V7.R

# SINGLE LINKAGE
mst_g <- create_mst(J=J, D=D)
invisible(capture.output( merge_singlelink <- hierarchical_clustering_v2(D=D, J=J, mst_g=mst_g) ) )

# GREEDY DETERMINISTIC
thres <- c(1, 0.5, 0) # 1 greedy deterministico, 0 = greedy random
invisible(capture.output( merge_greedy_detm <- hierarchical_clustering_probabilistic_greedy(D = D, J = J, trh= thres[1] ) ) )

# MODULARITY
N <- length(colnames(J))
g <- graph_from_adjacency_matrix(D, weighted=TRUE, mode="undirected", diag=FALSE)
fc <- cluster_fast_greedy(g)
# is_hierarchical(fc) # TRUE
# plot_dendrogram(fc)
result <- get_each_node_of_each_cluster(merge = fc$merges)
nodes_of_clusters <- result$nodes_of_clusters
merge <- result$merge
merge_modularity <- report_merge_from_modularity(merge = result$merge, node_of_clusters = node_of_clusters)
merge_modularity$cluster <- seq(from=1, to=N-1)


# formacion de data frame para ggplot
iteration <- rep(seq(from = 1, to = N-1), 3)
algo <- c(rep("MST", N-1), rep("Modularity", N-1), rep("Algorithm 1", N-1))
dc <- c(merge_singlelink$dc, merge_modularity$dc, merge_greedy_detm$dc)
df <- data.frame(Iteration = iteration, Algorithm = algo, CouplingDistance = dc)


# # # # PREPARING FOR figure_paths_energies_for_IEEE_eps.eps
# doing the ggplot
f1 <- ggplot(df, aes(x = iteration, y = log(CouplingDistance), colour = Algorithm)) + 
  geom_point(shape=20, size=1.5) + geom_line() + 
  #geom_smooth(show.legend = FALSE, fill = "gray", span = 0.3) +
  #geom_errorbar(aes(ymin=log(mean-sd/sqrt(I)), ymax=log(mean+sd/sqrt(I))), width=.2, position=position_dodge(.01)) +
  #geom_errorbar(aes(ymin=log(mean-sd/sqrt(1)), ymax=log(mean+sd/sqrt(1))), width=.2, position=position_dodge(.01)) +
  labs(x ="Iteration", y = "log of Coupling Distance") +
  theme_linedraw() + theme(legend.position="bottom") +
  scale_x_continuous(breaks = seq(from=1, to=N-1, by=2) ) +
  #scale_y_continuous(breaks = c(-1, 0, 1, 2, 3,4,5,6,7,8)) +
  #scale_y_continuous(breaks = seq(from=1, to=8) ) +
  theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "white"))
f1 <- f1 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#cairo_ps("figure.eps")
#f + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#dev.off()
#ggsave("figure2eps_10samplings.eps")
#save(f5, file = "figure5_from_simulation_hc_V6.Rdata")


# # # # # PREPARING FOR figure_dendograms_for_IEEE_eps.eps
# de L247 de simulation_hc_V7.R
# graficas de los dendogramas
hc_sl <- to_dendo2(merge_singlelink[, c(1,2,3,6) ], enames = c(1:ncol(D)))
hc_sl$labels <- real_names
hc_det <- to_dendo2(merge_greedy_detm[, c(1,2,3,6) ], enames = c(1:ncol(D)))
hc_det$labels <- real_names
hc_mod <- to_dendo2(merge_modularity[, c(1,2,3,4) ], enames = c(1:ncol(D)))
hc_mod$labels <- real_names
par(mfrow=c(3,1))    # set the plotting area into a 1*2 array
plot(hc_det, main="Greedy Algorithm 1 Dendrogram", cex.main=1, ylab="Coupling Distance")
plot(hc_sl, main="MST Dendrogram", cex.main=1, ylab="Coupling Distance")
plot(hc_mod, main="Greedy Modularity Optimization Dendrogram", cex.main=1, ylab="Coupling Distance")
library(cluster)
# https://stat.ethz.ch/R-manual/R-patched/library/cluster/html/coef.hclust.html
# https://bradleyboehmke.github.io/HOML/hierarchical.html
# agglomerative coefficient (AC), which measures the amount of clustering structure found.
# Generally speaking, the AC describes the strength of the clustering structure. Values closer to 1 
# suggest a more balanced clustering structure such as the complete linkage and Ward’s method dendrograms 
# in Figure 21.3. Values closer to 0 suggest less well-formed clusters such as the single linkage
# dendrogram in Figure 21.3. However, the AC tends to become larger as   n
# increases, so it should not be used to compare across data sets of very different sizes.
# https://www.umass.edu/landeco/teaching/multivariate/schedule/cluster1.pdf
# Agglomerative coeeficients:
# For each observation i, denote
# by m(i) its dissimilarity to the
# first cluster it is merged with,
# divided by the dissimilarity of
# the merger in the final step of
# the algorithm. The AC is the
# average of all 1 - m(i).

mean(1-merge_singlelink$dc/merge_singlelink$dc[nrow(merge_singlelink)]) # 0.8186595
mean(1-merge_greedy_detm$dc/merge_greedy_detm$dc[nrow(merge_greedy_detm)]) # 0.9222798
mean(1-merge_modularity$dc/merge_modularity$dc[nrow(merge_modularity)]) # 0.6481349









# jan 14, 2020
# # # # # # # # # # # # # # # # # # # # # # # SIMULACION # # # # # # # # # # # # # # # # # # # # # # # # #
# # SIMULATING PROBABILISTIC GREEDY WITH REAL DATA BASE
# making figure realcase_paths_from_simulation_hc_V8_150120.pdf

# execute lines L22 a L73
thres <- c(1, 0.75, 0.5, 0) # 1 greedy deterministico, 0 = greedy random
N <- length(colnames(J))

# # # # # # GREEDY
# dic, 18, 2019
invisible(capture.output( merge_greedy_detm <- hierarchical_clustering_probabilistic_greedy(D = D, J = J, trh= thres[1] ) ) )
# Begining the simulations for probabilistic greedy
I = 300
ware <- list()
ware[["greedy_semi05"]] = list()
ware[["greedy_semi075"]] = list()
ware[["greedy_prob"]] = list()
for (i in 1:I){
  
  invisible(capture.output( merge_greedy_semi075 <- hierarchical_clustering_probabilistic_greedy(D = D, J = J, trh= thres[2] ) ) )
  invisible(capture.output( merge_greedy_prob <- hierarchical_clustering_probabilistic_greedy(D = D, J = J, trh= thres[4] ) ) )
  invisible(capture.output( merge_greedy_semi05 <- hierarchical_clustering_probabilistic_greedy(D = D, J = J, trh= thres[3] ) ) )
  
  # # # # # # SAVING merges IN LISTS
  ware[["greedy_semi075"]][[i]] = merge_greedy_semi075[, c(1,2,6)]
  ware[["greedy_semi05"]][[i]] = merge_greedy_semi05[, c(1,2,6)]
  ware[["greedy_prob"]][[i]] = merge_greedy_prob[, c(1,2,6)]

  
  #
  cat("\n Iteration number :", i) 
  #flush.console()
}

# getting the means and creating the final dataframe to plot
# https://stackoverflow.com/questions/31839996/mean-of-each-row-of-data-frames-column-nested-in-list
Gsemi05 <- do.call(cbind, lapply(ware[["greedy_semi05"]], "[", "dc"))
dc_mean_Gsemi05 <- rowMeans(Gsemi05)
dc_sdev_Gsemi05 <- apply(Gsemi05,1,sd)

Gsemi075 <- do.call(cbind, lapply(ware[["greedy_semi075"]], "[", "dc"))
dc_mean_Gsemi075 <- rowMeans(Gsemi075)
dc_sdev_Gsemi075 <- apply(Gsemi075,1,sd)

Gprob <- do.call(cbind, lapply(ware[["greedy_prob"]], "[", "dc"))
dc_mean_Gprob <- rowMeans(Gprob)
dc_sdev_Gprob <- apply(Gprob,1,sd)
# 

# Calculating mean and sd of clustering coefficients
lista <- lapply(ware[["greedy_semi05"]], "[", "dc")
clustcoef_greedy_semi05 <-  (do.call(cbind,  lapply( lista, FUN=function(x)  mean(1-x$dc/x$dc[nrow(x)] ) ) ) )
clustcoef_sd_greedy_semi05 <- sd(do.call(cbind,  lapply( lista, FUN=function(x)  mean(1-x$dc/x$dc[nrow(x)] ) ) ) )

lista <- lapply(ware[["greedy_semi075"]], "[", "dc")
clustcoef_greedy_semi075 <-  (do.call(cbind,  lapply( lista, FUN=function(x)  mean(1-x$dc/x$dc[nrow(x)] ) ) ) )
clustcoef_sd_greedy_semi075 <- sd(do.call(cbind,  lapply( lista, FUN=function(x)  mean(1-x$dc/x$dc[nrow(x)] ) ) ) )

lista <- lapply(ware[["greedy_prob"]], "[", "dc")
clustcoef_greedy_prob <- (do.call(cbind,  lapply( lista, FUN=function(x)  mean(1-x$dc/x$dc[nrow(x)] ) ) ) )
clustcoef_sd_greedy_prob <- sd (do.call(cbind,  lapply( lista, FUN=function(x)  mean(1-x$dc/x$dc[nrow(x)] ) ) ) )
clustcoef_mean_greedy_det <- mean(1-merge_greedy_detm $dc/merge_greedy_detm $dc[nrow(merge_greedy_detm )])
# preparing dataframe to plot boxplot for cluster coefficients
algo <- c(rep("Gsemi0.5", I), rep("Gsemi0.75", I), rep("Gprob", I))
clust_coefs_values <- c(clustcoef_greedy_semi05,  clustcoef_greedy_semi075, clustcoef_greedy_prob )
df_coef <- data.frame(Algorithm = algo, ClusteringCoefficients = clust_coefs_values)
# BOXPLOT figure realcase_boxplot_from_simulation_hc_V8_150120.pdf boxplot of cluster coefficients
boxplot((df_coef$ClusteringCoefficients) ~ df_coef$Algorithm, col=rgb(0.3, 0.5, 0.4), 
        ylab="Clustering Coefficient", xlab="Threshold",
        ylim=c(0.8, 0.95)) 
        #yaxs="i", 
        #main="Clustering Coefficient means")
abline(h=clustcoef_mean_greedy_det, col="red")

   
# Preparing dataframe to ggplot realcase_paths_from_simulation_hc_V8_150120.pdf
iteration <- rep(seq(from = 1, to = N-1),4)
algo <- c(rep("Gdet", N-1), rep("Gsemi0.5", N-1), rep("Gsemi0.75", N-1), rep("Gprob", N-1))
medias <- c(merge_greedy_detm$dc, dc_mean_Gsemi05, dc_mean_Gsemi075, dc_mean_Gprob)
desvs <- c(rep(NA, N-1), dc_sdev_Gsemi05, dc_sdev_Gsemi075, dc_sdev_Gprob) 

df <- data.frame(iteration = iteration, Algorithm = algo, 
                 mean = medias, sd = desvs)
df$min <- ifelse(df$mean - df$sd <=0, 0, log(df$mean - df$sd) )
df$max <- ifelse(df$mean + df$sd <=0, 0, log(df$mean + df$sd) )


# doing the ggplot
#los boxplot vienen de simulation_hc_v7 L344 en adelante
# https://stats.stackexchange.com/questions/103393/how-can-i-dodge-the-position-of-geom-point-in-ggplot2
# https://stackoverflow.com/questions/33113322/dodging-intervals-in-ggplot
f <- ggplot(df, aes(x = iteration, y = log(mean), colour = Algorithm)) + 
  geom_point(shape=20, size=2, na.rm=TRUE, position=position_dodge(width=0.7) )  + 
  #geom_smooth(show.legend = FALSE, fill = "gray", span = 0.3) +
  #geom_errorbar(aes(ymin=min, ymax=max ), width=.2, position=position_dodge(.01)) +
  geom_errorbar(aes(x=iteration, ymin=log(mean-sd/sqrt(3)), ymax=log(mean+sd/sqrt(3))), width=.7, na.rm=TRUE, position=position_dodge(width=0.7)) +  
  labs(x ="Iteration", y = "mean of log(U*)") +
  theme_linedraw() + theme(legend.position="bottom") +
  #scale_x_continuous(breaks = seq(from=1, to=N-1, by=2) ) +
  scale_y_continuous(breaks = c(-1, 0, 1, 2, 3,4,5,6,7,8)) +
  #scale_y_continuous(breaks = seq(from=1, to=8) ) +
  theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "white"))

#cairo_ps("figure.eps")
f + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# # # # # # # # # # # # # # # # # # # # # # # SIMULACION # # # # # # # # # # # # # # # # # # # # # # # # #








# jan 17, 2020
# # # # # # # # # # # # # # # # # # # # # # # SIMULACION # # # # # # # # # # # # # # # # # # # # # # # # #
# # SIMULATING PROBABILISTIC GREEDY WITH REAL DATA BASE
# making figure realcase_dendoprobs_from_simulation_hc_V8_170220_pdf.pdf
load("workdata_170220_from_simulation_hc_V8.RData")
# execute lines L22 a L73
thres <- c(1, 0.75, 0.5, 0) # 1 greedy deterministico, 0 = greedy random
N <- length(colnames(J))
invisible(capture.output( merge_greedy_detm <- hierarchical_clustering_probabilistic_greedy(D = D, J = J, trh= thres[1] ) ) )
invisible(capture.output( merge_greedy_semi075 <- hierarchical_clustering_probabilistic_greedy(D = D, J = J, trh= thres[2] ) ) )
invisible(capture.output( merge_greedy_semi05 <- hierarchical_clustering_probabilistic_greedy(D = D, J = J, trh= thres[3] ) ) )
invisible(capture.output( merge_greedy_prob <- hierarchical_clustering_probabilistic_greedy(D = D, J = J, trh= thres[4] ) ) )
mst_g <- create_mst(J=J, D=D)
# dic, 18, 2019
invisible(capture.output( merge_singlelink <- hierarchical_clustering_v2(D=D, J=J, mst_g=mst_g) ) )
# original

hc_det <- to_dendo2(merge_greedy_detm[, c(1,2,3,6) ], enames = c(1:ncol(D)))
hc_det$labels <- real_names
hc_semi075 <- to_dendo2(merge_greedy_semi075[, c(1,2,3,6) ], enames = c(1:ncol(D)))
hc_semi075$labels <- real_names
hc_semi05 <- to_dendo2(merge_greedy_semi05[, c(1,2,3,6) ], enames = c(1:ncol(D)))
hc_semi05$labels <- real_names
hc_prob <- to_dendo2(merge_greedy_prob[, c(1,2,3,6) ], enames = c(1:ncol(D)))
hc_prob$labels <- real_names
hc_SL <- to_dendo2(merge_singlelink[, c(1,2,3,6) ], enames = c(1:ncol(D)))
hc_SL$labels <- real_names

par(mfrow=c(3,1))    # set the plotting area into a 1*2 array
plot(hc_semi075, main="Greedy Algorithm 2 Dendrogram with gamma=0.75", cex.main=1, ylab="Coupling Distance")
plot(hc_semi05, main="Greedy Algorithm 2 Dendrogram with gamma=0.5", cex.main=1, ylab="Coupling Distance")
plot(hc_prob, main="Greedy Algorithm 2 Dendrogram with gamma=0", cex.main=1, ylab="Coupling Distance")


# Para hacer la table de clusters o canastas tipos.

# https://r.789695.n4.nabble.com/Cutting-hierarchical-cluster-tree-at-specific-height-fails-td4693737.html
cutree.h <- function(tree,h) {
  # this line adapted from cutree(...) code
  k <- nrow(tree$merge) + 2L - apply(outer(c(tree$height, Inf), h, ">"), 2, which.max)
  return(cutree(tree,k=k))
}

members_det <- cutree.h(hc_det, h = 20)
members_075 <- cutree.h(hc_semi075, h = 20)
members_05 <- cutree.h(hc_semi05, h = 20)
members_prob <- cutree.h(hc_prob, h = 20)
members_SL <- cutree.h(hc_SL, h = 20)
# # # # # # # # # # # # # # # # # # # # # # # SIMULACION # # # # # # # # # # # # # # # # # # # # # # # # #


















# jan 18, 2020
# # # # # # # # # # # # # # # # # # # # # # # SIMULACION # # # # # # # # # # # # # # # # # # # # # # # # #
# # SIMULATING PROBABILISTIC GREEDY TO GET THE COPHONETIC CORRELATIONS
# https://people.revoledu.com/kardi/tutorial/Clustering/Cophenetic.htm
# https://stackoverflow.com/questions/5639794/in-r-how-can-i-plot-a-similarity-matrix-like-a-block-graph-after-clustering-d
# https://rdrr.io/cran/dendextend/man/cor_bakers_gamma.html
I = 100

thres <- c(1, 0.75, 0.5, 0) # 1 greedy deterministico, 0 = greedy random
N <- length(colnames(J))
# create the container coph_SL, coph_det, coph_075, coph_05, coph_prob
M <- matrix(NA, ncol=5, nrow=I)
invisible(capture.output( merge_greedy_detm <- hierarchical_clustering_probabilistic_greedy(D = D, J = J, trh= thres[1] ) ) )
mst_g <- create_mst(J=J, D=D)
invisible(capture.output( merge_singlelink <- hierarchical_clustering_v2(D=D, J=J, mst_g=mst_g) ) )
hc_det <- to_dendo2(merge_greedy_detm[, c(1,2,3,6) ], enames = c(1:ncol(D)))
hc_det$labels <- real_names
hc_SL <- to_dendo2(merge_singlelink[, c(1,2,3,6) ], enames = c(1:ncol(D)))
hc_SL$labels <- real_names
coph_det <- cophonetic_cor(D = D, ultr = as.matrix(cophenetic(hc_det)), hc_clust = hc_det)
coph_SL <- cophonetic_cor(D = D, ultr = as.matrix(cophenetic( hc_SL)), hc_clust = hc_SL)
for (i in 1:I) {
  # get the merge matrix 
  invisible(capture.output( merge_greedy_semi075 <- hierarchical_clustering_probabilistic_greedy(D = D, J = J, trh= thres[2] ) ) )
  invisible(capture.output( merge_greedy_semi05 <- hierarchical_clustering_probabilistic_greedy(D = D, J = J, trh= thres[3] ) ) )
  invisible(capture.output( merge_greedy_prob <- hierarchical_clustering_probabilistic_greedy(D = D, J = J, trh= thres[4] ) ) )

  # get the dendograms opbjects
  hc_semi075 <- to_dendo2(merge_greedy_semi075[, c(1,2,3,6) ], enames = c(1:ncol(D)))
  hc_semi075$labels <- real_names
  hc_semi05 <- to_dendo2(merge_greedy_semi05[, c(1,2,3,6) ], enames = c(1:ncol(D)))
  hc_semi05$labels <- real_names
  hc_prob <- to_dendo2(merge_greedy_prob[, c(1,2,3,6) ], enames = c(1:ncol(D)))
  hc_prob$labels <- real_names
 
  #get the cophonetic correlations
  coph_075 <- cophonetic_cor(D = D, ultr = as.matrix(cophenetic( hc_semi075)), hc_clust = hc_semi075)
  coph_05 <- cophonetic_cor(D = D, ultr = as.matrix(cophenetic( hc_semi05)), hc_clust = hc_semi05)
  coph_prob <- cophonetic_cor(D = D, ultr = as.matrix(cophenetic( hc_prob)), hc_clust = hc_prob)
  
  # put results in container
  corrs <- c(coph_SL, coph_det, coph_075, coph_05, coph_prob)
  M[i, ] <- corrs
  
  cat("\n Iteration number :", i) 
}
colMeans(M)

# preparing dataframe to plot boxplot for cophonetics correlations
# nota: hay que ejecutar L293 hacia abajo.
algo <- c(rep("Gsemi0.5", I), rep("Gsemi0.75", I), rep("Gprob", I))
clust_coefs_values <- c(clustcoef_greedy_semi05,  clustcoef_greedy_semi075, clustcoef_greedy_prob )
df_coef <- data.frame(Algorithm = algo, ClusteringCoefficients = clust_coefs_values)
# BOXPLOT figure realcase_boxplot_from_simulation_hc_V8_150120.pdf boxplot of cluster coefficients
boxplot((df_coef$ClusteringCoefficients) ~ df_coef$Algorithm, col="lightgrey", 
        ylab="Clustering Coefficient", 
        #xlab="Threshold",
        ylim=c(0.8, 0.95)) 
#yaxs="i", 
#main="Clustering Coefficient means")
abline(h=clustcoef_mean_greedy_det, col="red")
# # # # # # # # # # # # # # # # # # # # # # # SIMULACION # # # # # # # # # # # # # # # # # # # # # # # # #
















members_det <- (as.data.frame(members_det))
table(members_det$members_det) #6 clusters

members_075 <- as.data.frame(members_075)
table(members_075$members_075) # 8 clusters

members_05 <- as.data.frame(members_05)
table(members_05$members_05) # 9 clusters

members_prob <- as.data.frame(members_prob)
table(members_prob$members_prob) # 13 clusters

members_SL <- as.data.frame(members_SL)
table(members_SL$members_SL) # 19 clusters

#https://stackoverflow.com/questions/3789549/display-a-matrix-including-the-values-as-a-heatmap


heatmap( ultr, Rowv=NA, Colv=NA, col = heat.colors(256),  margins=c(5,10))
heatmap( ultr, Rowv=NA, Colv=NA, col = heat.colors(256) )

image(1:ncol(ultr), 1:nrow(ultr), t(ultr), col = terrain.colors(60), axes = FALSE)
axis(1, 1:ncol(ultr), colnames(ultr))
axis(2, 1:nrow(ultr), rownames(ultr))

library(ggplot2)
library(gplots)
heatmap.2( ultr, Rowv=FALSE, Colv=FALSE, dendrogram='none', 
           cellnote=ultr, notecol="black", trace='none', 
           key=FALSE, lwid = c(.01,.99), lhei = c(.01,.99) )


#https://stackoverflow.com/questions/34803842/ggplot2-geom-pointrange-facet-grid-with-coord-flip-and-free-scales










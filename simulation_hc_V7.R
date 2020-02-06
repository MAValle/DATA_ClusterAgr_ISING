# En este script simulamos un ejemplo con N=25 nodos, para mostar
# como caso puntual y comprar entre todos los algoritmos.


# Prcedure:
# 1 create a coupling matrix of size N
# 2 get the merge dataframe with coupling distances of each iteration for:
#     * single linkage
#     * modularity method
#     * deterministic greedy
#     * probabilistic greedy with gamma = 0.5
#     * probabilistic greedy with gamma = 0 (100% probabilistic)
# 3 save coupling distances
# Hacer los dendogramas y comparar grafica de merge: iteracion vs distancia de acople
# y con medida de cluster woigthed coefficient.





# actual name: simulation_hc_V6.R
# creation: 21.dic.19


# Notes:
# 21-dic-19: creation
# 03-ene-20: simulacion para figure2eps_10samplings.eps
# 06-ene-20: simulamos figura figure3eps_n25_from_simulation_hc_V7_060120.eps, en donde
#           simulamso una matriz de acople J, y luego comparamos greedy approach 
#           deterministico con greedy probabilistico con distintos niveles de umbrales.
# 07-ene-20: simulacion en que graficamos distancias de acople de iteracion
#             23 vs distancias de acople con la iteracion 1, segun threshold.


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




# # # # # # # # # # # # # # # # # # # # # # # SIMULACION # # # # # # # # # # # # # # # # # # # # # # # # #
# dic 25, 2019
# Using J distributed with N(0,1)
# tomando como ejemplo simulation_hc_V6.R
set.seed(123)
#N <- c(10, 25, 50, 100)
N <- 50
thres <- c(1, 0.5, 0) # 1 greedy deterministico, 0 = greedy random
Jmean = 1.5
sj = 1

# CREATING COUPLINGS
ot <- create_coupling(Nn = N, media = Jmean, sj=sj, lw=-3, up=3, kind="normal") # para J normal
#ot <- create_coupling(Nn = N, media = Jmean, sj=sj, lw=-3, up=3, kind="uniform") # para J uniforme
J <- ot$J
D <- ot$D
#J
# # # # # # GREEDY
# dic, 18, 2019
invisible(capture.output( merge_greedy_detm <- hierarchical_clustering_probabilistic_greedy(D = D, J = J, trh= thres[1] ) ) )
# Begining the simulations for probabilistic greedy
I = 300
ware <- list()
ware[["greedy_semi"]] = list()
ware[["greedy_prob"]] = list()
for (i in 1:I){
 
  invisible(capture.output( merge_greedy_semi <- hierarchical_clustering_probabilistic_greedy(D = D, J = J, trh= thres[2] ) ) )
  invisible(capture.output( merge_greedy_prob <- hierarchical_clustering_probabilistic_greedy(D = D, J = J, trh= thres[3] ) ) )
  
  # # # # # # SAVING merges IN LISTS
  ware[["greedy_semi"]][[i]] = merge_greedy_semi[, c(1,2,6)]
  ware[["greedy_prob"]][[i]] = merge_greedy_prob[, c(1,2,6)]
  
  #
  cat("\n Iteration number :", i) 
  #flush.console()
}

# getting the means and creating the final dataframe to plot
# https://stackoverflow.com/questions/31839996/mean-of-each-row-of-data-frames-column-nested-in-list
Gsemi <- do.call(cbind, lapply(ware[["greedy_semi"]], "[", "dc"))
dc_mean_Gsemi <- rowMeans(Gsemi)
dc_sdev_Gsemi <- apply(Gsemi,1,sd)
Gprob <- do.call(cbind, lapply(ware[["greedy_prob"]], "[", "dc"))
dc_mean_Gprob <- rowMeans(Gprob)
dc_sdev_Gprob <- apply(Gprob,1,sd)
# 

# Preparing dataframe to ggplot
iteration <- rep(seq(from = 1, to = N-1),3)
algo <- c(rep("Gdet", N-1), rep("Gsemi", N-1), rep("Gprob", N-1))
medias <- c(merge_greedy_detm$dc, dc_mean_Gsemi, dc_mean_Gprob)
desvs <- c(rep(NA, N-1), dc_sdev_Gsemi, dc_sdev_Gprob) 

df <- data.frame(iteration = iteration, Algorithm = algo, 
                 mean = medias, sd = desvs)
df$min <- ifelse(df$mean - df$sd <=0, 0, log(df$mean - df$sd) )
df$max <- ifelse(df$mean + df$sd <=0, 0, log(df$mean + df$sd) )


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

cairo_ps("figure.eps")
f + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
ggsave("figure2eps_10samplings.eps")
save(f5, file = "figure5_from_simulation_hc_V6.Rdata")
# Lo que vemos es que mientras le agregamos mas randomness al 
# proceso greedy, la *media* de la distancia de acoples es cada
# ves mas mala en relacion al greedy deterministico, lo cual 
# tiene sentido puesto que al ser mas random, es como ir con 
# los ojos cerrados, y tendremos a veces mala suerte o buena
# suerte en mejorar la distancia de acople en relacion al greedy
# deterministico.

# Grafica de la desviacion estandar de las distancias de acople para
# los algoritmos greedy.
fsd <-  ggplot(df, aes(x = iteration, y = sd, colour = Algorithm)) + 
  geom_point(shape=20, size=1.5) 
fsd
# Esto es importante, porque vemos que mientras mas random le 
# agregamos al proceso, habra mas varianza, lo que significa que
# podemos tener muy buenos resultados o muy malos resultados 
# en referencia al grreedy deterministico. 
# # # # # # # # # # # # # # # # # # # # # # # SIMULACION # # # # # # # # # # # # # # # # # # # # # # # # #






# NUEVA SIMULACION PARA GRAFICA PAPER figure2eps_10samplings.eps
# # # # # # # # # # # # # # # # # # # # # # # SIMULACION # # # # # # # # # # # # # # # # # # # # # # # # #
# ene 03, 2020
# Using J distributed with N(0,1)
# Grafica para el paper Figura 2:
# una grafica en donde tomamos un J particular con N=50 nodos, y lo simulamos con:
#   * single linkage SL, modularity MOS y greddy algo1 deterministico
set.seed(123)
#N <- c(10, 25, 50, 100)
N <- 25
thres <- c(1, 0.5, 0) # 1 greedy deterministico, 0 = greedy random
Jmean = 1.5
sj = 1

# CREATING COUPLINGS
ot <- create_coupling(Nn = N, media = Jmean, sj=sj, lw=-3, up=3, kind="normal") # para J normal
#ot <- create_coupling(Nn = N, media = Jmean, sj=sj, lw=-3, up=3, kind="uniform") # para J uniforme
J <- ot$J
D <- ot$D
#J


# SINGLE LINKAGE
mst_g <- create_mst(J=J, D=D)
invisible(capture.output( merge_singlelink <- hierarchical_clustering_v2(D=D, J=J, mst_g=mst_g) ) )

# GREEDY DETERMINISTIC
invisible(capture.output( merge_greedy_detm <- hierarchical_clustering_probabilistic_greedy(D = D, J = J, trh= thres[1] ) ) )

# MODULARITY
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
iteration <- rep(seq(from = 1, to = N-1),3)
algo <- c(rep("MST", N-1), rep("Modularity", N-1), rep("Algorithm 1", N-1))
dc <- c(merge_singlelink$dc, merge_modularity$dc, merge_greedy_detm$dc)
df <- data.frame(Iteration = iteration, Algorithm = algo, CouplingDistance = dc)

# doing the ggplot
f1 <- ggplot(df, aes(x = iteration, y = (CouplingDistance), colour = Algorithm)) + 
  geom_point(shape=20, size=1.5) + geom_line() + 
  #geom_smooth(show.legend = FALSE, fill = "gray", span = 0.3) +
  #geom_errorbar(aes(ymin=log(mean-sd/sqrt(I)), ymax=log(mean+sd/sqrt(I))), width=.2, position=position_dodge(.01)) +
  #geom_errorbar(aes(ymin=log(mean-sd/sqrt(1)), ymax=log(mean+sd/sqrt(1))), width=.2, position=position_dodge(.01)) +
  labs(x ="Iteration", y = "Coupling Distance") +
  theme_linedraw() + theme(legend.position="bottom") +
  scale_x_continuous(breaks = seq(from=1, to=N-1, by=1) ) +
  #scale_y_continuous(breaks = c(-1, 0, 1, 2, 3,4,5,6,7,8)) +
  #scale_y_continuous(breaks = seq(from=1, to=8) ) +
  theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "white"))
f1 <- f1 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#cairo_ps("figure.eps")
#f + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#dev.off()
#ggsave("figure2eps_10samplings.eps")
#save(f5, file = "figure5_from_simulation_hc_V6.Rdata")


# graficas de los dendogramas para dendograms_n25_from_simulation_hc_V7_030120_eps.eps
hc_sl <- to_dendo2(merge_singlelink[, c(1,2,3,6) ], enames = c(1:ncol(D)))
hc_det <- to_dendo2(merge_greedy_detm[, c(1,2,3,6) ], enames = c(1:ncol(D)))
hc_mod <- to_dendo2(merge_modularity[, c(1,2,3,4) ], enames = c(1:ncol(D)))
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

mean(1-merge_singlelink$dc/merge_singlelink$dc[nrow(merge_singlelink)]) # 0.8033335
mean(1-merge_greedy_detm$dc/merge_greedy_detm$dc[nrow(merge_greedy_detm)]) # 0.9254003
mean(1-merge_modularity$dc/merge_modularity$dc[nrow(merge_modularity)]) # 0.9022149
# # # # # # # # # # # # # # # # # # # # # # # SIMULACION # # # # # # # # # # # # # # # # # # # # # # # # #






# NUEVA SIMULACION PARA GRAFICA PAPER figure3eps_n25_from_simulation_hc_V7_060120.eps
# # # # # # # # # # # # # # # # # # # # # # # SIMULACION # # # # # # # # # # # # # # # # # # # # # # # # #
# ene 06, 2020
# Using J distributed with N(0,1)
# Grafica para el paper Figura 2:
# una grafica en donde tomamos un J particular con N=50 nodos, y lo simulamos con:
#   * single linkage SL, modularity MOS y greddy algo1 deterministico
set.seed(123)
#N <- c(10, 25, 50, 100)
N <- 25
thres <- c(1, 0.7, 0.3, 0) # 1 greedy deterministico, 0 = greedy random
Jmean = 1.5
sj = 1

# CREATING COUPLINGS
ot <- create_coupling(Nn = N, media = Jmean, sj=sj, lw=-3, up=3, kind="normal") # para J normal
#ot <- create_coupling(Nn = N, media = Jmean, sj=sj, lw=-3, up=3, kind="uniform") # para J uniforme
J <- ot$J
D <- ot$D
#J


# GREEDY DETERMINISTIC thres=1
invisible(capture.output( merge_greedy_detm <- hierarchical_clustering_probabilistic_greedy(D = D, J = J, trh= thres[1] ) ) )

# Begining the simulations for probabilistic greedy
I = 200
ware <- list()
ware[["greedy_semi07"]] = list()
ware[["greedy_semi03"]] = list()
ware[["greedy_prob"]] = list()
for (i in 1:I){
  
  invisible(capture.output( merge_greedy_semi07 <- hierarchical_clustering_probabilistic_greedy(D = D, J = J, trh= thres[2] ) ) )
  invisible(capture.output( merge_greedy_semi03 <- hierarchical_clustering_probabilistic_greedy(D = D, J = J, trh= thres[3] ) ) )
  invisible(capture.output( merge_greedy_prob <-   hierarchical_clustering_probabilistic_greedy(D = D, J = J, trh= thres[4] ) ) )
  
  # # # # # # SAVING merges IN LISTS
  ware[["greedy_semi07"]][[i]] = merge_greedy_semi07[, c(1,2,6)]
  ware[["greedy_semi03"]][[i]] = merge_greedy_semi03[, c(1,2,6)]
  ware[["greedy_prob"]][[i]] =   merge_greedy_prob[, c(1,2,6)]
  
  #
  cat("\n Iteration number :", i) 
  #flush.console()
}

# getting the means and creating the final dataframe to plot
# https://stackoverflow.com/questions/31839996/mean-of-each-row-of-data-frames-column-nested-in-list
Gsemi07 <- do.call(cbind, lapply(ware[["greedy_semi07"]], "[", "dc"))
dc_mean_Gsemi07 <- rowMeans(Gsemi07)
dc_sdev_Gsemi07 <- apply(Gsemi07,1,sd)
Gsemi03 <- do.call(cbind, lapply(ware[["greedy_semi03"]], "[", "dc"))
dc_mean_Gsemi03 <- rowMeans(Gsemi03)
dc_sdev_Gsemi03 <- apply(Gsemi03,1,sd)
Gprob <- do.call(cbind, lapply(ware[["greedy_prob"]], "[", "dc"))
dc_mean_Gprob <- rowMeans(Gprob)
dc_sdev_Gprob <- apply(Gprob,1,sd)
# 

# What are the clustering coefficients?
cof07 <- lapply(ware[["greedy_semi07"]],
       function(x) { mean(1-x[,"dc"]/x[nrow(x),"dc"]) }    )
(do.call(mean, cof07) ) # 0.9223314
sd(unlist(cof07)) #0.004995167
cof03 <- lapply(ware[["greedy_semi03"]],
                function(x) { mean(1-x[,"dc"]/x[nrow(x),"dc"]) }    )
(do.call(mean, cof03) ) # 0.916
sd(unlist(cof03)) #0.011
cof00 <- lapply(ware[["greedy_prob"]],
                function(x) { mean(1-x[,"dc"]/x[nrow(x),"dc"]) }    )
(do.call(mean, cof00) ) # 0.8994
sd(unlist(cof00)) #0.017


# what iterations do I want to compare?
#its <- c(1, 16, 23)
its <- c(1, 8, 17, 23)


# Preparing dataframe to ggplot
# cluster 1
base_value1 <- merge_greedy_detm[its[1], "dc"]
algo <- c(rep("Threshold=0.7", I), rep("Threshold=0.3", I), rep("Threshold=0", I) )
values <- c(as.numeric(Gsemi07[its[1], ]), as.numeric(Gsemi03[its[1], ]), as.numeric(Gprob[its[1], ]) )
df1 <- data.frame(Algorithm = algo, CuplingDistance = values)
df1$Algorithm <- droplevels(df1$Algorithm)


# cluster 8
base_value8 <- merge_greedy_detm[its[2], "dc"]
algo <- c(rep("Threshold=0.7", I), rep("Threshold=0.3", I), rep("Threshold=0", I) )
values <- c(as.numeric(Gsemi07[its[2], ] ), as.numeric(Gsemi03[its[2], ]), as.numeric(Gprob[its[2], ]) )
df8 <- data.frame(Algorithm = algo, CuplingDistance = values)
df8$Algorithm <- droplevels(df8$Algorithm)


# cluster 17
base_value17 <- merge_greedy_detm[its[3], "dc"]
algo <- c(rep("Threshold=0.7", I), rep("Threshold=0.3", I), rep("Threshold=0", I) )
values <- c(as.numeric(Gsemi07[its[3], ] ), as.numeric(Gsemi03[its[3], ]), as.numeric(Gprob[its[3], ]) )
df17 <- data.frame(Algorithm = algo, CuplingDistance = values)
df17$Algorithm <- droplevels(df17$Algorithm)


# cluster 23
base_value23 <- merge_greedy_detm[its[4], "dc"]
algo <- c(rep("Threshold=0.7", I), rep("Threshold=0.3", I), rep("Threshold=0", I) )
values <- c(as.numeric(Gsemi07[its[4], ]), as.numeric(Gsemi03[its[4], ]), as.numeric(Gprob[its[4], ]) )
df23 <- data.frame(Algorithm = algo, CuplingDistance = values)
df23$Algorithm <- droplevels(df23$Algorithm)



#cairo_ps("figure3eps_n25_from_simulation_hc_V7_060120.eps")
# figure3eps_n25_from_simulation_hc_V7_270120-eps
par(mfrow=c(2,2))
boxplot((df1$CuplingDistance) ~ df1$Algorithm, col="lightgrey", 
        ylab="Coupling Distance", 
        #xlab="Threshold",
        ylim=c(0, 2.5), yaxs="i", main="Iteration 1",
        names=c("0", "03", "07"))
abline(h=base_value1, col="red")

#
boxplot((df8$CuplingDistance) ~ df8$Algorithm, col="lightgrey", 
        ylab="Coupling Distance", 
        #xlab="Threshold",
        #ylim=c(-2.5, 6), 
        #yaxs="i", 
        main="Iteration 8",
        names=c("0", "03", "07"))
abline(h=base_value8, col="red")

#
boxplot((df17$CuplingDistance) ~ df17$Algorithm, col="lightgrey", 
        ylab="Coupling Distance", 
        #xlab="Threshold",
        #ylim=c(-2.5, 6), 
        #yaxs="i", 
        main="Iteration 17",
        names=c("0", "03", "07"))
abline(h=base_value17, col="red")

#
boxplot((df23$CuplingDistance) ~ df23$Algorithm, col="lightgrey", 
        ylab="Coupling Distance", 
        #xlab="Threshold",
        #ylim=c(142, 300), 
        #yaxs="i", 
        main="Iteration 23",
        names=c("0", "03", "07"))
abline(h=base_value23, col="red")


# Jan 27, 2020
# We repeat the before boxplot, but with ggplot2
library(ggplot2)
library(ggpubr)
# https://smcclatchy.github.io/r-pub-quality-graphics/02-boxplot/
# http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/77-facilitating-exploratory-data-visualization-application-to-tcga-genomic-data/
# http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/
# # # # # # # # # # # # # # # # # # # # # # # SIMULACION # # # # # # # # # # # # # # # # # # # # # # # # #





# # # # # # # # # # # # # # # # # # # # # # # SIMULACION # # # # # # # # # # # # # # # # # # # # # # # # #
# 07-ene-20
# # CONCLUSION:
# se puede ver que no existe relacion alguna entre la distancia de acople
# en la it=1 con la de it=23.
# NUEVA SIMULACION PARA GRAFICA DE ULTIMA ITERACION CON LA PRIMERA
# para ver si cuando parte con alta energia de acople, termina con 
# baja o viceversa.
set.seed(123)
#N <- c(10, 25, 50, 100)
N <- 25
thres <- c(1, 0.7, 0.3, 0) # 1 greedy deterministico, 0 = greedy random
Jmean = 1.5
sj = 1

# CREATING COUPLINGS
ot <- create_coupling(Nn = N, media = Jmean, sj=sj, lw=-3, up=3, kind="normal") # para J normal
#ot <- create_coupling(Nn = N, media = Jmean, sj=sj, lw=-3, up=3, kind="uniform") # para J uniforme
J <- ot$J
D <- ot$D
#J

# GREEDY DETERMINISTIC thres=1
invisible(capture.output( merge_greedy_detm <- hierarchical_clustering_probabilistic_greedy(D = D, J = J, trh= thres[1] ) ) )

# Begining the simulations for probabilistic greedy
I = 200
ware <- list()
ware[["greedy_semi07"]] = list()
ware[["greedy_semi03"]] = list()
ware[["greedy_prob"]] = list()
for (i in 1:I){
  
  invisible(capture.output( merge_greedy_semi07 <- hierarchical_clustering_probabilistic_greedy(D = D, J = J, trh= thres[2] ) ) )
  invisible(capture.output( merge_greedy_semi03 <- hierarchical_clustering_probabilistic_greedy(D = D, J = J, trh= thres[3] ) ) )
  invisible(capture.output( merge_greedy_prob <-   hierarchical_clustering_probabilistic_greedy(D = D, J = J, trh= thres[4] ) ) )
  
  # # # # # # SAVING merges IN LISTS
  ware[["greedy_semi07"]][[i]] = merge_greedy_semi07[, c(1,2,6)]
  ware[["greedy_semi03"]][[i]] = merge_greedy_semi03[, c(1,2,6)]
  ware[["greedy_prob"]][[i]] =   merge_greedy_prob[, c(1,2,6)]
  
  #
  cat("\n Iteration number :", i) 
  #flush.console()
}


# getting the means and creating the final dataframe to plot
# https://stackoverflow.com/questions/31839996/mean-of-each-row-of-data-frames-column-nested-in-list
Gsemi07 <- do.call(cbind, lapply(ware[["greedy_semi07"]], "[", "dc"))
dc_mean_Gsemi07 <- rowMeans(Gsemi07)
dc_sdev_Gsemi07 <- apply(Gsemi07,1,sd)
Gsemi03 <- do.call(cbind, lapply(ware[["greedy_semi03"]], "[", "dc"))
dc_mean_Gsemi03 <- rowMeans(Gsemi03)
dc_sdev_Gsemi03 <- apply(Gsemi03,1,sd)
Gprob <- do.call(cbind, lapply(ware[["greedy_prob"]], "[", "dc"))
dc_mean_Gprob <- rowMeans(Gprob)
dc_sdev_Gprob <- apply(Gprob,1,sd)
# 
# what iterations do I want to compare?
#its <- c(1, 16, 23)
its <- c(1, 8, 17, 23)

# Quiero extraer solo distancia de acople de la it 1 y 23 solamente.
scatter07 <- lapply(ware[["greedy_semi07"]], 
       function(x) { c(x[1,"dc"], x[23,"dc"]) }   )
scatter07 <- as.data.frame(do.call(rbind, scatter07) )
scatter03 <- lapply(ware[["greedy_semi03"]], 
                    function(x) { c(x[1,"dc"], x[23,"dc"]) }   )
scatter03 <- as.data.frame(do.call(rbind, scatter03) )
scatter00 <- lapply(ware[["greedy_prob"]], 
                    function(x) { c(x[1,"dc"], x[23,"dc"]) }   )
scatter00 <- as.data.frame(do.call(rbind, scatter00) )

# Preparing dataframe to ggplot
algo <- c(rep("Threshold=0.7", I), rep("Threshold=0.3", I), rep("Threshold=0", I) )
value1 <- c(scatter07$V1, scatter03$V1, scatter00$V1)
value2 <- c(scatter07$V2, scatter03$V2, scatter00$V2)
df <- data.frame(Algorithm = algo, CouplingDistance1 = value1, CouplingDistance2 = value2)

# ggplot
p <- ggplot(df, aes(CouplingDistance1, CouplingDistance2, colour = Algorithm))
p + geom_point(size = 3) + scale_shape_manual(values = 21:25)


ggplot(df, aes(x = CouplingDistance1, y = CouplingDistance2 )) +
  geom_point(aes(color = Algorithm ) )

# # # # # # # # # # # # # # # # # # # # # # # SIMULACION # # # # # # # # # # # # # # # # # # # # # # # # #

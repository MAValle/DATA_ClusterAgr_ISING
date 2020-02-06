# Dada una solucion a una distancia de acople h, encontramos consistencia interna
# utilizando la funcion my_own_silhouette en  my_own_silhouette_function.R.
# Ver pag. 275-277 apuntes.




# actual name: computing_silhouettes_V2.R
# creation: 04-feb-20


# Notes:
# 04-feb-20: creation
#

# recursos:
# https://rpkgs.datanovia.com/factoextra/reference/fviz_silhouette.html

# # # # # # # # # # # # # # # # # # # preLOADING ENVIRONMENT # # # # # # # # # # # # # # # # # # # # # # # #
# Carga de datos, 
# igual que en clusteragrr_v1.R
rm(list = ls())
load("workdata_170220_from_simulation_hc_V8.RData")
source("silhouette_byhand_function.R")

library(igraph)
library(Matrix)
library(ggplot2)
# https://r.789695.n4.nabble.com/Cutting-hierarchical-cluster-tree-at-specific-height-fails-td4693737.html
cutree.h <- function(tree,h) {
  # this line adapted from cutree(...) code
  k <- nrow(tree$merge) + 2L - apply(outer(c(tree$height, Inf), h, ">"), 2, which.max)
  return(cutree(tree,k=k))
}
# # # # # # # # # # # # # # # # # # # preLOADING ENVIRONMENT # # # # # # # # # # # # # # # # # # # # # # # #





aquiiiiiii

# # # # # # # # # # # # SILHOUETTE FOR  DETERMINISTIC GREEDY # # # # # # # # # # # # # # # # # # # # # # # #
# start
invisible(capture.output( merge_greedy_detm <- hierarchical_clustering_probabilistic_greedy(D = D, J = J, trh=1 ) ) )
hc_det <- to_dendo2(merge_greedy_detm[, c(1,2,3,6) ], enames = c(1:ncol(D)))
hc_det$labels <- real_names
low = 5 # low level of cut distance
high = 300 # max level of cut distance
cuts <- c(low:high)
salida <- matrix(NA, ncol=3, nrow=length(cuts))
for(i in 1:length(cuts)) {
  members_det <- cutree.h(hc_det, h = cuts[i])
  sil <- silhouette_byhand(members = members_det, D=D)
  sil <- as.data.frame(sil)
  salida[i,] <- c(cuts[i], mean(sil$a), mean(sil$sil))
}
which(salida[,3] == max(salida[,3])) #[1] 1 2 3 4
colnames(salida) <- c("hcut", "a", "sil")
plot(salida[,1], salida[,3]) #cut con sil
plot(salida[,1], salida[,2]) #cut con a
# vemos que a una distancia de k= entre 5 y 8, se logra maximo de silhouette
# igual a 0.02835263.
# # # # # # # # # # # # SILHOUETTE FOR  DETERMINISTIC GREEDY # # # # # # # # # # # # # # # # # # # # # # # #


# Jan 30, 2020
# # # # # # # # # # # # # # # #  # # #   SILHOUETTE FOR  MST # # # # # # # # # # # # # # # # # # # # # # # #
mst_g <- create_mst(J=J, D=D)
invisible(capture.output( merge_singlelink <- hierarchical_clustering_v2(D=D, J=J, mst_g=mst_g) ) )
hc_SL <- to_dendo2(merge_singlelink[, c(1,2,3,6) ], enames = c(1:ncol(D)))
hc_SL$labels <- real_names
low = 5 # low level of cut distance
high = 300 # max level of cut distance
cuts <- c(low:high)
salida_SL <- matrix(NA, ncol=3, nrow=length(cuts))
for(i in 1:length(cuts)) {
  members_det <- cutree.h(hc_SL, h = cuts[i])
  sil <- silhouette_byhand(members = members_det, D=D)
  sil <- as.data.frame(sil)
  salida_SL[i,] <- c(cuts[i], mean(sil$a), mean(sil$sil))
}
which(salida_SL[,3] == max(salida_SL[,3])) #[1] 1 2 3 4
plot(salida_SL[,1], salida_SL[,3])
plot(salida_SL[,1], salida_SL[,2])
# vemos que a una distancia de k= entre 84 y 283, se logra maximo de silhouette
# igual a 0.06.
# # # # # # # # # # # # # # # #  # # #   SILHOUETTE FOR  MST # # # # # # # # # # # # # # # # # # # # # # # #







# PARA MST
members_SL <- cutree.h(hc_SL, h = 84)
assg <- data.frame(node_names = as.numeric(names(members_SL)), cl = members_SL)
#source("post_hoc_analysis_function.R")
values <- get_mean_clusters(assg = assg)
M <- values$averages
coupl_energy_average <- M[,2]/H[,3]
M <- cbind(M, coupl_energy_average)
colnames(M) <- c("cluster", "sum_coupl_energy", "cardinality", "mean_energy")
M
# # # # # # # # # # # # # # # #  # # #   ALTERNATIVE SILHOUETTE # # # # # # # # # # # # # # #  # # # # # # #

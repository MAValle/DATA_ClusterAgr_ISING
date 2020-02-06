
# En este script nos enfocamos en hacer las graficas de 
# heapmap, que originalmente iniciamos en 
# simulation_hc_V8.R para el caso real de 25 nodos (datos de verdad)

# Este script reproduce los heatmaps que iran en el paper




# actual name: heatmaps_figures.R
# creation: 28.ene.20


# Notes:
# 28-ene-20: creation
#




# # # # # # # # # # # # # # # # # # # preLOADING ENVIRONMENT # # # # # # # # # # # # # # # # # # # # # # # #
# Carga de datos, 
# igual que en clusteragrr_v1.R
rm(list = ls())
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






# # # # # # # # # # # # # # # # # #  PLOTTING HEATMAP OF COUPLING DISTANCES AND COPHONETIC DISTANCES# # # # # 
#22.01.20
# PLOTTING HEATMAP OF COUPLING DISTANCES AND COPHONETIC DISTANCES
# https://stackoverflow.com/questions/5639794/in-r-how-can-i-plot-a-similarity-matrix-like-a-block-graph-after-clustering-d
# http://mres.uni-potsdam.de/index.php/2017/05/08/using-the-new-function-heatmap-to-display-a-distance-matrix-in-cluster-analysis/
colnames(D) <- rownames(D) <- real_names
diag(D) <- NA
cut_value = 8 # determinado en computing_silhouettes_V1.R el 290120.

# # # # # heatmap para greedy deterministico
ord <- cutree.h(hc_det, h = cut_value)
coph <- cophenetic(hc_det)
colnames(D) <- rownames(D) <- real_names
diag(D) <- 0
# nececario para poner bien los productos donde corresponde
names(ord)
nuevo_ord <- ord
names(nuevo_ord) <- 1:ncol(D)
nuevo_ord <- sort(nuevo_ord)
order_wanted <- as.numeric(names(nuevo_ord))
# fin
#layout(matrix(1:2, ncol = 1))
cairo_ps("heatmapk8_originaldist_deterministic_from_heatmaps_figures_290120_eps.eps")
par(mar = c(2, 2, 2, 2)) # Set the margin on all sides to 2
image((log(D) )[order_wanted, order_wanted], col = heat.colors(50, alpha = 1), xaxt="n", yaxt="n")
mtext(side=1, text=colnames(D)[order_wanted], las=3, at=seq(0, 1, 0.0416))
mtext(side=2, text=colnames(D)[order_wanted], las=1, at=seq(0, 1, 0.0416))
dev.off()
coph <- as.matrix(coph)
diag(coph) <- 0
cairo_ps("heatmapk8_cophdist_deterministic_from_heatmaps_figures_290120_eps.eps")
par(mar = c(2, 2, 2, 2)) # Set the margin on all sides to 2
image( log(coph)[order_wanted, order_wanted], col=  heat.colors(50, alpha = 1), xaxt="n", yaxt="n")
mtext(side=1, text=colnames(D)[order_wanted], las=3, at=seq(0, 1, 0.0416))
mtext(side=2, text=colnames(D)[order_wanted], las=1, at=seq(0, 1, 0.0416))
dev.off()
#coph correlation
modif_D <- D[order_wanted, order_wanted]
modif_coph <- coph[order_wanted, order_wanted]
d <- modif_D[upper.tri(modif_D, diag=FALSE)]
c <- modif_coph[upper.tri(modif_coph, diag=FALSE)]
cor.test(d,c, alternative="two.sided", method="pearson" )
# p-value = 0.0004889   # rho = 0.2001042 



# # # # #  heatmap para greedy con gamma=0.5
ord <- cutree.h(hc_semi05, h = cut_value)
coph <- cophenetic(hc_semi05)
colnames(D) <- rownames(D) <- real_names
diag(D) <- 0
# nececario para poner bien los productos donde corresponde
names(ord)
nuevo_ord <- ord
names(nuevo_ord) <- 1:ncol(D)
nuevo_ord <- sort(nuevo_ord)
order_wanted <- as.numeric(names(nuevo_ord))
# fin
#layout(matrix(1:2, ncol = 1))
cairo_ps("heatmapk8_originaldist_semi05_from_heatmaps_figures_290120_eps.eps")
par(mar = c(2, 2, 2, 2)) # Set the margin on all sides to 2
image((log(D) )[order_wanted, order_wanted], col = heat.colors(50, alpha = 1), xaxt="n", yaxt="n")
mtext(side=1, text=colnames(D)[order_wanted], las=3, at=seq(0, 1, 0.0416))
mtext(side=2, text=colnames(D)[order_wanted], las=1, at=seq(0, 1, 0.0416))
dev.off()
coph <- as.matrix(coph)
diag(coph) <- 0
cairo_ps("heatmapk8_cophdist_semi05_from_heatmaps_figures_290120_eps.eps")
par(mar = c(2, 2, 2, 2)) # Set the margin on all sides to 2
image( log(coph)[order_wanted, order_wanted], col=  heat.colors(50, alpha = 1), xaxt="n", yaxt="n")
mtext(side=1, text=colnames(D)[order_wanted], las=3, at=seq(0, 1, 0.0416))
mtext(side=2, text=colnames(D)[order_wanted], las=1, at=seq(0, 1, 0.0416))
dev.off()
#coph correlation
modif_D <- D[order_wanted, order_wanted]
modif_coph <- coph[order_wanted, order_wanted]
d <- modif_D[upper.tri(modif_D, diag=FALSE)]
c <- modif_coph[upper.tri(modif_coph, diag=FALSE)]
cor.test(d,c, alternative="two.sided", method="pearson" )
#  p-value = 0.4155   # rho = 0.0471789 

# # # # #  heatmap para SL
ord <- cutree.h(hc_SL, h = cut_value)
coph <- cophenetic(hc_SL)
colnames(D) <- rownames(D) <- real_names
diag(D) <- 0
# nececario para poner bien los productos donde corresponde
names(ord)
nuevo_ord <- ord
names(nuevo_ord) <- 1:ncol(D)
nuevo_ord <- sort(nuevo_ord)
order_wanted <- as.numeric(names(nuevo_ord))
# fin
#layout(matrix(1:2, ncol = 1))
#image((log(D) )[ord, ord], col = heat.colors(50, alpha = 1), xaxt="n", yaxt="n", main = "Original distances")
cairo_ps("heatmapk8_originaldist_SL_from_heatmaps_figures_290120_eps.eps")
par(mar = c(2, 2, 2, 2)) # Set the margin on all sides to 2
image((log(D) )[order_wanted, order_wanted], col = heat.colors(50, alpha = 1), xaxt="n", yaxt="n")
mtext(side=1, text=colnames(D)[order_wanted], las=3, at=seq(0, 1, 0.0416))
mtext(side=2, text=colnames(D)[order_wanted], las=1, at=seq(0, 1, 0.0416))
dev.off()
coph <- as.matrix(coph)
diag(coph) <- 0
#image( log(coph)[ord, ord], col=  heat.colors(50, alpha = 1), xaxt="n", yaxt="n", main = "Cophenetic distances")
cairo_ps("heatmapk8_cophdist_SL_from_heatmaps_figures_290120_eps.eps")
par(mar = c(2, 2, 2, 2)) # Set the margin on all sides to 2
image((log(coph) )[order_wanted, order_wanted], col = heat.colors(50, alpha = 1), xaxt="n", yaxt="n")
mtext(side=1, text=colnames(D)[order_wanted], las=3, at=seq(0, 1, 0.0416))
mtext(side=2, text=colnames(D)[order_wanted], las=1, at=seq(0, 1, 0.0416))
dev.off()
#coph correlation
modif_D <- D[order_wanted, order_wanted]
modif_coph <- coph[order_wanted, order_wanted]
d <- modif_D[upper.tri(modif_D, diag=FALSE)]
c <- modif_coph[upper.tri(modif_coph, diag=FALSE)]
cor.test(d,c, alternative="two.sided", method="pearson" )
# p-value = 2.2e-16  # rho = 0.6390883 

# # # # # # # # # # # # # # # # # #  PLOTTING HEATMAP OF COUPLING DISTANCES AND COPHONETIC DISTANCES# # # # # 

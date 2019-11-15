

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
# 02.nov.19 FUNCION
# funcion creada en script simulation_hc_V2.R.
# para hacerla funcionar necesita las siguiente funciones:
#        create_coupling, create_mst, hierarchical_clustering_v2 y hierarchical_clustering_v5

# Funcion para simular un cluster jerarquico.
# Se simulan acoples J con distribucion normal con media mu y desv estandar sj
# todo viene de simulation_hc_v1.R y simulation_hc_v2.R
# inputs:
# N = numero de nodos
# mu
# vj
# output
# merg2 y merge4  matrices merge de single linkage normal y modificado de hierarchical_clustering_v2(D)
# y hierarchical_clustering_v4(D) respectivamente.
# grafos g y mst_g
# Objeto hc2 y hc4 para ploteo de dendograma.

# 12-nov-19: incorporamos la funcion hierarchical_clustering_greedy, para considerar
# el clustering utilizando el greedy algortiym.

simulation_hc_v1 <- function(Nn, mu=Jmean, sj=1) {
  # # # # # # # # # # # # # # # # # # # # CREATING A J MATRIX # # # # # # # # # # # # # # # # # # # # #
  # creating artifitial matrix of couplings
  # set.seed(134)
  ot <- create_coupling(Nn=Nn, media=Jmean, sj=sj, lw=-3, up=3)
  J <- ot$J
  D <- ot$D
  # # # # # # # # # # # # # # # # # # # # CREATING A J MATRIX # # # # # # # # # # # # # # # # # # # # #
  
  # # # # # # # # # # # # # # # # # # # # # FIND MST # # # # # # # # # # # # # # # # # # # # # #
  # Conforming the COUPLING Network
  mst_g <- create_mst(J=J, D=D)
  # # # # # # # # # # # # # # # # # # # # # FIND MST # # # # # # # # # # # # # # # # # # # # # #
  
  
  # # # # # # # # # # # # # # # # # # # # # # # HC # # # # # # # # # # # # # # # # # # # # # # #
  merge2 <- hierarchical_clustering_v2(D=D, J=J, mst_g=mst_g)
  hc2 <- to_dendo(D, merge2[,c(1:4)], enames=c(1:ncol(D)) )
  #plot(hc2)
  # ultrametric distancias: http://www.bioinfo.org.cn/lectures/index-68.html
  # https://stats.stackexchange.com/questions/103767/extract-ultrametric-distances-from-hclust-or-dendrogram
  #cophon2 <- cophenetic(hc2)
  
  merge5 <- hierarchical_clustering_v5(D=D, J=J, mst_g=mst_g)
  hc5 <- to_dendo(D, merge5[,c(1:4)], enames=c(1:ncol(D)) )
  #plot(hc4)
  #cophon4 <- cophenetic(hc4)
  
  mergeG <-  hierarchical_clustering_greedy(D = D, J = J)
  hcG <- to_dendo(D, merge5[,c(1:4)], enames=c(1:ncol(D)) )
  # # # # # # # # # # # # # # # # # # # # # # # HC # # # # # # # # # # # # # # # # # # # # # # #
  
  return(list(mst_g=mst_g, merge2=merge2, hc2=hc2, merge5=merge5, hc5=hc5, mergeG=mergeG, hcG=hcG, J=J, D=D))
}
# ejemplo:
# obj <- simulation_hc_v1(Nn=6, mu=0, sj=1)
# m2 <- obj$merge2
# m5 <- obj$merge5
# mG <- obj$mergeG
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  

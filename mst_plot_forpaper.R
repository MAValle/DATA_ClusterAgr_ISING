
# Replica de las figuras de red MST 
# Esto se basa en paper_congreso_figure_replica3.R.


# Creation date: 31-july-2019
# Filename: mst_plot_forpaper.R

# Notas:




# # # # # # # # # # # # # # # # # # # # CARGA de datos # # # # # # # # # # # # # # # # # # # # 
# 040319
rm(list = ls())
library(dplyr)
library(ggplot2)
library(latex2exp)
library(gdata)
library(scales)  # makes pretty labels on the x-axis
library(colorspace)
library(purrr)
library(igraph)
source("entropies_functions.R")
source("network_functions_v1.R")

load('new_inferring_parameters_environment270219.RData')
# # # # # # # # # # # # # # # # # # # # CARGA de datos # # # # # # # # # # # # # # # # # # # # 


# Following L615 from mstdistance_and_energies.R
# # # # # # # # # # # # # # # # # # # PLOTTING NETWORKS # # # # # # # # # # # # # # # # # ## # 

# Conforming the COUPLING Network
colnames(J) <- colnames(wb)
rownames(J) <- colnames(wb)
#rownames(h) <- c(1:dim(J)[1])
#http://www.shizukalab.com/toolkits/sna/plotting-networks-pt-2
net <- graph_from_adjacency_matrix(as.matrix(J), mode="upper", weighted=TRUE, diag=FALSE)
signos_de_J <- as.matrix(sign(upperTriangle(as.matrix(J), diag = FALSE, byrow=TRUE))) #obtenemos los signos de los acoples
E(net)$signos <- signos_de_J
E(net)$coupling <- E(net)$weight
E(net)$color[E(net)$signos > 0] <- 'orange' #acople positivo.
E(net)$color[E(net)$signos < 0] <- 'blue' #acople negativo
V(net)$magn <- abs(t(H))
V(net)$signs <- sign(H)
V(net)$color <- 'lightgrey'
#V(net)$color[V(net)$signs > 0] <- 'green'


# Conforming the MST Network
dis <- sqrt(-J + 3) 
# Convertir la matriz de distancia en objeto igraph
g <- graph_from_adjacency_matrix(dis, mode = "undirected", weighted = TRUE, diag = FALSE, add.colnames = NULL, add.rownames = NA)
E(g)$signos <- signos_de_J
E(g)$coupling <- E(net)$weight
E(g)$color[E(net)$signos > 0] <- 'orange' #acople positivo.
E(g)$color[E(net)$signos < 0] <- 'blue' #acople negativo
V(g)$magn <- abs(t(H))
V(g)$signs <- sign(H)
V(g)$color <- 'lightgrey'
#V(g)$color[V(g)$signs > 0] <- 'green'
mst_g <- minimum.spanning.tree(g, algorithm="prim")



# # # # # # COMPLETE COUPLING NETWORK without removal of edges
library(GGally)
library(ggnet)
library(network)
library(sna)
library(intergraph) #correspondencia entre igraph y network #http://intergraph.r-forge.r-project.org/howto.html

p1 <- ggnet2(net, alpha = 0.85, # Red completa
             mode = "circle",
             #node.size = V(net)$magn,
             #node.color = V(net)$color, #red is positive
             node.color = "grey87", #para paper congreso
             size = 10, 
             label=TRUE, 
             label.size=5, 
             edge.alpha = 0.8,
             edge.color = E(net)$color,
             ggtitle("Net of couplings"),
             edge.size = 5*abs(E(net)$weight))


# # # # MST
#https://www.rdocumentation.org/packages/GGally/versions/1.3.2/topics/ggnet
#https://briatte.github.io/ggnet/
mst_plot <- ggnet2(mst_g, alpha = 0.85, # MST
             #mode = "circle",
             mode = "fruchtermanreingold",
             #node.size = V(net)$magn,
             #node.color = V(mst_g)$color, #red is positive
             node.color = "grey87", #para paper congreso
             size = 10, 
             label=TRUE, 
             label.size=5, 
             edge.alpha = 0.8,
             edge.color = E(mst_g)$color,
             ggtitle("MST"),
             edge.size = 9/((E(mst_g)$weight-min(E(mst_g)$weight))/sd(E(mst_g)$weight)+1))
#edge.size = 5*abs(E(mst_g)$weight))

multiplot(p1, p2, cols=1)


# Cargamos grafica de umbrales de:
# plots_from_network_plot_forpaper310719.RData que fueron creadas en network_plot_forpaper.R
# para ponerla junto con el MST.
load('plots_from_network_plot_forpaper310719.RData')
multiplot(plt, mst_plot, cols=2)




# # # MAXIMUM spanning tree
# Conforming the MST Network
dis_ <- 1/sqrt(-J + 3) 
# Convertir la matriz de distancia en objeto igraph
g_ <- graph_from_adjacency_matrix(dis_, mode = "undirected", weighted = TRUE, diag = FALSE, add.colnames = NULL, add.rownames = NA)
E(g_)$signos <- signos_de_J
E(g_)$coupling <- E(net)$weight
E(g_)$color[E(net)$signos > 0] <- 'orange' #acople positivo.
E(g_)$color[E(net)$signos < 0] <- 'blue' #acople negativo
V(g_)$magn <- abs(t(H))
V(g_)$signs <- sign(H)
V(g_)$color <- 'lightgrey'
#V(g)$color[V(g)$signs > 0] <- 'green'
mst_g_ <- minimum.spanning.tree(g_, algorithm="prim")

mst_plot_ <- ggnet2(mst_g_, alpha = 0.85, # MST
                   #mode = "circle",
                   mode = "fruchtermanreingold",
                   #node.size = V(net)$magn,
                   #node.color = V(mst_g)$color, #red is positive
                   node.color = "grey87", #para paper congreso
                   size = 10, 
                   label=TRUE, 
                   label.size=5, 
                   edge.alpha = 0.8,
                   edge.color = E(mst_g_)$color,
                   ggtitle("MST"),
                   edge.size = 9/((E(mst_g_)$weight-min(E(mst_g_)$weight))/sd(E(mst_g_)$weight)+1))




# Replica de las figuras de redes de acoples con distintos umbrales 
# a la distirbución de los acoples. 
# Esto se basa en paper_congreso_figure_replica3.R.


# Creation date: 26-july-2019
# Filename: networks_plot_forpaper.R

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


# # # # # # REMOVAL OF EDGES between percentil 37.5 and 62.5
ic = 0.5
hist(E(net)$weight)
qq <- as.numeric(quantile(E(net)$weight, c(0.5-ic/2, 0.5+ic/2))) # dejando fuera el 25% alrededor de la media
quantile(E(net)$weight, c(0.5)) # 0.059 es la media

drops_j <- which(E(net)$weight > qq[1]  & E(net)$weight < qq[2]  )

net25 <- delete_edges(net,  drops_j ) 
plot(net25)
p2 <- ggnet2(net25, alpha = 0.85, # Red completa
             mode = "circle",
             #node.size = V(net)$magn,
             #node.color = V(net)$color, #red is positive
             node.color = "grey87", #para paper congreso
             size = 10, 
             label=TRUE, 
             label.size=5, 
             edge.alpha = 0.8,
             edge.color = E(net25)$color,
             ggtitle("Net of couplings"),
             edge.size = 5*abs(E(net25)$weight))
# # # # # # REMOVAL OF EDGES between percentil 37.5 and 62.5



# # # # # # REMOVAL OF EDGES between percentil 25 and 75
ic = 0.7
qq <- as.numeric(quantile(E(net)$weight, c(0.5-ic/2, 0.5+ic/2))) # dejando fuera el 25% alrededor de la media
drops_j <- which(E(net)$weight > qq[1]  & E(net)$weight < qq[2]  )

net50 <- delete_edges(net,  drops_j ) 
plot(net50)
V(net50)[igraph::degree(net50) == 0] 
p3 <- ggnet2(net50, alpha = 0.85, # Red completa
             mode = "circle",
             #node.size = V(net)$magn,
             #node.color = V(net)$color, #red is positive
             node.color = "grey87", #para paper congreso
             size = 10, 
             label=TRUE, 
             label.size=5, 
             edge.alpha = 0.8,
             edge.color = E(net50)$color,
             ggtitle("Net of couplings"),
             edge.size = 5*abs(E(net50)$weight))
# # # # # # REMOVAL OF EDGES between percentil 25 and 75



# # # # # # REMOVAL OF EDGES between percentil 12.5 and 87.5
ic = 0.9
qq <- as.numeric(quantile(E(net)$weight, c(0.5-ic/2, 0.5+ic/2))) # dejando fuera el 25% alrededor de la media
drops_j <- which(E(net)$weight > qq[1]  & E(net)$weight < qq[2]  )

net75 <- delete_edges(net,  drops_j ) 
plot(net75)
# nos dice los vertices que quedaron desconectados
V(net75)[igraph::degree(net75) == 0]  # le puse igraph::degree para eviat el conflicto con el paquete sna

p4 <- ggnet2(net75, alpha = 0.85, # Red completa
             mode = "circle",
             #node.size = V(net)$magn,
             #node.color = V(net)$color, #red is positive
             node.color = "grey87", #para paper congreso
             size = 10, 
             label=TRUE, 
             label.size=5, 
             edge.alpha = 0.8,
             edge.color = E(net75)$color,
             ggtitle("Net of couplings"),
             edge.size = 5*abs(E(net75)$weight))
# # # # # # REMOVAL OF EDGES between percentil 12.5 and 87.5

multiplot(p1, p3, p2, p4, cols=2)


# # # # # # # # # # # # # # # # # # # THRESHOLDS PLOTTINGS # # # # # # # # # # # # # # # # # ## # 
# 26-JUL-19
# aqui vamos variando el valor del intervalo de remosion de couplings y luego vemos el 
# porcentaje de nodos que quedan aislados.
intervals <- seq(from=0.1, to=0.99, by=0.001)
salida <- matrix(NA, ncol=3, nrow=length(intervals)) # intervalo quantil a ser removido,
for (i in 1:length(intervals)) {
  ic <- intervals[i]
  qq <- as.numeric(quantile(E(net)$weight, c(0.5-ic/2, 0.5+ic/2))) # dejando fuera el 25% alrededor de la media
  drops_j <- which(E(net)$weight > qq[1]  & E(net)$weight < qq[2]  )
  new_net <- delete_edges(net,  drops_j ) 
  #plot(new_net)
  drops <- V(new_net)[igraph::degree(new_net) == 0]
  rmoved_ratio <- length(drops)/length(V(net))
  con <- as.numeric(is_connected(new_net)) # 0 DESCONECTADO, 1 = SIGUE CONECTADO
  salida[i,] <- c(ic, rmoved_ratio, con)
}
colnames(salida) <- c("ic", "removed", "connected")
df <- data.frame(salida)
plot((df$ic), df$removed)
plot(df$ic, df$removed, type="s",
     xlab="% of smallest couplingsremoved", ylab="% of nodes with null degree")
plot(df$ic, df$con, type="s",
     xlab="% of smallest couplings removed", ylab="is.connected")
plt <- ggplot(df, aes(x=ic, y=removed)) +  geom_step(size=1.5) + scale_x_continuous(breaks=c(seq(from=0.1, to=0.99, by=0.1),1) ) +
  xlab("% of smallest couplings removed") + geom_vline(xintercept = 0.387, colour="red") + 
  #geom_text(aes(x=0.44, label="cut=0.386", y=0.2), colour="black", vjust = -1) +
  ylab("% of nodes with null degree") + theme_bw()

# Vemos que cuando el ic es del 40%, comienzan a aparecer nodos aislados.
save.image("~/Dropbox/Research/Project_cluster_aggregation/DATA_ClusterAgr_ISING/plots_from_network_plot_forpaper310719.RData")

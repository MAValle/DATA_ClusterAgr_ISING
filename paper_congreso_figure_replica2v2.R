# Hace exactamente lo mismo que paper_congreso_figure_replica2.R, pero lo hace
# con datos de acoples inferidos con 25 productos: (ver new_inferring_parameters.R)
# Las inferencias se encuentran en: new_inferring_parameters_environment270219.RData

# Replica de la figura de la distribucion de los acoples
# Seguimos el mismo procedimiento establecido en correlations_for_physicaA.R L280

# Creation date: 040319
# Filename: paper_congreso_figure_replica2v2.R

# Notas:


# # # # # # # # # # # # # # # # # # # # CARGA de datos # # # # # # # # # # # # # # # # # # # # 
# 040319
rm(list = ls())
library(dplyr)
library(ggplot2)
library(latex2exp)
library(gdata)
library(scales)  # makes pretty labels on the x-axis
source("ising_functions_v3.R")

load('new_inferring_parameters_environment270219.RData')
# # # # # # # # # # # # # # # # # # # # CARGA de datos # # # # # # # # # # # # # # # # # # # # 


# # # # # # # # # # # # # # # # # # # # PDF # # # # # # # # # # # # # #  # # # # # # # #  # # 
js <- upperTriangle(J, diag=FALSE, byrow=FALSE)  
js <- as.data.frame(js) 

h <- hist(js$js, plot=FALSE)  
h$counts <- h$counts/1e6  
plot(h)  

pl <- ggplot(js, aes(x = js)) + geom_freqpoly(binwidth = 0.3) #+ coord_trans(y = "log10") 
pg <- ggplot_build(pl) 
inf <- pg$data[[1]] 
inf$x 
inf$count 
df <- data.frame(X = inf$x, Y = inf$count) 
df <- df[-c(1,12),] 
df 
df$logy <- log10(df$Y) 
plot(df$X, df$logy) 
plot(df$X, df$Y) 
df$p <- df$Y/sum(df$Y) 
plot(df$X, df$p) 

ggplot(df, aes(x = X, y=log10(df$p))) + geom_point(shape=21, fill="black", color="darkred", size=2) +
  labs(x = "Jik") + labs(y = "log f(Jik)") + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 

ggplot(df, aes(X, p)) + geom_bar(stat = "identity", fill="darkgrey", width = 0.25) + 
  scale_x_continuous(breaks = round(df$X, digits=2) ) + 
  #annotate("text", x = -0.5, y = 0.25, label = "-0.46") + 
  labs(x = "") + labs(y = "") + 
  geom_vline(xintercept=mean(js$js), colour="red", linetype="dashed") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     panel.background = element_blank(), axis.line = element_line(colour = "black")) 
# # # # # # # # # # # # # # # # # # # # PDF # # # # # # # # # # # # # #  # # # # # # # #  # # 




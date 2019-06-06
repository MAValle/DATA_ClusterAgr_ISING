# Replica de la figura de la distribucion de los acoples
# Desde mstdistance_and_energies., utilizamos muestra = 19 
# para comparar los observables reales y simulados con metropolis hasting.
# Seguimos el mismo procedimiento establecido en correlations_for_physicaA.R L280

# Creation date: 170219
# Filename: paper_congreso_figure_replica2.R

# Notas:


# # # # # # # # # # # # # # # # # # # # CARGA de datos # # # # # # # # # # # # # # # # # # # # 
# 170219
rm(list = ls())
library(dplyr)
library(ggplot2)
library(latex2exp)
library(gdata)
source("ising_functions_v3.R")
#carga de las 250 muestras de weide_baskets
all_wbaskets <- readRDS("wbaskets250_doing_parallel_for_physicaA_210618.rds")
# carga de las 250 acoples y magnetizaciones de boltzmann learning
couplings <- readRDS("fieldsandcouplings250_doing_parallel_for_physicaA_090718.rds")

# Ejemplo para verificar que la muestra i tiene los mismos productos que la muestra i de acoples y magnetizaciones
wb <- all_wbaskets[[80]]
temp <- wb[,colSums(wb^2) !=0]
colnames(temp)
#rownames(a[["h_80"]])
#colnames(a[["j_80"]])
rownames(couplings[[80]]$H)
colnames(couplings[[80]]$J)
rm(temp, wb)
# # # # # # # # # # # # # # # # # # # # CARGA de datos # # # # # # # # # # # # # # # # # # # # 


# # # # # # # # # # # # # # # # # # # # PDF # # # # # # # # # # # # # #  # # # # # # # #  # # 
# 240219
# sample
muestra <- 19
H <- couplings[[muestra]]$H
J <- couplings[[muestra]]$J


hist(J)
hist(J, freq=FALSE)
library(ggplot2)
library(scales)  # makes pretty labels on the x-axis

# primera alternativa
js <- upperTriangle(J, diag=FALSE, byrow=FALSE)
js <- as.data.frame(js)
hist_data = hist(js$js, plot=F)
hist_data$counts = log10(hist_data$counts)
plot(hist_data, ylab='log10(Frequency)') # o que es lo mismo:
plot(hist_data$mids, hist_data$counts)


# algo parecido
h <- hist(js$js, plot=FALSE)
h$counts <- h$counts/1e6
plot(h)


# Mejor alternativa
ggplot(js, aes(x = js)) + geom_histogram(binwidth = 0.1)
pl <- ggplot(js, aes(x = js)) + geom_freqpoly(binwidth = 0.3) #+ coord_trans(y = "log10")
pg <- ggplot_build(pl)
inf <- pg$data[[1]]
inf$x
inf$count
df <- data.frame(X = inf$x, Y = inf$count)
df <- df[-c(1,9),]
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




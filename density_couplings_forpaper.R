# Este script calcula el histigrama y densidad y sus respectiva grafica 
# de la pdf de los acoples para efectos de publicacion en el paper.


# Para detalles ver el paper de conferencia a ICANN, y pag. 225-228 de los apuntes.



# Nombre del acutual archivo: density_couplings_forpaper.R
# ubicacion: dropbox->Research->Project_cluster_aggregation->DATA_ClusterAgr_ISING
# creation date: julio 10, 2019


#Notas:
# 10-jul-19: creation
# 20 jul: creamos la pdf de los acoples
# 21 jul: creamso la pdf de las correlaciones y scatterplot entre J y rho


# # # # # # # # # # # # # # # # # # # # CARGA de datos # # # # # # # # # # # # # # # # # # # # 
# 040319
rm(list = ls())
library(dplyr)
library(ggplot2)
library(latex2exp)
library(gdata)
library(scales)  # makes pretty labels on the x-axis
source("ising_functions_v3.R")
load('new_inferring_parameters_environment270219.RData') # load H and J parameters ising inferred
# # # # # # # # # # # # # # # # # # # # CARGA de datos # # # # # # # # # # # # # # # # # # # # 


isSymmetric(J) # TRUE
rho <- cor(wb)
isSymmetric(rho) # TRUE




* density de la pdf de los acoples con su media y varianza
* Q-Q plot de los acoples y correlaciones https://www.sheffield.ac.uk/polopoly_fs/1.579191!/file/stcp-karadimitriou-normalR.pdf
* density de la pdf de las correlacion y sus significancias
* plot de la correlacion vs couplings coloreando segun si la correlacion es significativa o no.


# Nota: algunas cosas que vienen de paper_congreso_figure_replica2v2.R

# # # # # # # # # # # # # # # # # # # # PDF # # # # # # # # # # # # # #  # # # # # # # #  # # 
js <- upperTriangle(J, diag=FALSE, byrow=FALSE)  
js <- as.data.frame(js) 
media <- mean(js$js)
varia <- var(js$js)


# histograma
h <- hist(js$js, plot=FALSE)  
h$counts <- h$counts/nrow(js) 
plot(h)  # es lo mismo que plot(h$mids, h$counts)
d <- density(js$js, kernel="gaussian", adjust=2.5)
lines(d$x, 0.25*d$y, col=2) # aqui le ponemos la densidad ajustada para que se vea con el histograma

# una distribucion normal de media y varianza igual a la de los acoples
xfit <- seq(min(js$js),max(js$js),length=40) 
#xfit <- seq(min(js$js),max(js$js),length=40) 
yfit <- dnorm(xfit, mean=media, sd=sqrt(varia) ) 
#yfit <- yfit*diff(h$mids[1:2])*length(h$counts)
yfit <- yfit*diff(h$mids[1:2])
#yfit <- yfit/sum(yfit)
lines(xfit, yfit, col="blue", lwd=2)

# histograma y densidad al mismo tiempo
hist(js$js, probability=T, main="Histogram of Couplings", xlab="Approximately normally distributed data")
d <- density(js$js, kernel="gaussian", adjust=2.5)
lines(d$x, d$y) #o
lines(d, col=2)

# solo grafico de densidad
d <- density(js$js, kernel="gaussian", adjust=2.5)
plot(d, main="Kernel Density of Couplings")
polygon(d, col="grey", border="black")

# grafico q-q
qqnorm(js$js, main="QQ plot of couplings",pch=19)
qqline(js$js)

# test de shapiro 
shapiro.test(js$js)


# grafica de histograma y density en ggplot2
df <- data.frame(h$mids, h$counts)
den <- data.frame(x= d$x, y = 0.25*d$y)
theme_set(
  theme_bw() + 
    theme(legend.position = "none") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
)
gg <- ggplot(df, aes(x = h.mids, y = h.counts)) +
  geom_segment(
    aes(x = h.mids, xend = h.mids, y = 0, yend = h.counts), 
    color = "lightgray"
  ) + geom_point(color='darkblue', size=3)
gg <- gg + geom_line(data = den, aes(x = x , y = y), color="red") + ylab("") + 
   scale_x_continuous("", breaks = h$mids)
# qq plot en ggplot2  
p <- ggplot(js, aes(sample = js)) + stat_qq() + stat_qq_line() + 
 stat_qq_line() + xlab("Theoretical quantiles") + ylab("Sample quantiles") 

gg2 <- ggplot(df, aes(x = h.mids, y = h.counts)) + geom_bar(stat = "identity",  color="white", fill="grey")
gg2 <- gg2 + geom_line(data = den, aes(x = x , y = y), color="red") + ylab("") + 
  scale_x_continuous("", breaks = h$mids)

# Inset de qq plot p, dentro de gg
# https://www.youtube.com/watch?v=Dkew_9bSCwE
p$layers <- rev(p$layers)
#A <- gg + annotation_custom(ggplotGrob(p), xmin=0.4, xmax=1.7, ymin=0.1, ymax=0.3)
A <- gg2 + annotation_custom(ggplotGrob(p), xmin=0.4, xmax=1.7, ymin=0.1, ymax=0.3)
A
# # # # # # # # # # # # # # # # # # # # PDF # # # # # # # # # # # # # #  # # # # # # # #  # # 


# # # # # # # # # # # # # # # # # # CORRELACIONES # # # # # # # # # # #  # # # # # # # #  # # 
# Vamos a hacer un plor de rho vs J, con un inset de la pdf de las correlaciones
#rho <- cor(wb)
library(Hmisc)
rho <- rcorr(wb, type = "pearson")
pv <- rho$P
pv <- upperTriangle(pv, diag=FALSE, byrow=FALSE)   # p-values
sig <- ifelse(pv<0.05, "yes", "no") # is significative the correlation?
rho <- rho$r
rho <- upperTriangle(rho, diag=FALSE, byrow=FALSE)  
js <- upperTriangle(J, diag=FALSE, byrow=FALSE)  
js <- as.data.frame(js) 
rho <- as.data.frame(rho) 
df <- data.frame(J = js, rho=rho, sig=sig)
plot(df$js, df$rho)
plot( df$rho, df$js, col= ifelse(df$sig =="no", "red", "black"))

# scatterplot rho versus couplings
theme_set(
  theme_bw() + 
    theme(legend.position = "none") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
)
sc <- ggplot(df, aes(x=rho, y=js)) + 
  geom_point(color="red", fill="black", shape=21, alpha=0.5, size=3, stroke = 0.3) + 
  xlab("") + ylab("")
# histograma de las correlaciones
h <- hist(df$rho, 20, plot=FALSE)  
h$counts <- h$counts/sum(h$counts)
hh <- data.frame(mids = h$mids, pr = h$counts)
ht <- ggplot(hh, aes(x = mids, y = pr)) + geom_bar(stat = "identity",  color="white", fill="grey") +
  scale_x_continuous("", breaks = seq(from = -0.07, to =0.39, by=0.04)) + 
  theme(axis.text.x = element_text(angle = 90)) + ylab("") 


# Inset de qq plot ht, dentro de sc
# https://www.youtube.com/watch?v=Dkew_9bSCwE
ht$layers <- rev(ht$layers)
B <- sc + annotation_custom(ggplotGrob(ht), xmin=0.16, xmax=0.36, ymin=-1.3, ymax=0.3)
B

# # # # # # # # # # # # # # # # # # CORRELACIONES # # # # # # # # # # #  # # # # # # # #  # # 




# Este script es parte del proyecto que viene de mstdistances_and_energies.R
# y lo que hace es utilizar la maquina de boltzmann para inferir parametros, 
# pero esta vez con mayor numero de subcatewgorias de productos.
# En los proyectos anteriores hemos obtenidos 250 muestras de wide_basket 
# para llevar a acabo 250 procesos de inferencia.
# esta vez pretendemos utilizar solo una inferencia pero tomando en cuenta
# al menos el 50% de los productos mas vendidos en terminos de volumen

# Creation date: 25-feb-19
# Nombre del actual archivo: new_inferring_parameters_.R
# Notas:
# 25-feb-19: creation

# Nota:
# De donde vienen los datos?
# La data original viene de al_R_from_Lampa_081216_csv.csv,
# Este archivo ha sifo utilizado en todos los proyectos para crear elist_basket (1722170 X 3)
# y wide_basket (179610 X 221). Estos componentes se encuentran en:
# workdata280417.RData
# En el primer proyecto Fondecyt estos componentes se generan en analysisV2.R
# Estos componentes tambien se encuentran en workdata260717.RData del proyecto Fondecyt2 (Ising)


# # # # # # # # # # # # # # # CARGA DE DATA # # # # # # # # # # # # # # # 
# 25-feb-19
rm(list = ls())
library(dplyr)
library(ggplot2)
library(latex2exp)
library(gdata)
library(scales) #para plot semilog
library(progress)
library(colorspace)
library(purrr)
library(igraph)
load("workdata260717.RData") # recordemos que ultima columna de wide_basket posee newid
source("ising_functions_v3.R")
source("entropies_functions.R")
rm(matrix_rho, node_importance, g.rho, mst.rho, mst.rho.pruned, mst.rho.pruned_vremoved, vnames)
#dm <- wide_basket*2 - 1
rm(elist_basket)
# # # # # # # # # # # # # # # CARGA DE DATA # # # # # # # # # # # # # # # 


# # # # # # # # # # # # ## VOLUMEN DE COMPRAS # # # # # # # # # # # # # # 
wb <- wide_basket[, -c(ncol(wide_basket))]
volume <- as.data.frame(colSums(wb))
volume <- colSums(wb)
volume <- sort(volume, decreasing = TRUE)
volume <- as.data.frame(volume)
#colnames(volume) <- "purchasings"
volume$subcat <- rownames(volume)
volume$porcen <- 100*volume$volume/sum(volume$volume)
volume$acum <- cumsum(volume$porcen)
plot(volume$acum)
volume$fila <- c(1:nrow(volume))
# vemos que con 50 productos ya tenemos app el 80% de todo el volumen.
# con 18 productos tenemos el 50% de todo el volumen.
# Conclusion: utilizaremos 24 productos que representa el 60% del volumen
# # # # # # # # # # # # ## VOLUMEN DE COMPRAS # # # # # # # # # # # # # # 


# # # # # # # # # # # SELECCION DE SUBCATEGORIAS # # # # # # ## #  # # # 
n <- 25
selected_subcat <- volume[c(1:n), "subcat"]
wb <- wb[, selected_subcat]
dim(wb)
rm(volume)
# # # # # # # # # # # SELECCION DE SUBCATEGORIAS # # # # # # ## #  # # # 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# PARAMETROS DE LA SIMULACION
LearningRate=0.8 #parametro de aprendizaje 0.1 el original funcionando
decay = 0.0001 #0.015 el original funcionando   0.0007
N = ncol(data) #numero de nodos
M = 40000 #num iteraciones metropolis hasting
E = 400 #numero de epocas (gradient descent) de la maquina
temp = matrix(NA, ncol=1, nrow=E)
temp[1,]=LearningRate #learning parameter
for (i in c(2:E)) {
  temp[i,] = temp[i-1,] *(1/(1 + decay*i))
}
plot(c(1:E), temp[,1], pch=20, ylab="Learning Rate", xlab="iterations")
rm(temp)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 



# # # # # # # # # # # # # # BOLTZMANN LEARNING # # # # # # # # ## # # # 
#BM: doing_boltzmannlearning_for_physicaA.R 
# L50
# Tenemos que borrar las columas con puros ceros
#borrar <- as.numeric(which(colMeans(wb) == 0))
#wbasket <- wb[, -c(borrar)]
#rm(borrar)
start_time <- Sys.time()
# kmax (num it metropolis = 20000), Kmax = numero it de la maquina
machine_results <- bmachine(condInitial = "a2", muestra = wb, Kmax = 400, metrop_its=40000,
                            LearningRate = 0.8, decay = 0.0001, N_sa = 10, umbral = .7, every = 5) #el kmax=metrop_its numero de iteraciones del metropolis es de 20000
H <- machine_results[[4]]
J <- machine_results[[5]]
end_time <- Sys.time()
end_time - start_time
save.image(file='new_inferring_parameters_environment270219.RData')
# 11.09032 hours Kmax = 20000
# condInitial = "a2", muestra = wb, Kmax = 200, LearningRate = 0.8, decay = 0.0005, N_sa = 10, umbral = .7, every = 10
# # # # # # # # # # # # # # BOLTZMANN LEARNING # # # # # # # # ## # # # 




# # # # # # # # # # # # BOLTZMANN LEARNING RESULTS # # # # # # #  # # # 
library(gdata)
library(latex2exp)
library(ggplot2)
# Recuperacion de la informacion
just_for_track <- machine_results[[1]]
diver <- machine_results[[2]]
cd <- machine_results[[3]]
inferences <- machine_results[[6]]
last_mh_sampling <- machine_results[[7]] # equivale a sampling[[1]] del ultimo metropolis hasting efectuado.
rownames(J) <- rownames(H)
colnames(J) <- rownames(H)

# Ploteos de ising-performance -----
# # # # # # # # # # # # # # # # PLOTEOS # # # # # # # # # # # # # # # # 
# Ploteos importantes 1) True versus simulated correlations, 2) True versus simulated <sisj>, 3) True versus simulated <si>
## PLoteo de True versus simulated correlations
wbasket <- wb
cw <- cor(wbasket)
cwr <- upperTriangle(cw, diag = FALSE, byrow=TRUE)
#cs <- cor(sampling[[1]])
cs <- cor(last_mh_sampling ) # en caso de que quiera utilizar el ultimo metropolis hasting de la maquina.
csr <- upperTriangle(cs, diag = FALSE, byrow=TRUE)
df <- data.frame(rho_sample = cwr, rho_mhasting = csr)
df$pc <- predict(prcomp(~cwr+csr, df))[,1]
q1 <- ggplot(df, aes(cwr, csr, color = pc)) +
  geom_point(shape = 16, size = 2, show.legend = FALSE, alpha = .999) +
  theme_minimal()  + 
  geom_abline(slope=1, intercept=0, alpha=0.7) + 
  scale_color_gradient(low = "#f0650e", high = "#0091ff") +
  labs(x=TeX(' $\\rho_{ij}$  '), y=TeX( '$ \\hat{\\rho_{ij} } $' ) ) +
  ggtitle('True versus estimated correlations' )

### Ploteo de los  <SiSj> real  vs  <Si><Sj> muestra simulada con metropolis. 15-abr-18
pp_real <- upperTriangle(pairproducts(wbasket))
#pp_sample <- upperTriangle(pairproducts(sampling[[1]]))
pp_sample <- upperTriangle(pairproducts(last_mh_sampling ) ) # en caso de que quiera utilizar el ultimo metropolis hasting de la maquina.
df <- data.frame(pp_real = pp_real, pp_sample = pp_sample)
df$pc <- predict(prcomp(~pp_real+pp_sample, df))[,1]
q2 <- ggplot(df, aes(pp_real, pp_sample, color = pc)) +
  geom_point(shape = 16, size = 2, show.legend = FALSE, alpha = .999) +
  theme_minimal() + 
  geom_abline(slope=1, intercept=0, alpha=0.7) + 
  scale_color_gradient(low = "#f0650e", high = "#0091ff") +
  labs(x=TeX(' $\\S_{i}S_{j}$  '), y=TeX( '$ \\hat{\\S_{i}S_{j} } $' ) ) +
  ggtitle('True versus estimated pairwise connections' )
q2
q2_rmse <- rmse(lm(pp_sample ~ pp_real, data=df) )

### Ploteo de los  <Si> real  vs  <Si> muestra simulada con metropolis. 17-abr-18
si_simu <- colMeans(last_mh_sampling )
si_real <-  colMeans(wbasket)
df <- data.frame(si_simu = si_simu, si_real = si_real)
df$pc <- predict(prcomp(~si_simu+si_real, df))[,1]
q3 <- ggplot(df, aes(si_real, si_real, color = pc)) +
  geom_point(shape = 16, size = 2, show.legend = FALSE, alpha = .999) +
  theme_minimal() + 
  geom_abline(slope=1, intercept=0, alpha=0.7) + 
  scale_color_gradient(low = "#f0650e", high = "#0091ff") +
  labs(x=TeX(' $<s_{i}>$  '), y=TeX( '$ \\hat{<s_{i}> } $' ) ) +
  ggtitle('True versus estimated one-body averages' )
q3
q3_rmse <- rmse(lm(si_simu ~ si_real, data=df) )

# # # # # # # # # # # ## # # # COMPUTATION OF Cij# ## # # # # ## # # # # ## # # # # ##  # # # # ## 
# Modification date: 20-12-17
# Calculo de los Cij = <SiSj> - <Si><Sj> y plot. La funcion esta en ising_functions_v3.R
d <- pairwise_cij(wbasket, last_mh_sampling)
# # # # # # # # # # # ## # # # COMPUTATION OF Cij# ## # # # # ## # # # # ## # # # # ##  # # # # ## 
## Plot of Cij verdadero vs Cij estimado
d$pc <- predict(prcomp(~V4+V3, d))[,1]
q4 <- ggplot(d, aes(V3, V4, color = pc)) +
  geom_point(shape = 16, size = 2, show.legend = FALSE, alpha = .999) +
  theme_minimal() + 
  geom_abline(slope=1, intercept=0, alpha=0.7) + 
  scale_color_gradient(low = "#f0650e", high = "#0091ff") +
  labs(x=TeX(' $<C_{ij}>$  '), y=TeX( '$ \\hat{<C_{ij}> } $' ) ) +
  ggtitle('True versus estimated two body averages' )
q4
q4_rmse <- rmse(lm(V3 ~ V3, data=d) )

### Ploteo de la divergencia 
diver <-  diver[complete.cases(diver), ]
diver <- as.data.frame(diver)
q5 <- ggplot(data = diver, aes(x = V1, y = V2)) + geom_line(color = "red", size = 1) +
  labs(x="Epochs", y=TeX( '$ KL $' ) ) +
  ggtitle('Kullback Leibler Divergence' )
multiplot(q1, q2, q3, q4, q5, cols=1)
# # # # # # # # # # # # # # # # PLOTEOS # # # # # # # # # # # # # # # # 

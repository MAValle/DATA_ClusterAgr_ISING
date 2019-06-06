# Replica de la figura de recovered <si> , <sisj> y Cij para el paper congreso
# Desde mstdistance_and_energies., utilizamos muestra = 19 
# para comparar los observables reales y simulados con metropolis hasting.
# Seguimos el mismo procedimiento establecido en one_and_two_body_abalysisV2.R

# Creation date: 170219
# Filename: paper_congreso_figure_replica.R

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


# # # # # # # # # # # # # # #  # # # METROPOLIS-HASTING # # # # # # # # # # #  # # # # # # # # 
# 170219
# generacion de los datos
muestra = 19  # con 1,3 malos, con 2 mas o menos
wbasket <- all_wbaskets[[muestra]]
wbasket <- wbasket[,colSums(wbasket^2) !=0] # dejamos las columnas de productos en que no hay ceros.
wbasket <- wbasket[,-1] #quitamos columna de newid
#impar <- 2*(muestra-1) + 1
#par <- 2*muestra
#chequeo de que productos en wbasket = productos en h
if(  all.equal(colnames(wbasket), rownames(couplings[[muestra]]$H)) != TRUE ) message("Products between wbasket and magnetizations did not match") 
#chequeo de que productos en wbasket = productos en J
if( all.equal(colnames(wbasket), rownames(couplings[[muestra]]$J)) != TRUE ) message("Products between wbasket and couplings did not match") 

H <- couplings[[muestra]]$H
J <- couplings[[muestra]]$J
# Inicio de sampleo Metropolis-Hastings
initial_spins=rep(1,ncol(wbasket))
sampling <- metropolis.sampling(magnetizaciones=H, acoples = J, initial_spin = initial_spins, kmax=1.5*nrow(wbasket), ks=0, To=NULL, alpha=NULL)
output <- pairwise_cij(wbasket, sampling[[1]]) # i - j - Cij de wbasket (V3), Cij de sampling[[1]] (V4)
colnames(output) <- c("V1", "V2", "real_Cij", "est_Cij")
# # # # # # # # # # # # # # #  # # # METROPOLIS-HASTING # # # # # # # # # # #  # # # # # # # # 


# # # # # # # # # #  # ## # # grafica de <si> reales y estimados# # # # # # # # ## # #  # # # # 
si_simu <- colMeans(sampling[[1]] )
si_real <-  colMeans(wbasket)
df <- data.frame(si_simu = si_simu, si_real = si_real)
df$pc <- predict(prcomp(~si_simu+si_real, df))[,1]
si <- ggplot(df, aes(si_real, si_real, color = pc)) +
  #geom_point(shape = 16, size = 2, show.legend = FALSE, alpha = .999) +
  geom_point(shape = 16, size = 3, show.legend = FALSE) +
  #geom_point(shape = 21, size = 3, fill="white", stroke=1, show.legend = FALSE, alpha = 1.0) +
  #geom_point(shape = 21, size = 3, stroke=1, show.legend = FALSE, alpha = 1.0) +
  theme_bw()  + 
  geom_abline(slope=1, intercept=0, alpha=0.7) + 
  scale_color_gradient(low = "#f0650e", high = "#0091ff") +
  labs(x=TeX(' $<s_{i}>$  '), y=TeX( '$ \\hat{<s_{i}> } $' ) ) +
  #ggtitle('True versus estimated one-body averages' )
  ggtitle('A')
rmse(lm(si_simu ~ si_real, data=df) )
cor(df$si_real, df$si_simu)
si
# # # # # # # # # #  # ## # # grafica de <si> reales y estimados# # # # # # # # ## # #  # # # # 


# # # # # # # # # # # # # # grafica de <sisj>  reales y estimados # # # # # # # ## # #  # # # # 
pp_real <- upperTriangle(pairproducts(wbasket))
pp_sample <- upperTriangle(pairproducts(sampling[[1]]))
#pp_sample <- upperTriangle(pairproducts(last_mh_sampling ) ) # en caso de que quiera utilizar el ultimo metropolis hasting de la maquina.
df <- data.frame(pp_real = pp_real, pp_sample = pp_sample)
df$pc <- predict(prcomp(~pp_real+pp_sample, df))[,1]
sisj <- ggplot(df, aes(pp_real, pp_sample, color = pc)) +
  #geom_point(shape = 16, size = 2, show.legend = FALSE, alpha = .999) +
  geom_point(shape = 16, size = 2, show.legend = FALSE) +
  #geom_point(shape = 21, size = 3, fill="white", stroke=1, show.legend = FALSE, alpha = 1.0) +
  theme_bw()  + 
  geom_abline(slope=1, intercept=0, alpha=0.7) + 
  scale_color_gradient(low = "#f0650e", high = "#0091ff") +
  labs(x=TeX(' $<s_{i}s_{j}>$  '), y=TeX( '$ \\hat{<s_{i}s_{j}> } $' ) ) +
  #labs(x=TeX(' $\\<s_{i}s_{j}>$  '), y=TeX( '$ \\hat{\\<S_{i}S_{j}> } $' ) ) +
  #ggtitle('True versus estimated pairwise connections' )
  ggtitle('B' )
rmse(lm(pp_sample ~ pp_real, data=df) )
cor(df$pp_real, df$pp_sample)
sisj
# # # # # # # # # # # # # # grafica de <sisj>  reales y estimados # # # # # # # ## # #  # # # # 


# # # # # # # # # # # # # # # grafica de Cij  reales y estimados # # # # # # # # # # #  # # # # 
output$pc <- predict(prcomp(~ est_Cij + real_Cij, output))[,1]
cij <- ggplot(output, aes(real_Cij, est_Cij, color = pc)) +
  geom_point(shape = 16, size = 2, show.legend = FALSE) +
  #geom_errorbar(data=final, mapping=aes(ymin = mean_cij_sampling - sdev, ymax = mean_cij_sampling + sdev), size=0.3, color="black") + 
  #theme_minimal() + 
  theme_bw() + 
  geom_abline(slope=1, intercept=0, alpha=0.7) + 
  scale_color_gradient(low = "#f0650e", high = "#0091ff") +
  labs(x=TeX(' $<C_{ij}>$  '), y=TeX( '$ \\hat{<C_{ij}> } $' ) ) +
  #ggtitle('True versus estimated two body connections' ) 
  ggtitle('C' ) 
#annotate("text", x = -0.015, y = 0.025, label = "RMSE==.0075", parse = TRUE) + 
#annotate("text", x = -0.019, y = 0.020, label = "rho==.789", parse = TRUE) 
cij  
rmse(lm(real_Cij ~ est_Cij, data=output) )
cor(output$real_Cij, output$est_Cij)
# # # # # # # # # # # # # # # grafica de Cij  reales y estimados # # # # # # # # # # #  # # # # 



# # # # # # # # # # # # # # # # # # # GRAFICA FINAL # # # # # # # # # # # # # # # # # # # # # # 
multiplot(si, sisj, cij, cols=1)
# # # # # # # # # # # # # # # # # # # GRAFICA FINAL # # # # # # # # # # # # # # # # # # # # # # 

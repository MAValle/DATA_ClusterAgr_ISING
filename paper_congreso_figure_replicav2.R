# Hace exactamente lo mismo que paper_congreso_figure_replica.R, pero lo hace
# con datos de acoples inferidos con 25 productos: (ver new_inferring_parameters.R)
# Las inferencias se encuentran en: new_inferring_parameters_environment270219.RData

# Replica de la figura de recovered <si> , <sisj> y Cij para el paper congreso
# Seguimos el mismo procedimiento establecido en one_and_two_body_abalysisV2.R

# Creation date: 040319
# Filename: paper_congreso_figure_replicav2.R

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




# # # # # # # # # # # # # # # # Recuperacion de informacion # # # # # # ## # # # # # # # # # # 
just_for_track <- machine_results[[1]]
diver <- machine_results[[2]]
cd <- machine_results[[3]]
inferences <- machine_results[[6]]
last_mh_sampling <- machine_results[[7]] # equivale a sampling[[1]] del ultimo metropolis hasting efectuado.
rownames(J) <- rownames(H)
colnames(J) <- rownames(H)
output <- pairwise_cij(wb, last_mh_sampling ) # i - j - Cij de wbasket (V3), Cij de sampling[[1]] (V4)
colnames(output) <- c("V1", "V2", "real_Cij", "est_Cij")
# # # # # # # # # # # # # # # # Recuperacion de informacion # # # # # # ## # # # # # # # # # # 



# # # # # # # # # #  # ## # # grafica de <si> reales y estimados# # # # # # # # ## # #  # # # # 
si_simu <- colMeans( last_mh_sampling )
si_real <-  colMeans(wb)
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
pp_real <- upperTriangle(pairproducts(wb))
pp_sample <- upperTriangle(pairproducts( last_mh_sampling ) )
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

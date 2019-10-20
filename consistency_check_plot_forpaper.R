# Graficas de consistencia de <si> vs <si*>, y <Cik> vs <Cik*>-

# Replica de la figura de recovered <si> , <sisj> y Cij para el paper 
# Seguimos el mismo procedimiento establecido en one_and_two_body_abalysisV2.R
# y exactamente el mismo para el paper de congreso en paper_congreso_figure_replica2-R

# Creation date: 26-jul-19
# Filename: consistency_ckeck_plot_forpaper.R

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
theme_set(
  theme_bw() + 
    theme(legend.position = "none") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
)
si_simu <- colMeans( last_mh_sampling )
si_real <-  colMeans(wb)
df <- data.frame(si_simu = si_simu, si_real = si_real)
df$pc <- predict(prcomp(~si_simu+si_real, df))[,1]
si <- ggplot(df, aes(si_real, si_real, color = pc)) +
  #geom_point(shape = 16, size = 2, show.legend = FALSE, alpha = .999) +
  geom_point(shape = 16, size = 7, alpha=0.9, show.legend = FALSE) +
  #geom_point(shape = 21, size = 3, fill="white", stroke=1, show.legend = FALSE, alpha = 1.0) +
  #geom_point(shape = 21, size = 3, stroke=1, show.legend = FALSE, alpha = 1.0) +
  geom_abline(slope=1, intercept=0, alpha=0.7) + 
  scale_color_gradient(low = "#f0650e", high = "#0091ff") +
  #labs(x=TeX(' $<s_{i}>$  '), y=TeX( '$ \\hat{<s_{i}> } $' ) ) +
  labs(x="", y="" ) +
  #ggtitle('True versus estimated one-body averages' )
  ggtitle('A') +
  theme_bw() + 
    theme(
      plot.background = element_blank()
      ,panel.grid.major = element_blank()
      ,panel.grid.minor = element_blank()
    ) 
rmse(lm(si_simu ~ si_real, data=df) )
cor(df$si_real, df$si_simu)
si
# # # # # # # # # #  # ## # # grafica de <si> reales y estimados# # # # # # # # ## # #  # # # # 





# # # # # # # # # # # # # # # grafica de Cij  reales y estimados # # # # # # # # # # #  # # # # 
output$pc <- predict(prcomp(~ est_Cij + real_Cij, output))[,1]
cij <- ggplot(output, aes(real_Cij, est_Cij, color = pc)) +
  geom_point(shape = 16, size = 7, alpha = 0.9, show.legend = FALSE) +
  #geom_errorbar(data=final, mapping=aes(ymin = mean_cij_sampling - sdev, ymax = mean_cij_sampling + sdev), size=0.3, color="black") + 
  #theme_minimal() + 
  geom_abline(slope=1, intercept=0, alpha=0.7) + 
  scale_color_gradient(low = "#f0650e", high = "#0091ff") +
  labs(x="", y="" ) +
  #ggtitle('True versus estimated two body connections' ) 
  ggtitle('B' ) +
  theme_bw() + 
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
  ) 
#annotate("text", x = -0.015, y = 0.025, label = "RMSE==.0075", parse = TRUE) + 
#annotate("text", x = -0.019, y = 0.020, label = "rho==.789", parse = TRUE) 
cij  
rmse(lm(real_Cij ~ est_Cij, data=output) )
cor(output$real_Cij, output$est_Cij)
# # # # # # # # # # # # # # # grafica de Cij  reales y estimados # # # # # # # # # # #  # # # # 



# # # # # # # # # # # # # # # # # # # GRAFICA FINAL # # # # # # # # # # # # # # # # # # # # # # 
multiplot(si, cij, cols=2)
# # # # # # # # # # # # # # # # # # # GRAFICA FINAL # # # # # # # # # # # # # # # # # # # # # # 

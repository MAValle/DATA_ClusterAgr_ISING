# Funcion para calcular la correlacion cofonetica
# recursos:
# https://www.sciencedirect.com/science/article/abs/pii/0031320378900389
# https://journalofinequalitiesandapplications.springeropen.com/articles/10.1186/1029-242X-2013-203
# https://rdrr.io/cran/dendextend/man/cor_cophenetic.html
# https://cran.r-project.org/web/packages/mdendro/mdendro.pdf
# https://stats.stackexchange.com/questions/149852/validate-dendrogram-in-cluster-analysis-what-is-the-meaning-of-cophenetic-corre

# Jan 17, 2020
# name: cophonetic_cor_function.R

# inputs:
# matriz D de distancias (pueden ser las distancias de acople)
# matriz ultr, que son las distancias ultrametricas o de acople

cophonetic_cor <- function(D, ultr, hc_clust) {
  ultr <- as.matrix(cophenetic(hc_clust))
  ultr <- ultr[upper.tri(ultr)]
  dcopl <- D[upper.tri(D)]
  value <- cor(dcopl, ultr, method="pearson" )
  return(value)
}


# ejemplo
# coph <- cophonetic_cor(D = D, ultr = as.matrix(cophenetic(hc_clust)), hc_clust=hc_semi075 )
# donde hc_ckust es un dendogram object

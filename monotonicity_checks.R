# 14-ago-19
# el problema de la monotonocidad
# http://mathworld.wolfram.com/MonotonicFunction.html
# https://en.wikipedia.org/wiki/Monotonic_function


library(Matrix)
set.seed(13)
x <- Matrix(rnorm(400),20)  # los J vienen con media 0 y desv estandar 1
J <- forceSymmetric(x)
J
D <- sqrt(3-J)
D

j <- J[upper.tri(J)]
d <- D[upper.tri(D)]
plot(j,d)

n=2
tt <- choose(20,n)
Sa <- matrix(NA, ncol = 2, nrow = tt)
for (i in 1:tt){
  id <- sample(1:length(j),n)
  jsample <- j[id]
  dsample <- d[id]
  jSum <- sum(jsample)
  dSum <- sum(dsample)
  Sa[i,] <- c(jSum, dSum)
  plot(jsample,dsample)
}
plot(Sa[,1], Sa[,2])  #### vemos que parece no se presente la monotononicity of the function transformation sqrt(3-J)


# y si estandarizamos los J entre -1 y 1?
Je <- 2*(J- min(J) )/(max(J) - min(J)) -1
je <- Je[upper.tri(Je)]
D <- sqrt(3-Je)
#D <- sqrt(2*(1-Je)) # con esta distancia resulta peor
D
d <- D[upper.tri(D)]
plot(je,d)

tt <- choose(20,n)
Sa <- matrix(NA, ncol = 2, nrow = tt)
for (i in 1:tt){
  id <- sample(1:length(j),n)
  jsample <- je[id]
  dsample <- d[id]
  jSum <- sum(jsample)
  dSum <- sum(dsample)
  Sa[i,] <- c(jSum, dSum)
  plot(jsample,dsample)
}
plot(Sa[,1], Sa[,2])  # #### vemos que mejora , ya no se ven tantas no-monotonocidades
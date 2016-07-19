library(mixpack)
library(normalmultinomial)
load('data/Pigs.rda')

X = as.matrix(Pigs)

fit = normalmultinomial_fitting(X)
(MU <- fit[[1]])
(SIGMA <- fit[[2]])

A = stepE(X = X, mu = MU, sigma = SIGMA)

Prob = cbind(exp(A), 1) / apply(cbind(exp(A), 1), 1, sum)
X.replaced = as.data.frame(apply(X, 1, sum) * Prob)


library(fpc)

#RAW = kmeansruns(scale(X, scale = FALSE),krange=1:10,critout=TRUE,runs=1000,criterion="ch")
#plot(RAW$crit)

pka <- kmeansruns(ilr_coordinates(X.replaced),krange=1:10,critout=TRUE,runs=1000,criterion="ch")
plot(pka$crit)

KM = lapply(2:15, function(k) kmeans(ilr_coordinates(X.replaced), centers = k, nstart = 100))

plot(sapply(KM, function(km) km$tot.withinss))

biplot(princomp(clr_coordinates(X.replaced)), asp=1, scale=0, xlim = c(-5, 5), ylim = c(-5,5))

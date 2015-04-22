library(Rcpp)
sourceCpp('src/fitting.cpp')

library(compositions)
library(dplyr)
d <- read.csv("data/ContZero2Marc.csv", stringsAsFactors=FALSE)

set.seed(1)
ind = c(5,1,3)

X = d[,3+ind]
cluster = kmeans(d[,3+ind], 3)$cluster

X1 = unname(as.matrix(split(X, cluster)[[1]]))
fit1 = adjustNormalMultinomial(X1, prop = 0.4)
plot(ilr(fit1[[3]]))

fit = adjustNormalMultinomial(unname(as.matrix(X)))

plot(compositions::ilr(X))
plot(compositions::ilr(fit[[3]]), add=T, col='red')

compositions::plot.rcomp(X)
compositions::plot.acomp(fit[[3]], add=T, col='red')



library(Rcpp)
sourceCpp('src/fitting.cpp')

library(compositions)
library(dplyr)
d <- read.csv("data/ContZero2Marc.csv", sep=';', stringsAsFactors=FALSE)

set.seed(1)
ind = c(5,1,3)

X = d[,3+ind]

fit1 = normalmultinomial_fitting(unname(as.matrix(X)), niter= 20, nsim = 1000)

library(ggplot2)
library(ggtern)
library(mixpack)

X.fit1 = data.frame(fit1[[3]])
names(X.fit1) = names(X)
ggtern() + 
  geom_point(data=X, aes(x = FEEDER, y=BED, z = PASSAGE), size=3) +
  geom_point(data=X.fit1, aes(x = FEEDER, y=BED, z = PASSAGE), col='red', size=3)


df = ilr_coordinates(X)
df.fit1 = ilr_coordinates(fit1[[3]])
ggplot() + 
  geom_point(data=df, aes(x = coord.1, y=coord.2), size=3) +
  geom_point(data=df.fit1, aes(x = coord.1, y=coord.2), col='red', size=3) + xlim(-5, 4.5) + ylim(-3, 2.5)



cluster = as.factor(kmeans(X, 3)$cluster)
ggtern() + 
  geom_point(data=X, aes(x = FEEDER, y=BED, z = PASSAGE, shape=cluster), size=3)

fit2 = lapply(split(X, cluster), function(.X){
  normalmultinomial_fitting(unname(as.matrix(.X)), niter= 20, nsim = 1000)
})

X.fit2 = data.frame(do.call('rbind', lapply(fit2, function(.fit) .fit[[3]])))
names(X.fit2) = names(X)
X.fit2$cluster = rep(names(fit2), sapply(fit2, function(.fit) nrow(.fit[[3]])))

ggtern() + 
  geom_point(data=X, aes(x = FEEDER, y=BED, z = PASSAGE, shape=cluster), size=3) +
  geom_point(data=X.fit2, aes(x = FEEDER, y=BED, z = PASSAGE, shape=cluster), col='red', size=3)


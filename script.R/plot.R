library(Rcpp)
sourceCpp('src/fitting.cpp')

#Paramatres 
N = 1000
Mu = c(2, 0)
Sigma = matrix(c(1, -0.25,
                 -0.25, 1), nrow=2) * 0.71
Size = sample(50:100, N, replace=TRUE)

## Plot de la distribució normal
library(compositions)
X = rnorm.multinom(N, Mu, Sigma, Size = Size)
sum(X==0)

cc = adjustNormalMultinomial(X = X, 1e-1)

par(mfrow=c(1,2))
plot.rcomp(X)
plot.rcomp(alrInv(cc[[4]]), add=TRUE, col='red')
plot(ilr(X))
points(ilr(alrInv(cc[[4]])), col='red')
points(ilr( alrInv(cc[[1]]) ), col='blue', cex=1)
par(mfrow=c(1,1))

library(Rcpp)
sourceCpp('src/fitting.cpp')

d <- read.csv("data/ContZero2Marc.csv", stringsAsFactors=FALSE)
S = 2
X1 = c(S,0,0)
X2 = c(1,S,1)
X3 = c(1,1,S)
X = rbind(X1, X2, X3)
X.nonZero  = X
X.nonZero[X == 0] = 0.3

A = log(X.nonZero[,1:2]/X.nonZero[,3])


fit = adjustNormalMultinomial_internal(X = X, A=A, iter=10, eps=1e-10, minSigma=1e-10)

library(compositions)
plot(ilr(X))
plot(ilr(fit[[3]]), add=T, col='red')

plot(rcomp(X))
plot(rcomp(fit[[3]]), add=T, col='red')


a = matrix(alr(fit[[3]][1,]), nrow=1)
x = X[1,]
inv_sigma = solve(fit[[2]])
mu = fit[[1]]

normpart = (2*pi)^(-k/2) * det(inv_sigma)^(-0.5) * exp(-0.5 * (matrix(a-mu, ncol=2) %*% inv_sigma %*% t(matrix(a-mu, ncol=2)))[1,1] )

multpart = factorial(sum(x)) / prod(factorial(x)) * prod( as.numeric(compositions::alrInv(a))^x )

log(normpart) + log(multpart)

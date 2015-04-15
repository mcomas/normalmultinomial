library(Rcpp)
sourceCpp('src/fitting.cpp')

d <- read.csv("data/ContZero2Marc.csv", stringsAsFactors=FALSE)
S = 100
X1 = c(S,1,1)
X2 = c(1,S,1)
X3 = c(1,1,S)
X = rbind(X1, X2, X3)
X.nonZero  = X
X.nonZero[X == 0] = 0.3

A = log(X.nonZero[,1:2]/X.nonZero[,3])


a = A[1,]
mu = apply(A, 2, mean)
inv_sigma = solve(cov(A))
x = X[1,]

mvf(a = a, mu = mu, inv_sigma = inv_sigma, x = x)

k = length(a)

normpart = (2*pi)^(-k/2) * det(inv_sigma)^(-0.5) * exp(-0.5 * (matrix(a-mu, ncol=2) %*% inv_sigma %*% t(matrix(a-mu, ncol=2)))[1,1] )
multpart = factorial(sum(x)) / prod(factorial(x)) * prod( as.numeric(compositions::alrInv(a))^x )

log(normpart) + log(multpart )

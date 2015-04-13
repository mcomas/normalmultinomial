# library(Rcpp)
# sourceCpp('src/fitting.cpp')
# source('R/basic.R')
library(normalmultinomial)

set.seed(1)
#Paramatres 
N = 100
Mu = c(0, 0)
Sigma = matrix(c(1, -0.75,
                 -0.75, 1), nrow=2) 
Size = sample(1*100:500, N, replace=TRUE)

## Plot de la distribució normal
library(compositions)
X = rnorm.multinom(N, Mu, Sigma, Size = Size)
sum(apply(X!=0, 1, prod) == 0)

A = log(X[,1:2] / X[,3])

mu = apply(A, 2, mean)
inv_sigma = solve(cov(A))

#Mstep(A = A)

cc = adjustNormalMultinomial(X = X)

df = data.frame(rbind(matrix(alr(X), nrow=nrow(X)), cc[[3]]), 'src' = rep(c(1 , 2), each=nrow(X)))
df$src[c(rep(F, nrow(X)), apply(X!=0, 1, prod) == 0)] = 3

par(mfrow=c(1,2))
plot.rcomp(rcomp(X))
plot(alrInv(cc[[3]]), add=TRUE, col='red')
plot(df[,1:2], col=df[,3])
points( cc[[1]][1,1], cc[[1]][1,2], col='blue', cex=1)
par(mfrow=c(1,1))

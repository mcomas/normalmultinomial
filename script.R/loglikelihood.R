library(Rcpp)
sourceCpp('src/fitting.cpp')

set.seed(1)
N = 50
Mu = c(1, 0)
Sigma = matrix(c(1, +0.25,
                 +0.25, 1), nrow=2)
Size = sample(50:100, N, replace=TRUE)

X = rnormalmultinomial(mu = Mu, sigma = Sigma, size = Size)
logLikelihood(X = X, mu = Mu, sigma = Sigma, N = 1000)

plot.rcomp(X)
fit = adjustNormalMultinomial(X = X, iter=1)
plot.rcomp(fit[[3]], col='red', add=T)
fit = adjustNormalMultinomial(X = X, iter=5)
plot.rcomp(fit[[3]], col='blue', add=T)
fit = adjustNormalMultinomial(X = X, iter=10)
plot.rcomp(fit[[3]], col='green', add=T)
fit = adjustNormalMultinomial(X = X, iter=15)
plot.rcomp(fit[[3]], col='orange', add=T)


set.seed(1)
N = 300
Mu = c(3, 0)
Sigma = matrix(c(1, +0.25,
                 +0.25, 1), nrow=2) * 3
Size = sample(50:100, N, replace=TRUE)

X = rnormalmultinomial(mu = Mu, sigma = Sigma, size = Size)
logLikelihood(X = X, mu = Mu, sigma = Sigma, N = 1000)

plot.rcomp(X)
fit = adjustNormalMultinomial(X = X, nmont = 0)













logLikelihood(X = X, mu = fit[[1]], sigma = fit[[2]], N = 100000)


0K = ncol(X)
k = K -1
sum(X == 0)

X.nonZero  = X
X.nonZero[X == 0] = 0.3

A = log(X.nonZero[,1:k]/X.nonZero[,K])

x = seq(1,8,0.05)
y = seq(-6,3,0.01)

Mu
Sigma

# 975:  56  1  0
# 990:  36   30   10
iobs = 975
X[iobs,]

df = data.frame( expand.grid(x,y) )
df$z =  apply(df, 1, mvf, mu = Mu, sigma = Sigma, x = X[iobs,]) 
library(reshape2)
z = as.matrix( dcast(df, Var1~Var2, value.var = 'z') )[,-1]

contour(x,y,z, nlevels = 50)
points(A[iobs,1], A[iobs,2])

opt = Mstep(A = matrix(A[iobs,], nrow=1), mu = Mu, sigma = Sigma, X = matrix(X[iobs,], nrow=1))
points(opt, col='red')

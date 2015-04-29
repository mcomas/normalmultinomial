library(normalmultinomial)
library(microbenchmark)

mu = c(1,1,1,4)
Sigma = diag(c(1,2,1,1))

X = rnormalmultinomial(mu, Sigma, size = rep(100, 1000))

microbenchmark(
  normalmultinomial_fitting(X, niter = 2, nsim = 100),
  normalmultinomial_fitting(X, niter = 2, nsim = 100, version = 2),
  normalmultinomial_fitting(X, niter = 2, nsim = 100, version = 3, nthreads = 1),
  normalmultinomial_fitting(X, niter = 2, nsim = 100, version = 3, nthreads = 2),
  normalmultinomial_fitting(X, niter = 2, nsim = 100, version = 3, nthreads = 3),
  normalmultinomial_fitting(X, niter = 2, nsim = 100, version = 3, nthreads = 4), times = 10)


Xfit = normalmultinomial_fitting(X, niter = 20, nsim =200)
Xfit[[1]]
Xfit[[2]]

system.time( Xfit <- normalmultinomial_fitting(X, niter = 20, nsim = 1000, version = 2) )
Xfit[[1]]
Xfit[[2]]

system.time( Xfit <- normalmultinomial_fitting(X, niter = 20, nsim = 1000, version = 3, nthreads = 8) )
Xfit[[1]]
Xfit[[2]]

head(X)
head(Xfit[[3]] * apply(X, 1, sum))


library(mixpack)

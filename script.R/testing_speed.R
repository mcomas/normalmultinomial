library(normalmultinomial)
library(microbenchmark)

mu = c(1,1,1,1)
Sigma = diag(c(1,1,1,1))

X = rnormalmultinomial(mu, Sigma, size = rep(100, 1000))

microbenchmark(
  normalmultinomial_fitting(X, niter = 2, nsim = 100),
  normalmultinomial_fitting(X, niter = 2, nsim = 100, version = 2),
  normalmultinomial_fitting(X, niter = 2, nsim = 100, version = 3, nthreads = 1),
  normalmultinomial_fitting(X, niter = 2, nsim = 100, version = 3, nthreads = 2),
  normalmultinomial_fitting(X, niter = 2, nsim = 100, version = 3, nthreads = 3),
  normalmultinomial_fitting(X, niter = 2, nsim = 100, version = 3, nthreads = 4), times = 10)



library(mixpack)

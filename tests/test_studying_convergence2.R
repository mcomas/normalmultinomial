library(parallel)

stepEM = function(fit){
  list(
    matrix(apply(fit[[1]], 2, mean), nrow=1),
    matrix(apply(fit[[2]], 2, mean), nrow=1),
    apply(fit[[3]], c(1,2), mean))
}


mu.steps = function(expected_function, X, MU, SIGMA, NSIM = 1000, NITER = 100){
  MU.steps = matrix(0, nrow=NITER,ncol=k)
  for(iter in 1:NITER){
    FIT = lapply(split(X, 1:NROW(X)), expected_function, MU, SIGMA, nsim = NSIM)
    MOM = lapply(FIT, stepEM)
    MU = apply(sapply(MOM, function(mom) mom[[2]]), 1, mean)
    SIGMA = matrix(apply(sapply(MOM, function(mom) mom[[3]]), 1, mean), ncol = k)
    SIGMA = SIGMA - MU %*% t(MU)
    MU.steps[iter,] = MU
  }
  MU.steps
}
ma = function(x,n=10){
  filter(x,rep(1/n,n), sides=2)
}


X = as.matrix(Pigs)
k = NCOL(X) -1
MU = rep(0, k)
SIGMA = diag(k)


MU3.steps = mu.steps(expectedA3_all, X, MU, SIGMA, NSIM = 100, NITER = 500)
MU3.steps.smooth = apply(MU3.steps, 2, ma)

par(mfrow=c(3,2))
plot(MU3.steps[,1], type = 'l')
lines(MU3.steps.smooth[,1], col='red')
plot(MU3.steps[,2], type = 'l')
lines(MU3.steps.smooth[,2], col='red')
plot(MU3.steps[,3], type = 'l')
lines(MU3.steps.smooth[,3], col='red')
plot(MU3.steps[,4], type = 'l')
lines(MU3.steps.smooth[,4], col='red')
plot(MU3.steps[,5], type = 'l')
lines(MU3.steps.smooth[,5], col='red')
par(mfrow=c(1,1))

MU4.steps = mu.steps(expectedA4_all, MU, SIGMA, NSIM = 100, NITER = 500)
MU4.steps.smooth = apply(MU4.steps, 2, ma)

par(mfrow=c(3,2))
plot(MU4.steps[,1], type = 'l')
lines(MU4.steps.smooth[,1], col='red')
plot(MU4.steps[,2], type = 'l')
lines(MU4.steps.smooth[,2], col='red')
plot(MU4.steps[,3], type = 'l')
lines(MU4.steps.smooth[,3], col='red')
plot(MU4.steps[,4], type = 'l')
lines(MU4.steps.smooth[,4], col='red')
plot(MU4.steps[,5], type = 'l')
lines(MU4.steps.smooth[,5], col='red')
par(mfrow=c(1,1))

stepEM = function(fit){
  list(
    matrix(apply(fit[[1]], 2, mean), nrow=1),
    matrix(apply(fit[[2]], 2, mean), nrow=1),
    apply(fit[[3]], c(1,2), mean))
}

X = as.matrix(Pigs)
x1 = as.numeric(apply(X, 1, which.max) == 3)
x2 = as.numeric(apply(X, 1, which.max) == 5)
A = cbind(1,x1,x2)

expected_function = expectedA2_all



k = NCOL(X) -1
MU = matrix(0, nrow=NROW(X), ncol=k)
SIGMA = diag(k)
NSIM = 100
NITER = 100

mu.steps = function(expected_function, X, MU, SIGMA, NSIM = 1000, NITER = 500){
  coefs = list()
  for(iter in 1:NITER){
    FIT = lapply(1:NROW(X), function(i){ expected_function(x = X[i,], mu = MU[i,], sigma = SIGMA, nsim=NSIM) })
    MOM = lapply(FIT, stepEM)

    M1 = t(sapply(MOM, function(mom) mom[[2]]))
    M2 = matrix(apply(sapply(MOM, function(mom) mom[[3]]), 1, mean), ncol = k)

    B <- solve(t(A) %*% A) %*% t(A) %*% M1
    coefs[[iter]] <- B

    MU = A %*% B
    SIGMA = nearestPSD(SIGMA - (t(MU) %*% MU) / (NROW(M1)))
  }
  coefs
}

COEFS = mu.steps(expectedA3_all, X, MU, SIGMA, NSIM = 100, NITER = 2000)
par(mfrow=c(5,3))
plot(sapply(COEFS, function(m) m[1,1]), type = 'l')
plot(sapply(COEFS, function(m) m[1,2]), type = 'l')
plot(sapply(COEFS, function(m) m[1,3]), type = 'l')
plot(sapply(COEFS, function(m) m[1,4]), type = 'l')
plot(sapply(COEFS, function(m) m[1,5]), type = 'l')

plot(sapply(COEFS, function(m) m[2,1]), type = 'l')
plot(sapply(COEFS, function(m) m[2,2]), type = 'l')
plot(sapply(COEFS, function(m) m[2,3]), type = 'l')
plot(sapply(COEFS, function(m) m[2,4]), type = 'l')
plot(sapply(COEFS, function(m) m[2,5]), type = 'l')

plot(sapply(COEFS, function(m) m[3,1]), type = 'l')
plot(sapply(COEFS, function(m) m[3,2]), type = 'l')
plot(sapply(COEFS, function(m) m[3,3]), type = 'l')
plot(sapply(COEFS, function(m) m[3,4]), type = 'l')
plot(sapply(COEFS, function(m) m[3,5]), type = 'l')
par(mfrow=c(1,1))

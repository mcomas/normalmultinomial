library(normalmultinomial)
x = c(1,0)
MU = 0
SIGMA = matrix(1)

NSIM = 1000
# Using a Montecarlo algorithm with mu=MU and sigma=SIGMA
fit1 <- expectedA1_all(x, MU, SIGMA, nsim = NSIM)
# Using a Montecarlo algorithm with mu = maximum of f and sigma=SIGMA
fit2 <- expectedA2_all(x, MU, SIGMA, nsim = NSIM)
# Using a Montecarlo algorithm with mu = maximum of f and sigma=hessian of f
fit3 <- expectedA3_all(x, MU, SIGMA, nsim = NSIM)
# Using a Montecarlo algorithm with mu = maximum of f and sigma=hessian of f
fit4 <- expectedA4_all(x, MU, SIGMA, nsim = NSIM)

par(mfrow=c(1,3))

a1 = sample(fit1[[1]][,1])
a2 = sample(fit2[[1]][,1])
a3 = sample(fit3[[1]][,1])
a4 = rep(fit4[[1]][,1], each=2)

plot(cumsum(a1)/seq_along(a1), type='l', xlab='Iterations', ylab = 'Estimation', main='Zero moment', ylim = c(0.35,0.65))
lines(cumsum(a2)/seq_along(a2), col='red')
lines(cumsum(a3)/seq_along(a3), col='blue')
lines(cumsum(a4)/seq_along(a4), col='green')

a1 = sample(fit1[[2]])
a2 = sample(fit2[[2]])
a3 = sample(fit3[[2]])
a4 = rep(fit4[[2]], each=2)

plot(cumsum(a1)/seq_along(a1), type='l', xlab='Iterations', ylab = 'Estimation', main='First moment')
lines(cumsum(a2)/seq_along(a2), col='red')
lines(cumsum(a3)/seq_along(a3), col='blue')
lines(cumsum(a4)/seq_along(a4), col='green')


a1 = sample(fit1[[3]][1,1,])
a2 = sample(fit2[[3]][1,1,])
a3 = sample(fit3[[3]][1,1,])
a4 = rep(fit4[[3]][1,1,], each=2)

plot(cumsum(a1)/seq_along(a1), type='l', xlab='Iterations', ylab = 'Estimation', main='Second moment')
lines(cumsum(a2)/seq_along(a2), col='red')
lines(cumsum(a3)/seq_along(a3), col='blue')
lines(cumsum(a4)/seq_along(a4), col='green')
par(mfrow=c(1,1))

##########################
##########################3
NSIM = 100
# Using a Montecarlo algorithm with mu=MU and sigma=SIGMA
var(fit1 <- replicate(1000, apply(expectedA1_all(x, MU, SIGMA, nsim = NSIM)[[1]], 2, mean)))
# Using a Montecarlo algorithm with mu = maximum of f and sigma=SIGMA
var(fit2 <- replicate(1000, apply(expectedA2_all(x, MU, SIGMA, nsim = NSIM)[[1]], 2, mean)))
# Using a Montecarlo algorithm with mu = maximum of f and sigma=hessian of f
var(fit3 <- replicate(1000, apply(expectedA3_all(x, MU, SIGMA, nsim = NSIM)[[1]], 2, mean)))
# Using a Montecarlo algorithm with mu = maximum of f and sigma=hessian of f
var(fit4 <- replicate(1000, apply(expectedA4_all(x, MU, SIGMA, nsim = NSIM)[[1]], 2, mean)))

###########################
###########################
###########################

data(Pigs)
x = as.matrix(Pigs[1,])
k = NCOL(x) -1
MU = rep(0, k)
SIGMA = diag(k)

NSIM = 10000
# Using a Montecarlo algorithm with mu=MU and sigma=SIGMA
fit1 <- expectedA1_all(x, MU, SIGMA, nsim = NSIM)
# Using a Montecarlo algorithm with mu = maximum of f and sigma=SIGMA
fit2 <- expectedA2_all(x, MU, SIGMA, nsim = NSIM)
# Using a Montecarlo algorithm with mu = maximum of f and sigma=hessian of f
fit3 <- expectedA3_all(x, MU, SIGMA, nsim = NSIM)
# Using a Montecarlo algorithm with mu = maximum of f and sigma=hessian of f
fit4 <- expectedA4_all(x, MU, SIGMA, nsim = NSIM)

par(mfrow=c(2,3))
for(i in 1:5){
  #a1 = sample(fit1[[2]][,i])
  a2 = sample(fit2[[1]][,1])
  a3 = sample(fit3[[1]][,1])
  a4 = rep(fit4[[1]][,1], each=2)

  plot(cumsum(a2)/seq_along(a2), type='l', xlab='Iterations', ylab = 'Estimation', main='Zero moment')
  lines(cumsum(a3)/seq_along(a3), col='red')
  lines(cumsum(a4)/seq_along(a4), col='blue')
}

par(mfrow=c(2,3))
for(i in 1:5){
  #a1 = sample(fit1[[2]][,i])
  a2 = sample(fit2[[2]][,i])
  a3 = sample(fit3[[2]][,i])
  a4 = rep(fit4[[2]][,i], each=2)

  plot(cumsum(a2)/seq_along(a2), type='l', xlab='Iterations', ylab = 'Estimation', main='First moment')
  lines(cumsum(a3)/seq_along(a3), col='red')
  lines(cumsum(a4)/seq_along(a4), col='blue')
}

par(mfrow=c(2,3))
for(i in 1:5){
  #a1 = sample(fit1[[2]][,i])
  a2 = sample(fit2[[3]][i,i,])
  a3 = sample(fit3[[3]][i,i,])
  a4 = rep(fit4[[3]][i,i,], each=2)

  plot(cumsum(a2)/seq_along(a2), type='l', xlab='Iterations', ylab = 'Estimation', main=sprintf('Second moment (%d,%d)',i,i))
  lines(cumsum(a3)/seq_along(a3), col='red')
  lines(cumsum(a4)/seq_along(a4), col='blue')
}


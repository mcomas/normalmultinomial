library(normalmultinomial)
set.seed(1)
x = c(1,0)
MU = 0
SIGMA = matrix(10)

NSIM = 1000
fit0 <- expectedMonteCarlo(x, MU, SIGMA, nsim = NSIM, antithetic_variates = F, importance_sampling_mu = F)
# Using a Montecarlo algorithm with mu=MU and sigma=SIGMA
fit1 <- expectedMonteCarlo(x, MU, SIGMA, nsim = NSIM, antithetic_variates = T, importance_sampling_mu = F)
fit2 <- expectedMonteCarlo(x, MU, SIGMA, nsim = NSIM, antithetic_variates = F, importance_sampling_mu = T)
fit3 <- expectedMonteCarlo(x, MU, SIGMA, nsim = NSIM, antithetic_variates = T, importance_sampling_mu = T)
fit4 <- expectedMetropolis(x, MU, SIGMA, nsim = NSIM)

par(mfrow=c(2,1))
a0 = fit0[[2]]
a1 = fit1[[2]]
a2 = fit2[[2]]
a3 = fit3[[2]]
a4 = fit4[[2]]

plot(cumsum(a0)/seq_along(a0), type='l', xlab='Iterations', ylab = 'Estimation', main='First moment')
lines(cumsum(a1)/seq_along(a1), col='blue')
lines(cumsum(a2)/seq_along(a2), col='green')
lines(cumsum(a3)/seq_along(a3), col='orange')
lines(cumsum(a4)/seq_along(a4), col='red')

a0 = c(fit0[[3]])
a1 = c(fit1[[3]])
a2 = c(fit2[[3]])
a3 = c(fit3[[3]])
a4 = c(fit4[[3]])

plot(cumsum(a0)/seq_along(a0), type='l', xlab='Iterations', ylab = 'Estimation', main='Second moment')
lines(cumsum(a1)/seq_along(a1), col='blue')
lines(cumsum(a2)/seq_along(a2), col='green')
lines(cumsum(a3)/seq_along(a3), col='orange')
lines(cumsum(a4)/seq_along(a4), col='red')
par(mfrow=c(1,1))

res = list()
res[['Montecarlo']] = replicate(100, sapply(expectedMonteCarlo(x, MU, SIGMA, nsim = NSIM, antithetic_variates = F, importance_sampling_mu = F), mean))
res[['M_AV']] = replicate(100, sapply(expectedMonteCarlo(x, MU, SIGMA, nsim = NSIM, antithetic_variates = T, importance_sampling_mu = F), mean))
res[['M_IS']] = replicate(100, sapply(expectedMonteCarlo(x, MU, SIGMA, nsim = NSIM, antithetic_variates = F, importance_sampling_mu = T), mean))
res[['M_AV_IS']] = replicate(100, sapply(expectedMonteCarlo(x, MU, SIGMA, nsim = NSIM, antithetic_variates = T, importance_sampling_mu = T), mean))
res[['Metropolis']] = replicate(100, sapply(expectedMetropolis(x, MU, SIGMA, nsim = NSIM), mean))

t(sapply(res, function(r) setNames(apply(r,1,mean)[2:3], c('M1','M2'))))
t(sapply(res, function(r) setNames(apply(r,1,sd)[2:3], c('M1','M2'))))



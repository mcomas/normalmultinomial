library(mixpack)
library(dplyr)
library(tidyr)

set.seed(1)
MU = 1
SIGMA = matrix(10)
x = c(10,0)

df = data_frame(
  a = seq(-3,10,0.1),
  f = exp(sapply(a, mvf, mu = MU, inv_sigma = solve(SIGMA), x = x))) %>%
  mutate(f = f / (sum(f)*(a[2]-a[1])))
plot(df, type='l')
abline(v = maximize_mvf(mu = MU, inv_sigma = solve(SIGMA), X = matrix(x, nrow=1)), col='red')
abline(v = mean(expectedA1(x = x, mu = MU, sigma = SIGMA, nsim=10000)), col = 'blue')


MU = c(1,1)
SIGMA = matrix(c(2,1,1,2), nrow=2)
x = c(10,0, 5)
a1 = seq(-3,3,0.1)
a2 = seq(-3,3,0.1)
A = expand.grid(a1,a2)
df = data_frame(
  a1 = A[,1],
  a2 = A[,2],
  f = exp(apply(A, 1, mvf, mu = MU, inv_sigma = solve(SIGMA), x = x))) %>%
  mutate(f = f / (sum(f)*0.1*0.1))
z = as.matrix(spread(df, key = a2, value = f) %>% select(-a1))

contour(a1, a2, z, xlim = c(-1,3), ylim = c(-3,1))
points(maximize_mvf(mu = MU, inv_sigma = solve(SIGMA), X = matrix(x, nrow=1)), col = 'red', pch = 19, cex=1)
points(matrix(apply(expectedA1(x = x, mu = MU, sigma = SIGMA, nsim=10000), 2, mean), nrow=1), col = 'blue', pch = 19, cex=1)
points(matrix(apply(expectedA1(x = x, mu = MU, sigma = SIGMA, nsim=10000), 2, mean), nrow=1), col = 'blue', pch = 19, cex=1)




vec_mvf_norm(matrix(c(0,0), nrow=1), mu = c(0,0), inv_sigma = solve(matrix(c(2,1,1,2), nrow=2)))








maximum = choose_starting(x = x, mu = MU, inv_sigma = solve(SIGMA), prop = 0.06)


library(mvtnorm)

omega = function(a) 1 + sum(exp(a))
omega(a[1,])^(-sum(x)) * exp(sum(A[1,] * x[1:ncol(A)]))

Ap = rmvnorm(10000, mean = MU, sigma = matrix(SIGMA))
A1 = factorial(sum(x))/prod(factorial(x)) *
  Ap *
  apply(Ap, 1, omega)^( -sum(x) ) *
  exp(apply(Ap * matrix(x[1:ncol(Ap)], ncol=2, nrow = nrow(Ap), byrow = T),1,sum))
A2 = expectedA(x = x, mu = MU, sigma = matrix(SIGMA), nsim=10000)
par(mfrow=c(1,2))
plot(A1, xlim = c(-1,1), ylim=c(-1, 3))
plot(A2, xlim = c(-1,1), ylim=c(-1, 3))
par(mfrow=c(1,1))



plot(A, asp=1)
points(t(maximum), col=2)

x = c(1,1,3)
A = matrix(c(0,0,1,2), ncol = 2, byrow = T)
exp(vec_mvf_multinom_mult(A,x))


omega = function(a) 1 + sum(exp(a))
omega(a[1,])^(-sum(x)) * exp(sum(A[1,] * x[1:ncol(A)]))

expectedA(x = c(0,20,0), mu = MU, sigma = SIGMA, nsim=10)

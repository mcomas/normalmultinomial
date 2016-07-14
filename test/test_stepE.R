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
abline(v = expectedA1(x = x, mu = MU, sigma = SIGMA, nsim=1000)[[1]], col = 'blue')
abline(v = expectedA2(x = x, mu = MU, sigma = SIGMA, nsim=1000)[[1]], col = 'green')
abline(v = expectedA3(x = x, mu = MU, sigma = SIGMA, nsim=1000)[[1]], col = 'brown')

sd(replicate(1000,expectedA1(x = x, mu = MU, sigma = SIGMA, nsim=1000)[[1]]))
sd(replicate(1000,expectedA2(x = x, mu = MU, sigma = SIGMA, nsim=1000)[[1]]))
sd(replicate(1000,expectedA3(x = x, mu = MU, sigma = SIGMA, nsim=1000)[[1]]))


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
points(t(expectedA1(x = x, mu = MU, sigma = SIGMA, nsim=1000)[[1]]), col = 'blue', pch = 19, cex=1)
points(t(expectedA2(x = x, mu = MU, sigma = SIGMA, nsim=1000)[[1]]), col = 'green', pch = 19, cex=1)
points(t(expectedA3(x = x, mu = MU, sigma = SIGMA, nsim=1000)[[1]]), col = 'brown', pch = 19, cex=1)

sqrt(det(cov(t(replicate(1000,c(expectedA1(x = x, mu = MU, sigma = SIGMA, nsim=1000)[[1]]))))))
sqrt(det(cov(t(replicate(1000,c(expectedA2(x = x, mu = MU, sigma = SIGMA, nsim=1000)[[1]]))))))
sqrt(det(cov(t(replicate(1000,c(expectedA3(x = x, mu = MU, sigma = SIGMA, nsim=1000)[[1]]))))))


###### E-step
set.seed(1)
library(mixpack)

(N <- 1000)
(MU <- c(0,0))

SIGMA.ilr <- matrix(c(1,0,
                      0,1), nrow = 2)
ilrB = t(clr_coordinates(t(matrix(unlist(ilr_basis(3)), ncol=2))))
ilr_to_alr = MASS::ginv(ilrB) %*% matrix(c(1,0,-1,
                                           0,1,-1), ncol=2)
SIGMA = t(ilr_to_alr) %*% SIGMA.ilr %*% ilr_to_alr

sample_nm = rnormalmultinomial(N, 5, mu = MU, sigma = SIGMA, probs = TRUE)

A = log(sample_nm$probs[,1:2]/sample_nm$probs[,3])
counts = sample_nm$counts

A.est = stepE(X = counts, mu = MU, sigma = SIGMA)

pars = list(MU, SIGMA)
( pars <- stepEM(X = counts, mu = pars[[1]], sigma = pars[[2]], nsim = 100))

P.est = as.data.frame(cbind(exp(stepE(X = counts, mu = pars[[1]], sigma = pars[[2]])), 1))
plot(ilr_coordinates(P.est), asp = 1)

ggtern() +
  geom_point(data = P.est, aes(x=V1,y=V2,z=V3)) +
  theme_bw()



library(randtoolbox)
library(normalmultinomial)
data("eleccatparl2015")

stepE = function(x, mu, sigma, nblock, err = 0.1){
  Z = matrix(torus(nblock, dim = nrow(sigma), normal = TRUE), ncol=nrow(sigma))
  fit0 = expectedMonteCarloFixed(x, mu, sigma, Z)
  E = list(n = nblock, x = apply(fit0[[2]], 2, sum), x2 = apply(fit0[[3]], 1:2, sum))
  eps = Inf
  while(eps > err){
    Z = matrix(torus(nblock, dim = nrow(sigma), normal = TRUE, init = FALSE), ncol=nrow(sigma))
    fit_new = expectedMonteCarloFixed(x, mu, sigma, Z)
    E_new = list(n = E[['n']] + nblock,
                 x = E[['x']] + apply(fit_new[[2]], 2, sum),
                 x2 = E[['x2']] + apply(fit_new[[3]], 1:2, sum))
    eps = base::norm(E[['x']] / E[['n']]- E_new[['x']] / E_new[['n']], '2')
    E = E_new
  }
  list('x' = E[['x']]/E[['n']], 'x2' = E[['x2']]/E[['n']], 'nsteps' = E[['n']])
}

MU = ilr_coordinates(matrix(apply(X, 2, mean), nrow=1))
SIGMA = diag(ncol(MU))

m_est = MU
s_est = SIGMA

res = apply(X, 1, stepE, m_est, s_est, 500, 0.1)

library(parallel)
cl <- makeCluster(getOption("cl.cores", 2))
clusterExport(cl, list('stepE', 'torus', 'expectedMonteCarloFixed'))

E1 = parApply(cl, X, 1, stepE, m_est, s_est, 500, 0.1)

EMsteps = 10
res.m = rep(0, EMsteps) #matrix(0, nrow=EMsteps, ncol = k)
res.s = rep(0, EMsteps) #matrix(0, nrow=EMsteps, ncol = k)
for(i in 1:EMsteps){
  E1 = parApply(cl, X, 1, stepE, m_est, s_est, 500, 0.1) #E2 = apply(X, 1, stepE, m_est, s_est, 2000, 0.01)

  m_est <- Reduce(`+`, lapply(E1, `[[`, 1)) / nrow(X)
  s_est <- nearestPSD(Reduce(`+`, lapply(E1, `[[`, 2)) / nrow(X) - t(t(m_est)) %*% t(m_est))

  res.m[i] = sqrt(sum((m_est-MU)^2))
  res.s[i] = norm(s_est-SIGMA, type = '2')

  print(i)
}
stopCluster(cl)


###############
#### TESTING
stepE(X[79,], MU, SIGMA, 500)


x = X[79,]
mu = fit$mu
sigma = fit$sigma
nblock = 500
err = 0.01

Z = matrix(torus(nblock, dim = nrow(sigma), normal = TRUE), ncol=nrow(sigma))
fit0 = expectedMonteCarloFixed(x, mu, sigma, Z)
E = list(n = nblock, x = apply(fit0[[2]], 2, sum), x2 = apply(fit0[[3]], 1:2, sum))
eps = Inf
while(eps > err){
  Z = matrix(torus(nblock, dim = nrow(sigma), normal = TRUE, init = FALSE), ncol=nrow(sigma))
  fit_new = expectedMonteCarloFixed(x, mu, sigma, Z)
  E_new = list(n = E[['n']] + nblock,
               x = E[['x']] + apply(fit_new[[2]], 2, sum),
               x2 = E[['x2']] + apply(fit_new[[3]], 1:2, sum))
  eps = base::norm(E[['x']] / E[['n']]- E_new[['x']] / E_new[['n']], '2')
  E = E_new
}
list('x' = E[['x']]/E[['n']], 'x2' = E[['x2']]/E[['n']], 'nsteps' = E[['n']])

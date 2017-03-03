SEED = 10
D = 10
NSIM = 10000
SIZE = 1000

set.seed(SEED)

library(normalmultinomial)
library(randtoolbox)

#x = as.numeric(rmultinomial(1, D*10, 1:(D+1)/sum(1:(D+1))))
MU = ilr_coordinates(1:(D+1))
SIGMA = diag(1, D)

SIM = rnormalmultinomial(1, size = SIZE, MU, SIGMA, probs = TRUE)

x = SIM$counts[1,]
p.real = SIM$probs[1,]

MU_EXP = MU

Z = halton(NSIM, dim = D, normal = TRUE)
Z = sobol(NSIM, dim = D, normal = TRUE)

runMetropolisUSING_mcmcPackage = function(){
  loglik = function(z, x, mu, sigma){
    p = inv_ilr_coordinates(z)
    dmultinom(x, sum(x), p, log = TRUE)  + dmvnorm(z, mean = mu, sigma = sigma, log = TRUE)
  }
  out = metrop(loglik, rep(1,D), nbatch = 5000, blen = 10, x = x, mu = MU, sigma = SIGMA)
  apply(out$batch, 2, mean)
}
runMetropolis = function(burnin = 0){
  RES = expectedMetropolis(x, MU, SIGMA,
                           Z = matrix(rnormal(NSIM+burnin, mu = rep(0, D), sigma = diag(1, D)),
                                      ncol=D), mu_exp = rep(-1, D))
  if(burnin> 0){
    RES[[1]] = RES[[1]][-(1:burnin)]
    RES[[2]] = RES[[2]][-(1:burnin),]
    RES[[3]] = RES[[3]][,,-(1:burnin)]
  }
  M1 = colMeans(RES[[2]])
  M2 = apply(RES[[3]], 1:2, mean)
  list(M1, M2)
}

runMontecarlo = function(type, seed = 4711){
  if(type == 'pseudo'){
    Z = matrix(rnormal(NSIM, mu = rep(0, D), sigma = diag(1, D)), ncol=D)
  }
  if(type == 'halton'){
    Z = matrix(halton(NSIM, dim = D, normal = TRUE, usetime = T, init = FALSE), ncol=D)
  }
  if(type == 'sobol'){
    Z = matrix(sobol(NSIM, dim = D, normal = TRUE, scrambling = 1, seed = seed), ncol=D)
  }
  RES = expectedMonteCarlo(x, MU, SIGMA,
                           Z = Z,
                           mu_exp = MU_EXP)
  M1 = colMeans(RES[[2]])
  M2 = apply(RES[[3]], 1:2, mean)
  list(M1, M2)
}


res = list()
set.seed(SEED)
res[['Metropolis']] = replicate(50, runMetropolis(burnin = 500), simplify = FALSE)

set.seed(SEED)
res[['Pseudo generation']] = replicate(50, runMontecarlo('pseudo'), simplify = FALSE)

set.seed(SEED)
res[['Quasi-MC (Halton)']] = replicate(50, runMontecarlo('halton'), simplify = FALSE)

set.seed(SEED)
res[['Quasi-MC (Sobol)']] = sapply(1:50, function(seed)
  runMontecarlo('sobol', seed), simplify = FALSE)


M1.est = sapply(res, function(ires){
  apply(sapply(ires, dplyr::first), 1, mean, na.rm=T)
})

M1.sd = sapply(res, function(ires){
  apply(sapply(ires, dplyr::first), 1, sd, na.rm=T)
})

M2.est = lapply(res, function(ires){
  apply(array(sapply(ires, dplyr::last), c(D, D, NSIM)), 1:2, mean, na.rm=TRUE)
})

M2.sd = lapply(res, function(ires){
  apply(array(sapply(ires, dplyr::last), c(D, D, NSIM)), 1:2, sd, na.rm=TRUE)
})

#M1.est
dist(t(cbind(M1.est, 'Real' = ilr_coordinates(p.real))))

#M1.est
apply(M1.sd, 2, range)

#dist(t(M1.est))

plot(ilr_coordinates(p.real), ylim=c(-2,2))
points(M1.est[,1], col = 2)
points(M1.est[,2], col = 3)
points(M1.est[,3], col = 4)
points(M1.est[,4], col = 5)

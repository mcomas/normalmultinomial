evaluate_center = function(results){
  sapply(results, function(res){
    sqrt(sum((ilr_coordinates(res$p) - res$meanilr)^2))
  })
}

# RES.NM0 = lapply(RES.NM, lapply, first)
#
# tab.GBM = sapply(RES.GBM[1:6], evaluate_center)
# tab.NM = sapply(RES.NM0[1:6], evaluate_center)
#
# mean(tab.GBM)
# mean(tab.NM)
#
# res = RES.NM0[[6]][[1]]

n = 50
p = c(0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
      0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02,
      0.05, 0.05, 0.05, 0.05, 0.5)
N = 1000

set.seed(1)
cat(sprintf("%s\nN=%d\nn=%d\np=(%s)\n----\n", date(), N, n, paste(p, collapse=', ')))
XZ = t(stats::rmultinom(N, n, p))

XnZ.GBM = zCompositions::cmultRepl(XZ, method = 'GBM')
View(XnZ.GBM)


sigma = diag(ncol(X)-1)
eps = 0.001
nsim = 10000
parallel.cluster = NULL
max.em.iter = 100

MU0 = MU = ilr_coordinates(apply(XZ/apply(XZ, 1, sum), 2, mean))
SIGMA = sigma

if(nrow(SIGMA) <= 6){
  Z = matrix(randtoolbox::halton(nsim, dim = nrow(SIGMA), normal = TRUE), ncol=nrow(SIGMA))
}else{
  Z = matrix(randtoolbox::sobol(nsim, dim = nrow(SIGMA), normal = TRUE), ncol=nrow(SIGMA))
}

X = XZ

err = eps + 1
iter = 0
while(err > eps & iter < max.em.iter){
  if(!is.null(parallel.cluster)){
    FIT = parallel::parApply(parallel.cluster, X, 1, expectedMonteCarlo, MU, SIGMA, Z)
  }else{
    FIT = apply(X, 1, expectedMonteCarlo, MU, SIGMA, Z)
  }
  E = t(sapply(FIT, function(fit) apply(fit[[2]], 2, sum) / nsim))

  delta = (apply(E, 2, mean)-MU)
  err = sqrt(sum(delta^2))
  MU = MU + delta

  SIGMA = Reduce(`+`, lapply(FIT, function(fit) apply(fit[[3]], 1:2, sum)/ nsim)) / nrow(X) - MU %*% t(MU)
  iter = iter + 1
}
nm_fit(X, nsim = 10, eps = 0.001)

XnZ.NM = inv_ilr_coordinates(E)


p.NM = apply(XnZ.NM / apply(XnZ.NM, 1, sum), 2, mean)
p.GBM = apply(XnZ.GBM / apply(XnZ.GBM, 1, sum), 2, mean)



set.seed(1)
cat(sprintf("%s\nN=%d\nn=%d\np=(%s)\n----\n", date(), N, n, paste(p, collapse=', ')))
X = t(stats::rmultinom(N, n, p))
apply(X / apply(X, 1, sum), 2, mean)

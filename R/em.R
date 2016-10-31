
#' @export
nm_fit_mean = function(X, sigma = diag(ncol(X)-1), eps = 0.001, nsim = 1000, parallel.cluster = NULL,
                       max.em.iter = 100){
  D = ncol(X)

  ILR.TO.ILR = t(ilr_to_alr(D))
  inv.ILR.TO.ILR = solve(ILR.TO.ILR)

  MU = ilr_coordinates(matrix(apply(X/apply(X, 1, sum), 2, sum), nrow=1)) %*% ILR.TO.ILR
  SIGMA = t(ILR.TO.ILR) %*% sigma %*% ILR.TO.ILR

  if(nrow(SIGMA) <= 6){
    Z = matrix(randtoolbox::halton(nsim, dim = nrow(SIGMA), normal = TRUE), ncol=nrow(SIGMA))
  }else{
    Z = matrix(randtoolbox::sobol(nsim, dim = nrow(SIGMA), normal = TRUE), ncol=nrow(SIGMA))
  }

  err = eps + 1
  iter = 0
  while(err > eps & iter < max.em.iter){
    if(!is.null(parallel.cluster)){
      FIT = parallel::parApply(parallel.cluster, X, 1, expectedMoment1_alr, MU, SIGMA, Z)
    }else{
      FIT = apply(X, 1, expectedMoment1_alr, MU, SIGMA, Z)
    }
    E = t(sapply(FIT, function(fit) colMeans(na.omit(fit[[2]]))))

    MU.new = apply(E, 2, mean)
    delta = (MU.new %*% inv.ILR.TO.ILR - MU %*% inv.ILR.TO.ILR)
    err = sqrt(sum(delta^2))
    MU = MU.new
    iter = iter + 1
  }

  list(mu = MU %*% inv.ILR.TO.ILR,
       sigma = t(inv.ILR.TO.ILR) %*% SIGMA %*% inv.ILR.TO.ILR,
       iter = iter)
}

#' @export
nm_fit_spherical = function(X, sigma = diag(ncol(X)-1), eps = 0.001, nsim = 1000, parallel.cluster = NULL,
                            max.em.iter = 100){
  D = ncol(X)

  MU = ilr_coordinates(matrix(apply(X/apply(X, 1, sum), 2, sum), nrow=1))[1,]
  SIGMA = sigma

  if(nrow(SIGMA) <= 6){
    Z = matrix(randtoolbox::halton(nsim, dim = nrow(SIGMA), normal = TRUE), ncol=nrow(SIGMA))
  }else{
    Z = matrix(randtoolbox::sobol(nsim, dim = nrow(SIGMA), normal = TRUE), ncol=nrow(SIGMA))
  }

  err = eps + 1
  iter = 0
  while(err > eps & iter < max.em.iter){
    if(!is.null(parallel.cluster)){
      FIT = parallel::parApply(parallel.cluster, X, 1, expectedMonteCarlo, MU, SIGMA, Z)
    }else{
      FIT = apply(X, 1, expectedMonteCarlo, MU, SIGMA, Z)
    }
    E = t(sapply(FIT, function(fit) colMeans(na.omit(fit[[2]]))))

    delta = (apply(E, 2, mean)-MU)
    err = sqrt(sum(delta^2))
    MU = MU + delta

    SIGMA = diag(sum(diag(Reduce(`+`, lapply(FIT, function(fit) apply(fit[[3]], 1:2, sum)/ nsim)) / nrow(X) - MU %*% t(MU))), D-1)
    iter = iter + 1
  }

  list(mu = MU,
       sigma = SIGMA,
       iter = iter)
}

#' @export
nm_fit = function(X, sigma = diag(ncol(X)-1), eps = 0.001, nsim = 1000, parallel.cluster = NULL,
                  max.em.iter = 100){
  MU = ilr_coordinates(matrix(apply(X/apply(X, 1, sum), 2, sum), nrow=1))[1,]
  SIGMA = sigma

  if(nrow(SIGMA) <= 6){
    Z = matrix(randtoolbox::halton(nsim, dim = nrow(SIGMA), normal = TRUE), ncol=nrow(SIGMA))
  }else{
    Z = matrix(randtoolbox::sobol(nsim, dim = nrow(SIGMA), normal = TRUE), ncol=nrow(SIGMA))
  }

  err = eps + 1
  iter = 0
  while(err > eps & iter < max.em.iter){
    if(!is.null(parallel.cluster)){
      FIT = parallel::parApply(parallel.cluster, X, 1, expectedMonteCarlo, MU, SIGMA, Z)
    }else{
      FIT = apply(X, 1, expectedMonteCarlo, MU, SIGMA, Z)
    }
    E = t(sapply(FIT, function(fit) colMeans(na.omit(fit[[2]]))))

    delta = (apply(E, 2, mean)-MU)
    err = sqrt(sum(delta^2))
    MU = MU + delta
    #SIGMA = Reduce(`+`, lapply(FIT, function(fit) apply(fit[[3]], 1:2, function(v) mean(na.omit(v), trim=0.001)))) / nrow(X) - MU %*% t(MU)
    #SIGMA = (SIGMA + t(SIGMA))/2
    SIGMA = Reduce(`+`, lapply(FIT, function(fit) apply(fit[[3]], 1:2, sum)/ nsim)) / nrow(X) - MU %*% t(MU)
    iter = iter + 1
  }

  list(mu = MU,
       sigma = SIGMA,
       iter = iter)
}

#' @export
nm_expected = function(X, mu, sigma, nsim = 1000, parallel.cluster = NULL){
  if(nrow(sigma) <= 6){
    Z = matrix(randtoolbox::halton(nsim, dim = nrow(sigma), normal = TRUE), ncol=nrow(sigma))
  }else{
    Z = matrix(randtoolbox::sobol(nsim, dim = nrow(sigma), normal = TRUE), ncol=nrow(sigma))
  }
  if(!is.null(parallel.cluster)){
    FIT = parallel::parApply(parallel.cluster, X, 1, expectedMoment1, mu, sigma, Z)
  }else{
    FIT = apply(X, 1, expectedMoment1, mu, sigma, Z)
  }
  inv_ilr_coordinates(t(sapply(FIT, function(fit) colMeans(na.omit(fit[[2]])))))
}

#' @export
nm_fit = function(X, sigma = diag(ncol(X)-1), eps = 0.001, nsim = 1000, parallel.cluster = NULL){
  MU = ilr_coordinates(matrix(apply(X/apply(X, 1, sum), 2, sum), nrow=1))[1,]
  SIGMA = sigma

  Z = matrix(randtoolbox::torus(nsim, dim = nrow(SIGMA), normal = TRUE), ncol=nrow(SIGMA))
  err = eps + 1
  while(err > eps){
    if(!is.null(parallel.cluster)){
      FIT = parallel::parApply(parallel.cluster, X, 1, expectedMoment1, MU, SIGMA, Z)
    }else{
      FIT = apply(X, 1, expectedMoment1, MU, SIGMA, Z)
    }
    E = t(sapply(FIT, function(fit) apply(fit[[2]], 2, sum) / nsim))

    delta = (apply(E, 2, mean)-MU)
    err = sqrt(sum(delta^2))
    MU = MU + delta
  }

  list(mu = MU,
       sigma = SIGMA)
}

#' @export
nm_expected = function(X, mu, sigma, nsim = 1000, parallel.cluster = NULL){
  Z = matrix(randtoolbox::torus(nsim, dim = nrow(sigma), normal = TRUE), ncol=nrow(sigma))
  if(!is.null(parallel.cluster)){
    FIT = paralell::parApply(parallel.cluster, X, 1, expectedMoment1, mu, sigma, Z)
  }else{
    FIT = apply(X, 1, expectedMoment1, mu, sigma, Z)
  }
  inv_ilr_coordinates(t(sapply(FIT, function(fit) apply(fit[[2]], 2, sum) / nsim)))
}

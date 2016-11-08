

###################################
### OTHER FUNCTIONS
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
    E = t(sapply(FIT, function(fit) colMeans(stats::na.omit(fit[[2]]))))

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
    E = t(sapply(FIT, function(fit) colMeans(stats::na.omit(fit[[2]]))))

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



nm_fit2 = function(X, sigma = diag(ncol(X)-1), eps = 0.001, nsim = 1000, parallel.cluster = NULL,
                  max.em.iter = 100, verbose = FALSE){
  MU = ilr_coordinates(colSums(X/rowSums(X)))
  SIGMA = sigma

  MULTINOM_CONS = apply(X, 1, mvf_multinom_const)

  if(nrow(SIGMA) <= 6){
    Z = matrix(randtoolbox::halton(nsim, dim = nrow(SIGMA), normal = TRUE), ncol=nrow(SIGMA))
  }else{
    Z = matrix(randtoolbox::sobol(nsim, dim = nrow(SIGMA), normal = TRUE), ncol=nrow(SIGMA))
  }

  if(!is.null(parallel.cluster)){
    FIT = parallel::parApply(parallel.cluster, X, 1, expectedMonteCarlo2, MU, SIGMA, Z)
  }else{
    FIT = apply(X, 1, expectedMonteCarlo2, MU, SIGMA, Z)
  }
  E = t(sapply(FIT, function(fit) colMeans(stats::na.omit(fit[[2]]))))

  err = eps + 1
  iter = 0
  while(err > eps & iter < max.em.iter){
    if(!is.null(parallel.cluster)){
      FIT = parallel::parSapply(parallel.cluster, 1:nrow(X),
                                function(i)
                                  expectedMonteCarlo3(X[i,], MU, SIGMA, Z, E[i,]), simplify = FALSE)
    }else{
      FIT = sapply(1:nrow(X), function(i) expectedMonteCarlo3(X[i,], MU, SIGMA, Z, E[i,]), simplify = FALSE)
    }
    E = t(sapply(FIT, function(fit) colMeans(stats::na.omit(fit[[2]]))))

    delta = (apply(E, 2, mean)-MU)
    err = sqrt(sum(delta^2))
    MU = MU + delta
    L = lapply(FIT, function(fit){
      sel = apply(fit[[3]], 3, function(m) all(is.finite(m)))
      apply(fit[[3]][,,sel], 1:2, sum) / length(sel)
    })
    SIGMA = Reduce(`+`, L) / nrow(X) - MU %*% t(MU)
    iter = iter + 1
    if(verbose){ cat(sprintf('Step %d, error %f\n', iter, err)) }
  }

  list(mu = MU,
       sigma = SIGMA,
       iter = iter)
}


nm_fit3 = function(X, sigma = diag(ncol(X)-1), eps = 0.001, nsim = 1000, parallel.cluster = NULL,
                   max.em.iter = 100, verbose = FALSE){
  MU = ilr_coordinates(colSums(X/rowSums(X)))
  SIGMA = sigma

  MULTINOM_CONS = apply(X, 1, mvf_multinom_const)

  if(nrow(SIGMA) <= 6){
    Z = matrix(randtoolbox::halton(nsim, dim = nrow(SIGMA), normal = TRUE), ncol=nrow(SIGMA))
  }else{
    Z = matrix(randtoolbox::sobol(nsim, dim = nrow(SIGMA), normal = TRUE), ncol=nrow(SIGMA))
  }

  if(!is.null(parallel.cluster)){
    FIT = parallel::parApply(parallel.cluster, X, 1, expectedMonteCarlo2, MU, SIGMA, Z)
  }else{
    FIT = apply(X, 1, expectedMonteCarlo2, MU, SIGMA, Z)
  }
  E1 = t(sapply(FIT, function(fit) colMeans(stats::na.omit(fit[[2]]))))
  E2 = sapply(FIT, function(fit) margin.table(fit[[3]], 1:2)/nsim, simplify = FALSE)

  err = eps + 1
  iter = 0
  while(err > eps & iter < max.em.iter){
    if(!is.null(parallel.cluster)){
      FIT = parallel::parSapply(parallel.cluster, 1:nrow(X),
                                function(i)
                                  expectedMonteCarlo4(X[i,], MU, SIGMA, Z, E1[i,], diag(mean(diag(E2[[i]] - E1[i,] %*% t(E1[i,])), na.rm =T), ncol(E1))), simplify = FALSE)
    }else{
      FIT = sapply(1:nrow(X),
                   function(i)
                     expectedMonteCarlo4(X[i,], MU, SIGMA, Z, E1[i,], diag(mean(diag(E2[[i]] - E1[i,] %*% t(E1[i,])), na.rm =T), ncol(E1))), simplify = FALSE)
    }
    E1 = t(sapply(FIT, function(fit) colMeans(stats::na.omit(fit[[2]]))))
    E2 = sapply(FIT, function(fit) margin.table(fit[[3]], 1:2)/nsim, simplify = FALSE)

    delta = (apply(E1, 2, mean)-MU)
    err = sqrt(sum(delta^2))
    MU = MU + delta
    L = lapply(FIT, function(fit){
      sel = apply(fit[[3]], 3, function(m) all(is.finite(m)))
      apply(fit[[3]][,,sel], 1:2, sum) / length(sel)
    })
    SIGMA = Reduce(`+`, L) / nrow(X) - MU %*% t(MU)
    iter = iter + 1
    if(verbose){ cat(sprintf('Step %d, error %f\n', iter, err)) }
  }

  list(mu = MU,
       sigma = SIGMA,
       iter = iter)
}




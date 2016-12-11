generate_mv_normal_rnd = function(n, dim){
  if(dim <= 6){
    Z = matrix(randtoolbox::halton(n, dim = dim, normal = TRUE), ncol=dim)
  }else{
    Z = matrix(randtoolbox::sobol(n, dim = dim, normal = TRUE), ncol=dim)
  }
  return(Z)
}

initialize_with_bootstrap2 = function(X){
  B = t(replicate(1000, colMeans(X[sample(1:nrow(X),nrow(X), replace=TRUE),])))

  MU = ilr_coordinates(colMeans(X))
  SIGMA = nrow(X) * cov(ilr_coordinates(B))

  H = matrix(MU, nrow = nrow(X), ncol = ncol(X)-1, byrow = TRUE)

  list(H = H, MU = MU, SIGMA = SIGMA)
}

initialize_with_bootstrap3 = function(X){
  B = t(replicate(1000, colMeans(X[sample(1:nrow(X),nrow(X), replace=TRUE),])))

  MU = ilr_coordinates(colMeans(X))
  SIGMA = diag(mean(eigen(nrow(X) * cov(ilr_coordinates(B)))$values), ncol(X) - 1)

  H = matrix(MU, nrow = nrow(X), ncol = ncol(X)-1, byrow = TRUE)

  list(H = H, MU = MU, SIGMA = SIGMA)
}

initialize_with_bootstrap4 = function(X){
  B = t(replicate(1000, colMeans(X[sample(1:nrow(X),nrow(X), replace=TRUE),])))

  MU = ilr_coordinates(colMeans(X))
  SIGMA = diag(exp(mean(log((eigen(nrow(X) * cov(ilr_coordinates(B)))$values)))), ncol(X) - 1)

  H = matrix(MU, nrow = nrow(X), ncol = ncol(X)-1, byrow = TRUE)

  list(H = H, MU = MU, SIGMA = SIGMA)
}

initialize_with_bootstrap = function(X, parallel.cluster){
  B = t(replicate(1000, colMeans(X[sample(1:nrow(X),nrow(X), replace=TRUE),])))

  MU = ilr_coordinates(colMeans(X))
  SIGMA = nrow(X) * cov(ilr_coordinates(B))

  if(nrow(SIGMA) <= 6){
    Z = matrix(randtoolbox::halton(1000, dim = nrow(SIGMA), normal = TRUE), ncol=nrow(SIGMA))
  }else{
    Z = matrix(randtoolbox::sobol(1000, dim = nrow(SIGMA), normal = TRUE), ncol=nrow(SIGMA))
  }
  if(!is.null(parallel.cluster)){
    FIT = parallel::parApply(parallel.cluster, X, 1, expectedMonteCarlo2, MU, SIGMA, Z)
  }else{
    FIT = apply(X, 1, expectedMonteCarlo2, MU, SIGMA, Z)
  }
  H = t(sapply(FIT, function(fit) colMeans(stats::na.omit(fit[[2]]))))

  list(H = H, MU = MU, SIGMA = SIGMA)
}

initialize_with_bootstrapstep = function(X, steps = 1, parallel.cluster){
  B = t(replicate(1000, colMeans(X[sample(1:nrow(X),nrow(X), replace=TRUE),])))

  MU = ilr_coordinates(colMeans(X))
  SIGMA = nrow(X) * cov(ilr_coordinates(B))

  if(nrow(SIGMA) <= 6){
    Z = matrix(randtoolbox::halton(1000, dim = nrow(SIGMA), normal = TRUE), ncol=nrow(SIGMA))
  }else{
    Z = matrix(randtoolbox::sobol(1000, dim = nrow(SIGMA), normal = TRUE), ncol=nrow(SIGMA))
  }
  for(s in 1:steps){
    if(!is.null(parallel.cluster)){
      FIT = parallel::parApply(parallel.cluster, X, 1, expectedMonteCarlo2, MU, SIGMA, Z)
    }else{
      FIT = apply(X, 1, expectedMonteCarlo2, MU, SIGMA, Z)
    }
    H = t(sapply(FIT, function(fit) colMeans(stats::na.omit(fit[[2]]))))
    SIGMA = cov(H)
  }
  list(H = H, MU = MU, SIGMA = SIGMA)
}


initialize_with_dm = function(X){
  E = dm_fit(X)$expected

  H = ilr_coordinates(E)
  MU = colMeans(H)
  SIGMA = cov(H)
  list(H = H, MU = MU, SIGMA = SIGMA)
}

initialize_with_dmean = function(X){
  E = X + matrix(colMeans(X), ncol=ncol(X), nrow=nrow(X), byrow = TRUE)

  H = ilr_coordinates(E)
  MU = colMeans(H)
  SIGMA = cov(H)
  list(H = H, MU = MU, SIGMA = SIGMA)
}

initialize_with_dsum = function(X){
  E = X + matrix(colSums(X), ncol=ncol(X), nrow=nrow(X), byrow = TRUE)

  H = ilr_coordinates(E)
  MU = colMeans(H)
  SIGMA = cov(H)
  list(H = H, MU = MU, SIGMA = SIGMA)
}

initialize_with_dapprox = function(X){
  X.closured = X / rowSums(X)
  m = colMeans(X.closured)
  v = apply(X.closured, 2, var)
  K = ncol(X)
  E = X + matrix(m * exp(1/(K-1) * sum(log(m*(1-m)/v-1))), ncol=ncol(X), nrow=nrow(X), byrow = TRUE)

  H = ilr_coordinates(E)
  MU = colMeans(H)
  SIGMA = cov(H)
  list(H = H, MU = MU, SIGMA = SIGMA)
}

initialize_with_dapprox2 = function(X){
  XNZ = X + matrix(colMeans(X), ncol=ncol(X), nrow=nrow(X), byrow = TRUE)
  X.closured = XNZ / rowSums(XNZ)
  m = colMeans(X.closured)
  v = apply(X.closured, 2, var)
  K = ncol(X)
  E = X + matrix(m * exp(1/(K-1) * sum(log(m*(1-m)/v-1))), ncol=ncol(X), nrow=nrow(X), byrow = TRUE)

  H = ilr_coordinates(E)
  MU = colMeans(H)
  SIGMA = cov(H)
  list(H = H, MU = MU, SIGMA = SIGMA)
}

initialize_with_aitchison = function(X){
  P = colMeans(X)
  P = P / sum(P)
  MU = ilr_coordinates(P)
  B = t(replicate(1000, colMeans(X[sample(1:nrow(X),nrow(X), replace=TRUE),])))
  SIGMA = nrow(X) * cov(ilr_coordinates(B))


  G = P
  K = cov(X/rowSums(X))

  TAU = matrix(0, nrow=nrow(K), ncol=ncol(K))
  tau = function(i,j){
    K[i,i]/G[i]^2 - 2 * K[i,j]/(G[i]*G[j]) + K[j,j]/G[j]^2 + 1/4 * (K[i,i]/G[i]^2 - K[j,j]/G[j]^2)^2
  }
  for(i in 1:nrow(K)){for(j in 1:ncol(K)){ TAU[i,j] = tau(i,j) }}
  Fn = cbind(diag(1,nrow(K)-1), -1)
  S = -0.5 * Fn %*% TAU %*% t(Fn)
  SIGMA = solve(ilr_to_alr(nrow(K))) %*% S %*% t(solve(ilr_to_alr(nrow(K))))

  H = matrix(MU, nrow = nrow(X), ncol = ncol(X)-1, byrow = TRUE)

  list(H = H, MU = MU, SIGMA = SIGMA)
}
#'
#' Log-ratio normal-multinomial parameters estimation.
#'
#' The parameters mu and sigma are expressed with respect
#' basis given by function @ilrBase.
#'
#' @param X count data set
#' @param eps maximum error accepted on the last EM-step
#' @param nsim number of simulations used in the E-step
#' @param parallel.cluster parallel Socket Cluster created with function @makeCluster
#' @param max.em.iter maximum number of steps allowed in the EM-algorithm
#' @param expected if TRUE the expected probabilities are returned (default:TRUE)
#' @param verbose show information during estimation
#' @return A list with parameters mu and sigma and the number of iterations before convergence
#' @examples
#' X = rnormalmultinomial(100, 100, rep(0,4), diag(4))
#' nm_fit(X, verbose = T)
#' @export
nm_fit = function(X, eps = 0.001, nsim = 1000, parallel.cluster = NULL,
                  max.em.iter = 100, expected = TRUE, verbose = FALSE,
                  delta = 0.65, threshold = 0.5, init.method = 'dm',
                  development = FALSE){

  if(ncol(X) == 2){
    return( nm_fit_1d(X, eps, parallel.cluster, max.em.iter, expected, verbose) )
  }

  if(init.method == 'dm'){
    init = initialize_with_dm(X)
  }
  if(init.method == 'mean'){
    init = initialize_with_dmean(X)
  }
  if(init.method == 'sum'){
    init = initialize_with_dsum(X)
  }
  if(init.method == 'approx'){
    init = initialize_with_dapprox(X)
  }
  if(init.method == 'approx2'){
    init = initialize_with_dapprox2(X)
  }
  if(init.method == 'bootstrap'){
    init = initialize_with_bootstrap(X, parallel.cluster)
  }
  if(init.method == 'bootstrap2'){
    init = initialize_with_bootstrap2(X)
  }
  if(init.method == 'bootstrap3'){
    init = initialize_with_bootstrap3(X)
  }
  if(init.method == 'bootstrap4'){
    init = initialize_with_bootstrap4(X)
  }
  if(init.method == 'bootstrapstep'){
    init = initialize_with_bootstrapstep(X, 1, parallel.cluster)
  }
  if(init.method == 'bootstrapstep5'){
    init = initialize_with_bootstrapstep(X, 5, parallel.cluster)
  }
  if(init.method == 'aitchison'){
    init = initialize_with_aitchison(X)
  }
  E = init$H
  MU = init$MU
  SIGMA = init$SIGMA


  ##
  ##
  Z = generate_mv_normal_rnd(nsim, nrow(SIGMA))

  err_prev = err = eps + 1
  iter = 0
  while(err > eps & iter < max.em.iter){
    if(!is.null(parallel.cluster)){
      FIT = parallel::parSapply(parallel.cluster, 1:nrow(X),
                                function(i)
                                  expectedMonteCarlo(X[i,], MU, SIGMA, Z, E[i,]), simplify = FALSE)
      E = t(parallel::parSapply(parallel.cluster, FIT, function(fit) colMeans(stats::na.omit(fit[[2]]))))
      # L = parallel::parLapply(parallel.cluster, FIT, function(fit){
      #   sel = apply(fit[[3]], 3, function(m) all(is.finite(m)))
      #   apply(fit[[3]][,,sel], 1:2, sum) / length(sel)
      # })
      L = parallel::parLapply(parallel.cluster, FIT, function(fit){
        apply(fit[[3]], 1:2, mean)
      })
    }else{
      FIT = sapply(1:nrow(X), function(i) expectedMonteCarlo(X[i,], MU, SIGMA, Z, E[i,]), simplify = FALSE)
      E = t(sapply(FIT, function(fit) colMeans(stats::na.omit(fit[[2]]))))
      # L = lapply(FIT, function(fit){
      #   sel = apply(fit[[3]], 3, function(m) all(is.finite(m)))
      #   apply(fit[[3]][,,sel], 1:2, sum) / length(sel)
      # })
      L = lapply(FIT, function(fit){
        apply(fit[[3]], 1:2, mean)
      })
    }

    delta = (apply(E, 2, mean)-MU)
    err = sqrt(sum(delta^2))

    iter = iter + 1
    if(verbose | development){ cat(sprintf('Step %d, error %f\n', iter, err)) }

    if(err_prev < err){
      if(verbose | development){ cat(sprintf('nsim: %d\n', 2*nsim)) }
      nsim = 2 * nsim
      Z = generate_mv_normal_rnd(nsim, nrow(SIGMA))
    }
    err_prev = err

    MU = MU + delta
    SIGMA = Reduce(`+`, L) / nrow(X) - MU %*% t(MU)

  }
  if(development){
    return(
      list(
        mu = MU,
        sigma = SIGMA,
        iter = iter,
        expected = inv_ilr_coordinates(E),
        fit = FIT
      )
    )
  }
  if(expected){
    list(mu = MU,
         sigma = SIGMA,
         iter = iter,
         expected = inv_ilr_coordinates(E))
  }else{
    list(mu = MU,
         sigma = SIGMA,
         iter = iter)
  }
}

f.prob = function(n, k, mu, rho){
  K = lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1)
  function(x) sqrt(rho/pi) * exp(-rho * mu^2/2) * exp(K - rho*x^2 + (sqrt(2) * mu * rho + 2 * k - n) * x - n * log(exp(x)+exp(-x)))
}
f.m1 = function(n, k, mu, rho){
  K = lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1)
  function(x) sqrt(2) * x * sqrt(rho/pi) * exp(-rho * mu^2/2) * exp(K - rho*x^2 + (sqrt(2) * mu * rho + 2 * k - n) * x - n * log(exp(x)+exp(-x)))
}
f.m2 = function(n, k, mu, rho){
  K = lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1)
  function(x) 2 * x^2 * sqrt(rho/pi) * exp(-rho * mu^2/2) * exp(K - rho*x^2 + (sqrt(2) * mu * rho + 2 * k - n) * x - n * log(exp(x)+exp(-x)))
}

nm_fit_1d = function(X, eps = 0.00001, parallel.cluster = NULL,
                     max.em.iter = 100, expected = TRUE, verbose = FALSE){

  MU = ilr_coordinates(matrix(apply(X/apply(X, 1, sum), 2, sum), nrow=1))[1,]
  SIGMA = 1

  N = rowSums(X)

  err = eps + 1
  iter = 0
  while(err > eps & iter < max.em.iter){
    probs = sapply(1:nrow(X), function(i) integrate(f.prob(N[i], X[i,1], MU, 1/SIGMA), -Inf, Inf, rel.tol = 10^-8)$value)
    E = sapply(1:nrow(X), function(i) integrate(f.m1(N[i], X[i,1], MU, 1/SIGMA),
                                                -Inf, Inf, rel.tol = 10^-8)$value) / probs
    M2 = sapply(1:nrow(X), function(i) integrate(f.m2(N[i], X[i,1], MU, 1/SIGMA), -Inf, Inf)$value) / probs

    delta = (mean(E)-MU)
    err = sqrt(sum(delta^2))

    MU = MU + delta


    SIGMA = mean(M2) - MU * MU
    iter = iter + 1
    print(c('MU' = MU, 'SGIMA' = SIGMA))
    if(verbose){ cat(sprintf('Step %d, error %f\n', iter, err)) }
  }
  if(expected){
    list(mu = unname(MU),
         sigma = unname(SIGMA),
         iter = iter,
         expected = inv_ilr_coordinates(matrix(E, ncol=1)))
  }else{
    list(mu = MU,
         sigma = SIGMA,
         iter = iter)
  }
}

#'
#' Calculate the expected probabilities
#'
#' Calculate the probabilities conditioned on parameters mu and sigma
#' more likely to generate dataset X
#'
#' @param X count data set
#' @param mu paramater mu
#' @param sigma parameter sigma
#' @param nsim number of simulation used to calculated the expected values
#' @param parallel.cluster instance of a Socket Cluster (see function makeCluster from package parallel)
#' @return the expected probabilities
#' @examples
#' X = rnormalmultinomial(100, 100, rep(0,4), diag(4))
#' fit = nm_fit(X, verbose = T)
#' P = nm_expected(X, fit$mu, fit$sigma)
#' head(P)
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
  inv_ilr_coordinates(t(sapply(FIT, function(fit) colMeans(stats::na.omit(fit[[2]])))))
}

#'
#' Dirichlet-multinomial parameter estimation.
#'
#' @param X count data set
#' @param eps maximum error accepted for alpha parameters
#' @param max.iter maximum number of steps allowed
#' @return A list with parameter alpha and the number of iterations before convergence
#' @examples
#' X = rnormalmultinomial(100, 100, rep(0,4), diag(4))
#' dm_fit(X)
#' @export
dm_fit = function(X, expected = TRUE, eps = 10e-5, max.iter = 5000){
  res = c_dm_fit(X, eps, max.iter)
  res[[1]] = as.numeric(res[[1]])
  names(res) = c('alpha', 'iter')
  if(expected){
    res$expected = dm_expected(X, res$alpha)
  }
  res
}

#'
#' Calculate the expected probabilities
#'
#' Calculate the probabilities conditioned on parameters alpha
#' more likely to generate dataset X
#'
#' @param X count data set
#' @param alpha paramater alpha
#' @return the expected probabilities
#' @examples
#' X = rnormalmultinomial(100, 100, rep(0,4), diag(4))
#' fit = dm_fit(X)
#' P = dm_expected(X, fit$alpha)
#' head(P)
dm_expected = function(X, alpha){
  X = t(alpha+t(X))
  X/rowSums(X)
}

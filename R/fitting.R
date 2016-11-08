#'
#' Log-ratio normal-multinomial parameters estimation. The parameters mu and sigma are expressed with respect
#' basis given by function @ilrBase.
#'
#' @param X count data set
#' @param eps maximum error accepted on the last EM-step
#' @param nsim number of simulations used in the E-step
#' @param parallel.cluster parallel Socket Cluster created with function @makeCluster
#' @param max.em.iter maximum number of steps allowed in the EM-algorithm
#' @param verbose show information during estimation
#' @return A list with parameters mu and sigma and the number of iterations before convergence
#' @examples
#' X = rnormalmultinomial(100, 100, rep(0,4), diag(4))
#' nm_fit(X, verbose = T)
#' @export
nm_fit = function(X, eps = 0.001, nsim = 1000, parallel.cluster = NULL,
                  max.em.iter = 100, verbose = FALSE){

  MU = ilr_coordinates(matrix(apply(X/apply(X, 1, sum), 2, sum), nrow=1))[1,]
  SIGMA = diag(ncol(X)-1)

  if(nrow(SIGMA) <= 6){
    Z = matrix(randtoolbox::halton(nsim, dim = nrow(SIGMA), normal = TRUE), ncol=nrow(SIGMA))
  }else{
    Z = matrix(randtoolbox::sobol(nsim, dim = nrow(SIGMA), normal = TRUE), ncol=nrow(SIGMA))
  }

  err = eps + 1
  iter = 0
  while(err > eps & iter < max.em.iter){
    if(!is.null(parallel.cluster)){
      FIT = parallel::parApply(parallel.cluster, X, 1, expectedMonteCarlo2, MU, SIGMA, Z)
    }else{
      FIT = apply(X, 1, expectedMonteCarlo2, MU, SIGMA, Z)
    }
    E = t(sapply(FIT, function(fit) colMeans(na.omit(fit[[2]]))))

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

#'
#' Calculate the expected probabilities P to generate the dataset X conditioned on parameters
#' mu and sigma
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
#' @export
dm_fit = function(X, eps = 10e-5, max.iter = 5000, verbose = FALSE){
  res = c_dm_fit(X, eps, max.iter)
  res[[1]] = as.numeric(res[[1]])
  names(res) = c('alpha', 'iter')
  res
}

#'
#' Calculate the expected probabilities P to generate the dataset X conditioned on parameters
#' alpha
#' @param X count data set
#' @param alpha paramater alpha
#' @return the expected probabilities
#' @examples
#' X = rnormalmultinomial(100, 100, rep(0,4), diag(4))
#' fit = dm_fit(X, verbose = T)
#' P = dm_expected(X, fit$alpha)
#' head(P)
#' @export
#' @export
dm_expected = function(X, alpha){
  X = t(alpha+t(X))
  X/rowSums(X)
}

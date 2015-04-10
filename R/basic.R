#' Generate a sample from a normalmultinomial distribution
#'
#' This function generates a sample generated from a normal-multinomial (nm) distribution
#' @param N sample size
#' @param Mu alr mean indicating the center of the nm distribution
#' @param Sigma alt covariance indicating the variable of the nm distribution
#' @param Size Sampling size for the counting process
#' @return A sample comming from a normal-multinomial distribution
#'
#' @export
rnorm.multinom <- function (N, Mu, Sigma, Size) {
  M = cbind(exp(mvtnorm::rmvnorm(N, Mu, Sigma)), 1)
  M = M / apply(M, 1, sum)
  M = cbind(as.numeric(Size), M)
  t(apply(M, 1, function(p) as.numeric(rmultinom(1, as.numeric(p[1]), p[-1])) ))
}

#' Finds the mean and covariance of a normal multinomial distribution
#' 
#' @param X normal-multinomial sample
#' @export
adjustNormalMultinomial <- function(X) {
  .Call('adjustNormalMultinomial', PACKAGE = 'normalmultinomial', X)
}
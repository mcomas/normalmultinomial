#'
#' Simulate variables following a normal distribution.
#'
#' @param n sample size
#' @param mu mean parameter for the mean in a aln-normal distribution
#' @param sigma parameter for the sigma in a aln-normal distribution
#' @return Sample matrix
#' @export
rnormal <- function(n, mu, sigma) {
  .Call('normalmultinomial_rnormal', PACKAGE = 'normalmultinomial', n, mu, sigma)
}

#' Simulate variables following a multinomial distribution.
#'
#' @param A matrix with the probabilities used to generate the sample
#' @param size vector with the countings in each row of A
#' @return sample with the same dimension as A
#' @export
rmultinomial <- function(A, size, seed = 1L) {
  .Call('normalmultinomial_rmultinomial', PACKAGE = 'normalmultinomial', A, size, seed)
}

#' Simulate variables following a normal multinomial distribution.
#'
#' @param n sample size
#' @param mu mean parameter for the mean in a aln-normal distribution
#' @param sigma parameter for the sigma in a aln-normal distribution
#' @return Sample matrix
#' @export
rnormalmultinomial <- function(mu, sigma, size, seed = 1L) {
  .Call('normalmultinomial_rnormalmultinomial', PACKAGE = 'normalmultinomial', mu, sigma, size, seed)
}


#'
#' Simulate variables following a normal distribution.
#'
#' @param n sample size
#' @param mu mean parameter for the mean in a normal distribution
#' @param sigma parameter for the sigma in a normal distribution
#' @return Sample matrix
#' @examples
#' rnormal(10, mu = c(0,0), sigma = matrix(c(2,1,1,2), nrow=2))
#' @export
rnormal <- function(n, mu, sigma) {
  if( ! (is.vector(mu) & is.numeric(mu)) ){
    error("'mu' must be a numeric vector")
  }
  if( ! (is.matrix(sigma) & is.numeric(sigma)) ){
    error("'sigma' must be a numeric matrix")
  }
  if( det(sigma) <= 0 ){
    error("'sigma' must be a positive-definite matrix")
  }
  .Call('normalmultinomial_rnormal', PACKAGE = 'normalmultinomial', n, mu, sigma)
}

#' Simulate variables following a multinomial distribution.
#'
#' @param n matrix with the probabilities used to generate the sample
#' @param size number of trials en each observation
#' @param prob probability of observations (vector) or specific to each observation (matrix)
#' @return multinomial sample in a matrix
#' @examples
#' rmultinomial(10, c(10, 100), c(0.2,0.8))
#' @export
rmultinomial <- function(n, size, prob) {
  if( is.data.frame(prob) ){
    prob = as.matrix(prob)
  }
  if(! (is.vector(n) & is.numeric(n)) ){
    error("'n' must be numeric")
  }
  if(! (is.numeric(prob) & (is.vector(prob) | is.matrix(prob) )) ){
    error("'prob' must be numeric vector or a numeric matrix")
  }
  if( ! (is.vector(size) & is.numeric(size)) ){
    error("'size' must be a numeric vector")
  }
  if(length(n) == 1){
    N = n
  }else{
    N = length(n)
  }
  if(is.vector(prob)){
    A = matrix(rep(prob, N), nrow = N, byrow = TRUE)
  }else{
    if(nrow(prob) != N){
      stop("Number of rows for matrix prob different than n")
    }
    A = prob / apply(prob,1,sum)
  }
  SIZE = rep_len(size, N)
  .Call('normalmultinomial_rmultinomial', PACKAGE = 'normalmultinomial', A, SIZE, sample(.Random.seed, 1))
}

#' Simulate variables following a normal multinomial distribution.
#'
#' @param n number of observations
#' @param size number of trials
#' @param mu mean parameter for the mean in the aln-normal distribution
#' @param sigma parameter for the sigma in the aln-normal distribution
#' @return Sample matrix
#' @examples
#' rnormalmultinomial(10, 5, mu = c(0,0), sigma = matrix(c(2,1,1,2), nrow=2))
#' @export
rnormalmultinomial <- function(n, size, mu, sigma, probs = FALSE) {
  if(! (is.vector(n) & is.numeric(n)) ){
    error("'n' must be numeric")
  }
  if( ! (is.vector(size) & is.numeric(size)) ){
    error("'size' must be a numeric vector")
  }
  if(length(n) == 1){
    N = n
  }else{
    N = length(n)
  }
  if( ! (is.vector(mu) & is.numeric(mu)) ){
    error("'mu' must be a numeric vector")
  }
  if( ! (is.matrix(sigma) & is.numeric(sigma)) ){
    error("'sigma' must be a numeric matrix")
  }
  if( det(sigma) <= 0 ){
    error("'sigma' must be a positive-definite matrix")
  }
  SIZE = rep_len(size, N)
  res = .Call('normalmultinomial_rnormalmultinomial', PACKAGE = 'normalmultinomial', mu, sigma, SIZE, sample(.Random.seed, 1))
  if(probs){
    names(res) = c('probs', 'counts')
    res[['probs']] = res[['probs']] / apply(res[['probs']], 1, sum)
    res
  }else{
    res[[2]]
  }
}

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <random>
#include <vector>

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

const double log2pi = std::log(2.0 * M_PI);

// [[Rcpp::export]]
arma::mat rnormal(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

// [[Rcpp::export]]
arma::mat rmultinomial(arma::mat A, arma::vec size, int seed = 1) {
  int k = A.n_cols;
  int n = A.n_rows;

  arma::mat res = arma::zeros<arma::mat>(n, k);

  std::vector<double> p(k);
  std::random_device rd;
  std::mt19937 gen(seed);
  for(int i=0; i < n; i++){
    for(int j=0; j < k; j++) p[j] = A(i,j);
    std::discrete_distribution<> d(p.begin(),p.end());
    for(int l=0; l < size[i]; l++) ++res(i, d(gen));
  }
  return(res);
}

// [[Rcpp::export]]
arma::mat rnormalmultinomial(arma::vec mu, arma::mat sigma, arma::vec size, int seed = 1){
  arma::arma_rng::set_seed(seed);
  int n = size.n_elem;
  int k = mu.n_elem;
  arma::mat A = arma::ones<arma::mat>(n, k+1);
  A(arma::span::all, arma::span(0,k-1)) = exp(rnormal(n, mu, sigma));
  return(rmultinomial(A, size, seed));
}

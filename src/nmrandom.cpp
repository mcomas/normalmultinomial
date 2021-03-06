// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <random>
#include <vector>
#include "nmrandom.h"
#include "coda.h"
// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat c_rnormal(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

// [[Rcpp::export]]
arma::mat c_rmultinomial(arma::mat A, arma::vec size, int seed = 1) {
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
List c_rnormalmultinomial(arma::vec mu, arma::mat sigma, arma::vec size, int seed){
  int n = size.n_elem;
  int k = mu.n_elem;
  arma::mat ILR_TO_ALR = ilr_to_alr(k+1);
  arma::vec mu_alr = ILR_TO_ALR * mu;
  arma::mat sigma_alr = ILR_TO_ALR * sigma * ILR_TO_ALR.t();
  arma::mat A = arma::ones<arma::mat>(n, k+1);
  A(arma::span::all, arma::span(0,k-1)) = exp(c_rnormal(n, mu_alr, sigma_alr));
  return(List::create(A, c_rmultinomial(A, size, seed)));
}

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "stepE.h"
#include "loglik.h"
#include <random>
#include <vector>

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

//' @export
// [[Rcpp::export]]
Rcpp::List expectedA1_all(arma::vec x, arma::vec mu, arma::mat sigma, int nsim = 100){
  int K = x.size();
  int k = K - 1;
  int nsim2 = 2 * nsim;
  arma::mat Z1 = arma::randn(nsim, k);
  arma::mat Z2 = -Z1;
  arma::mat Z = join_cols(Z1, Z2);
  arma::mat Ap = arma::mat(nsim2, k);

  Ap = arma::repmat(mu.t(), nsim2, 1) + Z * arma::chol(sigma);

  arma::vec lik = exp( mvf_multinom_const(x) + mvf_multinom_mult(Ap, x) );
  arma::vec lik_st = lik / mean(lik);

  arma::mat expected_m0 = lik;
  arma::mat expected_m1 = arma::mat(nsim2, k);
  arma::cube expected_m2 = arma::cube(k,k, nsim2);
  for(int i = 0;i < k; i++){
    expected_m1.col(i) = Ap.col(i) % lik_st;
    for(int j = 0;j < k; j++){
      expected_m2.tube(i,j) = Ap.col(i) % Ap.col(j) % lik_st;
    }
  }
  return Rcpp::List::create(expected_m0, expected_m1, expected_m2);
}

//' @export
// [[Rcpp::export]]
Rcpp::List expectedA2_all(arma::vec x, arma::vec mu, arma::mat sigma, int nsim = 100){
  int K = x.size();
  int k = K - 1;
  int nsim2 = 2 * nsim;

  arma::mat Z1 = arma::randn(nsim, k);
  arma::mat Z2 = -Z1;
  arma::mat Z = join_cols(Z1, Z2);
  arma::mat Ap = arma::mat(nsim2, k);

  arma::mat inv_sigma = inv_sympd(sigma);
  arma::vec mu0 = mvf_maximum(x, mu, inv_sigma, 1e-8, 100, 0.66);

  Ap = arma::repmat(mu0.t(), nsim2, 1) + Z * arma::chol(sigma);

  arma::vec lik = exp( mvf_multinom_const(x) + mvf_multinom_mult(Ap, x) ) %
  mvf_norm(Ap, mu, inv_sigma)  / mvf_norm(Ap, mu0, inv_sigma);
  arma::vec lik_st = lik / mean(lik);

  arma::mat expected_m0 = lik;
  arma::mat expected_m1 = arma::mat(nsim2, k);
  arma::cube expected_m2 = arma::cube(k,k, nsim2);
  for(int i = 0;i < k; i++){
    expected_m1.col(i) = Ap.col(i) % lik_st;
    for(int j = 0;j < k; j++){
      expected_m2.tube(i,j) = Ap.col(i) % Ap.col(j) % lik_st;
    }
  }
  return Rcpp::List::create(expected_m0, expected_m1, expected_m2);
}

//' @export
// [[Rcpp::export]]
Rcpp::List expectedA3_all(arma::vec x, arma::vec mu, arma::mat sigma, int nsim = 100){
  int K = x.size();
  int k = K - 1;
  int nsim2 = 2 * nsim;
  arma::mat Z1 = arma::randn(nsim, k);
  arma::mat Z2 = -Z1;
  arma::mat Z = join_cols(Z1, Z2);
  arma::mat Ap = arma::mat(nsim2, k);

  arma::mat inv_sigma = inv_sympd(sigma);
  arma::vec mu0 = mvf_maximum(x, mu, inv_sigma, 1e-8, 100, 0.66);

  arma::mat H = arma::mat(k,k);
  for(int i=0;i<k;i++){
    for(int j=0;j<k;j++){
      H(i,j) = mvf_deriv2(i, j, mu0, mu, inv_sigma, x);
    }
  }
  arma::mat Hinv = inv_sympd(-H);
  //Rcout << -H << std::endl;

  Ap = arma::repmat(mu0.t(), nsim2, 1) + Z * arma::chol(Hinv);

  arma::vec lik = exp( mvf_multinom_const(x) + mvf_multinom_mult(Ap, x) ) %
  mvf_norm(Ap, mu, inv_sigma)  / mvf_norm(Ap, mu0, -H);
  arma::vec lik_st = lik / mean(lik);

  arma::mat expected_m0 = lik;
  arma::mat expected_m1 = arma::mat(nsim2, k);
  arma::cube expected_m2 = arma::cube(k,k, nsim2);
  for(int i = 0;i < k; i++){
    expected_m1.col(i) = Ap.col(i) % lik_st;
    for(int j = 0;j < k; j++){
      expected_m2.tube(i,j) = Ap.col(i) % Ap.col(j) % lik_st;
    }
  }
  return Rcpp::List::create(expected_m0, expected_m1, expected_m2);
}

//' @export
// [[Rcpp::export]]
Rcpp::List expectedA4_all(arma::vec x, arma::vec mu, arma::mat sigma, int nsim = 100){
  int K = x.size();
  int k = K - 1;

  arma::mat Z1 = arma::randn(nsim, k);
  arma::mat Z2 = -Z1;

  arma::mat Ap1 = arma::mat(nsim, k);
  arma::mat Ap2 = arma::mat(nsim, k);

  arma::mat inv_sigma = inv_sympd(sigma);
  arma::vec mu0 = mvf_maximum(x, mu, inv_sigma, 1e-8, 100, 0.66);

  arma::mat H = arma::mat(k,k);
  for(int i=0;i<k;i++){
    for(int j=0;j<k;j++){
      H(i,j) = mvf_deriv2(i, j, mu0, mu, inv_sigma, x);
    }
  }
  arma::mat Hinv = inv_sympd(-H);

  Ap1 = arma::repmat(mu0.t(), nsim, 1) + Z1 * arma::chol(Hinv);
  Ap2 = arma::repmat(mu0.t(), nsim, 1) + Z2 * arma::chol(Hinv);

  arma::vec lik1 = exp( mvf_multinom_const(x) + mvf_multinom_mult(Ap1, x) ) % mvf_norm(Ap1, mu, inv_sigma)  / mvf_norm(Ap1, mu0, -H);
  arma::vec lik2 = exp( mvf_multinom_const(x) + mvf_multinom_mult(Ap2, x) ) % mvf_norm(Ap2, mu, inv_sigma)  / mvf_norm(Ap2, mu0, -H);

  arma::vec lik1_st = lik1 / mean(lik1);
  arma::vec lik2_st = lik2 / mean(lik2);

  arma::vec lik = (lik1 + lik2) / 2;

  arma::mat expected_m0 = lik;
  arma::mat expected_m1 = arma::mat(nsim, k);
  arma::cube expected_m2 = arma::cube(k,k, nsim);
  for(int i = 0;i < k; i++){
    expected_m1.col(i) = (Ap1.col(i) % lik1_st + Ap2.col(i) % lik2_st) / 2;
    for(int j = 0;j < k; j++){
      expected_m2.tube(i,j) = (Ap1.col(i) % Ap1.col(j) % lik1_st + Ap2.col(i) % Ap2.col(j) % lik2_st) / 2;
    }
  }
  return Rcpp::List::create(expected_m0, expected_m1, expected_m2);
}

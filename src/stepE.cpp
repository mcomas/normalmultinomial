// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "stepE.h"
#include "loglik.h"
#include <random>
#include <vector>

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

Rcpp::List expectedA1(arma::vec x, arma::vec mu, arma::mat sigma, int nsim = 100){
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

  arma::vec expected_m1 = arma::vec(k);
  arma::mat expected_m2 = arma::mat(k,k);
  for(int i = 0;i < k; i++){
    expected_m1[i] = mean(Ap.col(i) % lik_st);
    for(int j = 0;j < k; j++){
      expected_m2(i,j) = mean(Ap.col(i) % Ap.col(j) % lik_st);
    }
  }
  return Rcpp::List::create(expected_m1, expected_m2);
}

Rcpp::List expectedA2(arma::vec x, arma::vec mu, arma::mat sigma, int nsim = 100){
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

  arma::vec expected_m1 = arma::vec(k);
  arma::mat expected_m2 = arma::mat(k,k);
  for(int i = 0;i < k; i++){
    expected_m1[i] = mean(Ap.col(i) % lik_st);
    for(int j = 0;j < k; j++){
      expected_m2(i,j) = mean(Ap.col(i) % Ap.col(j) % lik_st);
    }
  }
  return Rcpp::List::create(expected_m1, expected_m2);
}

Rcpp::List expectedA3(arma::vec x, arma::vec mu, arma::mat sigma, int nsim = 100){
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

  arma::vec expected_m1 = arma::vec(k);
  arma::mat expected_m2 = arma::mat(k,k);
  for(int i = 0;i < k; i++){
    expected_m1[i] = mean(Ap.col(i) % lik_st);
    for(int j = 0;j < k; j++){
      expected_m2(i,j) = mean(Ap.col(i) % Ap.col(j) % lik_st);
    }
  }
  return Rcpp::List::create(expected_m1, expected_m2);
}


Rcpp::List expectedA4(arma::vec a, arma::vec x, arma::vec mu, arma::mat sigma, int nsim = 100){
  int K = x.size();
  int k = K - 1;
  int nsim2 = 2 * nsim;
  arma::mat Z1 = arma::randn(nsim, k);
  arma::mat Z2 = -Z1;
  arma::mat Z = join_cols(Z1, Z2);
  arma::mat Ap = arma::mat(nsim2, k);

  arma::mat inv_sigma = inv_sympd(sigma);

  arma::mat Ha = arma::mat(k,k);
  for(int i=0;i<k;i++){
    for(int j=0;j<k;j++){
      Ha(i,j) = mvf_deriv2(i, j, a, mu, inv_sigma, x);
    }
  }
  arma::mat Hinv = inv_sympd(-Ha);
  //Rcout << -H << std::endl;

  Ap = arma::repmat(a.t(), nsim2, 1) + Z * arma::chol(Hinv);

  arma::vec lik = exp( mvf_multinom_const(x) + mvf_multinom_mult(Ap, x) ) %
  mvf_norm(Ap, mu, inv_sigma)  / mvf_norm(Ap, a, -Ha);
  arma::vec lik_st = lik / mean(lik);

  arma::vec expected_m1 = arma::vec(k);
  arma::mat expected_m2 = arma::mat(k,k);
  for(int i = 0;i < k; i++){
    expected_m1[i] = mean(Ap.col(i) % lik_st);
    for(int j = 0;j < k; j++){
      expected_m2(i,j) = mean(Ap.col(i) % Ap.col(j) % lik_st);
    }
  }
  return Rcpp::List::create(expected_m1, expected_m2);
}

//' @export
// [[Rcpp::export]]
arma::mat stepE(arma::mat X, arma::vec mu, arma::mat sigma, int nsim = 100){
  int n = X.n_rows;
  int K = X.n_cols;
  arma::mat A = arma::mat(n, K-1);
  for(int i = 0; i < n; i++){
    Rcpp::List moments = expectedA2(X.row(i).t(), mu, sigma, nsim);
    arma::mat m1 = moments[0];
    A.row(i) = m1.t();
  }
  return A;
}

//' @export
// [[Rcpp::export]]
Rcpp::List stepEM1(arma::mat X, arma::vec mu, arma::mat sigma, int nsim = 100){
  int n = X.n_rows;
  int K = X.n_cols;
  arma::vec mu_next = arma::zeros(K-1);
  arma::mat sigma_next = arma::zeros(K-1, K-1);
  for(int i = 0; i < n; i++){
    Rcpp::List moments = expectedA2(X.row(i).t(), mu, sigma, nsim);
    arma::mat m1 = moments[0];
    arma::mat m2 = moments[1];
    mu_next += m1;
    sigma_next += m2;
  }
  mu_next = mu_next / n;
  sigma_next = sigma_next / n - mu_next * mu_next.t();
  return Rcpp::List::create(mu_next, sigma_next);
}

//' @export
// [[Rcpp::export]]
Rcpp::List stepEM2(arma::mat A, arma::mat X, arma::vec mu, arma::mat sigma, int nsim = 100){
  int n = X.n_rows;
  int K = X.n_cols;
  arma::vec mu_next = arma::zeros(K-1);
  arma::mat sigma_next = arma::zeros(K-1, K-1);
  for(int i = 0; i < n; i++){
    Rcpp::List moments = expectedA4(A.row(i).t(), X.row(i).t(), mu, sigma, nsim);
    arma::mat m1 = moments[0];
    arma::mat m2 = moments[1];
    mu_next += m1;
    sigma_next += m2;
  }
  mu_next = mu_next / n;
  sigma_next = sigma_next / n - mu_next * mu_next.t();
  return Rcpp::List::create(mu_next, sigma_next);
}

// [[Rcpp::export]]
double c_dnormalmultinomial(arma::vec x, arma::vec mu, arma::mat sigma, int nsim = 100){
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
  return mean(lik);
}

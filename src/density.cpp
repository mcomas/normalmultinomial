// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "loglik.h"
#include "coda.h"

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

const double log2pi = std::log(2.0 * M_PI);

//' @export
// [[Rcpp::export]]
arma::mat log_dnormal(arma::mat A, arma::vec mu, arma::mat inv_sigma){
  int k = A.n_cols;
  int n = A.n_rows;

  double norm_const = -0.5 * k * log(2*PI) + 0.5 * log(det(inv_sigma));
  arma::mat log_norm = arma::mat(n, 1);

  for(int i = 0; i < n; i++){
    arma::mat cal = (A.row(i) - mu.t()) * inv_sigma * (A.row(i) - mu.t()).t();
    log_norm(i,0) = -0.5 * cal(0);
  }

  return(norm_const + log_norm);
}

//' @export
// [[Rcpp::export]]
arma::mat dnormal(arma::mat A, arma::vec mu, arma::mat inv_sigma){
  return arma::exp(log_dnormal(A,mu,inv_sigma));
}

double log_multinomial_constant(arma::rowvec x){
  int K = x.size();

  double x_total = 0;
  for(int i = 0; i < K; i++) x_total += x[i];

  double x_total_fact = log(x_total);
  for(int l = 1; l < x_total; l++) x_total_fact += log(l);

  double x_parts_fact = 0;
  for(int i = 0; i < K; i++){
    for(int l = 2; l <= x[i]; l++){
      x_parts_fact += log(l);
    }
  }
  double mult_const = 0;
  mult_const += (x_total_fact - x_parts_fact);

  return(mult_const);
}

//' @export
// [[Rcpp::export]]
arma::mat log_dmultinomial(arma::mat X, arma::mat A){

  arma::mat Mconst = arma::mat(X.n_rows, 1);
  for(int i = 0; i < X.n_rows; i++){
    Mconst(i,0) = log_multinomial_constant(X.row(i));
  }
  arma::mat Mprod = arma::sum(X % arma::log(inv_ilr_coordinates(A)), 1);

  return(Mconst + Mprod);
}

//' @export
// [[Rcpp::export]]
arma::mat dmultinomial(arma::mat X, arma::mat A){
  return(arma::exp(log_dmultinomial(X,A)));
}

//' @export
// [[Rcpp::export]]
arma::mat log_df_x_a(arma::mat X, arma::mat A, arma::vec mu, arma::mat inv_sigma){
  return( log_dmultinomial(X,A) + log_dnormal(A,mu,inv_sigma));
}

//' @export
// [[Rcpp::export]]
arma::mat df_x_a(arma::mat X, arma::mat A, arma::vec mu, arma::mat inv_sigma){
  return( arma::exp(log_df_x_a(X,A,mu,inv_sigma)));
}


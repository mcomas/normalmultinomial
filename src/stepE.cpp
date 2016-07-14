// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "stepE.h"
#include "loglik.h"
#include <random>
#include <vector>

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

bool DEBUG = true;

using namespace Rcpp;

const double log2pi = std::log(2.0 * M_PI);

//' @export
// [[Rcpp::export]]
arma::vec choose_starting(arma::vec x, arma::vec mu, arma::mat inv_sigma,
                    double eps = 1e-8, int max_iter = 100, double prop = 0.66) {
  int k = x.size() - 1;
  arma::vec a = arma::vec(k);

  for(int j = 0; j < k; j++){
    if(x[j] != 0 && x[k] != 0){
      a[j] = log(x[j]/x[k]);
    }else{
      if(x[j] == 0 && x[k] == 0){
        a[j] = 0;
      }else{
        if(x[j] == 0){
          a[j] = log(prop/x[k]);
        }else{
          a[j] = log(x[j]/prop);
        }
      }
    }
  }

  arma::vec out = arma::vec(a);
  arma::vec deriv = arma::zeros<arma::vec>(k);
  arma::mat deriv2 = arma::zeros<arma::mat>(k, k);
  arma::vec step = arma::zeros<arma::vec>(k);

  int current_iter = 0;
  do{
    current_iter++;
    for(int I=0; I<k; I++){
      deriv[I] =  mvf_deriv(I, out, mu, inv_sigma, x);
      for(int J=0; J<k; J++){
        deriv2(I,J) = mvf_deriv2(I, J, out, mu, inv_sigma, x);
      }
    }
    step = arma::solve(deriv2, deriv);
    out = out - step;
  }while( norm(step, 2) > eps && current_iter < max_iter);
  if(DEBUG){
    Rcout << "Iterations: " << current_iter << std::endl;
  }

  return out;
}

//' @export
// [[Rcpp::export]]
arma::vec vec_mvf_multinom_mult(arma::mat A, arma::vec x){
  int k = A.n_cols;
  int n = A.n_rows;
  arma::vec multinom = arma::vec(n);
  for(int i = 0; i < n; i++){
    double kappa = 1;
    for(int j = 0; j < k; j++) kappa += exp(A(i,j));
    multinom[i] = -x[k] * log(kappa);
    for(int j = 0; j < k; j++){
      multinom[i] += x[j] * ( A(i,j) - log(kappa));
    }
  }
  return(multinom);
}

//' @export
// [[Rcpp::export]]
arma::mat expectedA1(arma::vec x, arma::vec mu, arma::mat sigma, int nsim = 100){
  int K = x.size();
  int k = K - 1;
  int nsim2 = 2 * nsim;
  arma::mat Z1 = arma::randn(nsim, k);
  arma::mat Z2 = -Z1;
  arma::mat Z = join_cols(Z1, Z2);
  arma::mat Ap = arma::mat(nsim2, k);

  Ap = arma::repmat(mu.t(), nsim2, 1) + Z * arma::chol(sigma);

  arma::vec lik = exp( mvf_multinom_const(x) + vec_mvf_multinom_mult(Ap, x) );

  arma::mat expected = Ap % repmat(lik, 1, k) / mean(lik);

  return expected;
}

//' @export
// [[Rcpp::export]]
double vec_mvf_norm(arma::mat A, arma::vec mu, arma::mat inv_sigma){
  int k = A.n_cols;
  int n = A.n_rows;
  Rcout << mu << std::endl;
  arma::mat MU = repmat(mu, n, 1);
  Rcout << MU << std::endl;
  double norm_const = -0.5 * k * log(2*PI) - 0.5 * log(det(inv_sigma));
  arma::mat log_norm =  -0.5 * (A-MU).t() * inv_sigma * (A-MU);

  double norm = norm_const + log_norm(0);

  return(norm);
}

//' @export
// [[Rcpp::export]]
arma::mat expectedA2(arma::vec x, arma::vec mu, arma::mat sigma, int nsim = 100){
  int K = x.size();
  int k = K - 1;
  int nsim2 = 2 * nsim;
  arma::mat Z1 = arma::randn(nsim, k);
  arma::mat Z2 = -Z1;
  arma::mat Z = join_cols(Z1, Z2);
  arma::mat Ap = arma::mat(nsim2, k);

  arma::vec mu0 = choose_starting(x, mu, inv_sympd(sigma));

  Ap = arma::repmat(mu0.t(), nsim2, 1) + Z * arma::chol(sigma);

  arma::vec lik = exp( mvf_multinom_const(x) + vec_mvf_multinom_mult(Ap, x) );

  arma::mat expected = Ap % repmat(lik, 1, k) / mean(lik);

  return expected;
}

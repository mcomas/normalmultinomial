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
  // if(DEBUG){
  //   Rcout << "Iterations: " << current_iter << std::endl;
  // }

  return out;
}

//' @export
// [[Rcpp::export]]
arma::mat hessian(arma::vec a, arma::vec mu, arma::mat inv_sigma, arma::vec x){
  int k = x.size() - 1;
  arma::mat deriv2 = arma::mat(k,k);
  for(int I=0; I<k; I++){
    for(int J=0; J<k; J++){
      deriv2(I,J) = mvf_deriv2(I, J, a, mu, inv_sigma, x);
    }
  }
  return deriv2;
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
Rcpp::List expectedA1(arma::vec x, arma::vec mu, arma::mat sigma, int nsim = 100){
  int K = x.size();
  int k = K - 1;
  int nsim2 = 2 * nsim;
  arma::mat Z1 = arma::randn(nsim, k);
  arma::mat Z2 = -Z1;
  arma::mat Z = join_cols(Z1, Z2);
  arma::mat Ap = arma::mat(nsim2, k);

  Ap = arma::repmat(mu.t(), nsim2, 1) + Z * arma::chol(sigma);

  arma::vec lik = exp( mvf_multinom_const(x) + vec_mvf_multinom_mult(Ap, x) );
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
arma::vec vec_mvf_norm(arma::mat A, arma::vec mu, arma::mat inv_sigma){
  int k = A.n_cols;
  int n = A.n_rows;

  double norm_const = -0.5 * k * log(2*PI) + 0.5 * log(det(inv_sigma));
  arma::vec log_norm = arma::vec(n);

  for(int i = 0; i < n; i++){
    arma::mat cal = (A.row(i) - mu.t()) * inv_sigma * (A.row(i) - mu.t()).t();
    log_norm[i] = -0.5 * cal(0);
  }
  arma::vec norm = norm_const + log_norm;

  return(exp(norm));
}

//' @export
// [[Rcpp::export]]
Rcpp::List expectedA2(arma::vec x, arma::vec mu, arma::mat sigma, int nsim = 100){
  int K = x.size();
  int k = K - 1;
  int nsim2 = 2 * nsim;

  arma::mat Z1 = arma::randn(nsim, k);
  arma::mat Z2 = -Z1;
  arma::mat Z = join_cols(Z1, Z2);
  arma::mat Ap = arma::mat(nsim2, k);

  arma::mat inv_sigma = inv_sympd(sigma);
  arma::vec mu0 = choose_starting(x, mu, inv_sigma);

  Ap = arma::repmat(mu0.t(), nsim2, 1) + Z * arma::chol(sigma);

  arma::vec lik = exp( mvf_multinom_const(x) + vec_mvf_multinom_mult(Ap, x) ) %
    vec_mvf_norm(Ap, mu, inv_sigma)  / vec_mvf_norm(Ap, mu0, inv_sigma);
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
Rcpp::List expectedA3(arma::vec x, arma::vec mu, arma::mat sigma, int nsim = 100){
  int K = x.size();
  int k = K - 1;
  int nsim2 = 2 * nsim;
  arma::mat Z1 = arma::randn(nsim, k);
  arma::mat Z2 = -Z1;
  arma::mat Z = join_cols(Z1, Z2);
  arma::mat Ap = arma::mat(nsim2, k);

  arma::mat inv_sigma = inv_sympd(sigma);
  arma::vec mu0 = choose_starting(x, mu, inv_sigma);

  arma::mat H = arma::mat(k,k);
  for(int i=0;i<k;i++){
    for(int j=0;j<k;j++){
      H(i,j) = mvf_deriv2(i, j, mu0, mu, inv_sigma, x);
    }
  }
  arma::mat Hinv = inv_sympd(-H);
  //Rcout << -H << std::endl;

  Ap = arma::repmat(mu0.t(), nsim2, 1) + Z * arma::chol(Hinv);

  arma::vec lik = exp( mvf_multinom_const(x) + vec_mvf_multinom_mult(Ap, x) ) %
  vec_mvf_norm(Ap, mu, inv_sigma)  / vec_mvf_norm(Ap, mu0, -H);
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

  arma::vec lik = exp( mvf_multinom_const(x) + vec_mvf_multinom_mult(Ap, x) ) %
  vec_mvf_norm(Ap, mu, inv_sigma)  / vec_mvf_norm(Ap, a, -Ha);
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

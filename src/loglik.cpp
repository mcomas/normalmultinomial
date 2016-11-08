// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
# include "loglik.h"

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

const double log2pi = std::log(2.0 * M_PI);


double mvf_norm(arma::vec a, arma::vec mu, arma::mat inv_sigma){
  int k = a.size();

  double norm_const = -0.5 * k * log2pi - 0.5 * log(det(inv_sigma));
  arma::mat log_norm =  -0.5 * (a-mu).t() * inv_sigma * (a-mu);

  double norm = norm_const + log_norm(0);

  return(norm);
}

arma::vec mvf_norm(arma::mat A, arma::vec mu, arma::mat inv_sigma){
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

// [[Rcpp::export]]
double mvf_multinom_const(arma::vec x){
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

// [[Rcpp::export]]
double mvf_multinom_mult(arma::vec a, arma::vec x){
  int k = a.size();

  double kappa = 1;
  for(int i = 0; i < k; i++) kappa += exp(a[i]);
  double multinom = -x[k] * log(kappa);
  for(int i = 0; i < k; i++) multinom += x[i] * ( a[i] - log(kappa));

  return(multinom);
}

arma::vec mvf_multinom_mult(arma::mat A, arma::vec x){
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

double mvf_multinom(arma::vec a, arma::vec x){
  int k = a.size();
  int K = k +1;

  double x_total = 0;
  for(int i = 0; i < K; i++) x_total += x[i];

  double x_total_fact = log(x_total);
  for(int l = 1; l < x_total; l++) x_total_fact += log(l);

  double x_parts_fact = 0;
  for(int i = 0; i <= K; i++){
    for(int l = 1; l <= x[i]; l++){
      x_parts_fact += log(l);
    }
  }

  double mult_const = 0;
  mult_const += (x_total_fact - x_parts_fact);

  double kappa = 1;
  for(int i = 0; i < k; i++) kappa += exp(a[i]);
  double multinom = -x[k] * log(kappa);
  for(int i = 0; i < k; i++) multinom += x[i] * ( a[i] - log(kappa));

  double mult = mult_const + multinom;
  return(mult);
}

// [[Rcpp::export]]
double mvf(arma::vec a, arma::vec mu, arma::mat inv_sigma, arma::vec x){
  int k = a.size();
  int K = k +1;

  double norm_const = -0.5 * k * log2pi - 0.5 * log(det(inv_sigma));
  arma::mat log_norm =  -0.5 * (a-mu).t() * inv_sigma * (a-mu);

  double norm = norm_const + log_norm(0);

  double x_total = 0;
  for(int i = 0; i < K; i++) x_total += x[i];

  double x_total_fact = log(x_total);
  for(int l = 1; l < x_total; l++) x_total_fact += log(l);

  double x_parts_fact = 0;
  for(int i = 0; i <= K; i++){
    for(int l = 1; l <= x[i]; l++){
      x_parts_fact += log(l);
    }
  }

  double mult_const = 0;
  mult_const += (x_total_fact - x_parts_fact);

  double kappa = 1;
  for(int i = 0; i < k; i++) kappa += exp(a[i]);
  double multinom = -x[k] * log(kappa);
  for(int i = 0; i < k; i++) multinom += x[i] * ( a[i] - log(kappa));

  double mult = mult_const + multinom;
  return(norm + mult);
}

// [[Rcpp::export]]
double mvf_deriv(int I, arma::vec a, arma::vec mu, arma::mat inv_sigma, arma::vec x){
  int k = a.size();
  arma::mat log_norm =  -(a-mu).t() * inv_sigma(arma::span::all, I);
  double kappa = 1;
  for(int i = 0; i < k; i++) kappa += exp(a(i));
  double mult = 0;
  mult += x(I) * (kappa-exp(a(I))) / kappa;
  for(int i = 0; i < I; i++) mult += x(i) * (-exp(a(I))) / kappa;
  for(int i = I+1; i < k+1; i++) mult += x(i) * (-exp(a(I))) / kappa;
  return log_norm(0) + mult;
}

// [[Rcpp::export]]
double mvf_deriv2(int I, int J, arma::vec a, arma::vec mu, arma::mat inv_sigma, arma::vec x){
  int k = a.size();

  double kappa = 1;
  for(int i = 0; i < k; i++) kappa += exp(a(i));
  double mult = -inv_sigma(I, J);
  if(I == J){
    for(int i = 0; i < k + 1; i++) mult -= x(i) * exp(a(I)) * (kappa-exp(a(I))) / (kappa*kappa);
  }else{
    for(int i = 0; i < k + 1; i++) mult += x(i) * exp(a(I) + a(J)) / (kappa*kappa);
  }
  return mult;
}

// [[Rcpp::export]]
double logLike(arma::mat X, arma::mat A, arma::vec mu, arma::mat inv_sigma) {
  double loglik  = 0;
  for(int l = 0; l< A.n_rows; l++) loglik += mvf(A.row(l).t(), mu.row(0).t(), inv_sigma, X.row(l).t());
  return(loglik);
}

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


// [[Rcpp::export]]
arma::vec mvf_maximum(arma::vec x, arma::vec mu, arma::mat inv_sigma, double eps, int max_iter, double prop) {
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

  return out;
}

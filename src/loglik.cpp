// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

const double log2pi = std::log(2.0 * M_PI);

double mvf_norm(arma::vec a, arma::vec mu, arma::mat inv_sigma){
  int k = a.size();

  double norm_const = -0.5 * k * log(2*PI) - 0.5 * log(det(inv_sigma));
  arma::mat log_norm =  -0.5 * (a-mu).t() * inv_sigma * (a-mu);

  double norm = norm_const + log_norm(0);

  return(norm);
}

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

double mvf_multinom_mult(arma::vec a, arma::vec x){
  int k = a.size();

  double kappa = 1;
  for(int i = 0; i < k; i++) kappa += exp(a[i]);
  double multinom = -x[k] * log(kappa);
  for(int i = 0; i < k; i++) multinom += x[i] * ( a[i] - log(kappa));

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

//' @export
// [[Rcpp::export]]
double mvf(arma::vec a, arma::vec mu, arma::mat inv_sigma, arma::vec x){
  int k = a.size();
  int K = k +1;

  double norm_const = -0.5 * k * log(2*PI) - 0.5 * log(det(inv_sigma));
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


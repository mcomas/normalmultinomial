#include <RcppArmadillo.h>

double mvf_norm(arma::vec a, arma::vec mu, arma::mat inv_sigma);

double mvf_multinom_const(arma::vec x);

double mvf_multinom_mult(arma::vec a, arma::vec x);

double mvf_multinom(arma::vec a, arma::vec x);

double mvf(arma::vec a, arma::vec mu, arma::mat inv_sigma, arma::vec x);

double mvf_deriv(int I, arma::vec a, arma::vec mu, arma::mat inv_sigma, arma::vec x);

double mvf_deriv2(int I, int J, arma::vec a, arma::vec mu, arma::mat inv_sigma, arma::vec x);

double logLike(arma::mat X, arma::mat A, arma::vec mu, arma::mat inv_sigma);


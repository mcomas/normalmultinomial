#include <RcppArmadillo.h>


arma::vec choose_starting(arma::vec x, arma::vec mu, arma::mat inv_sigma, double eps, int max_iter, double prop);

Rcpp::List stepEM(arma::mat X, arma::vec mu, arma::mat sigma, int nsim);

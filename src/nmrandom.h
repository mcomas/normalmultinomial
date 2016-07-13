#include <RcppArmadillo.h>

arma::mat c_rnormal(int n, arma::vec mu, arma::mat sigma);
arma::mat c_rmultinomial(arma::mat A, arma::vec size, int seed);
Rcpp::List c_rnormalmultinomial(arma::vec mu, arma::mat sigma, arma::vec size, int seed);

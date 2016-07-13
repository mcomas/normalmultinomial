#include <RcppArmadillo.h>

using namespace Rcpp;

arma::mat c_rnormal(int n, arma::vec mu, arma::mat sigma);
arma::mat c_rmultinomial(arma::mat A, arma::vec size, int seed);
List c_rnormalmultinomial(arma::vec mu, arma::mat sigma, arma::vec size, int seed);

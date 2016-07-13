#include <RcppArmadillo.h>

using namespace Rcpp;

arma::mat rnormal(int n, arma::vec mu, arma::mat sigma);
arma::mat rmultinomial(arma::mat A, arma::vec size, int seed);
List rnormalmultinomial(arma::vec mu, arma::mat sigma, arma::vec size, int seed);

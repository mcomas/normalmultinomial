#include <RcppArmadillo.h>
using namespace Rcpp;

List expected_initial(arma::mat X, arma::vec mu, arma::mat sigma, double se_eps = 0.05);

List expected_guided(arma::mat X, arma::vec mu, arma::mat sigma, arma::mat mu_x, arma::cube sigma_x, double se_eps = 0.05);

#include <RcppArmadillo.h>


arma::vec choose_starting(arma::vec x, arma::vec mu, arma::mat inv_sigma, double eps, int max_iter, double prop);

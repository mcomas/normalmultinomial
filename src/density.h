#include <RcppArmadillo.h>

arma::mat df_x_a(arma::mat X, arma::mat A, arma::vec mu, arma::mat inv_sigma);

arma::mat dnormal(arma::mat A, arma::vec mu, arma::mat inv_sigma);

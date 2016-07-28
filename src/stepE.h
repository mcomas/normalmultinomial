#include <RcppArmadillo.h>


arma::vec choose_starting(arma::vec x, arma::vec mu, arma::mat inv_sigma, double eps, int max_iter, double prop);

Rcpp::List expectedA2(arma::vec x, arma::vec mu, arma::mat sigma, int nsim);

Rcpp::List expectedA4(arma::vec a, arma::vec x, arma::vec mu, arma::mat sigma, int nsim);

Rcpp::List stepEM1(arma::mat X, arma::vec mu, arma::mat sigma, int nsim);

Rcpp::List stepEM2(arma::mat A, arma::mat X, arma::vec mu, arma::mat sigma, int nsim);

double c_dnormalmultinomial(arma::vec x, arma::vec mu, arma::mat sigma, int nsim);

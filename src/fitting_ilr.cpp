// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "nmrandom.h"
#include "stepE.h"
#include "loglik.h"
#include "coda.h"
#include <random>
#include <vector>

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

//' Finds the mean and covariance of a normal multinomial distribution
//'
//' @param X normal-multinomial sample
//' @param nsim number of repetitions for the montecarlo integration process
//' @param niter number of iterations for the EM-algorithm
//' @param prop first 0 imputation
//' @export
// [[Rcpp::export]]
List nm_fit(arma::mat X, int nsim = 100, int niter = 20,
            double prop = 0.66, int version = 0){
  unsigned int K = X.n_cols;
  unsigned int k = K - 1;

  arma::mat mu = ilr_coordinates(arma::sum(X));
  arma::mat sigma = arma::eye<arma::mat>(k,k);

  List fit;
  for(int s = 0 ; s < niter; s++){
    fit = stepEM1(X, mu.t(), sigma, nsim);

    arma::mat mu_next = fit[0];
    arma::mat sigma_next = fit[1];
    mu = mu_next.t();
    sigma_next = sigma_next;
  }

  return Rcpp::List::create(mu, sigma);
}

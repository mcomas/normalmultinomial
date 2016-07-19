// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "nmrandom.h"
#include "stepE.h"
#include "loglik.h"
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
List normalmultinomial_fitting(arma::mat X, int nsim = 100, int niter = 20,
                               double prop = 0.66, int version = 0){
  int n = X.n_rows;
  int K = X.n_cols;
  int k = K - 1;

  // Initialize A, maybe better outsie C++ code
  arma::mat A = arma::mat(n,k);
  for(int i = 0; i < n; i++){
    for(int j = 0; j < k; j++){
      if(X(i,j) != 0 && X(i,k) != 0){
        A(i,j) = log(X(i,j)/X(i,k));
      }else{
        if(X(i,j) == 0 && X(i,k) == 0){
          A(i,j) = 0;
        }else{
          if(X(i,j) == 0){
            A(i,j) = log(prop/X(i,k));
          }else{
            A(i,j) = log(X(i,j)/prop);
          }
        }
      }
    }
  }
  arma::mat mu = mean(A);
  arma::mat sigma = cov(A);
  List fit;
  if(version == 0){
    for(int s = 0 ; s < niter; s++){
      fit = stepEM1(X, mu.t(), sigma, nsim);
      arma::mat mu_next = fit[0];
      arma::mat sigma_next = fit[1];
      mu = mu_next.t();
      sigma_next = sigma_next;
    }
  }
  if(version == 1){
    for(int s = 0 ; s < niter; s++){
      fit = stepEM2(A, X, mu.t(), sigma, nsim);
      arma::mat mu_next = fit[0];
      arma::mat sigma_next = fit[1];
      mu = mu_next.t();
      sigma_next = sigma_next;
    }
  }

  return fit;
}

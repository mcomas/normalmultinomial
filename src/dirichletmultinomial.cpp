// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "loglik.h"
#include <math.h>

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

const double log2pi = std::log(2.0 * M_PI);

//' @export
// [[Rcpp::export]]
double ddirmult(arma::vec x, arma::vec alpha){
  double left = R::gammafn( arma::sum(alpha) ) / R::gammafn( arma::sum(x+alpha) );
  double right = 1;
  for(int i=0;i<x.size();i++){
    right *= R::gammafn(x[i] + alpha[i])/R::gammafn(alpha[i]);
  }
  return( left * right );
}

//' @export
// [[Rcpp::export]]
double gradient_dirmult(arma::vec x, arma::vec alpha){
  double total = 0;
  for(int i=0;i<x.size();i++){
    right *= R::gammafn(x[i] + alpha[i])/R::gammafn(alpha[i]);
  }
  return( left * right );
}

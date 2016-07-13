// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "stepE.h"
#include "loglik.h"
#include <random>
#include <vector>

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

bool DEBUG = true;

using namespace Rcpp;

const double log2pi = std::log(2.0 * M_PI);

//' @export
// [[Rcpp::export]]
arma::vec choose_starting(arma::vec x, arma::vec mu, arma::mat inv_sigma,
                    double eps = 1e-8, int max_iter = 100, double prop = 0.66) {
  int k = x.size() - 1;
  arma::vec a = arma::vec(k);

  for(int j = 0; j < k; j++){
    if(x[j] != 0 && x[k] != 0){
      a[j] = log(a[j]/a[k]);
    }else{
      if(x[j] == 0 && x[k] == 0){
        a[j] = 0;
      }else{
        if(x[j] == 0){
          a[j] = log(prop/x[k]);
        }else{
          a[j] = log(x[j]/prop);
        }
      }
    }
  }

  arma::vec out = arma::vec(a);
  arma::vec deriv = arma::zeros<arma::vec>(k);
  arma::mat deriv2 = arma::zeros<arma::mat>(k, k);
  arma::vec step = arma::zeros<arma::vec>(k);

  int current_iter = 0;
  do{
    current_iter++;
    for(int I=0; I<k; I++){
      deriv[I] =  mvf_deriv(I, out, mu, inv_sigma, x);
      for(int J=0; J<k; J++){
        deriv2(I,J) = mvf_deriv2(I, J, out, mu, inv_sigma, x);
      }
    }
    step = arma::solve(deriv2, deriv);
    out = out - step;
  }while( norm(step, 2) > eps && current_iter < max_iter);
  if(DEBUG){
    Rcout << "Iterations: " << current_iter << std::endl;
  }

  return out;
}



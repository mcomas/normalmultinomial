// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "nmrandom.h"
#include "stepE.h"
#include "loglik.h"
#include "coda.h"
#include "expectations.h"
#include <random>
#include <vector>

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

//' @export
// [[Rcpp::export]]
arma::mat nearestPSD(arma::mat x, double eig_tol = 1e-06, int maxit = 100, double conv_tol = 1e-07, double posd_tol = 1e-03){
  int n = x.n_rows;

  arma::mat Ds = arma::zeros(n,n);

  arma::mat X = arma::mat(x);
  double iter = 0;
  bool converged = false;
  double conv = INFINITY;

  while(iter < maxit and !converged){
    arma::mat Y = arma::mat(X);

    arma::mat R = Y - Ds;
    arma::vec d;
    arma::mat Q;
    arma::eig_sym(d, Q, R);
    int npos = 0;
    for(int i =0; i<n;i++) if(d(i) > eig_tol * d(n-1)) npos++;

    arma::mat Qp = arma::mat(n,npos);

    for(int i =0, I=0; i<n;i++){
      if(d(i) > eig_tol * d(n-1)){
        Qp.col(I++) = Q.col(i) * sqrt(d(i));
      }
    }
    X = Qp * Qp.t();
    Ds = X - R;
    conv = arma::norm(Y-X,"inf") / arma::norm(Y,"inf");
    iter++;
    converged = conv <= conv_tol;
  }

  arma::vec d;
  arma::mat Q;
  arma::eig_sym(d, Q, X);
  double Eps = posd_tol * fabs(d(n-1));
  if(d(0) < Eps){
    for(int i=0; i<n;i++){
      if(d(i) < Eps) d(i) = Eps;
    }
    for(int i=0;i<n;i++){
      Q.col(i) = Q.col(i) * sqrt(d(i));
    }
    arma::vec orig_diag = X.diag();
    X = Q * Q.t();
    arma::vec D = arma::vec(n);
    for(int i=0;i<n;i++){
      D(i) = sqrt((Eps > orig_diag(i) ? Eps : orig_diag(i)) / X(i,i));
    }
    for(int i=0;i<n;i++){
      X.col(i) = X.col(i) * D(i);
      X.row(i) = X.row(i) * D(i);
    }
  }
  return X;
}

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

  unsigned int n = X.n_rows;
  unsigned int K = X.n_cols;
  unsigned int k = K - 1;

  arma::mat mu = ilr_coordinates(arma::sum(X));
  arma::mat sigma = arma::eye<arma::mat>(k,k);

  List fit = expected_initial(X, mu.t(), sigma);
  List fit2;
  arma::mat means;
  arma::mat vars;
  arma::mat MU;
  arma::mat SIGMA;

  arma::mat mu_x = arma::mat(k,n);
  arma::cube sigma_x = arma::cube(k,k,n);

  for(int s=0;s<5;s++){
    arma::mat means_temp = fit(0);
    arma::mat vars_temp = fit(1);
    means = means_temp;
    vars = vars_temp;
    for(int i =0; i < n;i++){
      int J = 1;
      for(int j1 = 0; j1 < k; j1++){
        mu_x(j1,i) = means(i,J);
        J++;
      }
      J = 1 + k;
      for(int j1 = 0; j1 < k; j1++){
        for(int j2 = 0; j2 <=j1; j2++){
          sigma_x(j2,j1,i) = sigma_x(j1,j2,i) = means(i,J);
          J++;
        }
      }
    }

    MU = arma::mean(mu_x, 1);
    SIGMA = arma::mean(sigma_x, 2);
    Rcout << means;
    Rcout << SIGMA;
    SIGMA = nearestPSD(SIGMA - MU * MU.t());

    for(int i =0; i < n;i++){
      sigma_x.slice(i) = arma::eye(k,k); //nearestPSD(sigma_x.slice(i) - mu_x.col(i) * mu_x.col(i).t());
    }
    fit = expected_guided(X, MU, SIGMA, mu_x, sigma_x);

  }

  return Rcpp::List::create(fit);
}

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "coda.h"
#include "loglik.h"
#include <random>
#include <vector>

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;



//' @export
// [[Rcpp::export]]
Rcpp::List expectedMonteCarlo(arma::vec x, arma::vec mu_ilr, arma::mat sigma_ilr, arma::mat Z, arma::vec mu_exp){

  int K = x.size();
  int k = K - 1;
  int nsim = Z.n_rows;

  arma::mat ILR_TO_ALR = ilr_to_alr(K);
  arma::mat ALR_TO_ILR = inv(ILR_TO_ALR);
  arma::mat ALR_TO_ILR_trans = ALR_TO_ILR.t();

  arma::vec mu = ILR_TO_ALR * mu_ilr;
  arma::mat sigma = ILR_TO_ALR * sigma_ilr * ILR_TO_ALR.t();
  arma::mat inv_sigma = inv_sympd(sigma);

  arma::mat sampling_mu =  (ILR_TO_ALR * mu_exp).t();

  arma::mat sampling_sigma = sigma;
  arma::mat inv_sampling_sigma = inv_sigma;

  arma::mat sampling_sigma_chol = arma::chol(sampling_sigma);

  arma::mat Z1 = Z;
  arma::mat Z2 = -Z1;

  arma::mat SAMPLING_MU = arma::repmat(sampling_mu, nsim, 1);
  arma::mat Ap1 = SAMPLING_MU + Z1 * sampling_sigma_chol;
  arma::mat Ap2 = SAMPLING_MU + Z2 * sampling_sigma_chol;

  arma::mat mu12 = (sampling_mu * inv_sigma * sampling_mu.t() - mu.t() * inv_sigma * mu)/2;
  double mult_const = mvf_multinom_const(x);
  arma::mat D = inv_sigma * (mu-sampling_mu.t());

  arma::vec loglik1 = mu12(0,0) + mult_const + mvf_multinom_mult(Ap1, x) + Ap1 * D;
  arma::vec loglik2 = mu12(0,0) + mult_const + mvf_multinom_mult(Ap2, x) + Ap2 * D;

  //Rcpp::Rcout << loglik1 << std::endl;
  //Rcpp::Rcout << loglik2 << std::endl;

  double cmax = std::max(max(loglik1), max(loglik2));

  arma::vec lik1 = exp(loglik1 - cmax);
  arma::vec lik2 = exp(loglik2 - cmax);

  //Rcpp::Rcout << lik1 << std::endl;
  //Rcpp::Rcout << lik2 << std::endl;

  //lik1.replace(arma::datum::nan, 0);
  //lik2.replace(arma::datum::nan, 0);


  // lik_st initialized to ones in case lik1 is too small
  arma::vec lik1_st = lik1 / mean(lik1);
  arma::vec lik2_st = lik2 / mean(lik2);

  //lik1_st.replace(arma::datum::nan, 1);
  //lik2_st.replace(arma::datum::nan, 1);

  //Rcpp::Rcout << lik1_st << std::endl;

  arma::vec M0 = arma::vec(nsim);
  arma::mat M1 = arma::mat(nsim, k);
  arma::cube M2 = arma::cube(k,k, nsim);

  M0 =  (lik1 + lik2) / 2;
  for(int i = 0;i < k; i++){
    arma::mat C1 = Ap1.col(i) % lik1_st;
    arma::mat C2 = Ap2.col(i) % lik2_st;
    M1.col(i) = (C1 + C2) / 2;
    for(int j = 0;j < k; j++){
      M2.tube(i,j) = (C1 % Ap1.col(j) + C2 % Ap2.col(j)) / 2;
    }
  }

  for(int s=0; s<nsim;s++){
    M1.row(s) = M1.row(s) * ALR_TO_ILR_trans;
    M2.slice(s) = ALR_TO_ILR * M2.slice(s) * ALR_TO_ILR_trans;
  }
  return Rcpp::List::create(M0, M1, M2);
}


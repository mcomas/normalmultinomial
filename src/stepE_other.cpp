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
Rcpp::List expectedMetropolis(arma::vec x, arma::vec mu_ilr, arma::mat sigma_ilr, int nsim = 100){
  int K = x.size();
  int k = K - 1;

  arma::mat ILR_TO_ALR = ilr_to_alr(K);
  arma::mat ALR_TO_ILR = inv(ILR_TO_ALR);
  arma::mat ALR_TO_ILR_trans = ALR_TO_ILR.t();

  arma::vec mu = ILR_TO_ALR * mu_ilr;
  arma::mat sigma = ILR_TO_ALR * sigma_ilr * ILR_TO_ALR.t();

  arma::mat Z1 = arma::randn(nsim, k).t();
  arma::vec U = arma::randu(nsim);

  arma::mat inv_sigma = inv_sympd(sigma);

  arma::mat Ap1 = arma::mat(k, nsim);
  Ap1.col(0) = mvf_maximum(x, mu, inv_sigma, 1e-8, 100, 0.66);

  arma::mat M0 = arma::mat(nsim, 1);
  for(int i=1;i<nsim;i++){
    arma::vec x_proposal = Ap1.col(i-1) + Z1.col(i-1);

    double f_prev = exp( mvf(Ap1.col(i-1), mu, inv_sigma, x) );
    double f_next = exp( mvf(x_proposal, mu, inv_sigma, x) );
    double alpha = f_next / f_prev;
    if(1 < alpha){
      Ap1.col(i) = x_proposal;
      M0(i,0) = f_next;
    }else{
      if(U(i) < alpha){
        Ap1.col(i) = x_proposal;
        M0(i,0) = f_next;
      }else{
        Ap1.col(i)=  Ap1.col(i-1);
        M0(i,0) = f_prev;
      }
    }
  }

  arma::mat M1 = arma::mat(nsim, k);
  arma::cube M2 = arma::cube(k,k, nsim);
  for(int i = 0;i < k; i++){
    M1.col(i) = Ap1.row(i).t();
    for(int j = 0;j < k; j++){
      M2.tube(i,j) = Ap1.row(i).t() % Ap1.row(j).t();
    }
  }
  for(int s=0; s<nsim;s++){
    M1.row(s) = M1.row(s) * ALR_TO_ILR_trans;
    M2.slice(s) = ALR_TO_ILR * M2.slice(s) * ALR_TO_ILR_trans;
  }
  return Rcpp::List::create(M0, M1, M2);
}

//' @export
// [[Rcpp::export]]
Rcpp::List expectedMoment1_alr(arma::vec x, arma::vec mu, arma::mat sigma, arma::mat Z){
  int K = x.size();
  int k = K - 1;
  int nsim = Z.n_rows;
  arma::vec x_noz = x.replace(0, 0.66);


  arma::mat inv_sigma = inv_sympd(sigma);

  arma::mat sampling_mu;
  arma::mat sampling_sigma;
  arma::mat inv_sampling_sigma;

  // Sampling from mode selection, if fails the sampling is done using vector x
  sampling_mu =  mvf_maximum(x, mu, inv_sigma, 1e-8, 100, 0.66).t();
  if( sampling_mu.has_nan() ){
    sampling_mu = log( x_noz(arma::span(0,k-1),0) / x_noz(k,0)).t();
  }

  sampling_sigma = sigma;
  inv_sampling_sigma = inv_sigma;

  arma::mat Z1 = Z;
  arma::mat Z2 = -Z1;

  arma::mat Ap1 = arma::repmat(sampling_mu, nsim, 1) + Z1 * arma::chol(sampling_sigma);
  arma::mat Ap2 = arma::repmat(sampling_mu, nsim, 1) + Z2 * arma::chol(sampling_sigma);

  arma::vec lik1;
  arma::vec lik2;

  lik1 = exp( mvf_multinom_const(x) + mvf_multinom_mult(Ap1, x) ) % mvf_norm(Ap1, mu, inv_sigma)  / mvf_norm(Ap1, sampling_mu.t(), inv_sampling_sigma);
  lik1.replace(arma::datum::nan, 0);

  lik2 = exp( mvf_multinom_const(x) + mvf_multinom_mult(Ap2, x) ) % mvf_norm(Ap2, mu, inv_sigma)  / mvf_norm(Ap2, sampling_mu.t(), inv_sampling_sigma);
  lik2.replace(arma::datum::nan, 0);

  arma::vec M0 = arma::vec(nsim);
  arma::mat M1 = arma::mat(nsim, k);

  // lik_st initialized to ones in case lik1 is too small
  arma::vec lik1_st = lik1 / mean(lik1);
  lik1_st.replace(arma::datum::nan, 1);

  arma::vec lik2_st = lik2 / mean(lik2);
  lik2_st.replace(arma::datum::nan, 1);

  M0 =  (lik1 + lik2) / 2;
  for(int i = 0;i < k; i++){
    M1.col(i) = (Ap1.col(i) % lik1_st + Ap2.col(i) % lik2_st) / 2;
  }

  return Rcpp::List::create(M0, M1);
}

//' @export
// [[Rcpp::export]]
Rcpp::List expectedMoment1(arma::vec x, arma::vec mu_ilr, arma::mat sigma_ilr, arma::mat Z){
  int K = x.size();
  int k = K - 1;
  int nsim = Z.n_rows;
  arma::vec x_noz = x.replace(0, 0.66);

  arma::mat ILR_TO_ALR = ilr_to_alr(K);
  arma::mat ALR_TO_ILR = inv(ILR_TO_ALR);
  arma::mat ALR_TO_ILR_trans = ALR_TO_ILR.t();

  arma::vec mu = ILR_TO_ALR * mu_ilr;
  arma::mat sigma = ILR_TO_ALR * sigma_ilr * ILR_TO_ALR.t();

  arma::mat inv_sigma = inv_sympd(sigma);

  arma::mat sampling_mu;
  arma::mat sampling_sigma;
  arma::mat inv_sampling_sigma;

  // Sampling from mode selection, if fails the sampling is done using vector x
  sampling_mu =  mvf_maximum(x, mu, inv_sigma, 1e-8, 100, 0.66).t();
  if( sampling_mu.has_nan() ){
    sampling_mu = log( x_noz(arma::span(0,k-1),0) / x_noz(k,0)).t();
  }

  sampling_sigma = sigma;
  inv_sampling_sigma = inv_sigma;

  arma::mat Z1 = Z;
  arma::mat Ap1 = arma::repmat(sampling_mu, nsim, 1) + Z1 * arma::chol(sampling_sigma);
  //Rcout << sampling_mu << arma::chol(sampling_sigma) << std::endl;
  arma::vec lik1;

  lik1 = exp( mvf_multinom_const(x) + mvf_multinom_mult(Ap1, x) ) % mvf_norm(Ap1, mu, inv_sigma)  / mvf_norm(Ap1, sampling_mu.t(), inv_sampling_sigma);
  lik1.replace(arma::datum::nan, 0);

  arma::vec M0 = arma::vec(nsim);
  arma::mat M1 = arma::mat(nsim, k);
  //arma::cube M2 = arma::cube(k,k, nsim);

  arma::mat Z2 = -Z1;
  arma::mat Ap2 = arma::repmat(sampling_mu, nsim, 1) + Z2 * arma::chol(sampling_sigma);
  arma::vec lik2;

  lik2 = exp( mvf_multinom_const(x) + mvf_multinom_mult(Ap2, x) ) % mvf_norm(Ap2, mu, inv_sigma)  / mvf_norm(Ap2, sampling_mu.t(), inv_sampling_sigma);
  lik2.replace(arma::datum::nan, 0);

  // lik_st initialized to ones in case lik1 is too small
  arma::vec lik1_st = lik1 / mean(lik1);
  lik1_st.replace(arma::datum::nan, 1);
  // arma::vec lik1_st = arma::ones(nsim);
  // if(mean(lik1) > 1e-10){
  //   lik1_st = lik1 / mean(lik1);
  // }
  arma::vec lik2_st = lik2 / mean(lik2);
  lik2_st.replace(arma::datum::nan, 1);
  // arma::vec lik2_st = arma::ones(nsim);
  // if(mean(lik2) > 1e-10){
  //   lik2_st = lik2 / mean(lik2);
  // }

  M0 =  (lik1 + lik2) / 2;
  for(int i = 0;i < k; i++){
    M1.col(i) = (Ap1.col(i) % lik1_st + Ap2.col(i) % lik2_st) / 2;
    // for(int j = 0;j < k; j++){
    //   M2.tube(i,j) = (Ap1.col(i) % Ap1.col(j) % lik1_st + Ap2.col(i) % Ap2.col(j) % lik2_st) / 2;
    // }
  }

  for(int s=0; s<nsim;s++){
    M1.row(s) = M1.row(s) * ALR_TO_ILR_trans;
    //M2.slice(s) = ALR_TO_ILR * M2.slice(s) * ALR_TO_ILR_trans;
  }
  return Rcpp::List::create(M0, M1);
}

//' @export
// [[Rcpp::export]]
Rcpp::List expectedMoment2spherical(arma::vec x, arma::vec mu_ilr, arma::mat sigma_ilr, arma::mat Z){
  int K = x.size();
  int k = K - 1;
  int nsim = Z.n_rows;
  arma::vec x_noz = x.replace(0, 0.66);

  arma::mat ILR_TO_ALR = ilr_to_alr(K);
  arma::mat ALR_TO_ILR = inv(ILR_TO_ALR);
  arma::mat ALR_TO_ILR_trans = ALR_TO_ILR.t();

  arma::vec mu = ILR_TO_ALR * mu_ilr;
  arma::mat sigma = ILR_TO_ALR * sigma_ilr * ILR_TO_ALR.t();

  arma::mat inv_sigma = inv_sympd(sigma);

  arma::mat sampling_mu;
  arma::mat sampling_sigma;
  arma::mat inv_sampling_sigma;

  // Sampling from mode selection, if fails the sampling is done using vector x
  sampling_mu =  mvf_maximum(x, mu, inv_sigma, 1e-8, 100, 0.66).t();
  if( sampling_mu.has_nan() ){
    sampling_mu = log( x_noz(arma::span(0,k-1),0) / x_noz(k,0)).t();
  }

  sampling_sigma = sigma;
  inv_sampling_sigma = inv_sigma;

  arma::mat Z1 = Z;
  arma::mat Ap1 = arma::repmat(sampling_mu, nsim, 1) + Z1 * arma::chol(sampling_sigma);
  //Rcout << sampling_mu << arma::chol(sampling_sigma) << std::endl;
  arma::vec lik1;

  lik1 = exp( mvf_multinom_const(x) + mvf_multinom_mult(Ap1, x) ) % mvf_norm(Ap1, mu, inv_sigma)  / mvf_norm(Ap1, sampling_mu.t(), inv_sampling_sigma);
  lik1.replace(arma::datum::nan, 0);

  arma::vec M0 = arma::vec(nsim);
  //arma::mat M1 = arma::mat(nsim, k);
  arma::cube M2 = arma::cube(k,k, nsim).zeros();

  arma::mat Z2 = -Z1;
  arma::mat Ap2 = arma::repmat(sampling_mu, nsim, 1) + Z2 * arma::chol(sampling_sigma);
  arma::vec lik2;

  lik2 = exp( mvf_multinom_const(x) + mvf_multinom_mult(Ap2, x) ) % mvf_norm(Ap2, mu, inv_sigma)  / mvf_norm(Ap2, sampling_mu.t(), inv_sampling_sigma);
  lik2.replace(arma::datum::nan, 0);

  // lik_st initialized to ones in case lik1 is too small
  arma::vec lik1_st = lik1 / mean(lik1);
  lik1_st.replace(arma::datum::nan, 1);

  arma::vec lik2_st = lik2 / mean(lik2);
  lik2_st.replace(arma::datum::nan, 1);

  M0 =  (lik1 + lik2) / 2;
  for(int i = 0;i < k; i++){
    //M1.col(i) = (Ap1.col(i) % lik1_st + Ap2.col(i) % lik2_st) / 2;
    M2.tube(i,i) = (Ap1.col(i) % Ap1.col(i) % lik1_st + Ap2.col(i) % Ap2.col(i) % lik2_st) / 2;
  }

  for(int s=0; s<nsim;s++){
    //M1.row(s) = M1.row(s) * ALR_TO_ILR_trans;
    M2.slice(s) = ALR_TO_ILR * M2.slice(s) * ALR_TO_ILR_trans;
  }
  return Rcpp::List::create(M0, M2);
}


//' @export
// [[Rcpp::export]]
Rcpp::List expectedMoment2(arma::vec x, arma::vec mu_ilr, arma::mat sigma_ilr, arma::mat Z){
  int K = x.size();
  int k = K - 1;
  int nsim = Z.n_rows;
  arma::vec x_noz = x.replace(0, 0.66);

  arma::mat ILR_TO_ALR = ilr_to_alr(K);
  arma::mat ALR_TO_ILR = inv(ILR_TO_ALR);
  arma::mat ALR_TO_ILR_trans = ALR_TO_ILR.t();

  arma::vec mu = ILR_TO_ALR * mu_ilr;
  arma::mat sigma = ILR_TO_ALR * sigma_ilr * ILR_TO_ALR.t();

  arma::mat inv_sigma = inv_sympd(sigma);

  arma::mat sampling_mu;
  arma::mat sampling_sigma;
  arma::mat inv_sampling_sigma;

  // Sampling from mode selection, if fails the sampling is done using vector x
  sampling_mu =  mvf_maximum(x, mu, inv_sigma, 1e-8, 100, 0.66).t();
  if( sampling_mu.has_nan() ){
    sampling_mu = log( x_noz(arma::span(0,k-1),0) / x_noz(k,0)).t();
  }

  sampling_sigma = sigma;
  inv_sampling_sigma = inv_sigma;

  arma::mat Z1 = Z;
  arma::mat Ap1 = arma::repmat(sampling_mu, nsim, 1) + Z1 * arma::chol(sampling_sigma);
  //Rcout << sampling_mu << arma::chol(sampling_sigma) << std::endl;
  arma::vec lik1;

  lik1 = exp( mvf_multinom_const(x) + mvf_multinom_mult(Ap1, x) ) % mvf_norm(Ap1, mu, inv_sigma)  / mvf_norm(Ap1, sampling_mu.t(), inv_sampling_sigma);
  lik1.replace(arma::datum::nan, 0);

  arma::vec M0 = arma::vec(nsim);
  //arma::mat M1 = arma::mat(nsim, k);
  arma::cube M2 = arma::cube(k,k, nsim);

  arma::mat Z2 = -Z1;
  arma::mat Ap2 = arma::repmat(sampling_mu, nsim, 1) + Z2 * arma::chol(sampling_sigma);
  arma::vec lik2;

  lik2 = exp( mvf_multinom_const(x) + mvf_multinom_mult(Ap2, x) ) % mvf_norm(Ap2, mu, inv_sigma)  / mvf_norm(Ap2, sampling_mu.t(), inv_sampling_sigma);
  lik2.replace(arma::datum::nan, 0);

  // lik_st initialized to ones in case lik1 is too small
  arma::vec lik1_st = lik1 / mean(lik1);
  lik1_st.replace(arma::datum::nan, 1);

  arma::vec lik2_st = lik2 / mean(lik2);
  lik2_st.replace(arma::datum::nan, 1);

  M0 =  (lik1 + lik2) / 2;
  for(int i = 0;i < k; i++){
    //M1.col(i) = (Ap1.col(i) % lik1_st + Ap2.col(i) % lik2_st) / 2;
    for(int j = 0;j < k; j++){
      M2.tube(i,j) = (Ap1.col(i) % Ap1.col(j) % lik1_st + Ap2.col(i) % Ap2.col(j) % lik2_st) / 2;
    }
  }

  for(int s=0; s<nsim;s++){
    //M1.row(s) = M1.row(s) * ALR_TO_ILR_trans;
    M2.slice(s) = ALR_TO_ILR * M2.slice(s) * ALR_TO_ILR_trans;
  }
  return Rcpp::List::create(M0, M2);
}

//' @export
// [[Rcpp::export]]
Rcpp::List expectedMonteCarlo_prev(arma::vec x, arma::vec mu_ilr, arma::mat sigma_ilr, arma::mat Z){

  int K = x.size();
  int k = K - 1;
  int nsim = Z.n_rows;
  arma::vec x_noz = x.replace(0, 0.66);

  arma::mat ILR_TO_ALR = ilr_to_alr(K);
  arma::mat ALR_TO_ILR = inv(ILR_TO_ALR);
  arma::mat ALR_TO_ILR_trans = ALR_TO_ILR.t();

  arma::vec mu = ILR_TO_ALR * mu_ilr;
  arma::mat sigma = ILR_TO_ALR * sigma_ilr * ILR_TO_ALR.t();

  arma::mat inv_sigma = inv_sympd(sigma);

  arma::mat sampling_mu;
  arma::mat sampling_sigma;
  arma::mat inv_sampling_sigma;

  sampling_mu =  mvf_maximum(x, mu, inv_sigma, 1e-8, 100, 0.66).t();
  if( sampling_mu.has_nan() ){
    sampling_mu = log( x_noz(arma::span(0,k-1),0) / x_noz(k,0)).t();
  }

  // Rcout << sampling_mu << std::endl;
  sampling_sigma = sigma;
  inv_sampling_sigma = inv_sigma;

  arma::mat Z1 = Z;
  arma::mat Ap1 = arma::repmat(sampling_mu, nsim, 1) + Z1 * arma::chol(sampling_sigma);
  //Rcout << sampling_mu << arma::chol(sampling_sigma) << std::endl;
  arma::vec lik1;

  lik1 = exp( mvf_multinom_const(x) + mvf_multinom_mult(Ap1, x) ) % mvf_norm(Ap1, mu, inv_sigma)  / mvf_norm(Ap1, sampling_mu.t(), inv_sampling_sigma);
  lik1.replace(arma::datum::nan, 0);
  arma::vec M0 = arma::vec(nsim);
  arma::mat M1 = arma::mat(nsim, k);
  arma::cube M2 = arma::cube(k,k, nsim);


  arma::mat Z2 = -Z1;
  arma::mat Ap2 = arma::repmat(sampling_mu, nsim, 1) + Z2 * arma::chol(sampling_sigma);
  arma::vec lik2;

  lik2 = exp( mvf_multinom_const(x) + mvf_multinom_mult(Ap2, x) ) % mvf_norm(Ap2, mu, inv_sigma)  / mvf_norm(Ap2, sampling_mu.t(), inv_sampling_sigma);
  lik2.replace(arma::datum::nan, 0);
  // lik_st initialized to ones in case lik1 is too small
  arma::vec lik1_st = lik1 / mean(lik1);
  lik1_st.replace(arma::datum::nan, 1);

  arma::vec lik2_st = lik2 / mean(lik2);
  lik2_st.replace(arma::datum::nan, 1);

  M0 =  (lik1 + lik2) / 2;
  for(int i = 0;i < k; i++){
    M1.col(i) = (Ap1.col(i) % lik1_st + Ap2.col(i) % lik2_st) / 2;
    for(int j = 0;j < k; j++){
      M2.tube(i,j) = (Ap1.col(i) % Ap1.col(j) % lik1_st + Ap2.col(i) % Ap2.col(j) % lik2_st) / 2;
    }
  }

  // Rcout << M0 << M1 << M2 << std::endl;
  for(int s=0; s<nsim;s++){
    M1.row(s) = M1.row(s) * ALR_TO_ILR_trans;
    M2.slice(s) = ALR_TO_ILR * M2.slice(s) * ALR_TO_ILR_trans;
  }
  return Rcpp::List::create(M0, M1, M2);
}

//' @export
// [[Rcpp::export]]
Rcpp::List expectedMonteCarlo2(arma::vec x, arma::vec mu_ilr, arma::mat sigma_ilr, arma::mat Z){

  int K = x.size();
  int k = K - 1;
  int nsim = Z.n_rows;
  arma::vec x_noz = x.replace(0, 0.66);

  arma::mat ILR_TO_ALR = ilr_to_alr(K);
  arma::mat ALR_TO_ILR = inv(ILR_TO_ALR);
  arma::mat ALR_TO_ILR_trans = ALR_TO_ILR.t();

  arma::vec mu = ILR_TO_ALR * mu_ilr;
  arma::mat sigma = ILR_TO_ALR * sigma_ilr * ILR_TO_ALR.t();
  arma::mat inv_sigma = inv_sympd(sigma);


  arma::mat sampling_mu =  mvf_maximum(x, mu, inv_sigma, 1e-8, 100, 0.66).t();
  if( sampling_mu.has_nan() ){
    sampling_mu = log( x_noz(arma::span(0,k-1),0) / x_noz(k,0)).t();
  }

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

  arma::vec lik1 = exp( mu12(0,0) + mult_const + mvf_multinom_mult(Ap1, x) + Ap1 * D);
  lik1.replace(arma::datum::nan, 0);
  arma::vec lik2 = exp( mu12(0,0) + mult_const + mvf_multinom_mult(Ap2, x) + Ap2 * D);
  lik2.replace(arma::datum::nan, 0);

  // lik_st initialized to ones in case lik1 is too small
  arma::vec lik1_st = lik1 / mean(lik1);
  lik1_st.replace(arma::datum::nan, 1);
  arma::vec lik2_st = lik2 / mean(lik2);
  lik2_st.replace(arma::datum::nan, 1);


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

//' @export
// [[Rcpp::export]]
Rcpp::List expectedMonteCarlo4(arma::vec x, arma::vec mu_ilr, arma::mat sigma_ilr, arma::mat Z,
                               arma::vec mu_exp, arma::mat sigma_exp){

  int K = x.size();
  int k = K - 1;
  int nsim = Z.n_rows;
  arma::vec x_noz = x.replace(0, 0.66);

  arma::mat ILR_TO_ALR = ilr_to_alr(K);
  arma::mat ALR_TO_ILR = inv(ILR_TO_ALR);
  arma::mat ALR_TO_ILR_trans = ALR_TO_ILR.t();

  arma::vec mu = ILR_TO_ALR * mu_ilr;
  arma::mat sigma = ILR_TO_ALR * sigma_ilr * ILR_TO_ALR.t();
  arma::mat inv_sigma = inv_sympd(sigma);

  arma::mat sampling_mu =  (ILR_TO_ALR * mu_exp).t();

  arma::mat sampling_sigma = ILR_TO_ALR * sigma_exp * ILR_TO_ALR.t();
  arma::mat inv_sampling_sigma = inv_sympd(sampling_sigma);

  arma::mat sampling_sigma_chol = arma::chol(sampling_sigma);

  arma::mat Z1 = Z;
  arma::mat Z2 = -Z1;

  arma::mat SAMPLING_MU = arma::repmat(sampling_mu, nsim, 1);
  arma::mat Ap1 = SAMPLING_MU + Z1 * sampling_sigma_chol;
  arma::mat Ap2 = SAMPLING_MU + Z2 * sampling_sigma_chol;

  double mult_const = mvf_multinom_const(x);

  arma::vec lik1 = exp( mult_const + mvf_multinom_mult(Ap1, x) ) % mvf_norm(Ap1, mu, inv_sigma)  / mvf_norm(Ap1, sampling_mu.t(), inv_sampling_sigma);
  lik1.replace(arma::datum::nan, 0);
  arma::vec lik2 = exp( mult_const + mvf_multinom_mult(Ap2, x) ) % mvf_norm(Ap1, mu, inv_sigma)  / mvf_norm(Ap2, sampling_mu.t(), inv_sampling_sigma);
  lik2.replace(arma::datum::nan, 0);

  // lik_st initialized to ones in case lik1 is too small
  arma::vec lik1_st = lik1 / mean(lik1);
  lik1_st.replace(arma::datum::nan, 1);
  arma::vec lik2_st = lik2 / mean(lik2);
  lik2_st.replace(arma::datum::nan, 1);


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
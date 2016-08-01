// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "stepE.h"
#include "density.h"
#include "loglik.h"
#include "coda.h"
#include <random>
#include <vector>

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

// Using defaulta mu and Sigma
//' @export
// [[Rcpp::export]]
List df_x_1(arma::mat X, arma::vec mu, arma::mat sigma, int nsim = 100){
  int K = X.n_cols;
  int n = X.n_rows;
  int k = K - 1;
  int nsim2 = 2 * nsim;

  arma::mat Z1 = arma::randn(nsim, k);
  arma::mat Z2 = -Z1;
  arma::mat Z = join_cols(Z1, Z2);
  arma::mat Ap = arma::mat(nsim2, k);

  arma::mat inv_sigma = inv_sympd(sigma);

  Ap = arma::repmat(mu.t(), nsim2, 1) + Z * arma::chol(sigma);
  arma::mat M = arma::zeros(n, 1);
  arma::mat M2 = arma::zeros(n, 1);
  for(int i = 0; i < Ap.n_rows; i++){
    arma::mat Mt = df_x_a(X, arma::repmat(Ap.row(i),n,1), mu, inv_sigma) /
                   arma::repmat(dnormal(Ap.row(i), mu, inv_sigma),n,1);

    arma::mat DELTA = (Mt - M);

    M = M + DELTA / (i+1);
    M2 = DELTA  % (Mt - M);
  }
  M2 = M2 / (Ap.n_rows-1);

  return Rcpp::List::create(M, M2);
}

// Using defaulta mu and Sigma
//' @export
// [[Rcpp::export]]
List df_x_2(arma::mat X, arma::vec mu, arma::mat sigma, int nsim = 100){
  int K = X.n_cols;
  int n = X.n_rows;
  int k = K - 1;
  int nsim2 = 2 * nsim;

  arma::mat Z1 = arma::randn(nsim, k);
  arma::mat Z2 = -Z1;
  arma::mat Z = join_cols(Z1, Z2);
  arma::mat Ap = arma::mat(nsim2, k);

  arma::mat inv_sigma = inv_sympd(sigma);

  Ap = arma::repmat(mu.t(), nsim2, 1) + Z * arma::chol(sigma);
  arma::mat M = arma::zeros(n, 1);
  arma::mat M2 = arma::zeros(n, 1);

  for(int i = 0; i < n; i++){
    int NSIM = Ap.n_rows;
    for(int s = 0; s < NSIM; s++){
      arma::mat Mt = df_x_a(X.row(i), Ap.row(s), mu, inv_sigma) /
                     dnormal(Ap.row(s), mu, inv_sigma);
      double DELTA = (Mt(0,0) - M(i,0));
      M(i,0) = M(i,0) + DELTA / (s+1);
      M2(i,0) = DELTA  * (Mt(0,0) - M(i,0));
    }
    M2(i,0)  = M2(i,0) / (NSIM-1);
  }

  return Rcpp::List::create(M, M2);
}

// Using defaulta mu and Sigma
//' @export
// [[Rcpp::export]]
List df_x_3(arma::mat X, arma::vec mu, arma::mat sigma, int nsim = 100){
  int K = X.n_cols;
  int n = X.n_rows;
  int k = K - 1;
  int nsim2 = 2 * nsim;

  arma::mat Z1 = arma::randn(nsim, k);
  arma::mat Z2 = -Z1;
  arma::mat Z = join_cols(Z1, Z2);
  arma::mat Ap0 = arma::mat(nsim2, k);

  arma::mat ALR_TO_ILR = arma::inv(ilr_to_alr(3));
  arma::mat inv_sigma = inv_sympd(sigma);

  Ap0 = Z * arma::chol(sigma);
  arma::mat M = arma::zeros(n, 1);
  arma::mat M2 = arma::zeros(n, 1);

  for(int i = 0; i < n; i++){
    arma::vec mu_max = ALR_TO_ILR * mvf_maximum(X.row(i).t(), mu, inv_sigma, 0.0001, 1000, 0.66);
    arma::mat Ap = arma::repmat(mu_max.t(), nsim2, 1) + Ap0;

    int NSIM = Ap.n_rows;
    for(int s = 0; s < NSIM; s++){
      arma::mat Mt = df_x_a(X.row(i), Ap.row(s), mu, inv_sigma) /
                     dnormal(Ap.row(s), mu_max, inv_sigma);
      double DELTA = (Mt(0,0) - M(i,0));
      M(i,0) = M(i,0) + DELTA / (s+1);
      M2(i,0) = DELTA  * (Mt(0,0) - M(i,0));
    }
    M2(i,0)  = M2(i,0) / (NSIM-1);
  }

  return Rcpp::List::create(M, M2);
}

// Using defaulta mu and Sigma
//' @export
// [[Rcpp::export]]
List expected1(arma::mat X, arma::vec mu, arma::mat sigma, int nsim = 100){
  int K = X.n_cols;
  int n = X.n_rows;
  int k = K - 1;
  int nsim2 = 2 * nsim;

  arma::mat Z1 = arma::randn(nsim, k);
  arma::mat Z2 = -Z1;
  arma::mat Z = join_cols(Z1, Z2);
  arma::mat Ap0 = arma::mat(nsim2, k);

  arma::mat inv_sigma = inv_sympd(sigma);

  Ap0 = Z * arma::chol(sigma);
  arma::mat M = arma::zeros(n, 1 + k + k * (k+1) / 2);
  arma::mat M2 = arma::zeros(n, 1 + k + k * (k+1) / 2);

  for(int i = 0; i < n; i++){
    arma::mat Ap = arma::repmat(mu.t(), nsim2, 1) + Ap0;

    int NSIM = Ap.n_rows;
    for(int s = 0; s < NSIM; s++){
      arma::mat Mt = df_x_a(X.row(i), Ap.row(s), mu, inv_sigma) /
                     dnormal(Ap.row(s), mu, inv_sigma);

      double DELTA = (Mt(0,0) - M(i,0));
      M(i,0) = M(i,0) + DELTA / (s+1);
      M2(i,0) = M2(i,0) + DELTA  * (Mt(0,0) - M(i,0));

      for(int j1 = 0; j1 < k; j1++){
        int J = 1+j1;
        DELTA = (Ap(s,j1) * Mt(0,0) - M(i,J));
        M(i,J) = M(i,J) + DELTA / (s+1);
        M2(i,J) = M2(i,J) + DELTA * (Ap(s,j1) * Mt(0,0) - M(i,J));
      }
      int J = 1 + k;
      for(int j1 = 0; j1 < k; j1++){
        for(int j2 = 0; j2 <=j1; j2++){
          DELTA = (Ap(s,j1) * Ap(s,j2) * Mt(0,0) - M(i,J));
          M(i,J) = M(i,J) + DELTA / (s+1);
          M2(i,J) = M2(i,J) + DELTA * (Ap(s,j1) * Ap(s,j2) * Mt(0,0) - M(i,J));
          J++;
        }
      }
    }
    double temp = M(i,0);
    M.row(i) = M.row(i)/M(i,0);
    M(i,0) = temp;
    M2.row(i) = M2.row(i) / (NSIM-1);
  }

  return Rcpp::List::create(M, M2);
}

// Using defaulta mu and Sigma
//' @export
// [[Rcpp::export]]
List expected2(arma::mat X, arma::vec mu, arma::mat sigma, double se_eps = 0.001){
  int K = X.n_cols;
  int n = X.n_rows;
  int k = K - 1;

  arma::mat M = arma::zeros(n, 1 + k + k * (k+1) / 2);
  arma::mat M2 = arma::zeros(n, 1 + k + k * (k+1) / 2);

  arma::mat inv_sigma = inv_sympd(sigma);

  arma::vec converged = arma::zeros(n);

  arma::mat ALR_TO_ILR = arma::inv(ilr_to_alr(3));

  bool running = true;
  int max_steps = 0, steps = 0;
  do{
    steps++;
    int nsim = 1000000;
    int nsim2 = 2 * nsim;

    arma::mat Z1 = arma::randn(nsim, k);
    arma::mat Z2 = -Z1;
    arma::mat Z = join_cols(Z1, Z2);
    arma::mat Ap0 = Z * arma::chol(sigma);

    for(int i=0; i < n; i++){
      arma::vec mu_max = ALR_TO_ILR * mvf_maximum(X.row(i).t(), mu, inv_sigma, 0.0001, 1000, 0.66);
      unsigned int s=0;
      do{
        arma::mat Ap = mu_max.t() + Ap0.row(s);

        arma::mat Mt = df_x_a(X.row(i), Ap, mu, inv_sigma) /
        dnormal(Ap, mu_max, inv_sigma);
        double DELTA = (Mt(0,0) - M(i,0));
        M(i,0) = M(i,0) + DELTA / (s+1);
        M2(i,0) = M2(i,0) + DELTA  * (Mt(0,0) - M(i,0));

        for(int j1 = 0; j1 < k; j1++){
          int J = 1+j1;
          DELTA = (Ap(0,j1) * Mt(0,0) - M(i,J));
          M(i,J) = M(i,J) + DELTA / (s+1);
          M2(i,J) = M2(i,J) + DELTA * (Ap(0,j1) * Mt(0,0) - M(i,J));
        }

        int J = 1 + k;
        for(int j1 = 0; j1 < k; j1++){
          for(int j2 = 0; j2 <=j1; j2++){
            DELTA = (Ap(j1) * Ap(j2) * Mt(0,0) - M(i,J));
            M(i,J) = M(i,J) + DELTA / (s+1);
            M2(i,J) = M2(i,J) + DELTA * (Ap(j1) * Ap(j2) * Mt(0,0) - M(i,J));
            J++;
          }
        }

        s++;
      } while ( (s < 10) or (s < nsim2 and se_eps < sqrt(M2(i,0))/s) );
      //Rcout << s << std::endl;
      double temp = M(i,0);
      M.row(i) = M.row(i)/M(i,0);
      M(i,0) = temp;
      M2.row(i) = M2.row(i) / (s-1);
    }

  } while(running & (steps < max_steps));

  return Rcpp::List::create(M, M2);

}

// Using defaulta mu and Sigma
//' @export
// [[Rcpp::export]]
List expected3(arma::mat X, arma::vec mu, arma::mat sigma, double se_eps = 0.001){
  int K = X.n_cols;
  int n = X.n_rows;
  int k = K - 1;

  arma::mat M = arma::zeros(n, 1 + k + k * (k+1) / 2);
  arma::mat M2 = arma::zeros(n, 1 + k + k * (k+1) / 2);

  arma::mat inv_sigma = inv_sympd(sigma);

  arma::vec converged = arma::zeros(n);

  arma::mat ALR_TO_ILR = arma::inv(ilr_to_alr(3));

  int nsim = 10000000;
  int nsim2 = 2 * nsim;

  arma::mat Z1 = arma::randn(nsim, k);
  arma::mat Z2 = -Z1;
  arma::mat Z = join_cols(Z1, Z2);
  arma::mat Ap0 = Z * arma::chol(sigma);

  for(int i=0; i < n; i++){
    arma::vec mu_max = ALR_TO_ILR * mvf_maximum(X.row(i).t(), mu, inv_sigma, 0.0001, 1000, 0.66);
    unsigned int s=0;
    double DELTA = 0;
    arma::mat Ap, Mt;
    do{
      Ap = mu_max.t() + Ap0.row(s);
      Mt = df_x_a(X.row(i), Ap, mu, inv_sigma) / dnormal(Ap, mu_max, inv_sigma);
      DELTA = (Mt(0,0) - M(i,0));
      M(i,0) = M(i,0) + DELTA / (s+1);
      M2(i,0) = M2(i,0) + DELTA  * (Mt(0,0) - M(i,0));
      s++;
    } while ( (s < 10) or (s < nsim2 and se_eps < sqrt(M2(i,0))/s) );
    Rcout << s << std::endl;
    M2(i,0) = M2(i,0) / (s-1);

    for(int j1 = 0; j1 < k; j1++){
      int J = 1+j1;
      s = 0;
      do{
        Ap = mu_max.t() + Ap0.row(s);
        Mt = df_x_a(X.row(i), Ap, mu, inv_sigma) / dnormal(Ap, mu_max, inv_sigma);
        DELTA = (Ap(0,j1) * Mt(0,0) - M(i,J));
        M(i,J) = M(i,J) + DELTA / (s+1);
        M2(i,J) = M2(i,J) + DELTA * (Ap(0,j1) * Mt(0,0) - M(i,J));
        s++;
      } while ( (s < 10) or (s < nsim2 and se_eps < sqrt(M2(i,J))/s));
      Rcout << s << std::endl;
      M(i,J) = M(i,J) / M(i,0);
      M2(i,J) = M2(i,J) / (s-1);
    }
    int J = 1 + k;
    for(int j1 = 0; j1 < k; j1++){
      for(int j2 = 0; j2 <=j1; j2++){
        s = 0;
        do{
          Ap = mu_max.t() + Ap0.row(s);
          Mt = df_x_a(X.row(i), Ap, mu, inv_sigma) / dnormal(Ap, mu_max, inv_sigma);

          DELTA = (Ap(j1) * Ap(j2) * Mt(0,0) - M(i,J));
          M(i,J) = M(i,J) + DELTA / (s+1);
          M2(i,J) = M2(i,J) + DELTA * (Ap(j1) * Ap(j2) * Mt(0,0) - M(i,J));
          s++;
        } while ((s < 10) or (s < nsim2 and se_eps < sqrt(M2(i,J))/s));
        Rcout << s << std::endl;
        M(i,J) = M(i,J) / M(i,0);
        M2(i,J) = M2(i,J) / (s-1);
        J++;
      }
    }
  }


  return Rcpp::List::create(M, M2);
}

// New sigma as a parameter
//' @export
// [[Rcpp::export]]
List expected4(arma::mat X, arma::vec mu, arma::mat sigma,
               arma::mat mu_x, arma::cube sigma_x,
               double se_eps = 0.0001){

  int K = X.n_cols;
  int n = X.n_rows;
  int k = K - 1;

  arma::mat M = arma::zeros(n, 1 + k + k * (k+1) / 2);
  arma::mat M2 = arma::zeros(n, 1 + k + k * (k+1) / 2);

  arma::mat inv_sigma = inv_sympd(sigma);
  arma::vec converged = arma::zeros(n);

  arma::mat ALR_TO_ILR = arma::inv(ilr_to_alr(3));

  int nsim = 10000000;
  int nsim2 = 2 * nsim;

  arma::mat Z1 = arma::randn(nsim, k);
  arma::mat Z2 = -Z1;
  arma::mat Z = join_cols(Z1, Z2);

  for(int i=0; i < n; i++){
    arma::mat inv_sigma_x = inv_sympd(sigma_x.slice(i));
    arma::mat Ap0 = Z * arma::chol(sigma_x.slice(i));
    unsigned int s=0;
    double DELTA = 0;
    arma::mat Ap, Mt;
    do{
      Ap = mu_x.col(i).t() + Ap0.row(s);
      Mt = df_x_a(X.row(i), Ap, mu, inv_sigma) / dnormal(Ap, mu_x.col(i), inv_sigma_x);

      DELTA = (Mt(0,0) - M(i,0));
      M(i,0) = M(i,0) + DELTA / (s+1);
      M2(i,0) = M2(i,0) + DELTA  * (Mt(0,0) - M(i,0));
      s++;
    } while ( (s < 10) or (s < nsim2 and se_eps < sqrt(M2(i,0))/s) );
    Rcout << s << std::endl;
    M2(i,0) = M2(i,0) / (s-1);

    for(int j1 = 0; j1 < k; j1++){
      int J = 1+j1;
      s = 0;
      do{
        Ap = mu_x.col(i).t() + Ap0.row(s);
        Mt = df_x_a(X.row(i), Ap, mu, inv_sigma) / dnormal(Ap, mu_x.col(i), inv_sigma_x);
        DELTA = (Ap(0,j1) * Mt(0,0) - M(i,J));
        M(i,J) = M(i,J) + DELTA / (s+1);
        M2(i,J) = M2(i,J) + DELTA * (Ap(0,j1) * Mt(0,0) - M(i,J));
        s++;
      } while ( (s < 10) or (s < nsim2 and se_eps < sqrt(M2(i,J))/s));
      Rcout << s << std::endl;
      M(i,J) = M(i,J) / M(i,0);
      M2(i,J) = M2(i,J) / (s-1);
    }
    int J = 1 + k;
    for(int j1 = 0; j1 < k; j1++){
      for(int j2 = 0; j2 <=j1; j2++){
        s = 0;
        do{
          Ap = mu_x.col(i).t() + Ap0.row(s);
          Mt = df_x_a(X.row(i), Ap, mu, inv_sigma) / dnormal(Ap, mu_x.col(i), inv_sigma_x);

          DELTA = (Ap(j1) * Ap(j2) * Mt(0,0) - M(i,J));
          M(i,J) = M(i,J) + DELTA / (s+1);
          M2(i,J) = M2(i,J) + DELTA * (Ap(j1) * Ap(j2) * Mt(0,0) - M(i,J));
          s++;
        } while ((s < 10) or (s < nsim2 and se_eps < sqrt(M2(i,J))/s));
        Rcout << s << std::endl;
        M(i,J) = M(i,J) / M(i,0);
        M2(i,J) = M2(i,J) / (s-1);
        J++;
      }
    }
  }


  return Rcpp::List::create(M, M2);
}

// Using defaulta mu and Sigma
//' @export
// [[Rcpp::export]]
List expected5(arma::mat X, arma::vec mu, arma::mat sigma,
               arma::mat mu_x, arma::cube sigma_x,
               double se_eps = 0.001){
  int K = X.n_cols;
  int n = X.n_rows;
  int k = K - 1;

  arma::mat M = arma::zeros(n, 1 + k + k * (k+1) / 2);
  arma::mat M2 = arma::zeros(n, 1 + k + k * (k+1) / 2);

  arma::mat inv_sigma = inv_sympd(sigma);

  int nsim = 1000000;
  int nsim2 = 2 * nsim;

  arma::mat Z1 = arma::randn(nsim, k);
  arma::mat Z = arma::mat(nsim2, k);
  for(int i = 0; i < nsim;i++){
    Z.row(2*i) = Z1.row(i);
    Z.row(2*i+1) = -Z1.row(i);
  }

  for(int i=0; i < n; i++){
    arma::mat inv_sigma_x = inv_sympd(sigma_x.slice(i));
    arma::mat Ap0 = Z * arma::chol(sigma_x.slice(i));
    unsigned int s=0;
    double DELTA = 0;
    arma::mat Ap, Mt;
    do{
      Ap = mu_x.col(i).t() + Ap0.row(s);
      Mt = df_x_a(X.row(i), Ap, mu, inv_sigma) / dnormal(Ap, mu_x.col(i), inv_sigma_x);

      DELTA = (Mt(0,0) - M(i,0));
      M(i,0) = M(i,0) + DELTA / (s+1);
      M2(i,0) = M2(i,0) + DELTA  * (Mt(0,0) - M(i,0));
      s++;
    } while ( (s < 10) or (s < nsim2 and se_eps < sqrt(M2(i,0))/s) );
    Rcout << s << std::endl;
    M2(i,0) = M2(i,0) / (s-1);

    for(int j1 = 0; j1 < k; j1++){
      int J = 1+j1;
      s = 0;
      do{
        Ap = mu_x.col(i).t() + Ap0.row(s);
        Mt = df_x_a(X.row(i), Ap, mu, inv_sigma) / dnormal(Ap, mu_x.col(i), inv_sigma_x);

        DELTA = (Ap(j1) * Mt(0,0) / M(i,0) - M(i,J));
        M(i,J) = M(i,J) + DELTA / (s+1);
        M2(i,J) = M2(i,J) + DELTA * (Ap(j1) * Mt(0,0) / M(i,0) - M(i,J));
        s++;
      } while ( (s < 10) or (s < nsim2 and se_eps < sqrt(M2(i,J))/s));
      Rcout << s << std::endl;
      M2(i,J) = M2(i,J) / (s-1);
    }

    int J = 1 + k;
    for(int j1 = 0; j1 < k; j1++){
      for(int j2 = 0; j2 <=j1; j2++){
        s = 0;
        do{
          Ap = mu_x.col(i).t() + Ap0.row(s);
          Mt = df_x_a(X.row(i), Ap, mu, inv_sigma) / dnormal(Ap, mu_x.col(i), inv_sigma_x);

          DELTA = (Ap(j1) * Ap(j2) * Mt(0,0) / M(i,0) - M(i,J));
          M(i,J) = M(i,J) + DELTA / (s+1);
          M2(i,J) = M2(i,J) + DELTA * (Ap(j1) * Ap(j2) * Mt(0,0) / M(i,0) - M(i,J));
          s++;
        } while ((s < 10) or (s < nsim2 and se_eps < sqrt(M2(i,J))/s));
        Rcout << s << std::endl;
        M2(i,J) = M2(i,J) / (s-1);
        J++;
      }
    }
  }


  return Rcpp::List::create(M, M2);
}


// Using defaulta mu and Sigma
//' @export
// [[Rcpp::export]]
List expected6(arma::mat X, arma::vec mu, arma::mat sigma,
               arma::mat mu_x, arma::cube sigma_x,
               double se_eps = 0.001){

  int K = X.n_cols;
  int n = X.n_rows;
  int k = K - 1;

  int nsim = 1000000;
  int nsim2 = 2 * nsim;

  arma::mat Z1 = arma::randn(nsim, k);
  arma::mat Z = arma::mat(nsim2, k);
  for(int i = 0; i < nsim;i++){
    Z.row(2*i) = Z1.row(i);
    Z.row(2*i+1) = -Z1.row(i);
  }

  arma::mat M = arma::zeros(n, 1 + k + k * (k+1) / 2);
  arma::vec M2 = arma::zeros(n);

  arma::mat inv_sigma = inv_sympd(sigma);

  for(int i=0; i < n; i++){
    arma::mat inv_sigma_x = inv_sympd(sigma_x.slice(i));
    arma::mat Ap0 = Z * arma::chol(sigma_x.slice(i));

    unsigned int s=0;
    double DELTA = 0;
    arma::mat Ap, Mt;
    do{
      Ap = mu_x.col(i).t() + Ap0.row(s);
      Mt = df_x_a(X.row(i), Ap, mu, inv_sigma) / dnormal(Ap, mu_x.col(i), inv_sigma_x);

      DELTA = (Mt(0,0) - M(i,0));
      M(i,0) = M(i,0) + DELTA / (s+1);
      M2(i) = M2(i) + DELTA  * (Mt(0,0) - M(i,0));

      int J = 1;
      for(int j1 = 0; j1 < k; j1++){

        double EVAL = Ap(0,j1) * Mt(0,0) / M(i,0);
        DELTA = (EVAL - M(i,J));
        M(i,J) = M(i,J) + DELTA / (s+1);
        //M2(i,J) = M2(i,J) + DELTA * (EVAL - M(i,J));
        J++;
      }
      J = 1 + k;
      for(int j1 = 0; j1 < k; j1++){
        for(int j2 = 0; j2 <=j1; j2++){
          double EVAL = Ap(j1) * Ap(j2) * Mt(0,0) / M(i,0);
          DELTA = (EVAL - M(i,J));
          M(i,J) = M(i,J) + DELTA / (s+1);
          //M2(i,J) = M2(i,J) + DELTA * (EVAL - M(i,J));
          J++;
        }
      }
      s++;
    } while ( (s < 10) or (s < nsim2 and se_eps < sqrt(M2(i))/s) );
    M2(i) = M2(i) / (s-1);

  }


  return Rcpp::List::create(M, M2);
}

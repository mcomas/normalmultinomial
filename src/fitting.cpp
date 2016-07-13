// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "nmrandom.h"
#include <random>
#include <vector>

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

const double log2pi = std::log(2.0 * M_PI);

// [[Rcpp::export]]
double mvf_norm(arma::vec a, arma::vec mu, arma::mat inv_sigma){
  int k = a.size();

  double norm_const = -0.5 * k * log(2*PI) - 0.5 * log(det(inv_sigma));
  arma::mat log_norm =  -0.5 * (a-mu).t() * inv_sigma * (a-mu);

  double norm = norm_const + log_norm(0);

  return(norm);
}

// [[Rcpp::export]]
double mvf_multinom_const(arma::vec x){
  int K = x.size();

  double x_total = 0;
  for(int i = 0; i < K; i++) x_total += x[i];

  double x_total_fact = log(x_total);
  for(int l = 1; l < x_total; l++) x_total_fact += log(l);

  double x_parts_fact = 0;
  for(int i = 0; i < K; i++){
    for(int l = 2; l <= x[i]; l++){
      x_parts_fact += log(l);
    }
  }

  double mult_const = 0;
  mult_const += (x_total_fact - x_parts_fact);
  return(mult_const);
}
// [[Rcpp::export]]
double mvf_multinom_mult(arma::vec a, arma::vec x){
  int k = a.size();

  double kappa = 1;
  for(int i = 0; i < k; i++) kappa += exp(a[i]);
  double multinom = -x[k] * log(kappa);
  for(int i = 0; i < k; i++) multinom += x[i] * ( a[i] - log(kappa));

  return(multinom);
}
// [[Rcpp::export]]
double mvf_multinom(arma::vec a, arma::vec x){
  int k = a.size();
  int K = k +1;

  double x_total = 0;
  for(int i = 0; i < K; i++) x_total += x[i];

  double x_total_fact = log(x_total);
  for(int l = 1; l < x_total; l++) x_total_fact += log(l);

  double x_parts_fact = 0;
  for(int i = 0; i <= K; i++){
    for(int l = 1; l <= x[i]; l++){
      x_parts_fact += log(l);
    }
  }

  double mult_const = 0;
  mult_const += (x_total_fact - x_parts_fact);

  double kappa = 1;
  for(int i = 0; i < k; i++) kappa += exp(a[i]);
  double multinom = -x[k] * log(kappa);
  for(int i = 0; i < k; i++) multinom += x[i] * ( a[i] - log(kappa));

  double mult = mult_const + multinom;
  return(mult);
}

//' Finds the mean and covariance of a normal multinomial distribution
//'
//' @param a aln hidden obeservation
//' @param mu mean parameter for the mean in a aln-normal distribution
//' @param sigma parameter for the sigma in a aln-normal distribution
//' @param x normal-multinomial observation
//' @return Loglikelihood og oberserved data
//' @export
// [[Rcpp::export]]
double mvf(arma::vec a, arma::vec mu, arma::mat inv_sigma, arma::vec x){
  int k = a.size();
  int K = k +1;

  double norm_const = -0.5 * k * log(2*PI) - 0.5 * log(det(inv_sigma));
  arma::mat log_norm =  -0.5 * (a-mu).t() * inv_sigma * (a-mu);

  double norm = norm_const + log_norm(0);

  double x_total = 0;
  for(int i = 0; i < K; i++) x_total += x[i];

  double x_total_fact = log(x_total);
  for(int l = 1; l < x_total; l++) x_total_fact += log(l);

  double x_parts_fact = 0;
  for(int i = 0; i <= K; i++){
    for(int l = 1; l <= x[i]; l++){
      x_parts_fact += log(l);
    }
  }

  double mult_const = 0;
  mult_const += (x_total_fact - x_parts_fact);

  double kappa = 1;
  for(int i = 0; i < k; i++) kappa += exp(a[i]);
  double multinom = -x[k] * log(kappa);
  for(int i = 0; i < k; i++) multinom += x[i] * ( a[i] - log(kappa));

  double mult = mult_const + multinom;
  return(norm + mult);
}

//' Finds the mean and covariance of a normal multinomial distribution
//'
//' @param a aln hidden obeservation
//' @param mu mean parameter for the mean in a aln-normal distribution
//' @param sigma parameter for the sigma in a aln-normal distribution
//' @param x normal-multinomial observation
//' @return Loglikelihood og oberserved data
//' @export
// [[Rcpp::export]]
double mvf_deriv(int I, arma::vec a, arma::vec mu, arma::mat inv_sigma, arma::vec x){
  int k = a.size();
  arma::mat log_norm =  -(a-mu).t() * inv_sigma(arma::span::all, I);
  double kappa = 1;
  for(int i = 0; i < k; i++) kappa += exp(a(i));
  double mult = 0;
  mult += x(I) * (kappa-exp(a(I))) / kappa;
  for(int i = 0; i < I; i++) mult += x(i) * (-exp(a(I))) / kappa;
  for(int i = I+1; i < k+1; i++) mult += x(i) * (-exp(a(I))) / kappa;
  return log_norm(0) + mult;
}

//' Finds the mean and covariance of a normal multinomial distribution
//'
//' @param a aln hidden obeservation
//' @param mu mean parameter for the mean in a aln-normal distribution
//' @param sigma parameter for the sigma in a aln-normal distribution
//' @param x normal-multinomial observation
//' @return Loglikelihood og oberserved data
//' @export
// [[Rcpp::export]]
double mvf_deriv2(int I, int J, arma::vec a, arma::vec mu, arma::mat inv_sigma, arma::vec x){
  int k = a.size();

  double kappa = 1;
  for(int i = 0; i < k; i++) kappa += exp(a(i));
  double mult = -inv_sigma(I, J);
  if(I == J){
    for(int i = 0; i < k + 1; i++) mult -= x(i) * exp(a(I)) * (kappa-exp(a(I))) / (kappa*kappa);
  }else{
    for(int i = 0; i < k + 1; i++) mult += x(i) * exp(a(I) + a(J)) / (kappa*kappa);
  }
  return mult;
}

//' Evaluates the log-likelihood of thenormal multinomial distribution given the hidden observation
//'
//' @param X normal-multinomial observations
//' @param A aln non-observed obeservations
//' @param mu mean parameter for the mean in a aln-normal distribution
//' @param sigma parameter for the sigma in a aln-normal distribution
//' @return Loglikelihood of the oberserved and non-observed data
//' @export
// [[Rcpp::export]]
double logLike(arma::mat X, arma::mat A, arma::vec mu, arma::mat inv_sigma) {
  double loglik  = 0;
  for(int l = 0; l< A.n_rows; l++) loglik += mvf(A.row(l).t(), mu.row(0).t(), inv_sigma, X.row(l).t());
  return(loglik);
}

//' Returns the gaussian more like to generate a sample X given parameters mu and sigma
//'
//' @param A aln hidden obeservations
//' @param mu mean parameter for the mean in a aln-normal distribution
//' @param sigma parameter for the sigma in a aln-normal distribution
//' @param X normal-multinomial observations
//' @return Loglikelihood oberserved data
//' @export
// [[Rcpp::export]]
arma::mat maximize_mvf(arma::vec mu, arma::mat inv_sigma, arma::mat X,
                       double eps = 1e-8, int max_iter = 100, double prop = 0.66) {
  int n = X.n_rows;
  int k = X.n_cols - 1;
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

  arma::mat out = arma::mat(A);
  arma::vec deriv = arma::zeros<arma::vec>(k);
  arma::mat deriv2 = arma::zeros<arma::mat>(k, k);

  arma::vec step = arma::zeros<arma::vec>(k);
  for(int l=0; l < n; l++){
    int current_iter = 0;
    do{
      current_iter++;
      for(int I=0; I<k; I++){
        deriv[I] =  mvf_deriv(I, out.row(l).t(), mu, inv_sigma, X.row(l).t());
        for(int J=0; J<k; J++){
          deriv2(I,J) = mvf_deriv2(I, J, out.row(l).t(), mu, inv_sigma, X.row(l).t());
        }
      }
      step = arma::solve(deriv2, deriv);
      out.row(l) = out.row(l) - step.t();
    }while( norm(step, 2) > eps && current_iter < max_iter);
  }

  return out;

}

//' Compute the mean and covariance of a normal multinomial distribution after itering an EM-algorithm
//'
//' @param X normal-multinomial sample
//' @param mu initial meean estimate
//' @param sigma initial covariance estimate
//' @param nsim number of repetitions for the montecarlo integration process
//' @param niter number of iteration
//' @export
// [[Rcpp::export]]
List EM_step(arma::mat X, arma::mat mu, arma::mat sigma, int nsim = 1000, int niter = 1){
  int n = X.n_rows;
  int K = X.n_cols;
  int k = K - 1;
  int nsim2 = 2 * nsim;
  arma::mat Z1 = arma::randn(nsim, k);
  arma::mat Z2 = -Z1;
  arma::mat Z = join_cols(Z1, Z2);
  arma::mat Ap = arma::mat(nsim2, k);

  arma::mat A = arma::mat(n, k);
  for(int iter = 0; iter < niter; iter++){
    // E-step
    Ap = arma::repmat(mu, nsim2, 1) + Z * arma::chol(sigma);
    arma::vec lik = arma::zeros<arma::vec>(nsim2);

    arma::mat Ap_m1 = arma::zeros(1, k);
    arma::mat Ap_m2 = arma::zeros(k, k);

    arma::vec mult_const = arma::vec(n);
    for(int i=0; i < n; i++){
      mult_const[i] = mvf_multinom_const(X.row(i).t());
    }

    for(int i=0; i < n; i++){
      for(int l=0; l < nsim2; l++){
        lik[l] = mult_const[i] + mvf_multinom_mult(Ap.row(l).t(), X.row(i).t());
      }
      lik = exp(lik);
      for(int j1 = 0; j1 < k; j1++){
        Ap_m1(0,j1) += A(i, j1) = mean(Ap.col(j1) % lik) / mean(lik);
        for(int j2 = 0; j2 < k; j2++){
          Ap_m2(j1, j2) += mean(Ap.col(j1) % Ap.col(j2) % lik) / mean(lik);
        }
      }
    }
    // M-step
    mu = Ap_m1/n;
    sigma = Ap_m2/n - mu.t() * mu;
  }

  arma::mat P = arma::mat(n, k+1);
  for(int i = 0;i<n;i++){
    double accum = 1;
    P(i,k) = 1;
    for(int j = 0; j<k; j++) accum += P(i,j) = exp(A(i,j));
    for(int j = 0; j<=k; j++) P(i,j) /= accum;
  }
  return List::create(mu, sigma, P);
}

//' Compute the mean and covariance of a normal multinomial distribution after itering an EM-algorithm
//'
//' @param X normal-multinomial sample
//' @param mu initial meean estimate
//' @param sigma initial covariance estimate
//' @param nsim number of repetitions for the montecarlo integration process
//' @param niter number of iteration
//' @export
// [[Rcpp::export]]
List E_step1(arma::mat X, arma::mat mu, arma::mat sigma, int nsim = 1000){
  int n = X.n_rows;
  int K = X.n_cols;
  int k = K - 1;

  arma::mat Z = arma::randn(nsim, k);
  arma::mat Ap = arma::mat(nsim, k);

  arma::mat A = arma::mat(n, k);
  arma::mat Ap_m1 = arma::zeros(1, k);
  arma::mat Ap_m2 = arma::zeros(k, k);

  arma::vec lik = arma::vec(nsim);
  arma::vec mult_const = arma::vec(n);
  for(int i=0; i < n; i++){
    mult_const[i] = mvf_multinom_const(X.row(i).t());
  }

  arma::mat MU = arma::repmat(mu, nsim, 1);
  arma::mat SIGMA = arma::chol(sigma);

  // E-step
  Ap = MU + Z * SIGMA;

  lik.zeros();

  for(int i=0; i < n; i++){
    for(int l=0; l < nsim; l++){
      lik[l] = mult_const[i] + mvf_multinom_mult(Ap.row(l).t(), X.row(i).t());
    }
    lik = exp(lik);
  }
  return List::create(Ap, lik);
}

//' Compute the mean and covariance of a normal multinomial distribution after itering an EM-algorithm
//'
//' @param X normal-multinomial sample
//' @param mu initial meean estimate
//' @param sigma initial covariance estimate
//' @param nsim number of repetitions for the montecarlo integration process
//' @param niter number of iteration
//' @export
// [[Rcpp::export]]
List E_step2(arma::mat X, arma::mat mu, arma::mat sigma, int nsim = 1000){
  int n = X.n_rows;
  int K = X.n_cols;
  int k = K - 1;
  nsim = nsim / 2;
  int nsim2 = 2 * nsim;
  arma::mat Z1 = arma::randn(nsim, k);
  arma::mat Z2 = -Z1;
  arma::mat Z = join_cols(Z1, Z2);
  arma::mat Ap = arma::mat(nsim2, k);

  arma::mat A = arma::mat(n, k);
  arma::mat Ap_m1 = arma::zeros(1, k);
  arma::mat Ap_m2 = arma::zeros(k, k);

  arma::vec lik = arma::vec(nsim2);
  arma::vec mult_const = arma::vec(n);
  for(int i=0; i < n; i++){
    mult_const[i] = mvf_multinom_const(X.row(i).t());
  }

  arma::mat MU = arma::repmat(mu, nsim2, 1);
  arma::mat SIGMA = arma::chol(sigma);

  // E-step
  Ap = MU + Z * SIGMA;

  lik.zeros();

  for(int i=0; i < n; i++){
    for(int l=0; l < nsim2; l++){
      lik[l] = mult_const[i] + mvf_multinom_mult(Ap.row(l).t(), X.row(i).t());
    }
    lik = exp(lik);
  }
  return List::create(Ap, lik);
}

//' Compute the mean and covariance of a normal multinomial distribution after itering an EM-algorithm
//'
//' @param X normal-multinomial sample
//' @param mu initial meean estimate
//' @param sigma initial covariance estimate
//' @param nsim number of repetitions for the montecarlo integration process
//' @param niter number of iteration
//' @export
// [[Rcpp::export]]
List E_step3(arma::mat X, arma::mat mu, arma::mat sigma, int nsim = 1000){
  int n = X.n_rows;
  int K = X.n_cols;
  int k = K - 1;
  nsim = nsim / 2;
  int nsim2 = 2 * nsim;
  arma::mat Z1 = arma::randn(nsim, k);
  arma::mat Z2 = -Z1;
  arma::mat Z = join_cols(Z1, Z2);
  arma::mat Ap = arma::mat(nsim2, k);

  arma::mat A = arma::mat(n, k);
  arma::mat Ap_m1 = arma::zeros(1, k);
  arma::mat Ap_m2 = arma::zeros(k, k);

  arma::vec lik = arma::vec(nsim2);
  arma::vec mult_const = arma::vec(n);
  for(int i=0; i < n; i++){
    mult_const[i] = mvf_multinom_const(X.row(i).t());
  }
  arma::mat inv_sigma = sigma.i();
  arma::mat MU_n = maximize_mvf(mu.row(0).t(), inv_sigma, X);
  arma::mat SIGMA = arma::chol(sigma);

  // E-step
  //Ap = MU + Z * SIGMA;

  lik.zeros();

  for(int i=0; i < n; i++){
    arma::mat MU = arma::repmat(MU_n.row(i), nsim2, 1);
    Ap = MU + Z * SIGMA;
    for(int l=0; l < nsim2; l++){
      lik[l] = mult_const[i] + mvf_multinom_mult(Ap.row(l).t(), X.row(i).t()) +
        mvf_norm(Ap.row(l).t(), mu.row(0).t(), inv_sigma) -  mvf_norm(Ap.row(l).t(), MU_n.row(i).t(), inv_sigma);
    }
    lik = exp(lik);
  }
  return List::create(Ap, lik);
}


//' Compute the mean and covariance of a normal multinomial distribution after itering an EM-algorithm
//'
//' @param X normal-multinomial sample
//' @param mu initial meean estimate
//' @param sigma initial covariance estimate
//' @param nsim number of repetitions for the montecarlo integration process
//' @param niter number of iteration
//' @export
// [[Rcpp::export]]
List E_step4(arma::mat X, arma::mat mu, arma::mat sigma, int nsim = 1000){
  int n = X.n_rows;
  int K = X.n_cols;
  int k = K - 1;

  arma::mat inv_sigma = sigma.i();
  arma::mat A = maximize_mvf(mu.row(0).t(), inv_sigma, X);

  arma::mat SIGMA = arma::chol(sigma);

  arma::mat z0 = arma::randn(1, k);
  arma::mat Z0 = arma::repmat(z0, n, 1);
  arma::mat Anext = A + Z0 * SIGMA;
  //temporal
  Anext = maximize_mvf(mu.row(0).t(), inv_sigma, X);
  arma::mat sample = arma::mat(nsim, k);
  for(int iter = 0; iter < nsim; iter++){
    arma::mat z = arma::randn(1, k);
    arma::mat Z = arma::repmat(z, n, 1);

    arma::mat Aprop = A + Z * SIGMA;
    arma::vec alpha = arma::vec(n);
    arma::vec unif = arma::randu(n);

    for(int i = 0; i < n; i++){
      alpha(i) = exp(mvf(Aprop.row(i).t(), mu.row(0).t(), inv_sigma, X.row(i).t()) -
        mvf(Anext.row(i).t(), mu.row(0).t(), inv_sigma, X.row(i).t()));
      if( unif(i) < alpha(i) ){
        Anext(i, arma::span::all) = Aprop(i, arma::span::all);
      }
    }
    sample(iter,arma::span::all) = Anext(0,arma::span::all);
  }


  return List::create(sample);
}


//' Finds the mean and covariance of a normal multinomial distribution
//'
//' @param X normal-multinomial sample
//' @param nsim number of repetitions for the montecarlo integration process
//' @param niter number of iterations for the EM-algorithm
//' @param prop first 0 imputation
//' @export
// [[Rcpp::export]]
List normalmultinomial_fitting(arma::mat X, int nsim = 1000, int niter = 20,
                               double prop = 0.66, int version = 0){
  int n = X.n_rows;
  int K = X.n_cols;
  int k = K - 1;

  //Rcout << "Fitting a normal-multinomial distribution" << std::endl;

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
  List list;
  //if(version == 2){
  //  list =EM_step2(X, mu, sigma, nsim, niter);
  //}else{
    list = EM_step(X, mu, sigma, nsim, niter);
  //}

  return list;
}

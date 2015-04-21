// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <map>
#include <random>
#include <vector>

using namespace Rcpp;

const double log2pi = std::log(2.0 * M_PI);



// [[Rcpp::export]]
arma::mat rnormal(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

// [[Rcpp::export]]
arma::mat rmultinom(arma::mat A, arma::vec size) {
  int k = A.n_cols;
  int n = A.n_rows;
  
  arma::mat res = arma::zeros<arma::mat>(n, k);
  
  std::vector<double> p(k);
  std::random_device rd;
  std::mt19937 gen(1); //r1); //d());
  for(int i=0; i < n; i++){
    for(int j=0; j < k; j++) p[j] = A(i,j);
    std::discrete_distribution<> d(p.begin(),p.end());
    for(int l=0; l < size[i]; l++) ++res(i, d(gen));
  }
  return(res);
}

// [[Rcpp::export]]
arma::mat rnormalmultinomial(arma::vec mu, arma::mat sigma, arma::vec size){
  arma::arma_rng::set_seed(1);
  int n = size.n_elem;
  int k = mu.n_elem;
  arma::mat A = arma::ones<arma::mat>(n, k+1);
  A(arma::span::all, arma::span(0,k-1)) = exp(rnormal(n, mu, sigma));
  return(rmultinom(A, size));
}

// [[Rcpp::export]]
double logLikelihood(arma::mat X, arma::vec mu, arma::mat sigma, int N = 100){
  int n = X.n_rows;
  int K = X.n_cols;
  int k = K - 1;
  
  arma::mat P = arma::ones<arma::mat>(n, K);
  arma::vec loglik_obs = arma::vec(n);
  double kappa;
  arma::vec E = arma::zeros<arma::vec>(n);
  for(int l = 0; l < N; l++){
    arma::mat A = rnormal(n, mu, sigma);
    P(arma::span::all, arma::span(0,k-1)) = exp(A);
    
    for(int i = 0; i < n; i++){
      kappa = 0; 
      for(int j=0; j < K; j++) kappa += P(i, j);
      
      loglik_obs[i] = 0;
      for(int  j=0; j < K; j++) loglik_obs[i] += X(i,j) * log( P(i,j) / kappa);
    }
    E += exp(loglik_obs);
    // NEgative
    P(arma::span::all, arma::span(0,k-1)) = exp(-A);
    
    for(int i = 0; i < n; i++){
      kappa = 0; 
      for(int j=0; j < K; j++) kappa += P(i, j);
      
      loglik_obs[i] = 0;
      for(int  j=0; j < K; j++) loglik_obs[i] += X(i,j) * log( P(i,j) / kappa);
    }
    E += exp(loglik_obs);
    //Rcout << exp(loglik_obs) << std::endl;
  }
  
  E /= (2*N);
  return(sum(log(E)));
}

double factorial2(double x){
  double res = 1;
  do{
    res *= (x--);
  }while(x > 1);
  
  return(res);
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
  //return( constant + log_norm(0) + multinom );
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


double loglike(arma::mat A, arma::vec mu, arma::mat inv_sigma, arma::mat X) {
  double loglik  = 0;
  for(int l = 0; l< A.n_rows; l++) loglik += mvf(A.row(l).t(), mu.row(0).t(), inv_sigma, X.row(l).t());
  return(loglik);
}

//' Finds the mean and covariance of a normal multinomial distribution
//' 
//' @param A aln hidden obeservations
//' @param mu mean parameter for the mean in a aln-normal distribution
//' @param sigma parameter for the sigma in a aln-normal distribution
//' @param X normal-multinomial observations
//' @return Loglikelihood oberserved data
//' @export
//' 
// [[Rcpp::export]]
arma::mat Mstep(arma::mat A, arma::vec mu, arma::mat inv_sigma, arma::mat X, 
                double eps = 1e-8, int max_iter = 100) {
  int n = A.n_rows;
  int k = A.n_cols;
  
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

//' Finds the mean and covariance of a normal multinomial distribution
//' 
//' @param X normal-multinomial sample
//' @param A initial values
//' @export
// [[Rcpp::export]]
List adjustNormalMultinomial_internal(arma::mat X, arma::mat A,
                                      double eps = 1e-04, int iter = 100, double minSigma = 1e-06){
  int n = X.n_rows;
  int K = X.n_cols;
  int k = K - 1;
  
  double tol = pow(eps, 2);
  int cur_iter = 0;
  
  arma::mat mu = mean(A);
  arma::mat sigma = cov(A);
  arma::mat inv_sigma = sigma.i();
    
  double loglik_prev, loglik = 0;
  
  loglik  = 0;
  for(int l = 0; l< n; l++) loglik += mvf(A.row(l).t(), mu.row(0).t(), inv_sigma, X.row(l).t());
  Rcout << "First guess LogLik: " << loglik << std::endl;
  
  do{
    cur_iter++;
    loglik_prev = loglik;
    
    
    for(int s = 0; s < 20; s++){
      A = Mstep(A, mu.row(0).t(), inv_sigma, X);
      
      //loglik  = 0;
      //for(int l = 0; l< n; l++) loglik += mvf(A.row(l).t(), mu.row(0).t(), inv_sigma, X.row(l).t());
      //Rcout << "After Mult LogLik: " << loglik << std::endl;
    }
    
    mu = mean(A);
    sigma = cov(A);
    inv_sigma = sigma.i();
    
    loglik  = 0;
    for(int l = 0; l< n; l++) loglik += mvf(A.row(l).t(), mu.row(0).t(), inv_sigma, X.row(l).t());
    //Rcout << "After Norm LogLik: " << loglik << std::endl;
    
    if( det(sigma) < 1e-20){
      Rcout << "Stop determinant close to zero" << std::endl;
      break;
    }
    Rcout << "LogLik:" << logLikelihood(X, mu.row(0).t(), sigma = sigma, 100000) << std::endl;
  } while (pow(loglik_prev - loglik, 2) > tol && cur_iter < iter); // arma::norm(mu-tmu) > eps &&sigma > minSigma && (pow(tmu-mu, 2) > tol || pow(tsigma-sigma, 2) > tol ) &&
  
  arma::mat A_comp = arma::zeros<arma::mat>(n, K);
  for(int i = 0; i < n; i ++){
    double kappa = 1;
    for(int j = 0; j < k; j++){
      kappa += A_comp(i,j) = exp(A(i,j));
    }
    A_comp(i,k) = 1;
    for(int j = 0; j < K; j++){
      A_comp(i,j) /= kappa;
    }
  }
  return List::create(mu, sigma, A_comp, cur_iter, A, loglik, loglik_prev);
}


//' Finds the mean and covariance of a normal multinomial distribution
//' 
//' @param X normal-multinomial sample
//' @export
// [[Rcpp::export]]
List adjustNormalMultinomial(arma::mat X,
                             double eps = 1e-15, int iter = 100, double minSigma = 1e-06, double prop = 0.3){
  int n = X.n_rows;
  int K = X.n_cols;
  int k = K - 1;
  
  Rcout << "Fitting a normal-multinomial distribution" << std::endl;
  
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
  List list = adjustNormalMultinomial_internal(X, A, eps, iter, minSigma);
  return list;
}



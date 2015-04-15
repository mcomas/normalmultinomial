// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;


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
  arma::mat log_norm =  -0.5 * (a-mu).t() * inv_sigma * (a-mu);
  double kappa = 1;
  for(int i = 0; i < k; i++) kappa += exp(a[i]);
  double multinom = -x[k] * log(kappa);
  for(int i = 0; i < k; i++) multinom += x[i] * ( a[i] - log(kappa));
  
  return( log_norm(0) + multinom );
}


double mvf2(int I, arma::vec a, arma::vec mu, arma::mat inv_sigma, arma::vec x){
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


double mvf3(int I, int J, arma::vec a, arma::vec mu, arma::mat inv_sigma, arma::vec x){
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


arma::mat loglike2(arma::mat A, arma::vec mu, arma::mat sigma, arma::mat X,
                   double eps = 1e-08, int iter = 100) {
  //double tol = pow(eps, 2);
  
  int n = A.n_rows;
  int k = A.n_cols;
  arma::mat inv_sigma = sigma.i();
  arma::mat out = arma::mat(A);
  
  arma::vec deriv = arma::zeros<arma::vec>(k);
  arma::mat deriv2 = arma::zeros<arma::mat>(k, k);
  arma::mat err = arma::mat(1,1);
  
  
  for(int l=0; l < n; l++){
    int cur_iter = 0;
    arma::rowvec temp = arma::zeros<arma::rowvec>(k);
    do{
      cur_iter++;
      temp = arma::rowvec(out.row(l));
      
      for(int I=0; I<k; I++){
        deriv(I) =  mvf2(I, out.row(l).t(), mu, inv_sigma, X.row(l).t());
        for(int J=0; J<k; J++){
          deriv2(I,J) = mvf3(I, J, out.row(l).t(), mu, inv_sigma, X.row(l).t());
        }
      }
      
      out.row(l) = out.row(l) - deriv.t() * deriv2.i();
      
    } while(arma::norm(temp - out.row(l)) > eps && cur_iter < iter); //
  }
  return out;
  
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
arma::mat Mstep(arma::mat A, arma::vec mu, arma::mat inv_sigma, arma::mat X) {
  
  int n = A.n_rows;
  int k = A.n_cols;
  
  arma::mat out = arma::mat(A);
  arma::vec deriv = arma::zeros<arma::vec>(k);
  arma::mat deriv2 = arma::zeros<arma::mat>(k, k);
  
  arma::vec step = arma::zeros<arma::vec>(k);

  for(int l=0; l < n; l++){
    for(int I=0; I<k; I++){
      deriv[I] =  mvf2(I, out.row(l).t(), mu, inv_sigma, X.row(l).t());
      for(int J=0; J<k; J++){
        deriv2(I,J) = mvf3(I, J, out.row(l).t(), mu, inv_sigma, X.row(l).t());
      }
    }
    step = arma::solve(deriv2, deriv);
    out.row(l) = out.row(l) - step.t();
  }
  return out;
  
}

//' Finds the mean and covariance of a normal multinomial distribution
//' 
//' @param X normal-multinomial sample
//' @param A initial values
//' @export
// [[Rcpp::export]]
List adjustNormalMultinomial2(arma::mat X, arma::mat A,
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
  for(int l = 0; l< n; l++) loglik += mvf(A.row(l).t(), mu.row(0).t(), inv_sigma, X.row(l).t());
  do{
    cur_iter++;
    loglik_prev = loglik;
    
    A = Mstep(A, mu.row(0).t(), inv_sigma, X);
    Rcout << "A:" << std::endl << A << std::endl;
    mu = mean(A);
    sigma = cov(A);
    inv_sigma = sigma.i();
    
    loglik  = 0;
    for(int l = 0; l< n; l++) loglik += mvf(A.row(l).t(), mu.row(0).t(), inv_sigma, X.row(l).t());
    Rcout << "LogLik: " << loglik << std::endl;
    
    if( det(sigma) < 1e-20){
      Rcout << "Stop determinant close to zero" << std::endl;
      break;
    }
    
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
                             double eps = 1e-04, int iter = 100, double minSigma = 1e-06){
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
            A(i,j) = log(0.3/X(i,k));
          }else{
            A(i,j) = log(X(i,j)/0.3);
          }          
        }

      }
    }
  }
  List list = adjustNormalMultinomial2(X, A, eps, iter, minSigma);
  return list;
}



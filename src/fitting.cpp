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


double mvf(arma::vec a, arma::vec mu, arma::mat sigma, arma::vec x){
  int k = a.size();
  arma::mat log_norm = - 0.5 * log(det(sigma)) -0.5 * (a-mu).t() * inv(sigma) * (a-mu);
  //arma::mat log_mult =
  
  return( pow(2 * PI, -0.5 * k) * exp(log_norm(0)));
  
  //pow(a-mu, 2) / (2*sigma*sigma) + x * log( exp(a)/(1+exp(a)) ) + x * log( 1/(1+exp(a)) );
}


double mvf2(int I, arma::vec a, arma::vec mu, arma::mat sigma_inv, arma::vec x){
  int k = a.size();
  arma::mat log_norm =  -(a-mu).t() * sigma_inv(arma::span::all, I);
  double kappa = 1;
  for(int i = 0; i < k; i++) kappa += exp(a(i));
  double mult = 0;
  mult += x(I) * (kappa-exp(a(I))) / kappa;
  for(int i = 0; i < I; i++) mult += x(i) * (-exp(a(I))) / kappa;
  for(int i = I+1; i < k+1; i++) mult += x(i) * (-exp(a(I))) / kappa;
  return log_norm(0) + mult;
}

// [[Rcpp::export]]
double mvf3(int I, int J, arma::vec a, arma::vec mu, arma::mat sigma_inv, arma::vec x){
  int k = a.size();
  
  double kappa = 1;
  for(int i = 0; i < k; i++) kappa += exp(a(i));
  double mult = -sigma_inv(I, J);
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


arma::mat mvloglike(arma::mat A, arma::vec mu, arma::mat sigma, arma::mat X,
                    double eps = 1e-08, int iter = 100) {
  //double tol = pow(eps, 2);
  
  int n = A.n_rows;
  int k = A.n_cols;
  arma::mat inv_sigma = sigma.i();
  arma::mat out = arma::mat(A);
  
  arma::vec deriv = arma::zeros<arma::vec>(k);
  arma::mat deriv2 = arma::zeros<arma::mat>(k, k);
  arma::mat err = arma::mat(1,1);
  arma::vec step = arma::zeros<arma::vec>(k);
  for(int l=0; l < n; l++){
    int cur_iter = 0;
    arma::rowvec temp = arma::zeros<arma::rowvec>(k);
    do{
      cur_iter++;
      temp = out.row(l);
      
      for(int I=0; I<k; I++){
        deriv[I] =  mvf2(I, out.row(l).t(), mu, inv_sigma, X.row(l).t());
        for(int J=0; J<k; J++){
          deriv2(I,J) = mvf3(I, J, out.row(l).t(), mu, inv_sigma, X.row(l).t());
        }
      }
      step = arma::solve(deriv2, deriv);
      out.row(l) = out.row(l) - step.t();
      
    } while(arma::norm(step) > eps && cur_iter < iter); //err(0) > tol &&
  }
  return out;
  
}

// [[Rcpp::export]]
List adjustNormalMultinomial(arma::mat X,
                             double eps = 1e-04, int iter = 100, double minSigma = 1e-06){
  int n = X.n_rows;
  int K = X.n_cols;
  int k = K - 1;
  List pars(2);
  arma::mat A = arma::mat(n,k);
  for(int i = 0; i < n; i++){
    for(int j = 0; j < k; j++){
      if(X(i,j) != 0 & X(i,k) != 0){
        A(i,j) = log(X(i,j)/X(i,k));
      }else{
        if(X(i,j) == 0){
          A(i,j) = -2;
        }else{
          A(i,j) = 2;
        }
      }
    }
  }
  
  arma::mat mu = mean(A);
  arma::mat sigma = cov(A);
  arma::mat tmu, tsigma;
  int cur_iter = 0;
  do{
    cur_iter++;
    tmu = arma::mat(mu);
    tsigma = arma::mat(sigma);
    A = mvloglike(A, mu.row(0).t(), sigma, X, eps);
    //mvloglike(A = A, mu = Mu, sigma = Sigma, X = X)
    mu = mean(A);
    sigma = cov(A);
  } while (arma::norm(mu-tmu) > eps && cur_iter < iter); //sigma > minSigma && (pow(tmu-mu, 2) > tol || pow(tsigma-sigma, 2) > tol ) &&
  
  return List::create(mu, sigma, cur_iter);
}


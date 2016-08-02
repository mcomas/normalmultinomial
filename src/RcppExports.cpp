// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// ilr_basis
arma::mat ilr_basis(unsigned int dim);
RcppExport SEXP normalmultinomial_ilr_basis(SEXP dimSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< unsigned int >::type dim(dimSEXP);
    __result = Rcpp::wrap(ilr_basis(dim));
    return __result;
END_RCPP
}
// ilr_basis_simplex
arma::mat ilr_basis_simplex(unsigned int dim);
RcppExport SEXP normalmultinomial_ilr_basis_simplex(SEXP dimSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< unsigned int >::type dim(dimSEXP);
    __result = Rcpp::wrap(ilr_basis_simplex(dim));
    return __result;
END_RCPP
}
// ilr_to_alr
arma::mat ilr_to_alr(unsigned int dim);
RcppExport SEXP normalmultinomial_ilr_to_alr(SEXP dimSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< unsigned int >::type dim(dimSEXP);
    __result = Rcpp::wrap(ilr_to_alr(dim));
    return __result;
END_RCPP
}
// clr_coordinates
arma::mat clr_coordinates(arma::mat X);
RcppExport SEXP normalmultinomial_clr_coordinates(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    __result = Rcpp::wrap(clr_coordinates(X));
    return __result;
END_RCPP
}
// inv_clr_coordinates
arma::mat inv_clr_coordinates(arma::mat clrX);
RcppExport SEXP normalmultinomial_inv_clr_coordinates(SEXP clrXSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type clrX(clrXSEXP);
    __result = Rcpp::wrap(inv_clr_coordinates(clrX));
    return __result;
END_RCPP
}
// ilr_coordinates
arma::mat ilr_coordinates(arma::mat X);
RcppExport SEXP normalmultinomial_ilr_coordinates(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    __result = Rcpp::wrap(ilr_coordinates(X));
    return __result;
END_RCPP
}
// inv_ilr_coordinates
arma::mat inv_ilr_coordinates(arma::mat ilrX);
RcppExport SEXP normalmultinomial_inv_ilr_coordinates(SEXP ilrXSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type ilrX(ilrXSEXP);
    __result = Rcpp::wrap(inv_ilr_coordinates(ilrX));
    return __result;
END_RCPP
}
// ilr_coordinates_with_basis
arma::mat ilr_coordinates_with_basis(arma::mat X, arma::mat B);
RcppExport SEXP normalmultinomial_ilr_coordinates_with_basis(SEXP XSEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B(BSEXP);
    __result = Rcpp::wrap(ilr_coordinates_with_basis(X, B));
    return __result;
END_RCPP
}
// log_dnormal
arma::mat log_dnormal(arma::mat A, arma::vec mu, arma::mat inv_sigma);
RcppExport SEXP normalmultinomial_log_dnormal(SEXP ASEXP, SEXP muSEXP, SEXP inv_sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type inv_sigma(inv_sigmaSEXP);
    __result = Rcpp::wrap(log_dnormal(A, mu, inv_sigma));
    return __result;
END_RCPP
}
// dnormal
arma::mat dnormal(arma::mat A, arma::vec mu, arma::mat inv_sigma);
RcppExport SEXP normalmultinomial_dnormal(SEXP ASEXP, SEXP muSEXP, SEXP inv_sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type inv_sigma(inv_sigmaSEXP);
    __result = Rcpp::wrap(dnormal(A, mu, inv_sigma));
    return __result;
END_RCPP
}
// log_dmultinomial
arma::mat log_dmultinomial(arma::mat X, arma::mat A);
RcppExport SEXP normalmultinomial_log_dmultinomial(SEXP XSEXP, SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    __result = Rcpp::wrap(log_dmultinomial(X, A));
    return __result;
END_RCPP
}
// dmultinomial
arma::mat dmultinomial(arma::mat X, arma::mat A);
RcppExport SEXP normalmultinomial_dmultinomial(SEXP XSEXP, SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    __result = Rcpp::wrap(dmultinomial(X, A));
    return __result;
END_RCPP
}
// log_df_x_a
arma::mat log_df_x_a(arma::mat X, arma::mat A, arma::vec mu, arma::mat inv_sigma);
RcppExport SEXP normalmultinomial_log_df_x_a(SEXP XSEXP, SEXP ASEXP, SEXP muSEXP, SEXP inv_sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type inv_sigma(inv_sigmaSEXP);
    __result = Rcpp::wrap(log_df_x_a(X, A, mu, inv_sigma));
    return __result;
END_RCPP
}
// df_x_a
arma::mat df_x_a(arma::mat X, arma::mat A, arma::vec mu, arma::mat inv_sigma);
RcppExport SEXP normalmultinomial_df_x_a(SEXP XSEXP, SEXP ASEXP, SEXP muSEXP, SEXP inv_sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type inv_sigma(inv_sigmaSEXP);
    __result = Rcpp::wrap(df_x_a(X, A, mu, inv_sigma));
    return __result;
END_RCPP
}
// df_x_1
List df_x_1(arma::mat X, arma::vec mu, arma::mat sigma, int nsim);
RcppExport SEXP normalmultinomial_df_x_1(SEXP XSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP nsimSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    __result = Rcpp::wrap(df_x_1(X, mu, sigma, nsim));
    return __result;
END_RCPP
}
// df_x_2
List df_x_2(arma::mat X, arma::vec mu, arma::mat sigma, int nsim);
RcppExport SEXP normalmultinomial_df_x_2(SEXP XSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP nsimSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    __result = Rcpp::wrap(df_x_2(X, mu, sigma, nsim));
    return __result;
END_RCPP
}
// df_x_3
List df_x_3(arma::mat X, arma::vec mu, arma::mat sigma, int nsim);
RcppExport SEXP normalmultinomial_df_x_3(SEXP XSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP nsimSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    __result = Rcpp::wrap(df_x_3(X, mu, sigma, nsim));
    return __result;
END_RCPP
}
// expected1
List expected1(arma::mat X, arma::vec mu, arma::mat sigma, int nsim);
RcppExport SEXP normalmultinomial_expected1(SEXP XSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP nsimSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    __result = Rcpp::wrap(expected1(X, mu, sigma, nsim));
    return __result;
END_RCPP
}
// expected2
List expected2(arma::mat X, arma::vec mu, arma::mat sigma, double se_eps);
RcppExport SEXP normalmultinomial_expected2(SEXP XSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP se_epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type se_eps(se_epsSEXP);
    __result = Rcpp::wrap(expected2(X, mu, sigma, se_eps));
    return __result;
END_RCPP
}
// expected3
List expected3(arma::mat X, arma::vec mu, arma::mat sigma, double se_eps);
RcppExport SEXP normalmultinomial_expected3(SEXP XSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP se_epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type se_eps(se_epsSEXP);
    __result = Rcpp::wrap(expected3(X, mu, sigma, se_eps));
    return __result;
END_RCPP
}
// expected4
List expected4(arma::mat X, arma::vec mu, arma::mat sigma, arma::mat mu_x, arma::cube sigma_x, double se_eps);
RcppExport SEXP normalmultinomial_expected4(SEXP XSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP mu_xSEXP, SEXP sigma_xSEXP, SEXP se_epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mu_x(mu_xSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type sigma_x(sigma_xSEXP);
    Rcpp::traits::input_parameter< double >::type se_eps(se_epsSEXP);
    __result = Rcpp::wrap(expected4(X, mu, sigma, mu_x, sigma_x, se_eps));
    return __result;
END_RCPP
}
// expected5
List expected5(arma::mat X, arma::vec mu, arma::mat sigma, arma::mat mu_x, arma::cube sigma_x, double se_eps);
RcppExport SEXP normalmultinomial_expected5(SEXP XSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP mu_xSEXP, SEXP sigma_xSEXP, SEXP se_epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mu_x(mu_xSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type sigma_x(sigma_xSEXP);
    Rcpp::traits::input_parameter< double >::type se_eps(se_epsSEXP);
    __result = Rcpp::wrap(expected5(X, mu, sigma, mu_x, sigma_x, se_eps));
    return __result;
END_RCPP
}
// expected6
List expected6(arma::mat X, arma::vec mu, arma::mat sigma, arma::mat mu_x, arma::cube sigma_x, double se_eps);
RcppExport SEXP normalmultinomial_expected6(SEXP XSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP mu_xSEXP, SEXP sigma_xSEXP, SEXP se_epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mu_x(mu_xSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type sigma_x(sigma_xSEXP);
    Rcpp::traits::input_parameter< double >::type se_eps(se_epsSEXP);
    __result = Rcpp::wrap(expected6(X, mu, sigma, mu_x, sigma_x, se_eps));
    return __result;
END_RCPP
}
// expected_initial
List expected_initial(arma::mat X, arma::vec mu, arma::mat sigma, double se_eps);
RcppExport SEXP normalmultinomial_expected_initial(SEXP XSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP se_epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type se_eps(se_epsSEXP);
    __result = Rcpp::wrap(expected_initial(X, mu, sigma, se_eps));
    return __result;
END_RCPP
}
// expected_guided
List expected_guided(arma::mat X, arma::vec mu, arma::mat sigma, arma::mat mu_x, arma::cube sigma_x, double se_eps);
RcppExport SEXP normalmultinomial_expected_guided(SEXP XSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP mu_xSEXP, SEXP sigma_xSEXP, SEXP se_epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mu_x(mu_xSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type sigma_x(sigma_xSEXP);
    Rcpp::traits::input_parameter< double >::type se_eps(se_epsSEXP);
    __result = Rcpp::wrap(expected_guided(X, mu, sigma, mu_x, sigma_x, se_eps));
    return __result;
END_RCPP
}
// normalmultinomial_fitting
List normalmultinomial_fitting(arma::mat X, int nsim, int niter, double prop, int version);
RcppExport SEXP normalmultinomial_normalmultinomial_fitting(SEXP XSEXP, SEXP nsimSEXP, SEXP niterSEXP, SEXP propSEXP, SEXP versionSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< double >::type prop(propSEXP);
    Rcpp::traits::input_parameter< int >::type version(versionSEXP);
    __result = Rcpp::wrap(normalmultinomial_fitting(X, nsim, niter, prop, version));
    return __result;
END_RCPP
}
// nearestPSD
arma::mat nearestPSD(arma::mat x, double eig_tol, int maxit, double conv_tol, double posd_tol);
RcppExport SEXP normalmultinomial_nearestPSD(SEXP xSEXP, SEXP eig_tolSEXP, SEXP maxitSEXP, SEXP conv_tolSEXP, SEXP posd_tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type eig_tol(eig_tolSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< double >::type conv_tol(conv_tolSEXP);
    Rcpp::traits::input_parameter< double >::type posd_tol(posd_tolSEXP);
    __result = Rcpp::wrap(nearestPSD(x, eig_tol, maxit, conv_tol, posd_tol));
    return __result;
END_RCPP
}
// nm_fit
List nm_fit(arma::mat X, int nsim, int niter, double prop, int version);
RcppExport SEXP normalmultinomial_nm_fit(SEXP XSEXP, SEXP nsimSEXP, SEXP niterSEXP, SEXP propSEXP, SEXP versionSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< double >::type prop(propSEXP);
    Rcpp::traits::input_parameter< int >::type version(versionSEXP);
    __result = Rcpp::wrap(nm_fit(X, nsim, niter, prop, version));
    return __result;
END_RCPP
}
// mvf
double mvf(arma::vec a, arma::vec mu, arma::mat inv_sigma, arma::vec x);
RcppExport SEXP normalmultinomial_mvf(SEXP aSEXP, SEXP muSEXP, SEXP inv_sigmaSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type inv_sigma(inv_sigmaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    __result = Rcpp::wrap(mvf(a, mu, inv_sigma, x));
    return __result;
END_RCPP
}
// mvf_deriv
double mvf_deriv(int I, arma::vec a, arma::vec mu, arma::mat inv_sigma, arma::vec x);
RcppExport SEXP normalmultinomial_mvf_deriv(SEXP ISEXP, SEXP aSEXP, SEXP muSEXP, SEXP inv_sigmaSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type I(ISEXP);
    Rcpp::traits::input_parameter< arma::vec >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type inv_sigma(inv_sigmaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    __result = Rcpp::wrap(mvf_deriv(I, a, mu, inv_sigma, x));
    return __result;
END_RCPP
}
// mvf_deriv2
double mvf_deriv2(int I, int J, arma::vec a, arma::vec mu, arma::mat inv_sigma, arma::vec x);
RcppExport SEXP normalmultinomial_mvf_deriv2(SEXP ISEXP, SEXP JSEXP, SEXP aSEXP, SEXP muSEXP, SEXP inv_sigmaSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type I(ISEXP);
    Rcpp::traits::input_parameter< int >::type J(JSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type inv_sigma(inv_sigmaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    __result = Rcpp::wrap(mvf_deriv2(I, J, a, mu, inv_sigma, x));
    return __result;
END_RCPP
}
// logLike
double logLike(arma::mat X, arma::mat A, arma::vec mu, arma::mat inv_sigma);
RcppExport SEXP normalmultinomial_logLike(SEXP XSEXP, SEXP ASEXP, SEXP muSEXP, SEXP inv_sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type inv_sigma(inv_sigmaSEXP);
    __result = Rcpp::wrap(logLike(X, A, mu, inv_sigma));
    return __result;
END_RCPP
}
// hessian
arma::mat hessian(arma::vec a, arma::vec mu, arma::mat inv_sigma, arma::vec x);
RcppExport SEXP normalmultinomial_hessian(SEXP aSEXP, SEXP muSEXP, SEXP inv_sigmaSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type inv_sigma(inv_sigmaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    __result = Rcpp::wrap(hessian(a, mu, inv_sigma, x));
    return __result;
END_RCPP
}
// mvf_maximum
arma::vec mvf_maximum(arma::vec x, arma::vec mu, arma::mat inv_sigma, double eps, int max_iter, double prop);
RcppExport SEXP normalmultinomial_mvf_maximum(SEXP xSEXP, SEXP muSEXP, SEXP inv_sigmaSEXP, SEXP epsSEXP, SEXP max_iterSEXP, SEXP propSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type inv_sigma(inv_sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type prop(propSEXP);
    __result = Rcpp::wrap(mvf_maximum(x, mu, inv_sigma, eps, max_iter, prop));
    return __result;
END_RCPP
}
// c_rnormal
arma::mat c_rnormal(int n, arma::vec mu, arma::mat sigma);
RcppExport SEXP normalmultinomial_c_rnormal(SEXP nSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    __result = Rcpp::wrap(c_rnormal(n, mu, sigma));
    return __result;
END_RCPP
}
// c_rmultinomial
arma::mat c_rmultinomial(arma::mat A, arma::vec size, int seed);
RcppExport SEXP normalmultinomial_c_rmultinomial(SEXP ASEXP, SEXP sizeSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    __result = Rcpp::wrap(c_rmultinomial(A, size, seed));
    return __result;
END_RCPP
}
// c_rnormalmultinomial
List c_rnormalmultinomial(arma::vec mu, arma::mat sigma, arma::vec size, int seed);
RcppExport SEXP normalmultinomial_c_rnormalmultinomial(SEXP muSEXP, SEXP sigmaSEXP, SEXP sizeSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    __result = Rcpp::wrap(c_rnormalmultinomial(mu, sigma, size, seed));
    return __result;
END_RCPP
}
// expectedA1
Rcpp::List expectedA1(arma::vec x, arma::vec mu, arma::mat sigma, int nsim);
RcppExport SEXP normalmultinomial_expectedA1(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP nsimSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    __result = Rcpp::wrap(expectedA1(x, mu, sigma, nsim));
    return __result;
END_RCPP
}
// expectedA2
Rcpp::List expectedA2(arma::vec x, arma::vec mu, arma::mat sigma, int nsim);
RcppExport SEXP normalmultinomial_expectedA2(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP nsimSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    __result = Rcpp::wrap(expectedA2(x, mu, sigma, nsim));
    return __result;
END_RCPP
}
// expectedA5
Rcpp::List expectedA5(arma::vec x, arma::vec mu, arma::mat sigma, arma::vec mu_x, arma::mat sigma_x, int nsim);
RcppExport SEXP normalmultinomial_expectedA5(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP mu_xSEXP, SEXP sigma_xSEXP, SEXP nsimSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu_x(mu_xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma_x(sigma_xSEXP);
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    __result = Rcpp::wrap(expectedA5(x, mu, sigma, mu_x, sigma_x, nsim));
    return __result;
END_RCPP
}
// stepE
arma::mat stepE(arma::mat X, arma::vec mu, arma::mat sigma, int nsim);
RcppExport SEXP normalmultinomial_stepE(SEXP XSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP nsimSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    __result = Rcpp::wrap(stepE(X, mu, sigma, nsim));
    return __result;
END_RCPP
}
// stepEM1
Rcpp::List stepEM1(arma::mat X, arma::vec mu, arma::mat sigma, int nsim);
RcppExport SEXP normalmultinomial_stepEM1(SEXP XSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP nsimSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    __result = Rcpp::wrap(stepEM1(X, mu, sigma, nsim));
    return __result;
END_RCPP
}
// stepEM2
Rcpp::List stepEM2(arma::mat A, arma::mat X, arma::vec mu, arma::mat sigma, int nsim);
RcppExport SEXP normalmultinomial_stepEM2(SEXP ASEXP, SEXP XSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP nsimSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    __result = Rcpp::wrap(stepEM2(A, X, mu, sigma, nsim));
    return __result;
END_RCPP
}
// stepEM3
Rcpp::List stepEM3(arma::mat X, arma::vec mu, arma::mat sigma, int nsim);
RcppExport SEXP normalmultinomial_stepEM3(SEXP XSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP nsimSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    __result = Rcpp::wrap(stepEM3(X, mu, sigma, nsim));
    return __result;
END_RCPP
}
// stepEM4
Rcpp::List stepEM4(arma::mat X, arma::vec mu, arma::mat sigma, arma::mat mu_x, arma::cube sigma_x, int nsim);
RcppExport SEXP normalmultinomial_stepEM4(SEXP XSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP mu_xSEXP, SEXP sigma_xSEXP, SEXP nsimSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mu_x(mu_xSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type sigma_x(sigma_xSEXP);
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    __result = Rcpp::wrap(stepEM4(X, mu, sigma, mu_x, sigma_x, nsim));
    return __result;
END_RCPP
}
// c_dnormalmultinomial
double c_dnormalmultinomial(arma::vec x, arma::vec mu, arma::mat sigma, int nsim);
RcppExport SEXP normalmultinomial_c_dnormalmultinomial(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP nsimSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    __result = Rcpp::wrap(c_dnormalmultinomial(x, mu, sigma, nsim));
    return __result;
END_RCPP
}
// expectedA1_all
Rcpp::List expectedA1_all(arma::vec x, arma::vec mu, arma::mat sigma, int nsim);
RcppExport SEXP normalmultinomial_expectedA1_all(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP nsimSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    __result = Rcpp::wrap(expectedA1_all(x, mu, sigma, nsim));
    return __result;
END_RCPP
}
// expectedA2_all
Rcpp::List expectedA2_all(arma::vec x, arma::vec mu, arma::mat sigma, int nsim);
RcppExport SEXP normalmultinomial_expectedA2_all(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP nsimSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    __result = Rcpp::wrap(expectedA2_all(x, mu, sigma, nsim));
    return __result;
END_RCPP
}
// expectedA3_all
Rcpp::List expectedA3_all(arma::vec x, arma::vec mu, arma::mat sigma, int nsim);
RcppExport SEXP normalmultinomial_expectedA3_all(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP nsimSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    __result = Rcpp::wrap(expectedA3_all(x, mu, sigma, nsim));
    return __result;
END_RCPP
}
// expectedA4_all
Rcpp::List expectedA4_all(arma::vec x, arma::vec mu, arma::mat sigma, int nsim);
RcppExport SEXP normalmultinomial_expectedA4_all(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP nsimSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    __result = Rcpp::wrap(expectedA4_all(x, mu, sigma, nsim));
    return __result;
END_RCPP
}

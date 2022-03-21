#ifndef BASEFUNCTIONS
#define BASEFUNCTIONS

#include <Rcpp.h>

Rcpp::NumericVector log_bound(double mu, double nu, double p);

Rcpp::NumericVector CP_rejection(int n, double mu, double nu);

double CP_unnorm_pdf(Rcpp::NumericVector Y, double mu, double nu, bool out_log);

Rcpp::NumericVector CP_unnorm_lpdf_terms(Rcpp::NumericVector Y, double mu, double nu);

Rcpp::List ratios_IS_cpp(Rcpp::NumericVector lambda, Rcpp::NumericVector nu, double omega, int N_aux);

Rcpp::List est_logZinv_cpp(Rcpp::NumericVector lambda, Rcpp::NumericVector nu, int N_aux);

Rcpp::List rcompois(int n, double mu, double nu);

#endif

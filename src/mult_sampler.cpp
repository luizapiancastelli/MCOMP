#include <RcppArmadillo.h>
#include "basefunctions.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]


double Phi(int dim, NumericVector Y_i, NumericVector ratios, double omega){
  return exp(-omega*Y_i[dim]) - ratios[dim];
}


NumericVector Phi_vector(NumericVector Y_i, NumericVector ratios, double omega){
  
  int D = Y_i.length();
  NumericVector phi_vec(D);
  
  for(int d = 0; d < D; d++){
    phi_vec[d] = Phi(d, Y_i, ratios, omega);
  }
  return phi_vec;
}


double kernel_one_obs(int i_R, NumericMatrix Y, NumericVector ratios, double omega, NumericMatrix delta){
  
  int d = Y.ncol();
  NumericVector Y_vec = Y(i_R-1,_);
  NumericVector phi_vector = Phi_vector(Y_vec, ratios, omega);
  
  arma::mat delta_m = as<arma::mat>(delta); //transforming to arma matrix
  arma::vec phi_vector_v = as<arma::vec>(phi_vector);  //transforming to arma vector
  
  arma::vec kernel = (phi_vector_v.t() * delta_m)*phi_vector_v;
  arma::vec output = 1 + (1/R::choose(d, 2))*kernel;
  
  return output[0];
}


NumericVector CP_unnorm_lpdf(double Y, double mu, double nu){

  NumericVector Y_(1);
  Y_[0] = Y;
  NumericVector f1 = nu*(Y*log(mu) - log(factorial(Y_)));
  return(f1);
}


double lpdf_marginals_i(NumericVector Y, NumericVector lambda, NumericVector nu){
  
  int d = Y.length();
  double marginals = 0;
  
  int dim = 0;
  while(dim < d){
    marginals += CP_unnorm_lpdf(Y[dim], pow(lambda[dim], 1/nu[dim]), nu[dim])[0];
    dim +=1;
  }
  
  return marginals;
  
}


// [[Rcpp::export]]
NumericVector subset_nv(NumericVector vec, int stop){
  return vec[Rcpp::Range(0, stop)];
}


NumericMatrix subset_matrix(NumericMatrix m, int dim){
  return m(Rcpp::Range(0, dim), Rcpp::Range(0, dim));
}


// [[Rcpp::export]]
double asymptotic_z(double lambda, double nu){
  double pi = 3.141592653589793238463;
  double num = exp(nu*pow(lambda,1/nu));
  double den = pow(lambda, (nu-1)/(2*nu))*pow(2*pi,(nu-1)/2)*sqrt(nu);
  return num/den;
}

// [[Rcpp::export]]
double estimate_loglik(NumericMatrix Y, NumericVector lambda, NumericVector nu, NumericMatrix delta, double omega, int N_aux_z, int N_aux_r){

  int n = Y.nrow();
  
  List est_ratios =  ratios_IS_cpp(lambda, nu, omega, N_aux_r);
  NumericVector ratios = est_ratios["ratios"];
  List est_Z = est_logZinv_cpp(lambda, nu, N_aux_z);
  arma::vec log_inv_z = est_Z["loginvZ"];
  
  arma::vec kernel_by_obs(n);
  arma::vec marginal_by_obs(n);
  
  for(int i = 0; i <n; i++){
    kernel_by_obs[i] = kernel_one_obs(i+1, Y, ratios, omega, delta);
    
    NumericVector Y_i = Y(i,_);
    
    marginal_by_obs[i] = lpdf_marginals_i(Y_i, lambda, nu) + sum(log_inv_z);
  }
  
  //List output = List::create(_["kernel"] = kernel_by_obs, 
  //                           _["marginal"] = marginal_by_obs,
  //                           _["ratios"] = ratios,
  //                           _["loginvZ"] = log_inv_z,
  //                           _["test"] = normalisation);
  double output = sum(kernel_by_obs) + sum(marginal_by_obs);
  return output;
}


// [[Rcpp::export]]
arma::vec loglik_given_est(NumericMatrix Y, NumericVector lambda, NumericVector nu, NumericMatrix delta, double omega, NumericVector ratios, arma::vec log_inv_z){
  
  int n = Y.nrow();
  
  arma::vec log_kernel_by_obs(n);
  arma::vec marginal_by_obs(n);
  
  for(int i = 0; i <n; i++){
    log_kernel_by_obs[i] = log(kernel_one_obs(i+1, Y, ratios, omega, delta));
    
    NumericVector Y_i = Y(i,_);
    
    marginal_by_obs[i] = lpdf_marginals_i(Y_i, lambda, nu) + sum(log_inv_z);
  }
  
  //List output = List::create(_["kernel"] = kernel_by_obs, 
  //                           _["marginal"] = marginal_by_obs,
  //                           _["ratios"] = ratios,
  //                           _["loginvZ"] = log_inv_z,
  //                           _["test"] = normalisation);
  arma::vec output = log_kernel_by_obs + marginal_by_obs;
  return output;
}



// [[Rcpp::export]]
List loglik_reg_given(NumericMatrix Y, NumericVector X, NumericVector gamma, NumericVector nu, NumericMatrix delta, double omega, NumericVector ratios, NumericVector log_inv_z){
  
  int n = Y.nrow();
  arma::vec lambda1(2);
  lambda1[0] = exp(gamma[0] + gamma[1]);
  lambda1[1] = exp(gamma[0] + gamma[1] + gamma[2]);
  double lambda2 = exp(gamma[0]);

  arma::vec kernel_by_obs(n);
  arma::vec marginal_nu1(n);
  arma::vec marginal_nu2(n);
  arma::vec marginal_lambda1(n);
  arma::vec marginal_lambda2(n);
  arma::vec ll_i(n);
  NumericVector ratios_i(2);
  ratios_i[1] = ratios[2];
  double log_inv_z2 = log_inv_z[2];
  
  for(int i = 0; i <n; i++){
    
    int d1_index;
    if(X[i] == 0){
      d1_index = 0;
    } else {
      d1_index =1;
    }
    
    ratios_i[0] = ratios[d1_index];
    
    NumericVector Y_i = Y(i,_);
    kernel_by_obs[i] = kernel_one_obs(i+1, Y, ratios_i, omega, delta);
    
    marginal_nu1[i] = -nu[0]*lgamma(Y_i[0]+1) + log_inv_z[d1_index];
    marginal_nu2[i] = -nu[1]*lgamma(Y_i[1]+1) + log_inv_z2;
    
    marginal_lambda1[i] = Y_i[0]*log(lambda1[d1_index]) + log_inv_z[d1_index];
    marginal_lambda2[i] = Y_i[1]*log(lambda2) + log_inv_z2;
    
    ll_i[i] = marginal_lambda1[i] + marginal_lambda2[i] - nu[0]*lgamma(Y_i[0]+1) -nu[1]*lgamma(Y_i[1]+1) + log(kernel_by_obs[i] );
      
  }
  
  return List::create(_["kernel"] = prod(kernel_by_obs), 
                      _["log_kernel"] = sum(log(kernel_by_obs)),
                      _["lambda1"] =  sum(marginal_lambda1) + sum(log(kernel_by_obs)),
                      _["lambda2"] =  sum(marginal_lambda2) + sum(log(kernel_by_obs)),
                      _["nu1"] = sum(marginal_nu1) + sum(log(kernel_by_obs)),
                      _["nu2"] = sum(marginal_nu2) + sum(log(kernel_by_obs)),
                      _["both_lambda"] = sum(marginal_lambda1) +sum(marginal_lambda2)+
                         sum(log(kernel_by_obs)),
                      _["ll_by_obs"] =  ll_i);

  
}

List loglik_reg_est(NumericMatrix Y, NumericVector X, NumericVector gamma, NumericVector nu, NumericMatrix delta, double omega, int N_aux_z, int N_aux_r){
  
  int n = Y.nrow();
  arma::vec lambda1(2);
  lambda1[0] = exp(gamma[0] + gamma[1]);
  lambda1[1] = exp(gamma[0] + gamma[1] + gamma[2]);
  double lambda2 = exp(gamma[0]);
  
  //Auxiliary vectors for using ratios_IS_cpp function:
  NumericVector lambda_vec_aux(3); NumericVector nu_vec_aux(3);
  nu_vec_aux[0]= nu_vec_aux[1]= nu[0]; nu_vec_aux[2] = nu[1];
  lambda_vec_aux[0] = lambda1[0];lambda_vec_aux[1] = lambda1[1]; lambda_vec_aux[2] = lambda2;
  
  List est_ratios =  ratios_IS_cpp(lambda_vec_aux, nu_vec_aux, omega, N_aux_r);
  NumericVector ratios = est_ratios["ratios"];
  
  List est_Z = est_logZinv_cpp(lambda_vec_aux, nu_vec_aux, N_aux_z);
  arma::vec log_inv_z = est_Z["loginvZ"];
  
  arma::vec kernel_by_obs(n);
  arma::vec marginal1(n);
  arma::vec marginal2(n);
  NumericVector ratios_i(2);
  ratios_i[1] = ratios[2];
  double log_inv_z2 = log_inv_z[2];
  
  for(int i = 0; i <n; i++){
    
    int d1_index;
    if(X[i] == 0){
      d1_index = 0;
    } else {
      d1_index =1;
    }
    
    ratios_i[0] = ratios[d1_index];
    
    NumericVector Y_i = Y(i,_);
    kernel_by_obs[i] = kernel_one_obs(i+1, Y, ratios_i, omega, delta);
    
    marginal1[i] = Y_i[0]*log(lambda1[d1_index]) + log_inv_z[d1_index] -nu[0]*lgamma(Y_i[0]+1);
    marginal2[i] = Y_i[1]*log(lambda2) + log_inv_z2 -nu[1]*lgamma(Y_i[1]+1);
  }
  
  return List::create(_["log_marginal1"] = sum(marginal1),
                     _["log_marginal2"] = sum(marginal2),
                     _["logkernel"] = sum(log(kernel_by_obs)));
}


// [[Rcpp::export]]
double kernel(NumericVector Y, NumericVector ratios, double omega, NumericMatrix delta){
  
  int d = Y.length();
  NumericVector phi_vector = Phi_vector(Y, ratios, omega);
  
  arma::mat delta_m = as<arma::mat>(delta); //transforming to arma matrix
  arma::vec phi_vector_v = as<arma::vec>(phi_vector);  //transforming to arma vector
  
  arma::vec kernel = (phi_vector_v.t() * delta_m)*phi_vector_v;
  arma::vec output = 1 + (1/R::choose(d, 2))*kernel;
  
  return output[0];
}


double unorm_pmf_cond_j(int value, NumericVector given, NumericVector lambda, NumericVector nu, NumericMatrix delta, double omega, NumericVector ratios){
  
  //j_R starts from 2
  int j_R = given.length() +1;
  
  NumericVector lambda_j = subset_nv(lambda, j_R-1);
  NumericVector nu_j = subset_nv(nu, j_R-1);
  NumericMatrix delta_j = subset_matrix(delta, j_R-1);
  NumericVector ratios_j = subset_nv(ratios, j_R-1);
  
  NumericVector Y(j_R);
  
  Y[Rcpp::Range(0, j_R-2)] = given;
  Y[j_R-1]= value;
  
  double k_j = kernel(Y, ratios_j, omega, delta_j);
  
  double k_j_minus1;
  if(j_R ==2){
    k_j_minus1 = 1.0;
  } else {
    
    NumericVector ratios_j_minus_1 = subset_nv(ratios, j_R-2);
    NumericMatrix delta_j_minus_1 = subset_matrix(delta, j_R-2);
    NumericVector Y_j_minus_1 = Y[Rcpp::Range(0, j_R-2)];
    
    k_j_minus1= kernel(Y_j_minus_1, ratios_j_minus_1, omega, delta_j_minus_1);
  }
  
  double q = pow(lambda[j_R-1],Y[j_R-1])/pow(gamma(Y[j_R-1]+1), nu[j_R-1]);
  double output =q*(k_j/k_j_minus1);
  
  return output;
  
} 


List Z_cpmf_j_tol(NumericVector given, NumericVector lambda, NumericVector nu, NumericMatrix delta, double omega, NumericVector ratios, double tol){
  
  // Normalising constant of conditional via truncation, for set T
  double Z = 0; int t = 0;
  double increment = 999; 
  double loop = true;
  int j = given.length()+1;
  
  NumericVector marginal_mean(1);
  marginal_mean[0] = pow(lambda[j-1], 1/nu[j-1]);
  
  NumericVector min_iter_vec = ceiling(marginal_mean);
  
  int min_iter= min_iter_vec[0] + 10;
  
  //Start checking after min_iter
  while(loop){
    
    increment =unorm_pmf_cond_j(t, given, lambda, nu, delta, omega, ratios);
    
    loop = (t< min_iter) | (increment > tol);
    
    Z = Z + increment;
    t = t +1;
  }
  
  return List::create(_["t_value"] = t, _["Z"] = Z);
  
}


NumericVector rmcomp_1_T(NumericVector lambda, NumericVector nu, NumericMatrix delta, double omega, NumericVector ratios, int max_it, double tol){
  
  int d = lambda.length(); 
  NumericVector Y = rep(NA_REAL, d);
  
  //Simulate from marginal:
  List Y1_sample = rcompois(1, pow(lambda[0], 1/nu[0]),  nu[0]);
  int y1 = Y1_sample["Y"];
  Y[0] = y1;
  
  NumericVector cpmf_Z = rep(NA_REAL, d);
  NumericVector tval_cpmf_Z = rep(NA_REAL, d);
  
  for(int dim = 2; dim <= d; dim++){ 
    
    //Estimate Z_trunc_dim: needs to be recomputed depending on Y[1:(dim-1)], so estimate n times to get n samples
    List est_cpmf_Z = Z_cpmf_j_tol(Y[Rcpp::Range(0, dim-2)], lambda, nu, delta, omega, ratios, tol);
    cpmf_Z[dim -1] =  est_cpmf_Z["Z"];
    tval_cpmf_Z[dim-1] = est_cpmf_Z["t_value"];
    
    //Inverse Method for Discrete R.V.
    double u = R::runif(0.0,1.0);
    // Rprintf( "\n Uniform sample %f :", u);
    int val = 0; 
    double lsup = unorm_pmf_cond_j(val, Y[Rcpp::Range(0, dim-2)], lambda, nu, delta, omega, ratios)/cpmf_Z[dim -1] ;
    
    // Rprintf( "\n Dimension is %i", dim, "\n");
    // 
    bool found = u < lsup;
    while(!found){
      
      val = val +1; 
      lsup = lsup + unorm_pmf_cond_j(val, Y[Rcpp::Range(0, dim-2)], lambda, nu, delta, omega, ratios)/cpmf_Z[dim-1];
      
      // Rprintf( "\n Lim Sup is %f", lsup, "\n",
      //          "\n Current value %i", val, "\n");
      
      found = u< lsup;
      
      if (val >= max_it) {
        stop("Reached max iter");
      }
      
      
    }
    Y[dim-1] = val;
  }
  return Y;
}


// [[Rcpp::export]]
NumericMatrix rmcomp_T(int n, NumericVector lambda, NumericVector nu, NumericMatrix delta, double omega, int N_r, int max_it = 1000, double tol = 0.001){
  
  int d = lambda.length();
  
  //Estimate the constants
  List est_ratios =  ratios_IS_cpp(lambda, nu, omega, N_r);
  NumericVector ratios = est_ratios["ratios"];
  NumericMatrix Y_sample(n,d);
  
  int sampled =0; 
  while(sampled < n){
    Y_sample(sampled,_)= rmcomp_1_T(lambda, nu, delta, omega, ratios, max_it, tol);
    sampled = sampled + 1;
  }
  
  return Y_sample;
}

// [[Rcpp::export]]
NumericMatrix rmcomp_T_given(int n, NumericVector lambda, NumericVector nu, NumericMatrix delta, double omega, NumericVector ratios, double tol){
  
  int max_it = 10000;
  int d = lambda.length();
  
  //Ratios provided
  NumericMatrix Y_sample(n,d);
  
  int sampled =0; 
  while(sampled < n){
    Y_sample(sampled,_)= rmcomp_1_T(lambda, nu, delta, omega, ratios, max_it, tol);
    sampled = sampled + 1;
  }
  
  return Y_sample;
}

// [[Rcpp::export]]
List un_loglik_given(int dim_R, NumericMatrix Y, NumericVector lambda, NumericVector nu, NumericMatrix delta, double omega, NumericVector ratios){
  
  int n = Y.nrow();
  int dim = dim_R -1;
  
  arma::vec kernel_by_obs(n);
  arma::vec marginal_lambda_dim(n);
  arma::vec marginal_nu_dim(n);
  
  for(int i = 0; i <n; i++){
    kernel_by_obs[i] = kernel_one_obs(i+1, Y, ratios, omega, delta);
    
    NumericVector Y_i = Y(i,_);
    marginal_nu_dim[i] = -nu[dim]*lgamma(Y_i[dim]+1);
    marginal_lambda_dim[i] = Y_i[dim]*log(lambda[dim]);
  }
  
  return List::create(_["kernel"] = prod(kernel_by_obs), 
                      _["log_kernel"] =sum(log(kernel_by_obs)),
                      _["lambda"] =  sum(marginal_lambda_dim) + sum(log(kernel_by_obs)),
                      _["nu"] = sum(marginal_nu_dim) + sum(log(kernel_by_obs)));
  
}


// [[Rcpp::export]]
List exchange_cpp(NumericMatrix Y, NumericVector X, NumericVector gamma, NumericVector nu, NumericMatrix delta, double omega, NumericVector ratios){
  
  int n = Y.nrow();
  arma::vec lambda1(2);
  lambda1[0] = exp(gamma[0] + gamma[1]);
  lambda1[1] = exp(gamma[0] + gamma[1] + gamma[2]);
  double lambda2 = exp(gamma[0]);
  
  arma::vec kernel_by_obs(n);
  arma::vec q_nu1(n);
  arma::vec q_nu2(n);
  arma::vec q_lambda1(n);
  arma::vec q_lambda2(n);
  NumericVector ratios_i(2);
  ratios_i[1] = ratios[2];
  
  for(int i = 0; i <n; i++){
    
    int d1_index;
    if(X[i] == 0){
      d1_index = 0;
    } else {
      d1_index =1;
    }
    
    ratios_i[0] = ratios[d1_index];
    
    NumericVector Y_i = Y(i,_);
    kernel_by_obs[i] = kernel_one_obs(i+1, Y, ratios_i, omega, delta);
    
    q_nu1[i] = -nu[0]*lgamma(Y_i[0]+1);
    q_nu2[i] = -nu[1]*lgamma(Y_i[1]+1);
    
    q_lambda1[i] = Y_i[0]*log(lambda1[d1_index]);
    q_lambda2[i] = Y_i[1]*log(lambda2);
    
  }
  
  return List::create(_["kernel"] = prod(kernel_by_obs), 
                      _["q_lambda1"] =  sum(q_lambda1) + sum(log(kernel_by_obs)),
                      _["q_lambda2"] =  sum(q_lambda2) + sum(log(kernel_by_obs)),
                      _["q_nu1"] = sum(q_nu1) + sum(log(kernel_by_obs)),
                      _["q_nu2"] = sum(q_nu2) + sum(log(kernel_by_obs)),
                      _["both_lambda"] = sum(q_lambda1) +sum(q_lambda2)+
                        sum(log(kernel_by_obs)));
  
  
}

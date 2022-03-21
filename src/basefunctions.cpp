#include "RcppArmadillo.h"
using namespace Rcpp;

// ----------------------------- Functions for rejection sampling --------------------- //
NumericVector log_bound(double mu, double nu, double p){
  
  NumericVector a;
  NumericVector b;
  NumericVector c;
  NumericVector denom;
  NumericVector B;
  
  if(nu < 1){
    a = pow((1-p),(1/nu));
    b= floor(mu/a);
    c = nu*b;
    denom = b*log(1-p) + nu*log(factorial(b));
    B = log(1/p) + c*log(mu) - denom;
    
  } else if(nu >=1){
    
    a = floor(mu);
    B = (nu-1)*( a*log(mu) - log(factorial(a)) );
  }
  
  return(B);
  
}

// [[Rcpp::export]]
NumericVector CP_rejection(int n, double mu, double nu){
  
  int n_accepts = 0;
  NumericVector x(n);
  NumericVector B;
  
  if(nu<1){
    
    double p = 2*nu/( 2*mu*nu +1+nu );
    B = log_bound(mu,nu,p);
    
    while(n_accepts < n) {
      
      //NumericVector u0 = runif(1);
      //NumericVector y= floor(log(u0)/log(1-p)); //Simulate from geometric envelope
      
      NumericVector y = rgeom(1, p);
      
      // NumericVector y =2;
      
      //Calculate acceptance
      NumericVector numer =nu*( y*log(mu) - log(factorial(y)) );
      NumericVector denom = B + y*log(1-p) + log(p);
      
      LogicalVector comp = log(runif(1)) <= numer - denom;
      
      if(comp[0]==true){
        double y_accepted = y[0];
        x[n_accepts] = y_accepted;
        n_accepts += 1;
      }
      
    }
    
  } else {
    
    B = log_bound(mu,nu,1);
    
    while(n_accepts < n) {
      
      NumericVector y= rpois(1, mu); //Simulate from geometric envelope
      NumericVector alpha = (nu-1)*( y*log(mu) - log(factorial(y)) ) - B;
      
      LogicalVector comp = log(runif(1)) <= alpha;
      
      if(comp[0]==true){
        double y_accepted = y[0];
        x[n_accepts] = y_accepted;
        n_accepts += 1;
      }
      
    }
    
  }
  
  return(x);
}


double CP_unnorm_pdf(NumericVector Y, double mu, double nu, bool out_log){
  
  NumericVector f1 = nu*(Y*log(mu) - log(factorial(Y)));
  
  double output;
  if(out_log == true){
    output= sum(f1);
  } else {
    output = sum(exp(f1));
  }
  return(output);
}



NumericVector CP_unnorm_lpdf_terms(NumericVector Y, double mu, double nu){
  
  NumericVector f1 = nu*(Y*log(mu) - log(factorial(Y)));
  
  return(f1);
}


// [[Rcpp::export]]
NumericVector delta_limits_cpp(NumericVector Z){
  double Phi1 = Z[0]; double Phi2 = Z[1];
  NumericVector linf_vec(2);
  NumericVector lsup_vec(2);
  NumericVector limits(2);
  
  linf_vec[0] = (1-Phi1)*(1-Phi2);
  linf_vec[1] = Phi1*Phi2;
  
  lsup_vec[0] = Phi1*(1-Phi2);
  lsup_vec[1] = Phi2*(1-Phi1);
  
  limits[0] = -1/max(linf_vec);
  limits[1] = 1/max(lsup_vec);

  return limits;    
}

// [[Rcpp::export]]
List rcompois(int n, double mu, double nu){
    
    int n_accepts = 0;
    int n_draws = 0;
    NumericVector x(n);
    NumericVector B;
    
    if(nu<1){
      
      double p = 2*nu/( 2*mu*nu +1+nu );
      B = log_bound(mu,nu,p);
      
      while(n_accepts < n) {
        
        NumericVector y = rgeom(1, p);
        n_draws += 1;
        
        //Calculate acceptance
        NumericVector numer =nu*( y*log(mu) - log(factorial(y)) );
        NumericVector denom = B + y*log(1-p) + log(p);
        LogicalVector comp = log(runif(1)) <= numer - denom;
        
        if(comp[0]==true){
          double y_accepted = y[0];
          x[n_accepts] = y_accepted;
          n_accepts += 1;
        }
        
      }
      
    } else {
      
      B = log_bound(mu,nu,1);
      
      while(n_accepts < n) {
        
        NumericVector y= rpois(1, mu); //Simulate from geometric envelope
        n_draws += 1;
        NumericVector alpha = (nu-1)*( y*log(mu) - log(factorial(y)) ) - B;
        
        LogicalVector comp = log(runif(1)) <= alpha;
        
        if(comp[0]==true){
          double y_accepted = y[0];
          x[n_accepts] = y_accepted;
          n_accepts += 1;
        }
        
      }
      
    }
    
    List output = List::create(_["Y"] = x, _["n_draws"] = n_draws);
    return(output);
    
  }

// [[Rcpp::export]]
List ratio_IS_d_cpp(int d_R, NumericVector lambda, NumericVector nu, double omega, int N_aux){
  
  if (d_R > lambda.length()) {
    stop("Exceeds parameter vector dimension.");
  }
  
  int d = d_R-1;
  double mu_d = pow(lambda[d], 1/nu[d]);
  
  List sampler = rcompois(N_aux, mu_d, nu[d]);
  NumericVector q_numerator = CP_unnorm_lpdf_terms(sampler["Y"], pow( exp(-omega)*lambda[d], 1/nu[d]), nu[d]);
  NumericVector q_den = CP_unnorm_lpdf_terms(sampler["Y"], pow(lambda[d],1/nu[d]), nu[d]);
  
  double is = mean(exp(q_numerator - q_den));
  List output = List::create(_["draws"] = sampler["Y"], _["ratio"] = is);
  return(output);
  
} 

// [[Rcpp::export]]
List ratios_IS_cpp(NumericVector lambda, NumericVector nu, double omega, int N_aux){
  
  //Estimate ratio for d=1, ..., D:
  int D = lambda.length();
  NumericVector estimates(D);
  NumericMatrix samples(N_aux, D);
  
  for(int d=1; d<= D; d++){
    List est_d = ratio_IS_d_cpp(d, lambda, nu, omega, N_aux);
    estimates[d-1] = est_d["ratio"];
    NumericVector draws_d = est_d["draws"];
    samples(_,d-1) = draws_d;
  }
  List output = List::create(_["ratios"] = estimates, _["draws"] = samples);
  return(output);
  
}

double logz_est_aux_cpp(double lambda, double nu, double M){
  //Compute estimate given the envelope efficiency M
  
  if (nu < 0) {
    stop("Nu > 0!");
  }
  
  double mu = pow(lambda,1/nu);
  double B;
  double Z_gamma;
  
  if(nu < 1){ //Geometric envelope:
    
    double p = (2*nu)/(2*mu*nu + 1 + nu);
    double exponent = floor( mu/ pow(1-p,1/nu) );
    B =(nu*exponent)*log(mu) - exponent*log(1-p) -nu*lgamma(exponent+1);
    Z_gamma = - log(p);
    
  } else { //Poisson envelope:
     B = (nu-1)*( floor(mu)*log(mu) - lgamma(floor(mu)+1)  );
     Z_gamma = mu;
  }
  double output = log(M) - B - Z_gamma;
  return(output);
}

// [[Rcpp::export]]
List est_logZinv_d_cpp(int d_R, NumericVector lambda, NumericVector nu, int N_aux){
  
  if (d_R > lambda.length()) {
    stop("Exceeds parameter vector dimension.");
  }
  
  int d = d_R-1;
  List run_sampler = rcompois(N_aux, pow(lambda[d],1/nu[d]), nu[d]);
  double n_draws = pow(run_sampler["n_draws"],1.0);
  double M = n_draws/pow(N_aux,1.0);
  
  double estimate = logz_est_aux_cpp(lambda[d], nu[d], M);
    
  List output = List::create(_["draws"] = run_sampler["Y"], _["log_inv_z"] = estimate,
                             _["M"] =M, _["lambda"] = lambda[d]);
  return(output);
}

// [[Rcpp::export]]
List est_logZinv_cpp(NumericVector lambda, NumericVector nu, int N_aux){
  
  //Estimate for d=1, ..., D:
  int D = lambda.length();
  NumericVector estimates(D);
  NumericMatrix samples(N_aux, D);
  
  for(int d=1; d<= D; d++){
    
    List est_d = est_logZinv_d_cpp(d, lambda, nu, N_aux);
    estimates[d-1] = est_d["log_inv_z"];
    NumericVector draws_d = est_d["draws"];
    samples(_,d-1) = draws_d;
  }
  List output = List::create(_["loginvZ"] = estimates, _["draws"] = samples);
  return(output);
  
}
  
// [[Rcpp::export]]
double update_propsd_double(double current_sd, double accept_rate, int nprops, double target_ac = 0.44){
  double adj = (1.0/(1.0*nprops))*(accept_rate - target_ac);
  double new_sd = current_sd*exp(adj);
  return new_sd;
}

// [[Rcpp::export]]
List ln_proposal(double current_val, double sigma){
  
  NumericVector proposed = rlnorm(1, log(current_val), sigma);
  double q_ratio = log(proposed[0]) - log(current_val);
  
  List output =  List::create(_["proposed"] = proposed, _["log_q_ratio"] = q_ratio);
  return(output);
}




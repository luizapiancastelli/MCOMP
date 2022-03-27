initial_reg_premier = function(Y, X){
  
  glm1 = COMPoissonReg::glm.cmp(formula.lambda = Y[,1] ~ X,
                                formula.nu = Y[,1] ~1)
  
  glm2 = COMPoissonReg::glm.cmp(formula.lambda = Y[,2] ~ 1,
                                formula.nu = Y[,2] ~1)
  
  nu1 = exp(coef(glm1)[3])
  nu2 = exp(coef(glm2)[2])
  gamma0 = coef(glm2)[1]
  gamma1 = coef(glm1)[1] -  coef(glm2)[1]
  gamma2 = coef(glm1)[2]
  
  inits = list();
  inits$gamma = as.numeric(c(gamma0, gamma1, gamma2))
  inits$nu = as.numeric(c(nu1, nu2))
  
  delta_init = sign(cor(Y))
  delta_init[lower.tri(delta_init, diag = TRUE)]<-0
  inits$delta = delta_init
  inits$omega =  1
  return(inits)
  
}

compute_lambda = function(gamma){
  l = c(exp(gamma[1] + gamma[2]),
        exp(sum(gamma)),
        exp(gamma[1]))
  return(l)
}

check_pmf_reg = function(parameters, ratios, data){
  Y = data$Y 
  X = data$X
  
  check_ratios = sum(is.na(ratios)) ==0 
  
  ll_new = loglik_reg_given(Y, X, parameters$gamma, parameters$nu, parameters$delta, parameters$omega, ratios, c(NA, NA, NA))
  check_no_na = !(is.na(ll_new$kernel) | is.na(ll_new$log_kernel))
  check_positive= ifelse(check_no_na, ll_new$kernel>0, FALSE)
  
  check = (check_positive + check_no_na + check_positive)== 3
  return(check)
}

exchange_gamma0 = function(Y, X, gamma_current, nu_current, delta_current, omega_current, sigma_gamma0, ratios, draws_r, Y_aux_list, tol, prior_gamma0){
  
  lambda_current = compute_lambda(gamma_current)
  
  d = ncol(Y); n = nrow(Y);
  N_aux_r = nrow(draws_r);
  parameters = list('gamma' = gamma_current,'nu' = nu_current,'delta' = delta_current, 'omega' = omega_current)
  
  gamma_prime = gamma_current #create a copy
  gamma_prime[1] =  rnorm(1, gamma_current[1], sigma_gamma0) #update the element
  
  #Estimate all 3 ratios and save draws:
  lambda_prime = compute_lambda(gamma_prime)
  ratio_estimates = ratios_IS_cpp(lambda_prime, c(nu_current[1], nu_current[1], nu_current[2]), omega_current, N_aux_r) 
  ratios_prime = ratio_estimates$ratios
  draws_prime_r = ratio_estimates$draws
  
  parameters$gamma = gamma_prime
  
  #Exchange draws:
  Y_prime1 = rmcomp_T_given(668, lambda_prime[c(1,3)], nu_current, delta_current, omega_current, ratios_prime[c(1,3)], tol)
  Y_prime2 = rmcomp_T_given(472, lambda_prime[c(2,3)], nu_current, delta_current, omega_current, ratios_prime[c(2,3)], tol)
  Y_aux_list_prime= list("Pre" = Y_prime1, "During" = Y_prime2)
  Y_prime = rbind(Y_aux_list_prime$Pre, Y_aux_list_prime$During)
  
  check = check_pmf_reg(parameters, ratios_prime, data=list("Y" = Y_prime, "X" = X))   #Check if kernel is computable and yields proper pmf
  
  while(check ==FALSE){
    gamma_prime[1] =  rnorm(1, gamma_current[1], sigma_gamma0)
    
    #Estimate all 3 ratios and save draws:
    lambda_prime = compute_lambda(gamma_prime)
    ratio_estimates = ratios_IS_cpp(lambda_prime, c(nu_current[1], nu_current[1], nu_current[2]), omega_current, N_aux_r) 
    ratios_prime = ratio_estimates$ratios
    draws_prime_r = ratio_estimates$draws
    
    parameters$gamma = gamma_prime
    
    #Exchange draws:
    Y_prime1 = rmcomp_T_given(668, lambda_prime[c(1,3)], nu_current, delta_current, omega_current, ratios_prime[c(1,3)], tol)
    Y_prime2 = rmcomp_T_given(472, lambda_prime[c(2,3)], nu_current, delta_current, omega_current, ratios_prime[c(2,3)], tol)
    Y_aux_list_prime= list("Pre" = Y_prime1, "During" = Y_prime2)
    Y_prime = rbind(Y_aux_list_prime$Pre, Y_aux_list_prime$During)
    
    check = check_pmf_reg(parameters, ratios_prime, data=list("Y" = Y_prime, "X" = X))   
  }
  
  #Priors 
  prior_current = dnorm(gamma_current[1], mean = 0, sd = prior_gamma0$sd, log = TRUE)
  prior_prime = dnorm(gamma_prime[1], mean = 0, sd = prior_gamma0$sd, log = TRUE)
  
  #Exchange Y for Yprime
  num = exchange_cpp(Y, X, gamma_prime, nu_current, delta_current, omega_current, ratios_prime)$both_lambda + #{q(Y|prime)}
    exchange_cpp(Y_prime, X, gamma_current, nu_current, delta_current, omega_current, ratios)$both_lambda +    #{q(Y_prime|current)}
    prior_prime
  
  den = exchange_cpp(Y, X, gamma_current, nu_current, delta_current, omega_current, ratios)$both_lambda +  #{q(Y|current)}
    exchange_cpp(Y_prime, X, gamma_prime, nu_current, delta_current, omega_current, ratios_prime)$both_lambda +  #{q(Y_prime|prime)}
    prior_current
  
  #Acceptance:
  alpha = num - den
  
  if(log(runif(1)) < alpha){
    gamma_out = gamma_prime
    ratios_out = ratios_prime
    Y_aux_out = Y_aux_list_prime
    draws_r_out = draws_prime_r
    accepted = 1;
  } else {
    gamma_out = gamma_current
    ratios_out = ratios
    Y_aux_out = Y_aux_list
    draws_r_out = draws_r
    accepted = 0;
  }
  
  output = list("gamma0" = gamma_out[1], "accepted" = accepted, 
                "ratios" = ratios_out, "Y_aux_list" = Y_aux_out, "draws_r" = draws_r_out)
  return(output)
}

exchange_gamma = function(gamma_index, Y, X, gamma_current, nu_current, omega_current, delta_current, ratios, draws_r, sigma_gamma, Y_aux_list,tol, prior_gamma){
  
  lambda_current = compute_lambda(gamma_current)
  
  d = ncol(Y); n = nrow(Y);
  N_aux_r = nrow(draws_r); 
  parameters = list('gamma' = gamma_current,
                    'nu' = nu_current,
                    'delta' = delta_current,
                    'omega' = omega_current)
  
  gamma_prime = gamma_current
  gamma_prime[gamma_index+1] = rnorm(1, mean = gamma_current[gamma_index+1], sd = sigma_gamma[gamma_index])
  lambda_prime = compute_lambda(gamma_prime)
  ratios_prime = ratios; draws_prime_r = draws_r;
  
  #Estimate ratios[1:2] and save draws:
  r_first_element = ratio_IS_d_cpp(1, lambda_prime, c(nu_current[1], nu_current[1], nu_current[2]), omega_current, N_aux_r)
  r_second_element = ratio_IS_d_cpp(2, lambda_prime, c(nu_current[1], nu_current[1], nu_current[2]), omega_current, N_aux_r)
  
  ratios_prime[1] = r_first_element$ratio; 
  draws_prime_r[,1]= r_first_element$draws
  
  ratios_prime[2] = r_second_element$ratio; 
  draws_prime_r[,2] = r_second_element$draws
  
  #Exchange draws
  if(gamma_index==1){
    Y_aux_list_prime = list("Pre" = rmcomp_T_given(668, lambda_prime[c(1,3)], nu_current, delta_current, omega_current, ratios_prime[c(1,3)], tol),
                            "During" = rmcomp_T_given(472, lambda_prime[c(2,3)], nu_current, delta_current, omega_current, ratios_prime[c(2,3)], tol))
  } else if(gamma_index ==2){
    Y_aux_list_prime= Y_aux_list;
    Y_aux_list_prime$During = rmcomp_T_given(472, lambda_prime[c(2,3)], nu_current, delta_current, omega_current, ratios_prime[c(2,3)], tol)
  }
  
  Y_prime = rbind(Y_aux_list_prime$Pre, Y_aux_list_prime$During)
  
  #PMF check
  parameters$gamma = gamma_prime
  check = check_pmf_reg(parameters, ratios_prime, data = list('Y' = Y_prime, 'X'=X))  
  
  while(check == FALSE){
    gamma_prime[gamma_index+1] = rnorm(1, mean = gamma_current[gamma_index+1], sd = sigma_gamma[gamma_index])
    lambda_prime = compute_lambda(gamma_prime)
    ratios_prime = ratios; draws_prime_r = draws_r;
    
    #Estimate ratios[1:2] and save draws:
    r_first_element = ratio_IS_d_cpp(1, lambda_prime, c(nu_current[1], nu_current[1], nu_current[2]), omega_current, N_aux_r)
    r_second_element = ratio_IS_d_cpp(2, lambda_prime, c(nu_current[1], nu_current[1], nu_current[2]), omega_current, N_aux_r)
    
    ratios_prime[1] = r_first_element$ratio; 
    draws_prime_r[,1]= r_first_element$draws
    
    ratios_prime[2] = r_second_element$ratio; 
    draws_prime_r[,2] = r_second_element$draws
    
    #Exchange draws
    if(gamma_index==1){
      Y_aux_list_prime = list("Pre" = rmcomp_T_given(668, lambda_prime[c(1,3)], nu_current, delta_current, omega_current, ratios_prime[c(1,3)], tol),
                              "During" = rmcomp_T_given(472, lambda_prime[c(2,3)], nu_current, delta_current, omega_current, ratios_prime[c(2,3)], tol))
    } else if(gamma_index ==2){
      Y_aux_list_prime= Y_aux_list;
      Y_aux_list_prime$During = rmcomp_T_given(472, lambda_prime[c(2,3)], nu_current, delta_current, omega_current, ratios_prime[c(2,3)], tol)
    }
    
    Y_prime = rbind(Y_aux_list_prime$Pre, Y_aux_list_prime$During)
    
    #PMF check
    parameters$gamma = gamma_prime
    check = check_pmf_reg(parameters, ratios_prime, data = list('Y' = Y_prime, 'X'=X))  
  }
  
  #Priors
  prior_current = dnorm(gamma_current[gamma_index], sd = prior_gamma$sd,log = TRUE)
  prior_prime = dnorm(gamma_prime[gamma_index], sd= prior_gamma$sd,log = TRUE)
  
  #Exchange Y for Yprime
  num = exchange_cpp(Y, X, gamma_prime, nu_current, delta_current, omega_current, ratios_prime)$q_lambda1 + #{q(Y|prime)}
    exchange_cpp(Y_prime, X, gamma_current, nu_current, delta_current, omega_current, ratios)$q_lambda1 +    #{q(Y_prime|current)}
    prior_prime
  
  den = exchange_cpp(Y, X, gamma_current, nu_current, delta_current, omega_current, ratios)$q_lambda1 +  #{q(Y|current)}
    exchange_cpp(Y_prime, X, gamma_prime, nu_current, delta_current, omega_current, ratios_prime)$q_lambda1 +  #{q(Y_prime|prime)}
    prior_current
  
  #Acceptance:
  alpha = num - den
  
  if(log(runif(1)) < alpha){
    gamma_out = gamma_prime
    ratios_out = ratios_prime
    Y_aux_out = Y_aux_list_prime
    draws_r_out = draws_prime_r
    accepted = 1;
  } else {
    gamma_out = gamma_current
    ratios_out = ratios
    Y_aux_out = Y_aux_list
    draws_r_out = draws_r
    accepted = 0;
  }
  
  output = list("gamma" = gamma_out[gamma_index+1], 
                "accepted" = accepted, 
                "ratios" = ratios_out,
                "Y_aux_list" = Y_aux_out,
                "draws_r" = draws_r_out)
  return(output)
}

exchange_nu_reg = function(dim, Y, X, gamma_current, nu_current, omega_current, delta_current, ratios, draws_r, sigma_nu, Y_aux_list, tol, prior_nu){
  
  lambda_current = compute_lambda(gamma_current)
  
  d = ncol(Y); n = nrow(Y)
  N_aux_r = nrow(draws_r); 
  nu_prime = nu_current
  parameters = list('gamma' = gamma_current,
                    'nu' = nu_current,
                    'delta' = delta_current,
                    'omega' = omega_current)
  
  ratios_prime = ratios; 
  draws_prime_r = draws_r;
  
  update_nu = ln_proposal(nu_current[dim], sigma_nu[dim])
  nu_prime[dim] = update_nu$proposed
  q_ratio_nu = update_nu$log_q_ratio
  
  if(dim == 2){
    
    #Estimate ratios[3] and save draws:
    r_third_element = ratio_IS_d_cpp(3, compute_lambda(gamma_current), c(nu_prime[1], nu_prime[1], nu_prime[2]), omega_current, N_aux_r)
    
    ratios_prime[3] = r_third_element$ratio; 
    draws_prime_r[,3]= r_third_element$draws
    
  }  
  if(dim == 1){
    
    #Estimate ratios[1:2] and save draws:
    r_first_element = ratio_IS_d_cpp(1, compute_lambda(gamma_current), c(nu_prime[1], nu_prime[1], nu_prime[2]), omega_current, N_aux_r)
    r_second_element = ratio_IS_d_cpp(2, compute_lambda(gamma_current), c(nu_prime[1], nu_prime[1], nu_prime[2]), omega_current, N_aux_r)
    
    ratios_prime[1] = r_first_element$ratio; 
    draws_prime_r[,1]= r_first_element$draws;
    
    ratios_prime[2] = r_second_element$ratio;
    draws_prime_r[,2] = r_second_element$draws;
    
  }
  
  #Exchange draws:
  Y_aux_list_prime= list("Pre" = rmcomp_T_given(668, lambda_current[c(1,3)], nu_prime, delta_current, omega_current, ratios_prime[c(1,3)], tol),
                         "During" = rmcomp_T_given(472, lambda_current[c(2,3)], nu_prime, delta_current, omega_current, ratios_prime[c(2,3)], tol))
  Y_prime = rbind(Y_aux_list_prime$Pre, Y_aux_list_prime$During)
  
  #Check if kernel is computable and yields proper pmf:
  parameters$nu = nu_prime
  check = check_pmf_reg(parameters, ratios_prime, data = list("Y" = Y_prime, "X"=X))  
  
  while(check == FALSE){
    update_nu = ln_proposal(nu_current[dim], sigma_nu[dim])
    nu_prime[dim] = update_nu$proposed
    q_ratio_nu = update_nu$log_q_ratio
    
    if(dim == 2){
      
      #Estimate ratios[3] and save draws:
      r_third_element = ratio_IS_d_cpp(3, compute_lambda(gamma_current), c(nu_prime[1], nu_prime[1], nu_prime[2]), omega_current, N_aux_r)
      
      ratios_prime[3] = r_third_element$ratio; 
      draws_prime_r[,3]= r_third_element$draws
      
    }  
    if(dim == 1){
      
      #Estimate ratios[1:2] and save draws:
      r_first_element = ratio_IS_d_cpp(1, compute_lambda(gamma_current), c(nu_prime[1], nu_prime[1], nu_prime[2]), omega_current, N_aux_r)
      r_second_element = ratio_IS_d_cpp(2, compute_lambda(gamma_current), c(nu_prime[1], nu_prime[1], nu_prime[2]), omega_current, N_aux_r)
      
      ratios_prime[1] = r_first_element$ratio; 
      draws_prime_r[,1]= r_first_element$draws;
      
      ratios_prime[2] = r_second_element$ratio;
      draws_prime_r[,2] = r_second_element$draws;
      
    }
    
    #Exchange draws:
    Y_aux_list_prime= list("Pre" = rmcomp_T_given(668, lambda_current[c(1,3)], nu_prime, delta_current, omega_current, ratios_prime[c(1,3)], tol),
                           "During" = rmcomp_T_given(472, lambda_current[c(2,3)], nu_prime, delta_current, omega_current, ratios_prime[c(2,3)], tol))
    Y_prime = rbind(Y_aux_list_prime$Pre, Y_aux_list_prime$During)
    
    #Check if kernel is computable and yields proper pmf:
    parameters$nu = nu_prime
    check = check_pmf_reg(parameters, ratios_prime, data = list("Y" = Y_prime, "X"=X))  
  }
  
  #Priors
  prior_current = log( dtruncnorm(nu_current[dim], mean =0, a =0, sd = prior_nu$sd) )
  prior_prime = log( dtruncnorm(nu_prime[dim], mean =0, a =0, sd = prior_nu$sd) )
  
  num = exchange_cpp(Y, X, gamma_current, nu_prime, delta_current, omega_current, ratios_prime)[[paste0("q_nu",dim)]] + #{q(Y|prime)}
    exchange_cpp(Y_prime, X, gamma_current, nu_current, delta_current, omega_current, ratios)[[paste0("q_nu",dim)]] +    #{q(Y_prime|current)}
    prior_prime
  
  den = exchange_cpp(Y, X, gamma_current, nu_current, delta_current, omega_current, ratios)[[paste0("q_nu",dim)]] +  #{q(Y|current)}
    exchange_cpp(Y_prime, X, gamma_current, nu_prime, delta_current, omega_current, ratios_prime)[[paste0("q_nu",dim)]] +  #{q(Y_prime|prime)}
    prior_current
  
  #Acceptance:
  alpha = num - den + q_ratio_nu
  
  if(log(runif(1)) < alpha){
    nu_out = nu_prime
    ratios_out = ratios_prime
    Y_aux_out = Y_aux_list_prime
    draws_r_out = draws_prime_r
    accepted = 1;
  } else {
    nu_out = nu_current
    ratios_out = ratios
    Y_aux_out = Y_aux_list
    draws_r_out = draws_r
    accepted = 0;
  }
  
  output = list("nu" = nu_out,"accepted" = accepted, 
                "ratios" = ratios_out, "Y_aux_list" = Y_aux_out,"draws_r" = draws_r_out)
  return(output)
}

exchange_omega_reg = function(Y, X, gamma_current, nu_current, omega_current, delta_current, ratios, draws_r, sigma_omega, Y_aux_list, tol, prior_omega){
  
  d= ncol(Y);  n = nrow(Y)
  lambda_current = compute_lambda(gamma_current)
  
  parameters = list('gamma' = gamma_current,
                    'nu' = nu_current,
                    'delta' = delta_current,
                    'omega' = omega_current)
  
  update = ln_proposal(omega_current, sigma_omega)
  omega_prime = update$proposed
  q_ratio = update$log_q_ratio
  
  #Ratios: with recycled draws
  r1 = mean(((exp(-omega_prime)*lambda_current[1])/lambda_current[1])^draws_r[,1])
  r2 = mean(((exp(-omega_prime)*lambda_current[2])/lambda_current[2])^draws_r[,2])
  r3 = mean(((exp(-omega_prime)*lambda_current[3])/lambda_current[3])^draws_r[,3])
  ratios_prime = c(r1,r2,r3)
  
  #Exchange draws:
  Y_aux_list_prime= list("Pre" = rmcomp_T_given(668, lambda_current[c(1,3)], nu_current, delta_current, omega_prime, ratios_prime[c(1,3)], tol),
                         "During" = rmcomp_T_given(472, lambda_current[c(2,3)], nu_current, delta_current, omega_prime, ratios_prime[c(2,3)], tol))
  Y_prime = rbind(Y_aux_list_prime$Pre, Y_aux_list_prime$During)
  
  #Check if kernel is computable and yields proper pmf
  parameters$omega = omega_prime
  check = check_pmf_reg(parameters, ratios_prime, data = list("Y" = Y_prime, "X" = X))
  
  while(check == FALSE){
    update = ln_proposal(omega_current, sigma_omega)
    omega_prime = update$proposed
    q_ratio = update$log_q_ratio
    
    r1 = mean(((exp(-omega_prime)*lambda_current[1])/lambda_current[1])^draws_r[,1])
    r2 = mean(((exp(-omega_prime)*lambda_current[2])/lambda_current[2])^draws_r[,2])
    r3 = mean(((exp(-omega_prime)*lambda_current[3])/lambda_current[3])^draws_r[,3])
    ratios_prime = c(r1,r2,r3)
    
    #Exchange draws:
    Y_aux_list_prime= list("Pre" = rmcomp_T_given(668, lambda_current[c(1,3)], nu_current, delta_current, omega_prime, ratios_prime[c(1,3)], tol),
                           "During" = rmcomp_T_given(472, lambda_current[c(2,3)], nu_current, delta_current, omega_prime, ratios_prime[c(2,3)], tol))
    Y_prime = rbind(Y_aux_list_prime$Pre, Y_aux_list_prime$During)
    
    #Check if kernel is computable and yields proper pmf
    parameters$omega = omega_prime
    check = check_pmf_reg(parameters, ratios_prime, data = list("Y" = Y_prime, "X" = X))
  }
  
  #Priors
  prior_current = log( dtruncnorm(omega_current, mean =0, a =0, sd = prior_omega$sd) )
  prior_prime = log( dtruncnorm(omega_prime, mean =0, a =0, sd = prior_omega$sd) )
  
  t1 = exchange_cpp(Y, X, gamma_current, nu_current, delta_current, omega_prime, ratios_prime)$kernel/exchange_cpp(Y, X, gamma_current, nu_current, delta_current, omega_current, ratios)$kernel    #{q(Y|prime)}/{q(Y_prime|current)}
  t2 = exchange_cpp(Y_prime, X, gamma_current, nu_current, delta_current, omega_current, ratios)$kernel/exchange_cpp(Y_prime, X, gamma_current, nu_current, delta_current, omega_prime, ratios_prime)$kernel    #{q(Y|prime)}/{q(Y_prime|current)}
  
  alpha = (t1*t2)*exp(q_ratio)*exp(prior_prime- prior_current)
  alpha = min(c(1, alpha))
  
  if(runif(1) < alpha){
    omega_out = omega_prime
    ratios_out = ratios_prime
    Y_aux_out = Y_aux_list_prime
    accepted = 1;
  } else {
    omega_out = omega_current
    ratios_out = ratios
    Y_aux_out = Y_aux_list
    accepted = 0;
  }
  
  return(list("omega" = omega_out, "accepted" = accepted, "ratios" = ratios_out, "Y_aux_list" = Y_aux_out))
  
}

exchange_delta_reg = function(Y, X, gamma_current, nu_current, omega_current, delta_current, ratios, draws_r, sigma_delta, Y_aux_list, tol){
  
  parameters = list('gamma' = gamma_current,
                    'nu' = nu_current,
                    'delta' = delta_current,
                    'omega' = omega_current)
  
  lambda_current = compute_lambda(gamma_current)
  
  lim1 =delta_limits_cpp(c(ratios[1], ratios[3]) )
  lim2= delta_limits_cpp(c(ratios[2], ratios[3]) )
  range_delta = c( min(c(lim1[1], lim2[1])),   max(c(lim1[2], lim2[2])))
  
  delta_prime = delta_current #create a copy
  delta_prime[1,2] = rtruncnorm(1, a = range_delta[1], b = range_delta[2], mean = delta_current[1,2], sd = sigma_delta[1,2])
  parameters$delta = delta_prime
  
  #Exchange draws
  Y_aux_list_prime= list("Pre" = rmcomp_T_given(668, lambda_current[c(1,3)], nu_current, delta_prime, omega_current, ratios[c(1,3)], tol),
                         "During" = rmcomp_T_given(472, lambda_current[c(2,3)], nu_current, delta_prime, omega_current, ratios[c(2,3)], tol))
  Y_prime = rbind(Y_aux_list_prime$Pre, Y_aux_list_prime$During)
  
  check = check_pmf_reg(parameters, ratios, data = list("Y" = Y_prime, "X" = X))
  
  while(check == FALSE){
    delta_prime = delta_current #create a copy
    delta_prime[1,2] = rtruncnorm(1, a = range_delta[1], b = range_delta[2], mean = delta_current[1,2], sd = sigma_delta[1,2])
    parameters$delta = delta_prime
    
    #Exchange draws
    Y_aux_list_prime= list("Pre" = rmcomp_T_given(668, lambda_current[c(1,3)], nu_current, delta_prime, omega_current, ratios[c(1,3)], tol),
                           "During" = rmcomp_T_given(472, lambda_current[c(2,3)], nu_current, delta_prime, omega_current, ratios[c(2,3)], tol))
    Y_prime = rbind(Y_aux_list_prime$Pre, Y_aux_list_prime$During)
    
    check = check_pmf_reg(parameters, ratios, data = list("Y" = Y_prime, "X" = X))
  }
  
  t1 = exchange_cpp(Y, X, gamma_current, nu_current, delta_prime, omega_current, ratios)$kernel/exchange_cpp(Y, X, gamma_current, nu_current, delta_current, omega_current, ratios)$kernel    #{q(Y|prime)}/{q(Y_prime|current)}
  t2 = exchange_cpp(Y_prime, X, gamma_current, nu_current, delta_current, omega_current, ratios)$kernel/exchange_cpp(Y_prime, X, gamma_current, nu_current, delta_prime, omega_current, ratios)$kernel    #{q(Y|prime)}/{q(Y_prime|current)}
  
  #Acceptance:
  alpha = min(c(1, (t1*t2)))
  
  if(runif(1) < alpha){
    delta_out = delta_prime
    Y_aux_out = Y_aux_list_prime
    accepted = 1;
  } else {
    delta_out = delta_current
    Y_aux_out = Y_aux_list
    accepted = 0;
  }
  
  return(list("delta" = delta_out, "accepted" = accepted, "Y_aux_list" = Y_aux_out))
  
}

Exchange_regression_wrapper = function(burn_in, n_iter, initialise, Y, X, N_aux_r, priors_list, tol){
  
  n = nrow(Y); d = ncol(Y);
  
  # Initial values =----------------------------------------------------------------------------
  if(initialise == 'random'){
    
    check = FALSE
    while(check == FALSE){
      gamma_init = runif(n=3, min =-1, max =1)
      nu_init = runif(n=2, min =0.5, max =2)
      omega_init = runif(n=1, min =, max = 3)
      delta_init = matrix(0, ncol = d, nrow = d)
      delta_current[upper.tri(delta_current)] = rnorm(1, mean = 0, sd = 0.5)
      
      ratios = ratios_IS_cpp(compute_lambda(gamma_init), c(nu_init[1], nu_init[1], nu_init[2]), omega_init, N_aux_r)$ratios
      
      inits = list('gamma' = gamma_init,'nu' = nu_init,'delta' = delta_init, 'omega' = omega_init)
      check = check_pmf_reg(parameters, ratios, data = list("Y" = Y, "X"=X))
    }
    
    gamma_current = inits$gamma
    nu_current = inits$nu
    delta_current = inits$delta
    omega_current = inits$omega
    
  } else {
    inits = initial_reg_premier(Y,X)
    gamma_current = inits$gamma
    nu_current = inits$nu
    delta_current = inits$delta
    omega_current = inits$omega
  }
  
  cat('Initial values: ------------------',"\n",
      "gamma:", round(gamma_current,3), '\n',
      'nu:', round(nu_current,3), '\n',
      'delta;' , round(delta_current[1,2],3), '\n', 
      'omega:', round(omega_current,3), '\n',
      '----------------------------------', '\n')
  
  # Initialise chains/Rates ---------------------------------------------------------------------
  target_ac = 0.44
  sigma_omega = 0.2;  #SD
  sigma_delta = matrix(0.1, ncol = d, nrow =d) 
  sigma_delta[lower.tri(sigma_delta,TRUE)]<- NA #SD
  
  sigma_gamma0 = 0.1 #For common intercept
  sigma_gamma = rep(0.1,2) #Gamma1, gamma2
  sigma_nu = c(0.1,0.1)
  
  accept_omega = 0;
  accept_delta =0
  accept_gamma = rep(0, 3);
  accept_nu = rep(0, 2);
  
  gamma_chain = matrix(NA, nrow = burn_in+n_iter, ncol = 3)
  nu_chain = matrix(NA, nrow = burn_in+n_iter, ncol = d)
  omega_chain = rep(NA, burn_in+n_iter); 
  delta_chain = rep(NA, burn_in+n_iter); 
  
  # MCMC -----------------------------------------------------------------------------------------------
  lambda_current = compute_lambda(gamma_current)
  
  # Compute ratios and save draws
  ratio_estimates = ratios_IS_cpp(lambda_current, c(nu_current[1], nu_current[1], nu_current[2]), omega_current, N_aux_r) 
  ratios = ratio_estimates$ratios
  draws_r = ratio_estimates$draws
  
  iter =1
  start_time = Sys.time()
  while(iter <= (burn_in+n_iter)){
    
    Y_aux_list = list("Pre" = rmcomp_T_given(668, lambda_current[c(1,3)], nu_current, delta_current, omega_current, ratios[c(1,3)],tol),
                      "During" = rmcomp_T_given(472, lambda_current[c(2,3)], nu_current, delta_current, omega_current, ratios[c(2,3)], tol))
    
    
    # ----------------- 1. Common Intercept Update  --------------------
    intercepts = exchange_gamma0(Y, X, gamma_current, nu_current, delta_current, omega_current, sigma_gamma0, ratios, draws_r, Y_aux_list, tol, priors_list$gamma0)
    
    gamma_current[1] = intercepts$gamma0
    ratios = intercepts$ratios
    Y_aux_list = intercepts$Y_aux_list; 
    draws_r = intercepts$draws_r
    
    accept_gamma[1] = accept_gamma[1] + intercepts$accepted
    sigma_gamma0 = update_propsd_double(sigma_gamma0,  accept_gamma[1]/iter,  iter, target_ac)
    
    # ----------------- 2. Gamma1 ------------------------------------
    gamma_1 = exchange_gamma(1, Y, X, gamma_current, nu_current, omega_current, delta_current, ratios, draws_r, sigma_gamma, Y_aux_list, tol, priors_list$gamma)
    
    gamma_current[2]= gamma_1$gamma
    ratios = gamma_1$ratios
    Y_aux_list = gamma_1$Y_aux_list; 
    draws_r = gamma_1$draws_r
    
    accept_gamma[2] =accept_gamma[2] + gamma_1$accepted
    
    sigma_gamma[1]=update_propsd_double(sigma_gamma[1], accept_gamma[2]/iter,  iter, target_ac)
    
    # ----------------- 3. Gamma2 ------------------------------------
    gamma_2 = exchange_gamma(2, Y, X, gamma_current, nu_current, omega_current, delta_current, ratios, draws_r, sigma_gamma, Y_aux_list, tol, priors_list$gamma)
    
    gamma_current[3]= gamma_2$gamma
    ratios = gamma_2$ratios
    Y_aux_list = gamma_2$Y_aux_list; 
    draws_r = gamma_2$draws_r
    
    accept_gamma[3] =accept_gamma[3] + gamma_2$accepted
    
    sigma_gamma[2]=update_propsd_double(sigma_gamma[2], accept_gamma[3]/iter,  iter, target_ac)
    
    # ----------------- 4. Nu1, 4. Nu2 update  ----------------------------
    for(dim in 1:2){
      nu_update = exchange_nu_reg(dim, Y, X, gamma_current, nu_current, omega_current, delta_current, ratios, draws_r, sigma_nu, Y_aux_list, tol, priors_list$nu)
      nu_current = nu_update$nu
      ratios = nu_update$ratios
      Y_aux_list = nu_update$Y_aux_list; 
      draws_r = nu_update$draws_r
      
      accept_nu[dim] = accept_nu[dim] + nu_update$accepted
      
      sigma_nu[dim] =update_propsd_double( sigma_nu[dim], accept_nu[dim]/iter,  iter, target_ac)
      
    }
    
    # ----------------- 5. Update omega ----------------------------
    omega_update =exchange_omega_reg(Y, X, gamma_current, nu_current, omega_current, delta_current, ratios, draws_r, sigma_omega, Y_aux_list, tol, priors_list$omega)
    
    omega_current = omega_update$omega
    ratios = omega_update$ratios; 
    Y_aux_list = omega_update$Y_aux_list; 
    
    #Update sd
    accept_omega = accept_omega + omega_update$accepted
    sigma_omega <- update_propsd_double(sigma_omega,  accept_omega/iter,  iter, target_ac)
    
    # ----------------- 6. Update delta ---------------------------
    update_delta = exchange_delta_reg(Y, X, gamma_current, nu_current, omega_current, delta_current, ratios, draws_r, sigma_delta, Y_aux_list,tol, priors_list$delta)
    
    delta_current = update_delta$delta
    accept_delta = accept_delta + update_delta$accepted
    Y_aux_list = update_delta$Y_aux_list; 
    
    sigma_delta[1,2] <- update_propsd_double(sigma_delta[1,2], accept_delta/iter, iter, target_ac)
    
    # Store -------------------------------------------------------
    gamma_chain[iter,] = gamma_current
    nu_chain[iter,] = nu_current
    delta_chain[iter] = delta_current[1,2]
    omega_chain[iter] = omega_current
  
    iter = iter +1
    
    if(iter %% 100==0){
      
      current_time = Sys.time()
      
      cat("ITERATION:", iter, "----------------------------", "\n")
      
      cat("REGRESSION  --------------------------------:","\n",
          "Accepted gamma", round(accept_gamma/iter,3),"\n",
          "Current lambda:", round(compute_lambda(gamma_current),3), "\n",
          
          "Nu report ---------------------------------------:","\n",
          "Accepted nu", round(accept_nu/iter,3),"\n",
          "Current nu", round(nu_current,3),"\n",
          
          "KERNEL --------------------------------------------:","\n",
          "Omega accepted", round(accept_omega/iter,3),"\n",
          "Current omega:", round(omega_current,3), "\n",
          "Sigma omega:", round(sigma_omega,3), "\n",
          
          "Delta accepted",  round(accept_delta/iter,3),"\n",
          "Current delta:", round(delta_current[1,2],3), "\n",
          
          "TIMING --------------------------------------------", "\n",
          current_time-start_time, "\n",
          "-------------------------------------------------------","\n",
          "-------------------------------------------------------")
      
    }
    
  }
  end_time = Sys.time();
  
  # Return arguments --------------------------------------------------------------------
  total_iter = burn_in + n_iter
  output = list('ac_rates' = list('gamma' = accept_gamma/total_iter,
                                  'nu' = accept_nu/total_iter,
                                  'omega' = accept_omega/total_iter,
                                  'delta' = accept_delta/total_iter),
                "time" = end_time-start_time,
                "Nr" = N_aux_r, "initial_values" = inits,
                'proposal_parameters' = list("gamma0" = sigma_gamma0,
                                             "gamma1_2" = sigma_gamma,
                                             "omega" = sigma_omega,
                                             "delta" = sigma_delta),
                "gamma_chain" = gamma_chain[(burn_in+1):total_iter, ],
                "nu_chain" = nu_chain[(burn_in+1):total_iter, ],
                "omega_chain" = omega_chain[(burn_in+1):total_iter],
                "delta_chain" = delta_chain[(burn_in+1):total_iter])
  return(output)
}


Exchange_regression = function(burn_in, n_iter, initialise, Y, X, N_aux_r, priors_list, tol =0.001, chains =1, ncores = 5){
  
  if(chains == 1){
    start = Sys.time()
    run = Exchange_regression_wrapper(burn_in, n_iter, initialise, Y, X, N_aux_r, priors_list, tol)
    end = Sys.time()
    timing = end-start
    chains = list(run)
    
  } else if(chains >1){
    
    start = Sys.time()
    cl = parallel::makeCluster(ncores, setup_strategy = "sequential")
    parallel::clusterEvalQ(cl, library("multcp"))
    parallel::clusterEvalQ(cl, library("truncnorm"))
    parallel::clusterEvalQ(cl, library("COMPoissonReg"))
    parallel::clusterEvalQ(cl, library("mvtnorm"))
    
    chains = parallel::parSapply(cl, 1:chains, 
                                 function(times, burn_in, n_iter, initialise, Y, X, N_aux_r, priors_list, tol){
                                   Exchange_regression_wrapper(burn_in, n_iter, initialise, Y, X, N_aux_r, priors_list, tol)
                                 },burn_in, n_iter, initialise, Y, X, N_aux_r, priors_list, tol,simplify = F)
    
    parallel::stopCluster(cl) 
    end = Sys.time()
    timing = end-start
  }
  
  mcmc = list("mcmc" = chains,"time" = timing)
  return(mcmc)
  
}

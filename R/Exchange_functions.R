inits_marginal = function(dim, Y){
  glm = COMPoissonReg::glm.cmp(formula.lambda = Y[,dim] ~ 1,
                formula.nu = Y[,dim] ~1)
  inits = exp(coef(glm)) 
  names(inits)=c("lambda", "nu")
  return(inits)
}

check_move = function(current_list, prime_list, Y){
  
  check_ratios = sum(is.na(prime_list$ratios)) ==0
  # Y'|theta 
  q1 = un_loglik_given(1,prime_list$Y, current_list$params$lambda, current_list$params$nu, current_list$params$delta, current_list$params$omega, current_list$ratios)
  # Y| theta'
  q2 = un_loglik_given(1, Y, prime_list$params$lambda, prime_list$params$nu, prime_list$params$delta, prime_list$params$omega, prime_list$ratios)
  # Y'|theta'
  q3 = un_loglik_given(1, prime_list$Y,  prime_list$params$lambda, prime_list$params$nu, prime_list$params$delta, prime_list$params$omega, prime_list$ratios)
  
  check_nan = !any(is.na( c(q1$log_kernel, q2$log_kernel, q3$log_kernel)) )
  check_positive= ifelse(check_nan,  prod(c(q1$kernel >0, q2$kernel >0, q3$kernel>0))==1, FALSE)
  
  check = check_ratios&check_nan&check_positive
  return(check)
}

ratio_IS_recycle_d = function(dim, lambda, omega, draws_r){
  
  q_log_num = draws_r[,dim]*log( exp(-omega)*lambda[dim] )
  q_log_den = draws_r[,dim]*log(lambda[dim])
  
  is = mean(exp(q_log_num-q_log_den))
  return(is)
}

ratios_IS_recycle = function(lambda, omega, draws_r){
  output = sapply(1:length(lambda), ratio_IS_recycle_d, lambda, omega, draws_r)
  return(output)
}

Exchange_wrapper = function(n_iter, burn_in, Y, N_aux_r, initialise, priors_list, tol){
  
  n = nrow(Y); d = ncol(Y);
  
  # Initial values ----------------------------------------------------------------------------
  delta_current = matrix(0, ncol =d, nrow = d)
  combs = t(combn(seq(1,d,1),2)) 
  
  if(initialise == 'random'){
    lambda_current = sample(seq(0.5, 2, 0.1),d)
    nu_current = sample(seq(0.5, 2, 0.1),d)
    
    while( any(lambda_current^(1/nu_current)>3) ){
      lambda_current = sample(seq(0.5, 2, 0.1),d)
      nu_current = sample(seq(0.5, 2, 0.1),d)
    }
    delta_current = matrix(0, ncol = d, nrow = d)
    delta_current[upper.tri(delta_current)] = rnorm(choose(d,2), mean = 0, sd = 0.5)
    omega_current = sample(seq(0.1, 2, 0.5),1)
    
    inits = list("lambda" = lambda_current,"nu" = nu_current,
                 "delta" = delta_current,"omega" = omega_current)
  } else {
    inits= sapply(1:d, inits_marginal, Y) 
    lambda_current = inits[1,]
    nu_current = inits[2,]
    omega_current = 1
    for(i in 1:nrow(combs)){
      delta_current[combs[i,1], combs[i,2]] = sign(cor(Y[,combs[i,1]], Y[,combs[i,2]]))
    }
    
  }
  
  cat('Initial values: ------------------',"\n",
      "lambda:", round(lambda_current,3), '\n',
      'nu:', round(nu_current,3), '\n',
      'delta:' , round(delta_current,3), '\n', 
      'omega:', round(omega_current,3), '\n',
      '----------------------------------', '\n')
  
  # Initialise chains/Rates ---------------------------------------------------------------------
  target_ac = 0.44
  sigma_omega = 0.05;  
  sigma_delta = matrix(0.1, ncol = d, nrow =d) 
  sigma_delta[lower.tri(sigma_delta,TRUE)]<- NA
  sigma_lambda = rep(0.1,d)
  sigma_nu = rep(0.1,d)
  
  accept_omega = 0;
  accept_delta = matrix(0, ncol = d, nrow =d)
  accept_delta[lower.tri(accept_delta,TRUE)]<- NA
  accept_lambda = rep(0, d);
  accept_nu = rep(0, d);
  
  lambda_chain = matrix(NA, nrow = burn_in+n_iter, ncol = d)
  nu_chain = matrix(NA, nrow = burn_in+n_iter, ncol = d)
  omega_chain = rep(NA, burn_in+n_iter); 
  delta_chain = rep(list(matrix(NA, nrow = d, ncol = d)), burn_in+n_iter)
  
  # MCMC -----------------------------------------------------------------------------------------------
  
  # Compute ratios and save draws
  ratio_estimates = ratios_IS_cpp(lambda_current, nu_current, omega_current, N_aux_r) 
  ratios = ratio_estimates$ratios
  draws_r = ratio_estimates$draws
  
  iter =1
  start_time = Sys.time()
  while(iter <= (burn_in+n_iter)){
    
    # ----------------- 1. Lambda Update  --------------------
    for(dim in 1:d){
      update_l_d = exchange_lambda(dim, Y, lambda_current, nu_current, omega_current, delta_current, ratios, draws_r, sigma_lambda, priors_list$lambda,tol)
      
      lambda_current[dim] =update_l_d$lambda
      ratios = update_l_d$ratios
      draws_r = update_l_d$draws_r
      
      accept_lambda[dim] =   accept_lambda[dim] + update_l_d$accepted
      sigma_lambda[dim] = update_propsd_double(sigma_lambda[dim],  accept_lambda[dim]/iter,  iter, target_ac)

    }
    
    # ----------------- 2. Nu Update ------------------------------------
    for(dim in 1:d){
      update_nu_d = exchange_nu(dim, Y, lambda_current, nu_current, omega_current, delta_current, ratios, draws_r, sigma_nu, priors_list$nu, tol)
      
      nu_current[dim] =update_nu_d$nu
      ratios = update_nu_d$ratios
      draws_r = update_nu_d$draws_r
      
      accept_nu[dim] =   accept_nu[dim] + update_nu_d$accepted
      sigma_nu[dim] = update_propsd_double(sigma_nu[dim],  accept_nu[dim]/iter,  iter, target_ac)

    }
    
    # ----------------- 3. Omega Update ------------------------------------
    omega_update = exchange_omega(Y, lambda_current, nu_current, omega_current, delta_current, ratios, draws_r, sigma_omega, priors_list$omega,tol)
    omega_current = omega_update$omega
    ratios = omega_update$ratios; 
    
    #Update sd
    accept_omega = accept_omega + omega_update$accepted
    sigma_omega <- update_propsd_double(sigma_omega,  accept_omega/iter,  iter, target_ac)
    
    # ----------------- 4. Update delta ---------------------------
    for(index in 1:nrow(combs)){
      update_d = exchange_delta(index, combs, Y, lambda_current, nu_current, omega_current, delta_current, ratios, draws_r, sigma_delta, priors_list$delta,tol)
      
      delta_current = update_d$delta
      accept_delta[combs[index,1],combs[index,2]]<- accept_delta[combs[index,1],combs[index,2]]+ update_d$accepted
      
      #Update sd
      sigma_delta[combs[index,1],combs[index,2]] <- update_propsd_double(sigma_delta[combs[index,1],combs[index,2]], 
                                                                         accept_delta[combs[index,1],combs[index,2]]/iter,
                                                                         iter, target_ac)
    }
    
    # Store -------------------------------------------------------
    lambda_chain[iter,] = lambda_current
    nu_chain[iter,] = nu_current
    delta_chain[[iter]] = delta_current
    omega_chain[iter] = omega_current
    
    iter = iter +1

    if(iter %% 100==0){
      
      cat("ITERATION:", iter, "----------------------------", "\n")
      
      cat("Lambda report  --------------------------------:","\n",
          "Accepted lambda", round(accept_lambda/iter,3),"\n",
          "Current lambda:", round(lambda_current,3), "\n",
          
          "Nu report ---------------------------------------:","\n",
          "Accepted nu", round(accept_nu/iter,3),"\n",
          "Current nu", round(nu_current,3),"\n",
          
          "KERNEL --------------------------------------------:","\n",
          "Omega accepted", round(accept_omega/iter,3),"\n",
          "Current omega:", round(omega_current,3), "\n",
          "Sigma omega:", round(sigma_omega,3), "\n",
          
          "Delta accepted",  round(accept_delta/iter,3),"\n",
          "Current delta:", round(delta_current[1,2],3), "\n",
          
          "-------------------------------------------------------","\n",
          "-------------------------------------------------------")
    }
    
  }
  end_time = Sys.time();
  
  # Return arguments --------------------------------------------------------------------
  total_iter = burn_in + n_iter
  output = list('ac_rates' = list('lambda' = accept_lambda/total_iter,
                                  'nu' = accept_nu/total_iter,
                                  'omega' = accept_omega/total_iter,
                                  'delta' = accept_delta/total_iter),
                "time" = end_time-start_time,
                "Nr" = N_aux_r, "initial_values" = inits,
                'proposal_parameters' = list("lambda" = sigma_lambda,
                                             'nu' = sigma_nu,
                                             "omega" = sigma_omega,
                                             "delta" = sigma_delta),
                "lambda_chain" = lambda_chain[(burn_in+1):total_iter, ],
                "nu_chain" = nu_chain[(burn_in+1):total_iter, ],
                "omega_chain" = omega_chain[(burn_in+1):total_iter],
                "delta_chain" = delta_chain[(burn_in+1):total_iter])
  return(output)
}

exchange_lambda = function(dim, Y, lambda_current, nu_current, omega_current, delta_current, ratios, draws_r, sigma_lambda, prior_lambda, tol){

  parameters = list('lambda' = lambda_current, 'nu' = nu_current, 'omega'= omega_current, 'delta' = delta_current)
  d = ncol(Y); n = nrow(Y);
  N_aux_r = nrow(draws_r); 
  
  lambda_prime = lambda_current
  update_ = ln_proposal(lambda_current[dim], sigma_lambda[dim])
  lambda_prime[dim] = update_$proposed
  q_ratio_lambda = update_$log_q_ratio
  
  #Estimate ratios[dim] and save draws:
  ratios_prime = ratios; draws_prime_r = draws_r;
  run_ratio = ratio_IS_d_cpp(dim, lambda_prime, nu_current, omega_current, N_aux_r)
  ratios_prime[dim] = run_ratio$ratio
  draws_prime_r[,dim]= run_ratio$draws
  
  #Exchange draws
  Y_prime = rmcomp_T_given(n, lambda_prime, nu_current, delta_current, omega_current, ratios_prime, tol)
  
  #PMF check
  current_list=list(); current_list$params = parameters; current_list$ratios = ratios; prime_list = current_list;
  prime_list$ratios = ratios_prime
  prime_list$params$lambda = lambda_prime
  prime_list$Y = Y_prime
  
  check = check_move(current_list, prime_list, Y)
  while(check == FALSE){
    update_ = ln_proposal(lambda_current[dim], sigma_lambda[dim])
    lambda_prime[dim] = update_$proposed
    q_ratio_lambda = update_$log_q_ratio
    
    #Estimate ratios[dim] and save draws:
    ratios_prime = ratios; draws_prime_r = draws_r;
    run_ratio = ratio_IS_d_cpp(dim, lambda_prime, nu_current, omega_current, N_aux_r)
    ratios_prime[dim] = run_ratio$ratio
    draws_prime_r[,dim]= run_ratio$draws
    
    #Exchange draws
    Y_prime = rmcomp_T_given(n, lambda_prime, nu_current, delta_current, omega_current, ratios_prime, tol)
    
    #PMF check
    prime_list$ratios = ratios_prime
    prime_list$params$lambda = lambda_prime
    prime_list$Y = Y_prime
    
    check = check_move(current_list, prime_list, Y)
  }
  
  #Priors
  prior_current = dgamma(lambda_current[dim], shape = prior_lambda$shape, rate = prior_lambda$rate, log=TRUE)
  prior_prime = dgamma(lambda_prime[dim], shape = prior_lambda$shape, rate = prior_lambda$rate, log=TRUE)
  
  #Exchange Y for Yprime
  num =  un_loglik_given(dim, Y, lambda_prime, nu_current, delta_current, omega_current, ratios_prime)$lambda + #{q(Y|prime)}
    un_loglik_given(dim, Y_prime, lambda_current, nu_current, delta_current, omega_current, ratios)$lambda + #{q(Y_prime|current)}+
    prior_prime
  
  den =  un_loglik_given(dim, Y, lambda_current, nu_current, delta_current, omega_current, ratios)$lambda + #{q(Y|current)}
    un_loglik_given(dim, Y_prime, lambda_prime, nu_current, delta_current, omega_current, ratios_prime)$lambda + #{q(Y_prime|prime)}
    prior_current
  
  #Acceptance:
  alpha = num - den + q_ratio_lambda;
  
  if(log(runif(1)) < alpha){
    lambda_out = lambda_prime
    ratios_out = ratios_prime
    draws_r_out = draws_prime_r
    accepted = 1;
  } else {
    lambda_out = lambda_current
    ratios_out = ratios
    draws_r_out = draws_r
    accepted = 0;
  }
  
  output = list("lambda" = lambda_out[dim], 
                "accepted" = accepted, 
                "ratios" = ratios_out,
                "draws_r" = draws_r_out)
  return(output)
}

exchange_lambda2 = function(dim, Y, lambda_current, nu_current, omega_current, delta_current, ratios, draws_r, sigma_lambda, prior_lambda, tol){
  
  parameters = list('lambda' = lambda_current, 'nu' = nu_current, 'omega'= omega_current, 'delta' = delta_current)
  d = ncol(Y); n = nrow(Y);
  N_aux_r = nrow(draws_r); 
  
  lambda_prime = lambda_current
  update_ = ln_proposal(lambda_current[dim], sigma_lambda[dim])
  lambda_prime[dim] = update_$proposed
  q_ratio_lambda = update_$log_q_ratio
  
  #Estimate ratios[dim] and save draws:
  ratios_prime = ratios; draws_prime_r = draws_r;
  run_ratio = ratio_IS_d_cpp(dim, lambda_prime, nu_current, omega_current, N_aux_r)
  ratios_prime[dim] = run_ratio$ratio
  draws_prime_r[,dim]= run_ratio$draws
  
  #Exchange draws
  Y_prime = rmcomp_T_given(n, lambda_prime, nu_current, delta_current, omega_current, ratios_prime, tol)
  
  #PMF check
  current_list=list(); current_list$params = parameters; current_list$ratios = ratios; prime_list = current_list;
  prime_list$ratios = ratios_prime
  prime_list$params$lambda = lambda_prime
  prime_list$Y = Y_prime
  
  check = check_move(current_list, prime_list, Y)
  while(check == FALSE){
    update_ = ln_proposal(lambda_current[dim], sigma_lambda[dim])
    lambda_prime[dim] = update_$proposed
    q_ratio_lambda = update_$log_q_ratio
    
    #Estimate ratios[dim] and save draws:
    ratios_prime = ratios; draws_prime_r = draws_r;
    run_ratio = ratio_IS_d_cpp(dim, lambda_prime, nu_current, omega_current, N_aux_r)
    ratios_prime[dim] = run_ratio$ratio
    draws_prime_r[,dim]= run_ratio$draws
    
    #Exchange draws
    Y_prime = rmcomp_T_given(n, lambda_prime, nu_current, delta_current, omega_current, ratios_prime, tol)
    
    #PMF check
    prime_list$ratios = ratios_prime
    prime_list$params$lambda = lambda_prime
    prime_list$Y = Y_prime
    
    check = check_move(current_list, prime_list, Y)
  }
  
  #Priors
  prior_current = log(dtruncnorm(lambda_current[dim], a = 0, mean = 0, sd = prior_lambda$sd))
  prior_prime = log(dtruncnorm(lambda_prime[dim], a = 0, mean = 0, sd = prior_lambda$sd))
  
  #Exchange Y for Yprime
  num =  un_loglik_given(dim, Y, lambda_prime, nu_current, delta_current, omega_current, ratios_prime)$lambda + #{q(Y|prime)}
    un_loglik_given(dim, Y_prime, lambda_current, nu_current, delta_current, omega_current, ratios)$lambda + #{q(Y_prime|current)}+
    prior_prime
  
  den =  un_loglik_given(dim, Y, lambda_current, nu_current, delta_current, omega_current, ratios)$lambda + #{q(Y|current)}
    un_loglik_given(dim, Y_prime, lambda_prime, nu_current, delta_current, omega_current, ratios_prime)$lambda + #{q(Y_prime|prime)}
    prior_current
  
  #Acceptance:
  alpha = num - den + q_ratio_lambda;
  
  if(log(runif(1)) < alpha){
    lambda_out = lambda_prime
    ratios_out = ratios_prime
    draws_r_out = draws_prime_r
    accepted = 1;
  } else {
    lambda_out = lambda_current
    ratios_out = ratios
    draws_r_out = draws_r
    accepted = 0;
  }
  
  output = list("lambda" = lambda_out[dim], 
                "accepted" = accepted, 
                "ratios" = ratios_out,
                "draws_r" = draws_r_out)
  return(output)
}


exchange_nu = function(dim, Y, lambda_current, nu_current, omega_current, delta_current, ratios, draws_r, sigma_nu, prior_nu, tol){
  
  parameters = list('lambda' = lambda_current, 'nu' = nu_current, 'omega'= omega_current, 'delta' = delta_current)
  d = ncol(Y); n = nrow(Y);
  N_aux_r = nrow(draws_r); 
  
  nu_prime = nu_current
  update_ = ln_proposal(nu_current[dim], sigma_nu[dim])
  nu_prime[dim] = update_$proposed
  q_ratio_nu = update_$log_q_ratio
  
  #Estimate ratios[dim] and save draws:
  ratios_prime = ratios; draws_prime_r = draws_r;
  run_ratio = ratio_IS_d_cpp(dim, lambda_current, nu_prime, omega_current, N_aux_r)
  ratios_prime[dim] = run_ratio$ratio
  draws_prime_r[,dim]= run_ratio$draws
  
  #Exchange draws
  Y_prime = rmcomp_T_given(n, lambda_current, nu_prime, delta_current, omega_current, ratios_prime, tol)
  
  #PMF check
  current_list=list(); current_list$params = parameters; current_list$ratios = ratios; prime_list = current_list;
  prime_list$ratios = ratios_prime
  prime_list$params$nu = nu_prime
  prime_list$Y = Y_prime
  
  check = check_move(current_list, prime_list, Y)
  while(check == FALSE){
    update_ = ln_proposal(nu_current[dim], sigma_nu[dim])
    nu_prime[dim] = update_$proposed
    q_ratio_nu = update_$log_q_ratio
    
    #Estimate ratios[dim] and save draws:
    ratios_prime = ratios; draws_prime_r = draws_r;
    run_ratio = ratio_IS_d_cpp(dim, lambda_current, nu_prime, omega_current, N_aux_r)
    ratios_prime[dim] = run_ratio$ratio
    draws_prime_r[,dim]= run_ratio$draws
    
    #Exchange draws
    Y_prime = rmcomp_T_given(n, lambda_current, nu_prime, delta_current, omega_current, ratios_prime, tol)
    
    #PMF check
    prime_list$ratios = ratios_prime
    prime_list$params$nu = nu_prime
    prime_list$Y = Y_prime
    
    check = check_move(current_list, prime_list, Y)
  }
  
  #Priors
  prior_current = dgamma(nu_current[dim], log = TRUE, shape = prior_nu$shape, rate =prior_nu$rate)
  prior_prime = dgamma(nu_prime[dim], log = TRUE, shape = prior_nu$shape, rate =prior_nu$rate)
  
  #Exchange Y for Yprime
  num =  un_loglik_given(dim, Y, lambda_current, nu_prime, delta_current, omega_current, ratios_prime)$nu + #{q(Y|prime)}
    un_loglik_given(dim, Y_prime, lambda_current, nu_current, delta_current, omega_current, ratios)$nu + #{q(Y_prime|current)}+
    prior_prime
  
  den =  un_loglik_given(dim, Y, lambda_current, nu_current, delta_current, omega_current, ratios)$nu + #{q(Y|current)}
    un_loglik_given(dim, Y_prime, lambda_current, nu_prime, delta_current, omega_current, ratios_prime)$nu + #{q(Y_prime|prime)}
    prior_current
  
  #Acceptance:
  alpha = num - den + q_ratio_nu;
  
  if(log(runif(1)) < alpha){
    nu_out = nu_prime
    ratios_out = ratios_prime
    draws_r_out = draws_prime_r
    accepted = 1;
  } else {
    nu_out = nu_current
    ratios_out = ratios
    draws_r_out = draws_r
    accepted = 0;
  }
  
  output = list("nu" = nu_out[dim], 
                "accepted" = accepted, 
                "ratios" = ratios_out,
                "draws_r" = draws_r_out)
  return(output)
}

exchange_nu2 = function(dim, Y, lambda_current, nu_current, omega_current, delta_current, ratios, draws_r, sigma_nu, prior_nu, tol){
  
  parameters = list('lambda' = lambda_current, 'nu' = nu_current, 'omega'= omega_current, 'delta' = delta_current)
  d = ncol(Y); n = nrow(Y);
  N_aux_r = nrow(draws_r); 
  
  nu_prime = nu_current
  update_ = ln_proposal(nu_current[dim], sigma_nu[dim])
  nu_prime[dim] = update_$proposed
  q_ratio_nu = update_$log_q_ratio
  
  #Estimate ratios[dim] and save draws:
  ratios_prime = ratios; draws_prime_r = draws_r;
  run_ratio = ratio_IS_d_cpp(dim, lambda_current, nu_prime, omega_current, N_aux_r)
  ratios_prime[dim] = run_ratio$ratio
  draws_prime_r[,dim]= run_ratio$draws
  
  #Exchange draws
  Y_prime = rmcomp_T_given(n, lambda_current, nu_prime, delta_current, omega_current, ratios_prime, tol)
  
  #PMF check
  current_list=list(); current_list$params = parameters; current_list$ratios = ratios; prime_list = current_list;
  prime_list$ratios = ratios_prime
  prime_list$params$nu = nu_prime
  prime_list$Y = Y_prime
  
  check = check_move(current_list, prime_list, Y)
  while(check == FALSE){
    update_ = ln_proposal(nu_current[dim], sigma_nu[dim])
    nu_prime[dim] = update_$proposed
    q_ratio_nu = update_$log_q_ratio
    
    #Estimate ratios[dim] and save draws:
    ratios_prime = ratios; draws_prime_r = draws_r;
    run_ratio = ratio_IS_d_cpp(dim, lambda_current, nu_prime, omega_current, N_aux_r)
    ratios_prime[dim] = run_ratio$ratio
    draws_prime_r[,dim]= run_ratio$draws
    
    #Exchange draws
    Y_prime = rmcomp_T_given(n, lambda_current, nu_prime, delta_current, omega_current, ratios_prime, tol)
    
    #PMF check
    prime_list$ratios = ratios_prime
    prime_list$params$nu = nu_prime
    prime_list$Y = Y_prime
    
    check = check_move(current_list, prime_list, Y)
  }
  
  #Priors
  prior_current = log(dtruncnorm(nu_current[dim], a = 0, mean = 0, sd = prior_nu$sd))
  prior_prime = log(dtruncnorm(nu_prime[dim], a = 0, mean = 0, sd = prior_nu$sd))
  
  #Exchange Y for Yprime
  num =  un_loglik_given(dim, Y, lambda_current, nu_prime, delta_current, omega_current, ratios_prime)$nu + #{q(Y|prime)}
    un_loglik_given(dim, Y_prime, lambda_current, nu_current, delta_current, omega_current, ratios)$nu + #{q(Y_prime|current)}+
    prior_prime
  
  den =  un_loglik_given(dim, Y, lambda_current, nu_current, delta_current, omega_current, ratios)$nu + #{q(Y|current)}
    un_loglik_given(dim, Y_prime, lambda_current, nu_prime, delta_current, omega_current, ratios_prime)$nu + #{q(Y_prime|prime)}
    prior_current
  
  #Acceptance:
  alpha = num - den + q_ratio_nu;
  
  if(log(runif(1)) < alpha){
    nu_out = nu_prime
    ratios_out = ratios_prime
    draws_r_out = draws_prime_r
    accepted = 1;
  } else {
    nu_out = nu_current
    ratios_out = ratios
    draws_r_out = draws_r
    accepted = 0;
  }
  
  output = list("nu" = nu_out[dim], 
                "accepted" = accepted, 
                "ratios" = ratios_out,
                "draws_r" = draws_r_out)
  return(output)
}

exchange_omega= function(Y, lambda_current, nu_current, omega_current, delta_current, ratios, draws_r, sigma_omega, prior_omega, tol){
  
  parameters = list('lambda' = lambda_current, 'nu' = nu_current, 'omega'= omega_current, 'delta' = delta_current)
  d= ncol(Y);  n = nrow(Y)
  
  update = ln_proposal(omega_current, sigma_omega)
  omega_prime = update$proposed
  q_ratio = update$log_q_ratio
  
  #Ratios: with recycled draws
  ratios_prime = ratios_IS_recycle(lambda_current, omega_prime, draws_r)
  
  #Exchange draws:
  Y_prime = rmcomp_T_given(n, lambda_current, nu_current, delta_current, omega_prime, ratios_prime, tol)
  
  #Check if kernel is computable and yields proper pmf
  current_list=list(); current_list$params = parameters; current_list$ratios = ratios; prime_list = current_list;
  prime_list$ratios = ratios_prime
  prime_list$params$omega = omega_prime
  prime_list$Y = Y_prime
  
  check = check_move(current_list, prime_list, Y)
  while(check == FALSE){
    update = ln_proposal(omega_current, sigma_omega)
    omega_prime = update$proposed
    q_ratio = update$log_q_ratio
    
    #Ratios: with recycled draws
    ratios_prime = ratios_IS_recycle(lambda_current, omega_prime, draws_r)
    
    #Exchange draws:
    Y_prime = rmcomp_T_given(n, lambda_current, nu_current, delta_current, omega_prime, ratios_prime, tol)
    
    #Check if kernel is computable and yields proper pmf
    prime_list$ratios = ratios_prime
    prime_list$params$omega = omega_prime
    prime_list$Y = Y_prime
    
    check = check_move(current_list, prime_list, Y)
  }
  
  #Priors
  prior_current =  dgamma(omega_current, log = TRUE, shape = prior_omega$shape, rate =prior_omega$rate)
  prior_prime =  dgamma(omega_prime, log = TRUE, shape = prior_omega$shape, rate =prior_omega$rate)
  
  y_log_ratio =  un_loglik_given(1, Y, lambda_current, nu_current, delta_current, omega_prime, ratios_prime)$log_kernel-
    un_loglik_given(1, Y, lambda_current, nu_current, delta_current, omega_current, ratios)$log_kernel 
  
  yprime_log_ratio =  un_loglik_given(1, Y_prime, lambda_current, nu_current, delta_current, omega_current, ratios)$log_kernel-
    un_loglik_given(1, Y_prime, lambda_current, nu_current, delta_current, omega_prime, ratios_prime)$log_kernel
  
  # cat("Omega current, proposed", round(omega_current,3), round(omega_prime,3), "\n",
  #     'Exchange Y:', round(Y_exchange,3),"\n",
  #     'Exchange Yprime:', round(Yprime_exchange,3),"\n",
  #     "Ratios prime:", round(ratios_prime,3))
  
  alpha = y_log_ratio + yprime_log_ratio + q_ratio + prior_prime - prior_current;
  
  if(log(runif(1)) < alpha){
    omega_out = omega_prime
    ratios_out = ratios_prime
    accepted = 1;
  } else {
    omega_out = omega_current
    ratios_out = ratios
    accepted = 0;
  }
  
  return(list("omega" = omega_out, "accepted" = accepted, "ratios" = ratios_out))
  
}

exchange_omega2= function(Y, lambda_current, nu_current, omega_current, delta_current, ratios, draws_r, sigma_omega, prior_omega, tol){
  
  parameters = list('lambda' = lambda_current, 'nu' = nu_current, 'omega'= omega_current, 'delta' = delta_current)
  d= ncol(Y);  n = nrow(Y)
  
  update = ln_proposal(omega_current, sigma_omega)
  omega_prime = update$proposed
  q_ratio = update$log_q_ratio
  
  #Ratios: with recycled draws
  ratios_prime = ratios_IS_recycle(lambda_current, omega_prime, draws_r)
  
  #Exchange draws:
  Y_prime = rmcomp_T_given(n, lambda_current, nu_current, delta_current, omega_prime, ratios_prime, tol)
  
  #Check if kernel is computable and yields proper pmf
  current_list=list(); current_list$params = parameters; current_list$ratios = ratios; prime_list = current_list;
  prime_list$ratios = ratios_prime
  prime_list$params$omega = omega_prime
  prime_list$Y = Y_prime
  
  check = check_move(current_list, prime_list, Y)
  while(check == FALSE){
    update = ln_proposal(omega_current, sigma_omega)
    omega_prime = update$proposed
    q_ratio = update$log_q_ratio
    
    #Ratios: with recycled draws
    ratios_prime = ratios_IS_recycle(lambda_current, omega_prime, draws_r)
    
    #Exchange draws:
    Y_prime = rmcomp_T_given(n, lambda_current, nu_current, delta_current, omega_prime, ratios_prime, tol)
    
    #Check if kernel is computable and yields proper pmf
    prime_list$ratios = ratios_prime
    prime_list$params$omega = omega_prime
    prime_list$Y = Y_prime
    
    check = check_move(current_list, prime_list, Y)
  }
  
  #Priors
  prior_current = log(dtruncnorm(omega_current, a = 0, mean = 0, sd = prior_omega$sd))
  prior_prime =  log(dtruncnorm(omega_prime, a = 0, mean = 0, sd = prior_omega$sd))
  
  y_log_ratio =  un_loglik_given(1, Y, lambda_current, nu_current, delta_current, omega_prime, ratios_prime)$log_kernel-
    un_loglik_given(1, Y, lambda_current, nu_current, delta_current, omega_current, ratios)$log_kernel 
  
  yprime_log_ratio =  un_loglik_given(1, Y_prime, lambda_current, nu_current, delta_current, omega_current, ratios)$log_kernel-
    un_loglik_given(1, Y_prime, lambda_current, nu_current, delta_current, omega_prime, ratios_prime)$log_kernel
  
  alpha = y_log_ratio + yprime_log_ratio + q_ratio + prior_prime - prior_current;
  
  if(log(runif(1)) < alpha){
    omega_out = omega_prime
    ratios_out = ratios_prime
    accepted = 1;
  } else {
    omega_out = omega_current
    ratios_out = ratios
    accepted = 0;
  }
  
  return(list("omega" = omega_out, "accepted" = accepted, "ratios" = ratios_out))
  
}

exchange_delta = function(index, combs, Y, lambda_current, nu_current, omega_current, delta_current, ratios, draws_r, sigma_delta, prior_delta, tol){
  
  n = nrow(Y)
  parameters = list('lambda' = lambda_current, 'nu' = nu_current, 'omega'= omega_current, 'delta' = delta_current)
  i = combs[index,1];j = combs[index,2]
  range_delta = delta_limits_cpp(ratios[c(i,j)])
  
  delta_prime = delta_current #create a copy
  delta_prime[i,j] = rtruncnorm(1, a = range_delta[1], b = range_delta[2], mean = delta_current[i,j], sd = sigma_delta[i,j])
  
  #Exchange draws
  Y_prime = rmcomp_T_given(n, lambda_current, nu_current, delta_prime, omega_current, ratios, tol)
  
  #PMF check
  current_list=list(); current_list$params = parameters; current_list$ratios = ratios; prime_list = current_list;
  prime_list$params$delta = delta_prime
  prime_list$Y = Y_prime
  check = check_move(current_list, prime_list, Y)
  
  while(check == FALSE){
    delta_prime[i,j] = rtruncnorm(1, a = range_delta[1], b = range_delta[2], mean = delta_current[i,j], sd = sigma_delta[i,j])
    
    #Exchange draws
    Y_prime = rmcomp_T_given(n, lambda_current, nu_current, delta_prime, omega_current, ratios, tol)
    
    #PMF check
    prime_list$params$delta = delta_prime
    prime_list$Y = Y_prime
    check = check_move(current_list, prime_list, Y)
  }
  
  #Priors
  limits = delta_limits_cpp(ratios[c(i,j)])
  prior_current = log(dtruncnorm(delta_current[i,j], mean = 0, sd = prior_delta$sd, a = limits[1], b = limits[2]))
  prior_prime = log(dtruncnorm(delta_prime[i,j], mean = 0, sd = prior_delta$sd, a = limits[1], b = limits[2]))
  
  Y_exchange =  un_loglik_given(1, Y, lambda_current, nu_current, delta_prime, omega_current, ratios)$kernel/ 
    un_loglik_given(1, Y, lambda_current, nu_current, delta_current, omega_current, ratios)$kernel 
  
  Yprime_exchange =  un_loglik_given(1, Y_prime, lambda_current, nu_current, delta_current, omega_current, ratios)$kernel/
    un_loglik_given(1, Y_prime, lambda_current, nu_current, delta_prime, omega_current, ratios)$kernel
  
  alpha =  min(c(1, (Y_exchange*Yprime_exchange)) )
  
  if(runif(1) < alpha){
    delta_out = delta_prime
    accepted = 1;
  } else {
    delta_out = delta_current
    accepted = 0;
  }
  
  return(list("delta" = delta_out, "accepted" = accepted))
  
}

exchange_delta2 = function(index, combs, Y, lambda_current, nu_current, omega_current, delta_current, ratios, draws_r, sigma_delta, tol){
  
  n = nrow(Y)
  parameters = list('lambda' = lambda_current, 'nu' = nu_current, 'omega'= omega_current, 'delta' = delta_current)
  i = combs[index,1];j = combs[index,2]
  range_delta = delta_limits_cpp(ratios[c(i,j)])
  
  delta_prime = delta_current #create a copy
  delta_prime[i,j] = rtruncnorm(1, a = range_delta[1], b = range_delta[2], mean = delta_current[i,j], sd = sigma_delta[i,j])
  
  #Exchange draws
  Y_prime = rmcomp_T_given(n, lambda_current, nu_current, delta_prime, omega_current, ratios, tol)
  
  #PMF check
  current_list=list(); current_list$params = parameters; current_list$ratios = ratios; prime_list = current_list;
  prime_list$params$delta = delta_prime
  prime_list$Y = Y_prime
  check = check_move(current_list, prime_list, Y)
  
  while(check == FALSE){
    delta_prime[i,j] = rtruncnorm(1, a = range_delta[1], b = range_delta[2], mean = delta_current[i,j], sd = sigma_delta[i,j])
    
    #Exchange draws
    Y_prime = rmcomp_T_given(n, lambda_current, nu_current, delta_prime, omega_current, ratios, tol)
    
    #PMF check
    prime_list$params$delta = delta_prime
    prime_list$Y = Y_prime
    check = check_move(current_list, prime_list, Y)
  }
  
  Y_exchange =  un_loglik_given(1, Y, lambda_current, nu_current, delta_prime, omega_current, ratios)$kernel/ 
    un_loglik_given(1, Y, lambda_current, nu_current, delta_current, omega_current, ratios)$kernel 
  
  Yprime_exchange =  un_loglik_given(1, Y_prime, lambda_current, nu_current, delta_current, omega_current, ratios)$kernel/
    un_loglik_given(1, Y_prime, lambda_current, nu_current, delta_prime, omega_current, ratios)$kernel
  
  alpha =  min(c(1, (Y_exchange*Yprime_exchange)) )
  
  if(runif(1) < alpha){
    delta_out = delta_prime
    accepted = 1;
  } else {
    delta_out = delta_current
    accepted = 0;
  }
  
  return(list("delta" = delta_out, "accepted" = accepted))
  
}



Exchange = function(n_iter, burn_in, Y, N_aux_r, priors_list, initialise = 'regression', tol =0.001, chains=1, ncores = 3){
  
 if(chains == 1 & length(priors_list$lambda)==2){
   start = Sys.time()
   run = Exchange_wrapper(n_iter, burn_in, Y, N_aux_r, initialise, priors_list, tol)
   end = Sys.time()
   timing = end-start
   
   chains = list(run)
 } else if(chains == 1 & length(priors_list$lambda)==1){
   
   start = Sys.time()
   run = Exchange_wrapper2(n_iter, burn_in, Y, N_aux_r, initialise, priors_list, tol)
   end = Sys.time()
   timing = end-start
   
   chains = list(run)
 } else if(chains>1 & length(priors_list$lambda)==2 ){
   
   start = Sys.time()
   cl = parallel::makeCluster(ncores, setup_strategy = "sequential")
   parallel::clusterEvalQ(cl, library("multcp"))
   parallel::clusterEvalQ(cl, library("truncnorm"))
   parallel::clusterEvalQ(cl, library("COMPoissonReg"))
   parallel::clusterEvalQ(cl, library("mvtnorm"))
   
   chains = parallel::parSapply(cl, 1:chains, 
                                function(times, n_iter, burn_in, Y, N_aux_r, initialise, priors_list, tol){
                                  Exchange_wrapper(n_iter, burn_in, Y, N_aux_r, initialise, priors_list, tol) },
                                n_iter, burn_in, Y, N_aux_r, initialise, priors_list, tol,
                                simplify = F)
   
   parallel::stopCluster(cl) 
   end = Sys.time()
   timing = end-start
   
 } else if(chains >1 & length(priors_list$lambda)==1){
   
   start = Sys.time()
   cl = parallel::makeCluster(ncores, setup_strategy = "sequential")
   parallel::clusterEvalQ(cl, library("multcp"))
   parallel::clusterEvalQ(cl, library("truncnorm"))
   parallel::clusterEvalQ(cl, library("COMPoissonReg"))
   parallel::clusterEvalQ(cl, library("mvtnorm"))
   
   chains = parallel::parSapply(cl, 1:chains, 
                                function(times, n_iter, burn_in, Y, N_aux_r, initialise, priors_list, tol){
                                  Exchange_wrapper2(n_iter, burn_in, Y, N_aux_r, initialise, priors_list, tol) },
                                n_iter, burn_in, Y, N_aux_r, initialise, priors_list, tol,
                                simplify = F)
   
   parallel::stopCluster(cl) 
   end = Sys.time()
   timing = end-start
   
 }
   
  mcmc = list("mcmc" = chains,"time" = timing)
  return(mcmc)
  
}


Exchange_wrapper2 = function(n_iter, burn_in, Y, N_aux_r, initialise, priors_list, tol){
  
  n = nrow(Y); d = ncol(Y);
  
  # Initial values ----------------------------------------------------------------------------
  delta_current = matrix(0, ncol =d, nrow = d)
  combs = t(combn(seq(1,d,1),2)) 
  
  if(initialise == 'random'){
    lambda_current = sample(seq(0.5, 2, 0.1),d)
    nu_current = sample(seq(0.5, 2, 0.1),d)
    
    while( any(lambda_current^(1/nu_current)>3) ){
      lambda_current = sample(seq(0.5, 2, 0.1),d)
      nu_current = sample(seq(0.5, 2, 0.1),d)
    }
    delta_current = matrix(0, ncol = d, nrow = d)
    delta_current[upper.tri(delta_current)] = rnorm(choose(d,2), mean = 0, sd = 0.5)
    omega_current = sample(seq(0.1, 2, 0.5),1)
    
    inits = list("lambda" = lambda_current,"nu" = nu_current,
                 "delta" = delta_current,"omega" = omega_current)
  } else {
    inits= sapply(1:d, inits_marginal, Y) 
    lambda_current = inits[1,]
    nu_current = inits[2,]
    omega_current = 1
    for(i in 1:nrow(combs)){
      delta_current[combs[i,1], combs[i,2]] = sign(cor(Y[,combs[i,1]], Y[,combs[i,2]]))
    }
    
  }
  
  cat('Initial values: ------------------',"\n",
      "lambda:", round(lambda_current,3), '\n',
      'nu:', round(nu_current,3), '\n',
      'delta:' , round(delta_current,3), '\n', 
      'omega:', round(omega_current,3), '\n',
      '----------------------------------', '\n')
  
  # Initialise chains/Rates ---------------------------------------------------------------------
  target_ac = 0.44
  sigma_omega = 0.05;  
  sigma_delta = matrix(0.1, ncol = d, nrow =d) 
  sigma_delta[lower.tri(sigma_delta,TRUE)]<- NA
  sigma_lambda = rep(0.1,d)
  sigma_nu = rep(0.1,d)
  
  accept_omega = 0;
  accept_delta = matrix(0, ncol = d, nrow =d)
  accept_delta[lower.tri(accept_delta,TRUE)]<- NA
  accept_lambda = rep(0, d);
  accept_nu = rep(0, d);
  
  lambda_chain = matrix(NA, nrow = burn_in+n_iter, ncol = d)
  nu_chain = matrix(NA, nrow = burn_in+n_iter, ncol = d)
  omega_chain = rep(NA, burn_in+n_iter); 
  delta_chain = rep(list(matrix(NA, nrow = d, ncol = d)), burn_in+n_iter)
  
  # MCMC -----------------------------------------------------------------------------------------------
  
  # Compute ratios and save draws
  ratio_estimates = ratios_IS_cpp(lambda_current, nu_current, omega_current, N_aux_r) 
  ratios = ratio_estimates$ratios
  draws_r = ratio_estimates$draws
  
  iter =1
  start_time = Sys.time()
  while(iter <= (burn_in+n_iter)){
    
    # ----------------- 1. Lambda Update  --------------------
    for(dim in 1:d){
      update_l_d = exchange_lambda2(dim, Y, lambda_current, nu_current, omega_current, delta_current, ratios, draws_r, sigma_lambda, priors_list$lambda,tol)
      
      lambda_current[dim] =update_l_d$lambda
      ratios = update_l_d$ratios
      draws_r = update_l_d$draws_r
      
      accept_lambda[dim] =   accept_lambda[dim] + update_l_d$accepted
      sigma_lambda[dim] = update_propsd_double(sigma_lambda[dim],  accept_lambda[dim]/iter,  iter, target_ac)
      
    }
    
    # ----------------- 2. Nu Update ------------------------------------
    for(dim in 1:d){
      update_nu_d = exchange_nu2(dim, Y, lambda_current, nu_current, omega_current, delta_current, ratios, draws_r, sigma_nu, priors_list$nu, tol)
      
      nu_current[dim] =update_nu_d$nu
      ratios = update_nu_d$ratios
      draws_r = update_nu_d$draws_r
      
      accept_nu[dim] =   accept_nu[dim] + update_nu_d$accepted
      sigma_nu[dim] = update_propsd_double(sigma_nu[dim],  accept_nu[dim]/iter,  iter, target_ac)
      
    }
    
    # ----------------- 3. Omega Update ------------------------------------
    omega_update = exchange_omega2(Y, lambda_current, nu_current, omega_current, delta_current, ratios, draws_r, sigma_omega, priors_list$omega,tol)
    omega_current = omega_update$omega
    ratios = omega_update$ratios; 
    
    #Update sd
    accept_omega = accept_omega + omega_update$accepted
    sigma_omega <- update_propsd_double(sigma_omega,  accept_omega/iter,  iter, target_ac)
    
    # ----------------- 4. Update delta ---------------------------
    for(index in 1:nrow(combs)){
      update_d = exchange_delta2(index, combs, Y, lambda_current, nu_current, omega_current, delta_current, ratios, draws_r, sigma_delta,tol)
      
      delta_current = update_d$delta
      accept_delta[combs[index,1],combs[index,2]]<- accept_delta[combs[index,1],combs[index,2]]+ update_d$accepted
      
      #Update sd
      sigma_delta[combs[index,1],combs[index,2]] <- update_propsd_double(sigma_delta[combs[index,1],combs[index,2]], 
                                                                         accept_delta[combs[index,1],combs[index,2]]/iter,
                                                                         iter, target_ac)
    }
    
    # Store -------------------------------------------------------
    lambda_chain[iter,] = lambda_current
    nu_chain[iter,] = nu_current
    delta_chain[[iter]] = delta_current
    omega_chain[iter] = omega_current
    
    iter = iter +1
    
    if(iter %% 100==0){
      
      cat("ITERATION:", iter, "----------------------------", "\n")
      
      cat("Lambda report  --------------------------------:","\n",
          "Accepted lambda", round(accept_lambda/iter,3),"\n",
          "Current lambda:", round(lambda_current,3), "\n",
          
          "Nu report ---------------------------------------:","\n",
          "Accepted nu", round(accept_nu/iter,3),"\n",
          "Current nu", round(nu_current,3),"\n",
          
          "KERNEL --------------------------------------------:","\n",
          "Omega accepted", round(accept_omega/iter,3),"\n",
          "Current omega:", round(omega_current,3), "\n",
          "Sigma omega:", round(sigma_omega,3), "\n",
          
          "Delta accepted",  round(accept_delta/iter,3),"\n",
          "Current delta:", round(delta_current[1,2],3), "\n",
          
          "-------------------------------------------------------","\n",
          "-------------------------------------------------------")
    }
    
  }
  end_time = Sys.time();
  
  # Return arguments --------------------------------------------------------------------
  total_iter = burn_in + n_iter
  output = list('ac_rates' = list('lambda' = accept_lambda/total_iter,
                                  'nu' = accept_nu/total_iter,
                                  'omega' = accept_omega/total_iter,
                                  'delta' = accept_delta/total_iter),
                "time" = end_time-start_time,
                "Nr" = N_aux_r, "initial_values" = inits,
                'proposal_parameters' = list("lambda" = sigma_lambda,
                                             'nu' = sigma_nu,
                                             "omega" = sigma_omega,
                                             "delta" = sigma_delta),
                "lambda_chain" = lambda_chain[(burn_in+1):total_iter, ],
                "nu_chain" = nu_chain[(burn_in+1):total_iter, ],
                "omega_chain" = omega_chain[(burn_in+1):total_iter],
                "delta_chain" = delta_chain[(burn_in+1):total_iter])
  return(output)
}

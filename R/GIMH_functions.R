
update_omega = function(Y, lambda_current, nu_current, delta_current, omega_current, sigma_omega, ratios, log_inv_z, draws_r, prior_omega){
  
  d= ncol(Y);  n = nrow(Y)
  update = ln_proposal(omega_current, sigma_omega)
  omega_prime = update$proposed
  q_ratio = update$log_q_ratio
  
  #Ratios: with recycled draws
  ratios_prime = ratios_IS_recycle(lambda_current, omega_prime, draws_r)
  
  check = sum(is.na(exp(loglik_given_est(Y, lambda_current, nu_current, delta_current, omega_prime, ratios_prime, log_inv_z))))+
    sum(is.na(ratios_prime)) + sum(is.infinite(ratios_prime))
  
  while(check != 0){
    update = ln_proposal(omega_current, sigma_omega)
    omega_prime = update$proposed
    q_ratio = update$log_q_ratio
    
    #Ratios: with recycled draws
    ratios_prime = ratios_IS_recycle(lambda_current, omega_prime, draws_r)
    
    check = sum(is.na(exp(loglik_given_est(Y, lambda_current, nu_current, delta_current, omega_prime, ratios_prime, log_inv_z))))+
      sum(is.na(ratios_prime)) + sum(is.infinite(ratios_prime))
  }
  
  #Priors
  prior_current = dgamma(omega_current, log = TRUE, shape = prior_omega$shape, rate = prior_omega$rate)
  prior_prime = dgamma(omega_prime, log = TRUE, shape = prior_omega$shape, rate = prior_omega$rate)
  
  #Loglik estimates
  loglik_prime = un_loglik_given(1, Y, lambda_current, nu_current, delta_current, omega_prime, ratios_prime)$log_kernel
  loglik_current = un_loglik_given(1, Y, lambda_current, nu_current, delta_current, omega_current, ratios)$log_kernel
  
  num = loglik_prime + prior_prime 
  den = loglik_current + prior_current
  alpha = num - den + q_ratio
  
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


update_delta = function(index, combs, Y, lambda_current, nu_current, delta_current, omega_current, sigma_delta, ratios, prior_delta){
  
  d = ncol(Y); n = nrow(Y)
  i = combs[index,1];j = combs[index,2]
  
  range_delta = delta_limits_cpp(ratios[c(i,j)])
  
  delta_prime = delta_current #create a copy
  delta_prime[i,j] = rtruncnorm(1, a = range_delta[1], b = range_delta[2], mean = delta_current[i,j], sd = sigma_delta[i,j])
  
  #Loglik
  loglik_prime = un_loglik_given(1, Y, lambda_current, nu_current, delta_prime, omega_current, ratios)$log_kernel
  loglik_current = un_loglik_given(1, Y, lambda_current, nu_current, delta_current, omega_current, ratios)$log_kernel
  
  #Acceptance:
  num = loglik_prime 
  den = loglik_current 
  alpha = num - den
  
  check = is.numeric(alpha) & !is.nan(alpha)
  while(check == FALSE){
    delta_prime[i,j] = rtruncnorm(1, a = range_delta[1], b = range_delta[2], mean = delta_current[i,j], sd = sigma_delta[i,j])
    
    #Loglik
    loglik_prime = un_loglik_given(1, Y, lambda_current, nu_current, delta_prime, omega_current, ratios)$log_kernel
    loglik_current = un_loglik_given(1, Y, lambda_current, nu_current, delta_current, omega_current, ratios)$log_kernel
    
    #Acceptance:
    num = loglik_prime 
    den = loglik_current 
    alpha = num - den
    
    check = is.numeric(alpha) & !is.nan(alpha)
  }
  
  if(log(runif(1)) < alpha){
    delta_out = delta_prime
    accepted = 1;
  } else {
    delta_out = delta_current
    accepted = 0;
  }
  
  return(list("delta" = delta_out, "accepted" = accepted))
  
}


update_marginal_d = function(dim, Y, lambda_current, nu_current, omega_current, delta_current, ratios, log_inv_z, draws_r, draws_z, Sigma_marginals, prior_marginal){
  
  d = ncol(Y); n = nrow(Y)
  N_aux_r = nrow(draws_r); 
  N_aux_z = nrow(draws_z)
  Sigma_prop = Sigma_marginals[[dim]]
  
  lambda_prime = lambda_current; nu_prime = nu_current
  ratios_prime = ratios; log_inv_z_prime = log_inv_z;
  draws_r_prime = draws_r; draws_z_prime = draws_z;
  
  marginal_current = c(lambda_current[dim], nu_current[dim])
  marginal_prime = exp(mvtnorm::rmvnorm(1, log(marginal_current), Sigma_prop) )  
  
  lambda_prime[dim] = marginal_prime[1]
  nu_prime[dim] = marginal_prime[2]
  
  ratios_run = ratio_IS_d_cpp(dim, lambda_prime, nu_prime, omega_current, N_aux_r)
  linvz_run = est_logZinv_d_cpp(dim, lambda_prime, nu_prime, N_aux_z)
  ratios_prime[dim] = ratios_run$ratio; draws_r_prime[,dim] = ratios_run$draws
  log_inv_z_prime[dim] = linvz_run$log_inv_z; draws_z_prime[,dim] = linvz_run$draws
  
  check =  sum(is.na(exp(loglik_given_est(Y, lambda_prime, nu_prime, delta_current, omega_current, ratios_prime, log_inv_z_prime))))
  
  while(check != 0){
    marginal_prime = exp(mvtnorm::rmvnorm(1, log(marginal_current), Sigma_prop) )  
    
    lambda_prime[dim] = marginal_prime[1]
    nu_prime[dim] = marginal_prime[2]
    
    ratios_run = ratio_IS_d_cpp(dim, lambda_prime, nu_prime, omega_current, N_aux_r)
    linvz_run = est_logZinv_d_cpp(dim, lambda_prime, nu_prime, N_aux_z)
    ratios_prime[dim] = ratios_run$ratio; draws_r_prime[,dim] = ratios_run$draws
    log_inv_z_prime[dim] = linvz_run$log_inv_z; draws_z_prime[,dim] = linvz_run$draws
    
    check =  sum(is.na(exp(loglik_given_est(Y, lambda_prime, nu_prime, delta_current, omega_current, ratios_prime, log_inv_z_prime))))
  }
  
  q_curr_to_prime = mvtnorm::dmvnorm(log(marginal_prime), log(marginal_current), Sigma_prop, log=TRUE) -sum(log(marginal_prime))
  q_prime_to_curr = mvtnorm::dmvnorm(log(marginal_current), log(marginal_prime), Sigma_prop,log=TRUE) - sum(log(marginal_current))
  q_ratio = q_prime_to_curr - q_curr_to_prime
  
  #Priors
  prior_current = dgamma(marginal_current[1], shape = prior_marginal$lambda$shape,rate = prior_marginal$lambda$rate, log=TRUE)+
                  dgamma(marginal_current[2], shape = prior_marginal$nu$shape,rate = prior_marginal$nu$rate, log=TRUE)
  
  prior_prime= dgamma(marginal_prime[1], shape = prior_marginal$lambda$shape,rate = prior_marginal$lambda$rate, log=TRUE)+
    dgamma(marginal_prime[2], shape = prior_marginal$nu$shape,rate = prior_marginal$nu$rate, log=TRUE)
  
  #Loglik estimates
  loglik_prime = un_loglik_given(dim, Y, lambda_prime, nu_prime, delta_current, omega_current, ratios_prime)
  loglik_prime = sum(unlist(loglik_prime[c('lambda', 'nu')])) - loglik_prime$log_kernel + n*log_inv_z_prime[dim]
  
  loglik_current = un_loglik_given(dim, Y, lambda_current, nu_current, delta_current, omega_current, ratios)
  loglik_current = sum(unlist(loglik_current[c('lambda', 'nu')])) - loglik_current$log_kernel + n*log_inv_z[dim]
  
  #Acceptance:
  num = loglik_prime + prior_prime 
  den = loglik_current + prior_current
  alpha = num - den + q_ratio
  
  if(log(runif(1)) < alpha){
    params_out = marginal_prime
    ratios_out = ratios_prime
    log_inv_z_out = log_inv_z_prime
    draws_r_out = draws_r_prime
    draws_z_out = draws_z_prime
    accepted = 1;
  } else {
    params_out = marginal_current
    ratios_out = ratios
    log_inv_z_out = log_inv_z
    draws_r_out = draws_r
    draws_z_out = draws_z
    accepted = 0;
  }
  
  output = list("lambda" = params_out[1],
                'nu' = params_out[2],
                "accepted" = accepted, 
                "ratios" = ratios_out, "log_inv_z" = log_inv_z_out, 
                "draws_z" = draws_z_out, "draws_r" = draws_r_out)
  return(output)
}

GIMH_wrapper = function(Y, burn_in, n_iter, N_aux_r, N_aux_z, initialise, priors_list){
  
  target_ac=0.44
  n = nrow(Y); d = ncol(Y)
  combs = t(combn(seq(1,d,1),2))
  
  # Starting values ------------------------------------------------------------------------
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
  }else if(initialise == "uni_regression"){
    
    inits = list()
    inits_ = sapply(1:ncol(Y), inits_marginal, Y) 
    inits$lambda = inits_[1,]
    inits$nu = inits_[2,]
    delta_init = sign(cor(Y))
    delta_init[lower.tri(delta_init, diag = TRUE)]<-0
    inits$delta = delta_init
    inits$omega =  1
    
  }
  lambda_current = inits$lambda
  nu_current = inits$nu
  omega_current = inits$omega
  delta_current= inits$delta
  
  cat('Initial values: ---------------',"\n",
      "lambda:", lambda_current, '\n',
      'nu:', nu_current, '\n',
      'omega', omega_current, '\n',
      'delta', c(delta_current[upper.tri(delta_current)]), '\n',
      '----------------------------------', '\n')
  
  # Nz ------------------------------------------------------------------------
  
  if(class(N_aux_z)=='list'){   #Run adaptation if the list is supplied
    target_sd =N_aux_z$target_sd
    N_aux_z = N_aux_z$init
    cur_sd = 99; 
    while(cur_sd > target_sd){
      
      loglik_vec = replicate(100, estimate_loglik(Y, lambda_current, nu_current, delta_current, omega_current, N_aux_z, N_aux_r))
      cur_sd = sd(loglik_vec)
      
      cat('"Running Adaptation of Nz: -----',"\n",
          "Current loglik sd:", round(cur_sd,3),"\n",
          "Target loglik sd:", target_sd, "\n",
          "Nz:", N_aux_z, '\n',
        "------------------------------", "\n")
      
      adapt = cur_sd/target_sd
      N_aux_z = ifelse(adapt > 1,  ceiling(adapt*N_aux_z), N_aux_z)
    }
  }else{
    N_aux_z= N_aux_z
  }
  
  cat('Auxiliary draws: --------------------', '\n',
      'Nz:', N_aux_z, "\n",
      'Nr:', N_aux_r, "\n",
      '-------------------------------------', '\n')
  
  # Compute ratios and save draws
  ratio_estimates = ratios_IS_cpp(lambda_current, nu_current, omega_current, N_aux_r) 
  ratios = ratio_estimates$ratios
  draws_r = ratio_estimates$draws
  
  # Compute inverse constants and save draws
  inv_Z_estimates = est_logZinv_cpp(lambda_current, nu_current, N_aux_z)
  log_inv_z = inv_Z_estimates$loginvZ
  draws_z = inv_Z_estimates$draws
  
  Sigma_marginals =  rep(list(matrix(c(0.1, 0, 0, 0.1), ncol=2, nrow =2)),d)
  sigma_omega = 0.1
  sigma_delta = matrix(0, nrow=d, ncol =d)
  sigma_delta[upper.tri(sigma_delta)]<- 0.1
  
  accept_marginal = rep(0,d);
  accept_omega = 0;
  accept_delta = matrix(0, ncol = d, nrow =d)
  accept_delta[lower.tri(accept_delta,TRUE)]<- NA
  
  lambda_chain = matrix(NA, nrow = burn_in+n_iter, ncol = d)
  nu_chain = matrix(NA, nrow = burn_in+n_iter, ncol = d)
  omega_chain = rep(NA, burn_in+n_iter); 
  delta_chain = rep(list(matrix(NA, nrow = d, ncol = d)), burn_in+n_iter)
  loglik_chain = matrix(NA, ncol = n, nrow = burn_in+n_iter)
  
  max_ll= sum(loglik_given_est(Y, lambda_current, nu_current, delta_current, omega_current, ratios, log_inv_z))
  
  iter = 1
  start_time = Sys.time()
  while(iter <= (n_iter+burn_in)){
    
    # ----------------- 1. Marginal Update  ----------------------------
    for(dim in 1:d){
      joint_d = update_marginal_d(dim, Y, lambda_current, nu_current, omega_current, delta_current,
                                  ratios, log_inv_z, draws_r, draws_z, Sigma_marginals, priors_list)
      lambda_current[dim] = joint_d$lambda
      nu_current[dim] = joint_d$nu
      
      ratios = joint_d$ratios
      log_inv_z = joint_d$log_inv_z; 
      draws_z = joint_d$draws_z; draws_r = joint_d$draws_r
      
      accept_marginal[dim] = accept_marginal[dim]  + joint_d$accepted
      Sigma_marginals[[dim]][1,1] = Sigma_marginals[[dim]][2,2]  = update_propsd_double( Sigma_marginals[[dim]][1,1],
                                                                                         accept_marginal[dim]/iter,  iter, target_ac)
      
    }
    
    # ----------------- 2. Update omega ----------------------------
    omega_update = update_omega(Y, lambda_current, nu_current, delta_current, omega_current, sigma_omega, ratios, log_inv_z, draws_r, priors_list$omega)
    omega_current = omega_update$omega
    ratios = omega_update$ratios; 
    
    accept_omega = accept_omega + omega_update$accepted
    sigma_omega <- update_propsd_double(sigma_omega,  accept_omega/iter,  iter, target_ac)
    
    # ----------------- 3. Update deltas ---------------------------
    for(index in 1:nrow(combs)){
      update_d = update_delta(index, combs, Y, lambda_current, nu_current, delta_current, omega_current, sigma_delta, ratios, priors_list$delta)
      
      delta_current = update_d$delta
      accept_delta[combs[index,1],combs[index,2]]<- accept_delta[combs[index,1],combs[index,2]]+ update_d$accepted
      sigma_delta[combs[index,1],combs[index,2]] <- update_propsd_double(sigma_delta[combs[index,1],combs[index,2]], 
                                                                         accept_delta[combs[index,1],combs[index,2]]/iter,
                                                                         iter, target_ac)
      
    }
    
    #Storing --------------------------------------------------------------
    
    lambda_chain[iter,] = lambda_current
    nu_chain[iter,] = nu_current
    delta_chain[[iter]] = delta_current
    omega_chain[iter] = omega_current
    loglik_chain[iter,] = loglik_given_est(Y, lambda_current, nu_current, delta_current, omega_current, ratios, log_inv_z)
    
    maxloglik = max(max_ll, sum(loglik_chain[iter,]))
    iter = iter +1
    
    if(iter %% 100==0){
       cat("ITERATION:", iter, "----------------------------", "\n",
           
           "Marginal report ---------:","\n",
           "Accepted", round(accept_marginal/iter,3),"\n",
           "Current lambda:", round(lambda_current,3), "\n",
           "Current nu:", round(nu_current,3), "\n",
           
           "Omega report ---------:","\n",
           "Accepted", round(accept_omega/iter,3),"\n",
           "Current:", round(omega_current,3), "\n",
           
           "Delta report ---------:","\n",
           "Accepted",  round(accept_delta/iter,3),"\n",
           "Current:", round(delta_current,3), "\n",
           "--------------------------------------------------")}
  }
  end_time = Sys.time()
  
  total_iter = burn_in+ n_iter
  output = list('ac_rates' = list('marginal'= accept_marginal/n_iter, 
                                  'omega' = accept_omega/total_iter,
                                  'delta' = accept_delta/total_iter),
                "time" = end_time-start_time,
                "Nz" = N_aux_z, 
                "Nr" = N_aux_r, "initial_values" = inits,
                'proposal_parameters' = list("marginals" = Sigma_marginals,
                                             "omega" = sigma_omega,
                                             "delta" = sigma_delta),
                "lambda_chain" = lambda_chain[(burn_in+1):total_iter, ],
                "nu_chain" = nu_chain[(burn_in+1):total_iter, ],
                "omega_chain" = omega_chain[(burn_in+1):total_iter],
                "delta_chain" = delta_chain[(burn_in+1):total_iter],
                "loglik" = loglik_chain[((burn_in+1):total_iter),] )
  
  return(output)
  
}

GIMH = function(Y, burn_in, n_iter, N_aux_r, N_aux_z, priors_list, initialise = "uni_regression", chains=1, ncores = 3){
  
  if(chains == 1){
    start = Sys.time()
    run = GIMH_wrapper(Y, burn_in, n_iter, N_aux_r, N_aux_z, initialise, priors_list)
    end = Sys.time()
    timing = end-start
    chains = list(run)
  } else {
    
    start = Sys.time()
    cl = parallel::makeCluster(ncores, setup_strategy = "sequential")
    parallel::clusterEvalQ(cl, library("multcp"))
    parallel::clusterEvalQ(cl, library("truncnorm"))
    parallel::clusterEvalQ(cl, library("COMPoissonReg"))
    parallel::clusterEvalQ(cl, library("mvtnorm"))
    
    chains = parallel::parSapply(cl, 1:chains, 
                                 function(times,Y, burn_in, n_iter, N_aux_r, N_aux_z, initialise, priors_list){
                                   GIMH_wrapper(Y, burn_in, n_iter, N_aux_r, N_aux_z, initialise, priors_list)},
                                 Y, burn_in, n_iter, N_aux_r, N_aux_z, initialise, priors_list,
                                 simplify = F)
    
    parallel::stopCluster(cl) 
    end = Sys.time()
    timing = end-start
    
  }
  mcmc = list("mcmc" = chains,"time" = timing)
  return(mcmc)
  
}
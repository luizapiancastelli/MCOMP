
s2_chain_param = function(param_index, phi_ij, phi_means){
  sum((phi_ij[,param_index] - phi_means[param_index])^2)/(nrow(phi_ij-1))
}

s2_chain = function(chain_index, chains){
  n_params = ncol(chains[[chain_index]]$parameters_wide) -1
  phi_ij = chains[[chain_index]]$parameters_wide[,1:n_params]
  phi_means = apply(phi_ij, 2, mean)
  return(sapply(1:n_params, s2_chain_param, phi_ij, phi_means))
}

param_means_fixed_j = function(chain_index, chains){
  n_params = ncol(chains[[chain_index]]$parameters_wide) -1
  means = apply(chains[[chain_index]]$parameters_wide[,1:n_params],2, mean)
  return(means)
}

B_param = function(param_index, means_by_chain){
  sum( ( means_by_chain[param_index,] - mean(means_by_chain[param_index,]) )^2 )
}

diagnostic = function(chains){
  
  n_iter = nrow(chains[[1]]$parameters_wide)
  n_params = ncol(chains[[1]]$parameters_wide) -1
  param_names = names(chains[[1]]$parameters_wide)[1:n_params]
  n_chains = length(chains)
  
  means_by_chain = sapply(1:n_chains,param_means_fixed_j ,chains)
  B = (n_iter/(n_chains-1))*sapply(1:n_params, B_param, means_by_chain)
  
  W = rowMeans(sapply(1:n_chains, s2_chain, chains))
  var_sup = ((n_iter-1)/n_iter)*W + (1/n_iter)*B
  
  Rhat = sqrt(var_sup/W)
  names(Rhat) = param_names

  
  return(Rhat)
  
}

process_mcmc = function(mcmc){
  
  for(chain in 1:length(mcmc)){
    delta_ = mcmc[[chain]]$delta
    
    post_chains = names(mcmc[[chain]])[stringr::str_detect(names(mcmc[[chain]]), "chain")]
    post_chains =  post_chains[!stringr::str_detect(post_chains, "delta")]
    col_names =  gsub("_chain", "", post_chains)
    
    df1 = data.frame(mcmc[[chain]][[post_chains[1]]]);
    start = ifelse( stringr::str_detect(col_names[1], 'gamma'),0,1)
    end = ifelse( stringr::str_detect(col_names[1], 'gamma'),ncol(df1)-1,ncol(df1))
    names(df1) = paste0(col_names[1], start:end)
    
    df2 = data.frame(mcmc[[chain]][[post_chains[2]]]);
    start = ifelse( stringr::str_detect(col_names[2], 'gamma'),0,1)
    end = ifelse( stringr::str_detect(col_names[2], 'gamma'),ncol(df2)-1,ncol(df2))
    names(df2) = paste0(col_names[2], start:end)
    
    df3 = data.frame(mcmc[[chain]][[post_chains[3]]]);
    start = ifelse( stringr::str_detect(col_names[3], 'gamma'),0,1)
    end = ifelse( stringr::str_detect(col_names[3], 'gamma'),ncol(df3)-1,ncol(df3))
    names(df3) = paste0(col_names[3], start:end)
    
    df = cbind(df1, df2, df3);
    names(df)[stringr::str_detect(names(df), "omega")]<- 'omega'
    
    #Handling delta (can be a list)
    if(class(mcmc[[chain]]$delta_chain)=='numeric'){
      df$delta = mcmc[[chain]]$delta_chain
    } else {
      d = ncol(mcmc[[chain]]$lambda_chain)
      combs = t(combn(seq(1,d,1),2))
      
      
      delta_list = list();
      for(i in 1:nrow(combs)){
        delta_list[[i]]= unlist(lapply(mcmc[[chain]]$delta_chain, function(x,i,j){x[[i,j]]}, i =combs[i,1], j=combs[i,2]))
      }
      
      df4 = data.frame(do.call(cbind, delta_list))
      names(df4) = apply(combs, 1, function(x){paste0("delta",paste0(x, collapse = ''))})
      df = cbind(df, df4)
    }
    
    df$iter = 1:nrow(df)
    mcmc[[chain]]$parameters_wide = df
  }
  return(mcmc)
}

posterior_summaries = function(mcmc){
  
  mcmc = process_mcmc(mcmc)
  Rhat = diagnostic(mcmc)
  
  by_chain = lapply(mcmc, "[[", "parameters_wide")
  for(i in 1:length(mcmc)){
    by_chain[[i]]$chain = i
  }
  parameters_wide_all_chains = do.call(rbind, by_chain)
  parameters_long = reshape2::melt(parameters_wide_all_chains, id.vars =c("iter", 'chain') )
  
  density_plots = ggplot2::ggplot(parameters_long, ggplot2::aes(x = value))+ggplot2::geom_density(adjust = 8)+
    ggplot2::facet_wrap(~variable, scales = 'free')+
    ggplot2::theme_bw()
  
  parameters_long$chain = as.factor(parameters_long$chain)
  trace_plots = ggplot2::ggplot(parameters_long, ggplot2::aes(x = iter, y = value, color = chain))+
    ggplot2::geom_line()+ggplot2::theme_bw()+  ggplot2::facet_wrap(~variable, scales = 'free')
  
  stats = parameters_long %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise('Q5' = quantile(value, 0.05),
              'Q10' = quantile(value, 0.10),
              'Median' = median(value), 
              'Q90' = quantile(value, 0.90),
              'Q95' = quantile(value, 0.95),
              'Mean' = mean(value),
              'SD' = sd(value))
  
  return(list('Stats' = stats,
              'Rhat' = Rhat,
              'Density_plots' = density_plots,
              'Trace_plots' = trace_plots))
}

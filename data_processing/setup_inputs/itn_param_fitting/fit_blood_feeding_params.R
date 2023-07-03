# fit_blood_feeding_params.R
# Contact: mambrose
# Created: 2021
# Contents: functions to fit the relationship between hut mortality and blood-feeding from extracted datasets

hbhi_dir = 'C:/Users/moniqueam/Dropbox (IDM)/Malaria Team Folder/projects/burundi_hbhi'
fit_blood_feeding = TRUE



# ======================================== #
# functions for fitting blood-feeding
# ======================================== #

exp_fed = function(k1,k2,hut_mort){
  #' Calculate expected blood-feeding given mortality and parameters
  #' @param k1 fit parameter corresponding to ITN extraction function
  #' @param k2 fit parameter corresponding to ITN extraction function
  #' @param hut_mort the hut mortality value
  return(k1 * exp(k2 * (1-hut_mort-0.5)))
}

get_dif_exp_obs = function(k1,k2,hut_mort,obs_fed){
  #' Calculate the difference between expected and observed blood-feeding
  #' @param k1 fit parameter corresponding to ITN extraction function
  #' @param k2 fit parameter corresponding to ITN extraction function
  #' @param hut_mort the hut mortality value
  #' @param obs_fed the observed blood-feeding value
  return(obs_fed - exp_fed(k1,k2,hut_mort))
}

get_sum_square_dif_dataset = function(k1,k2,all_hut_morts,all_obs_fed){
  #' Calculate the sum of squares for the full dataset given specific parameters
  #' @param k1 fit parameter corresponding to ITN extraction function
  #' @param k2 fit parameter corresponding to ITN extraction function
  #' @param all_hut_morts all hut mortality values
  #' @param all_obs_fed all matched observed blood-feeding values
  ss_each = rep(NA, length(all_hut_morts))
  for(ii in 1:length(all_hut_morts)){
    ss_each[ii] = (get_dif_exp_obs(k1,k2,all_hut_morts[ii], all_obs_fed[ii]))^2
  }
  return(sum(ss_each))
}






#########################################################################
# fit the blood-feeding rate as a function of hut-trial mortality
#########################################################################
if(fit_blood_feeding){
  # read in extracted dataset on hut-trial blood-feeding versus mortality
  dd = read.csv(paste0(hbhi_dir, "/ento/blood feeding and mortality/ITN_extract_BF_mort_hut.csv"))[-1,]
  dd = dd[!is.na(dd$Mortality_rate),]
  dd = dd[!is.na(dd$BF_rate),]
  fit_with_treated_nets_only = TRUE
  
  if(fit_with_treated_nets_only){
    ############ dataset = treated nets #######################
    dd_treat = dd[dd$treated_net == 1,]
    
    # save outputs across all param sets
    k1_vals = seq(0,1,0.005)
    k2_vals = seq(0,4,0.05)
    all_hut_morts = as.numeric(dd_treat$Mortality_rate)/100
    all_obs_fed = as.numeric(dd_treat$BF_rate)/100
    # all_obs_fed = 0.25 * exp(1.5 * (all_hut_morts-0.5)) + runif(n=length(all_hut_morts),-0.2,0.2)
    ss_mat = matrix(NA, nrow=length(k1_vals), ncol=length(k2_vals))
    for(i_k1 in 1:length(k1_vals)){
      for(i_k2 in 1:length(k2_vals)){
        ss_mat[i_k1, i_k2] = get_sum_square_dif_dataset(k1=k1_vals[i_k1],k2=k2_vals[i_k2],all_hut_morts, all_obs_fed)
      }
    }
    
    rownames(ss_mat) = k1_vals
    colnames(ss_mat) = k2_vals
    
    best = which(ss_mat <= (min(ss_mat, na.rm=TRUE)+0.000001), arr.ind=TRUE)
    k1_fit = k1_vals[best[1]]
    k2_fit = k2_vals[best[2]]
    plot((as.numeric(dd_treat$Mortality_rate)/100), (as.numeric(dd_treat$BF_rate)/100), bty='L', ylab='fraction blood fed', xlab='hut trial mortality rate', pch=c(21,20)[as.numeric(dd_treat$treated_net)+1])
    lines(seq(0,1,0.1),exp_fed(k1=k1_vals[best[1]],k2=k2_vals[best[2]],hut_mort=seq(0,1,0.1)), col='dodgerblue', lwd=3)
    lines(seq(0,1,0.1), (1 - (exp(0.04*(1-exp(4.66*(1-seq(0,1,0.1))))/4.66)))/(1-seq(0,1,0.1)), col='darkred', lwd=2)  # from Nash et al. 2021
    
    # text(x=0.6,y=0.9, paste0('k1=', k1_vals[best[1]],'; k2=', k2_vals[best[1]]), col='dodgerblue')
    text(x=0.5,y=0.98, paste0('BF = ', k1_vals[best[1]],' * exp(', k2_vals[best[2]], ' * (0.5 - hut_mort))'), col='dodgerblue', cex=0.8)
    
  } else{
    ################ dataset = all nets ########################
    # save outputs across all param sets
    k1_vals = seq(0,1,0.005)
    k2_vals = seq(0,4,0.05)
    all_hut_morts = as.numeric(dd$Mortality_rate)/100
    all_obs_fed = as.numeric(dd$BF_rate)/100
    ss_mat = matrix(NA, nrow=length(k1_vals), ncol=length(k2_vals))
    for(i_k1 in 1:length(k1_vals)){
      for(i_k2 in 1:length(k2_vals)){
        ss_mat[i_k1, i_k2] = get_sum_square_dif_dataset(k1=k1_vals[i_k1],k2=k2_vals[i_k2],all_hut_morts, all_obs_fed)
      }
    }
    rownames(ss_mat) = k1_vals
    colnames(ss_mat) = k2_vals
    
    best = which(ss_mat <= (min(ss_mat, na.rm=TRUE)+0.000001), arr.ind=TRUE)
    k1_fit = k1_vals[best[1]]
    k2_fit = k2_vals[best[2]]
    
    # plot(all_hut_morts, all_obs_fed)
    plot((as.numeric(dd$Mortality_rate)/100), as.numeric(dd$BF_rate)/100, pch=c(21,20)[as.numeric(dd$treated_net)+1], bty='L', ylab='fraction of mosquitoes blood fed', xlab='hut trial mortality rate')
    lines(seq(0,1,0.1),exp_fed(k1=k1_vals[best[1]],k2=k2_vals[best[2]],hut_mort=seq(0,1,0.1)), col='dodgerblue', lwd=3)
    # text(x=0.6,y=0.9, paste0('k1=', k1_vals[best[1]],'; k2=', k2_vals[best[1]]), col='dodgerblue')
    text(x=0.6,y=0.9, paste0('BF = ', k1_vals[best[1]],' * exp(', k2_vals[best[1]], ' * (0.5 - hut_mort))'), col='dodgerblue')
  }
}

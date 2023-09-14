# add_hut_trial_mort_block.R
# Contact: mambrose
# Created: 2021
# Contents: functions that calculate blood-feeding and killing rates for PBO and non-PBO nets given permethrin bioassay mortality,
#   translate the values into DTK-killing and DTK-blocking, and add the relevant values to ITN simulation input files.

# fixed parameters for calculations - obtained from fit_blood_feeding_params.R
k1_fit = 0.27
k2_fit = 1.3
ig2_kill_rate = 0.75
use_new_ig2 = TRUE  # if TRUE, calculates IG2 values based on expected hut mortality and BF instead of directly getting EMOD kill and block params

# ================================================================== #
# functions for getting dtk parameters from bioassay mortality
# ================================================================== #

get_hut_mortality_from_bioassay_mortality = function(bioassay_mortality, fit_version='loglogistic_Nash'){
  #' Given a bioassay mortality, get the hut mortality using one of several possible functions
  #' @param bioassay_mortality the bioassay mortality value (can be a scalar or vector)
  #' @param fit_version string corresponding to the relationship to use to calculate hut mortality from bioassay mortality
  
  if(fit_version == 'loglogistic_Nash'){
    return((1 / (1 + ((bioassay_mortality)/0.89)^(-0.47))))  # log-logistic fit from Nash et al. 2021
  } else if(fit_version == 'logistic_Nash'){
    return(1 / (1 + exp(-1*((bioassay_mortality) - 0.7)*3.57)))  # logistic fit from Nash et al. 2021
  } else if(fit_version == 'ave_loglog_log_Nash'){  # value in between the log-logistic and logistic fit from Nash et al. 2021
    return(((1 / (1 + ((bioassay_mortality)/0.89)^(-0.47))) + (1 / (1 + exp(-1*((bioassay_mortality) - 0.7)*3.57)))) / 2)  # log-logistic fit from Nash et al. 2021
  } else if(fit_version == 'linear_fit'){
    return(0.25+0.44*bioassay_mortality)  # linear fit data-grabbed from presentation
  } else if(fit_version == '2021_pres'){
    return(0.53-0.1*(1-bioassay_mortality)/(0.2+bioassay_mortality))  # rough guess of log logit from the presentation
  } else if(fit_version == '2016_Churcher'){
    return(1 / (1 + exp(-1*(0.63 + 4.0 * (bioassay_mortality - 0.5)))))  # parameters from Churcher et al 2016
  } else warning('Name for the bioassay- to hut-trial-mortality function not found.')
}

get_hut_BF_from_hut_mortality = function(k1,k2,hut_mortality, fit_version_BF='ITN_extraction'){
  #' Given the hut mortality value, get the hut bloodfeeding estimate using one of several possible functions
  #' @param k1 fit parameter corresponding to ITN extraction function
  #' @param k2 fit parameter corresponding to ITN extraction function
  #' @param hut_mortality the hut mortality value (can be a scalar or vector)
  #' @param fit_version_BF string corresponding to the relationship to use to calculate hut bloodfeeding from hut mortality
  
  if(fit_version_BF == 'ITN_extraction'){
    return(k1 * exp(k2 * (0.5-hut_mortality)))
  } else if(fit_version_BF == 'Nash_2021'){
    return((1 - (exp(0.04*(1-exp(4.66*(1-hut_mortality)))/4.66)))/(1-hut_mortality))  # in ms, they adjusted to get BF and survive by assuming independent mortality. Adjust back here by dividing by (1-P(die)). 
  } else warning('Name for the hut-trial-mortality to blood-feeding-fraction function not found.')
}

get_PBO_bioassay_from_bioassay = function(bioassay_mortality){
  #' Given a permethrin bioassay mortality, get estimate of PBO bioassay mortality
  #' @param bioassay_mortality the bioassay mortality value (can be a scalar or vector)
  
  1 / (1 + exp(-(3.41 + (5.88 * (bioassay_mortality - 0.5) / (1 + 0.78 * (bioassay_mortality - 0.5)))) ))
}

get_dtk_block_kill_from_hut_mort_BF = function(hut_mortality, hut_BF, frac_reassign_feed_survive){
  #' Given the hut trial mortality and bloodfeeding values, estimate the dtk block and kill rates
  #' @param hut_mortality the hut mortality value
  #' @param hut_BF the hut bloodfeeding value
  #' @param frac_reassign_feed_survive the fraction of mosquitoes that fall in the 'feed and die' category that are reassigned to 'feed and survive'
  #' Additional details: reassign the mosquitoes that eat and then die in the hut trials to either eat and survive or die without feeding for the dtk
  #'                     we assume that hut mortality and hut blood-feeding are independent events

  hut_frac_BF_survive = hut_BF * (1-hut_mortality)
  hut_frac_BF_die = hut_BF * hut_mortality
  hut_frac_noBF_survive = (1-hut_BF) * (1-hut_mortality)
  hut_frac_noBF_die = (1-hut_BF) * hut_mortality
  
  dtk_frac_BF_survive = sapply((hut_frac_BF_survive + hut_frac_BF_die * frac_reassign_feed_survive), min, 1)
  dtk_frac_noBF_survive = hut_frac_noBF_survive
  dtk_frac_noBF_die = sapply((hut_frac_noBF_die + hut_frac_BF_die * (1-frac_reassign_feed_survive)), min, 1)
  # renormalize, if necessary
  dtk_total = dtk_frac_BF_survive + dtk_frac_noBF_survive + dtk_frac_noBF_die
  dtk_frac_BF_survive = dtk_frac_BF_survive / dtk_total
  dtk_frac_noBF_survive = dtk_frac_noBF_survive / dtk_total
  dtk_frac_noBF_die = dtk_frac_noBF_die / dtk_total
  
  dtk_blocking_rate = 1 - dtk_frac_BF_survive
  dtk_killing_rate = dtk_frac_noBF_die / (dtk_frac_noBF_die + dtk_frac_noBF_survive)
  
  return(list(list(dtk_blocking_rate, dtk_killing_rate), list(dtk_frac_BF_survive, dtk_frac_noBF_survive, dtk_frac_noBF_die), list(hut_frac_BF_survive, hut_frac_BF_die, hut_frac_noBF_survive, hut_frac_noBF_die)))
}

get_reverse_hut_mort_BF_from_dtk_block_kill = function(dtk_blocking_rate, dtk_killing_rate, frac_reassign_feed_survive){
  #' Given the dtk block and kill rates, reverse calculation to get the hut trial mortality and bloodfeeding values 
  #' @param dtk_blocking_rate the dtk initial blocking rate parameter
  #' @param dtk_killing_rate the dtk initial killing rate parameter
  #' @param frac_reassign_feed_survive the fraction of mosquitoes that fall in the 'feed and die' category that are reassigned to 'feed and survive'
  #' Additional details: reassign the mosquitoes that eat and then die in the hut trials to either eat and survive or die without feeding for the dtk
  #'                     we assume that hut mortality and hut blood-feeding are independent events
  
  loss = function(vars){
    x_hut_mortality = vars[1]
    x_hut_BF = vars[2]
    return(sum(c(1 - (x_hut_BF * (1-x_hut_mortality) + x_hut_BF * x_hut_mortality * frac_reassign_feed_survive) - dtk_blocking_rate,
                 ((1-x_hut_BF) * x_hut_mortality + x_hut_BF * x_hut_mortality * (1-frac_reassign_feed_survive)) / (((1-x_hut_BF) * x_hut_mortality + x_hut_BF * x_hut_mortality * (1-frac_reassign_feed_survive)) + (1-x_hut_BF) * (1-x_hut_mortality)) - dtk_killing_rate
                 )^2))
  }
  results = nlm(loss, c(0.5,0.5))$estimate
  return(c(results[1], results[2]))
}

get_dtk_outcome_from_dtk_block_kill = function(dtk_killing_rate, dtk_blocking_rate){
  #' Given the dtk parameters, return the probability of each outcome of a mosquito attempting an indoor feed
  #' @param dtk_killing_rate the net's kill_initial parameter
  #' @param dtk_blocking_rate the net's block_initial parameter

  dtk_frac_BF_survive = 1 - dtk_blocking_rate
  dtk_frac_noBF_survive = dtk_blocking_rate * (1 - dtk_killing_rate)
  dtk_frac_noBF_die = dtk_blocking_rate * dtk_killing_rate
  
  return(list(dtk_frac_BF_survive, dtk_frac_noBF_survive, dtk_frac_noBF_die))
}

get_dtk_block_kill_from_bioassay = function(bioassay_mortality, k1, k2, frac_reassign_feed_survive, fit_version='loglogistic_Nash', fit_version_BF='ITN_extraction'){
  #' Given a bioassay mortality, get the blocking and killing rates for the dtk
  #' @param bioassay_mortality the bioassay mortality value (can be a scalar or vector)
  #' @param k1 fit parameter corresponding to ITN extraction function
  #' @param k2 fit parameter corresponding to ITN extraction function
  #' @param frac_reassign_feed_survive the fraction of mosquitoes that fall in the 'feed and die' category that are reassigned to 'feed and survive'
  #' @param fit_version string corresponding to the relationship to use to calculate hut mortality from bioassay mortality
  #' @param fit_version_BF string corresponding to the relationship to use to calculate hut bloodfeeding from hut mortality
  
  hut_mortality = get_hut_mortality_from_bioassay_mortality(bioassay_mortality, fit_version=fit_version)
  hut_BF = get_hut_BF_from_hut_mortality(k1,k2,hut_mortality, fit_version_BF=fit_version_BF)
  dtk_block_kill = get_dtk_block_kill_from_hut_mort_BF(hut_mortality, hut_BF, frac_reassign_feed_survive)[[1]]
  return(dtk_block_kill)
}





#########################################################################
# add dtk parameters to simulation input files from bioassay values
#########################################################################



add_kill_block_to_df = function(itn_df, frac_reassign_feed_survive=1, filename_cur = ''){
  #' add killing and blocking parameters to a data frame that currently contains bioassay mortality
  #' @param itn_df data frame containing ITN distribution information. Must contain bioassay mortality in column called 'bio_mortality'. 
  #'               May include different seeds (in column 'seed'), which will be assigned different functions. 
  #'               May include column indicating whether or not the net type is PBO (column name 'net_type_PBO').
  #' @param frac_reassign_feed_survive the fraction of mosquitoes that fall in the 'feed and die' category that are reassigned to 'feed and survive'
  #' @param filename_cur name of the ITN file, which is used to determine whether it is part of calibration (which doesn't have different ITN functions for each seed)
  #' @return the ITN data frame with new columns for the blocking and killing rates to use for the DTK simulations
  
  itn_df$block_initial = NA
  itn_df$kill_initial = NA
  
  # check whether there are multiple seeds in the input file. if so, use different functions to calculate the kill/block for seeds
  #    - seed %% 2 = 0 --> log-logistic Nash 2021 fit for mortality
  #    - seed %% 2 = 1 -->  logistic Nash 2021 fit for mortality
  #    - seed %% 4 <= 1 --> ITN_extraction for blood-feeding
  #    - seed %% 4 > 1 --> Nash_2021 for blood-feeding
  if('seed' %in% colnames(itn_df)){
    num_seeds = length(unique(itn_df$seed))
  } else num_seeds = 0
  if (('seed' %in% colnames(itn_df)) & !(grepl('calib', filename_cur)) & (num_seeds > 0)){
    # mortality_function_names = c('loglogistic_Nash','logistic_Nash','loglogistic_Nash','logistic_Nash')
    # blocking_function_names = c('ITN_extraction','ITN_extraction','Nash_2021','Nash_2021')  # needs to be the same length as mortality_function_names
    mortality_function_names = c('logistic_Nash')
    blocking_function_names = c('Nash_2021')
    if(length(blocking_function_names) != length(mortality_function_names)) warning('PROBLEM DETECTED: length of mortality and blocking function vectors is not the same. Need to fix this for code to work properly.')
    
    rows_each_function_type = list()
    for(ff in 1:length(mortality_function_names)){
      ff_seed_mod = ff %% length(mortality_function_names)
      rows_each_function_type[[ff]] = which((itn_df$seed %% length(mortality_function_names)) == ff_seed_mod)
    }
  } else{
    mortality_function_names = c('logistic_Nash')
    blocking_function_names = c('Nash_2021')
    rows_each_function_type = list(1:nrow(itn_df))
  }
  
  # check whether any nets are PBO
  if('net_type_PBO' %in% colnames(itn_df)){
    rows_PBO = which(itn_df$net_type_PBO == 1)
    rows_standard = which(itn_df$net_type_PBO == 0)
  } else{
    rows_PBO = c()
    rows_IG2 = c()
    rows_standard = 1:nrow(itn_df)
  }
  if('llin_type' %in% colnames(itn_df)){
    rows_PBO = grep('PBO', itn_df$llin_type)
    rows_IG2 = grep('IG2', itn_df$llin_type)
    if((length(rows_PBO) + length(rows_IG2))>0) {
      rows_standard = seq(1, nrow(itn_df))[-c(rows_PBO, rows_IG2)]
    } else{
      rows_standard = seq(1, nrow(itn_df))
    }
  } else{
    rows_PBO = c()
    rows_IG2 = c()
    rows_standard = 1:nrow(itn_df)
  }
  
  
  # calculate and add kill rates and blocking rates
  for(ff in 1:length(mortality_function_names)){
    # non-PBO nets
    rows_cur = intersect(rows_each_function_type[[ff]], rows_standard)
    if(length(rows_cur)>0){
      hut_mortality = get_hut_mortality_from_bioassay_mortality(itn_df$bio_mortality[rows_cur], fit_version=mortality_function_names[ff])
      hut_BF = get_hut_BF_from_hut_mortality(k1_fit, k2_fit, hut_mortality, fit_version_BF=blocking_function_names[ff])
      params_and_fractions = get_dtk_block_kill_from_hut_mort_BF(hut_mortality, hut_BF, frac_reassign_feed_survive)
      # add dtk parameters to files
      itn_df$block_initial[rows_cur] = params_and_fractions[[1]][[1]]
      itn_df$kill_initial[rows_cur] = params_and_fractions[[1]][[2]]
    }
    
    # PBO nets
    rows_cur_PBO = intersect(rows_each_function_type[[ff]], rows_PBO)
    if(length(rows_cur_PBO)>0){
      PBO_equivalent_bioassay = get_PBO_bioassay_from_bioassay(itn_df$bio_mortality[rows_cur_PBO])
      PBO_hut_mortality = get_hut_mortality_from_bioassay_mortality(PBO_equivalent_bioassay, fit_version=mortality_function_names[ff])
      PBO_hut_BF = get_hut_BF_from_hut_mortality(k1_fit, k2_fit, PBO_hut_mortality, fit_version_BF=blocking_function_names[ff])
      params_and_fractions = get_dtk_block_kill_from_hut_mort_BF(PBO_hut_mortality, PBO_hut_BF, frac_reassign_feed_survive)
      # add dtk parameters to files
      itn_df$block_initial[rows_cur_PBO] = params_and_fractions[[1]][[1]]
      itn_df$kill_initial[rows_cur_PBO] = params_and_fractions[[1]][[2]]
    }
    
    # IG2 nets
    rows_cur_IG2 = intersect(rows_each_function_type[[ff]], rows_IG2)
    if(length(rows_cur_IG2)>0){
      if(use_new_ig2){
        # new (4/19/2023) way of doing IG2 based on expected hut BF and mortality rather than dtk block and kill parameters
        # values for pyrethroid-only and PBO nets (used to calculate IG2 estimates)
        hut_mortality = get_hut_mortality_from_bioassay_mortality(itn_df$bio_mortality[rows_cur_IG2], fit_version=mortality_function_names[ff])
        hut_BF = get_hut_BF_from_hut_mortality(k1_fit, k2_fit, hut_mortality, fit_version_BF=blocking_function_names[ff])
        PBO_equivalent_bioassay = get_PBO_bioassay_from_bioassay(itn_df$bio_mortality[rows_cur_IG2])
        PBO_hut_mortality = get_hut_mortality_from_bioassay_mortality(PBO_equivalent_bioassay, fit_version=mortality_function_names[ff])
        PBO_hut_BF = get_hut_BF_from_hut_mortality(k1_fit, k2_fit, PBO_hut_mortality, fit_version_BF=blocking_function_names[ff])
        
        # from PBO and pyrethroid hut values, estimate expected hut outcomes for IG2
        IG2_hut_mortality = (PBO_hut_mortality+3*ig2_kill_rate)/4
        IG2_hut_BF = (hut_BF + 2*PBO_hut_BF)/3
        IG2_params_and_fractions = get_dtk_block_kill_from_hut_mort_BF(IG2_hut_mortality, IG2_hut_BF, frac_reassign_feed_survive=frac_reassign_feed_survive)
        # add dtk parameters to files
        itn_df$block_initial[rows_cur_IG2] = IG2_params_and_fractions[[1]][[1]]
        itn_df$kill_initial[rows_cur_IG2] = IG2_params_and_fractions[[1]][[2]]
        
      } else{
        # add dtk parameters to files - IG2 killing is set to a constant value and blocking is based on the values for other net types
        itn_df$kill_initial[rows_cur_IG2] = ig2_kill_rate
        # get pyrethroid-only parameters to inform IG2 blocking
        hut_mortality = get_hut_mortality_from_bioassay_mortality(itn_df$bio_mortality[rows_cur_IG2], fit_version=mortality_function_names[ff])
        hut_BF = get_hut_BF_from_hut_mortality(k1_fit, k2_fit, hut_mortality, fit_version_BF=blocking_function_names[ff])
        params_and_fractions = get_dtk_block_kill_from_hut_mort_BF(hut_mortality, hut_BF, frac_reassign_feed_survive)
        itn_df$block_initial[rows_cur_IG2] = params_and_fractions[[1]][[1]]
      }
    }
  }
  return(itn_df)
}




add_kill_block_params_to_itn_files = function(hbhi_dir, frac_reassign_feed_survive=1, indoor_net_protection=0.75, create_itn_kill_block_plots=TRUE){
  #' add killing and blocking parameters to all ITN data frames in a project that currently only contain bioassay mortality
  #' @param hbhi_dir project directory
  #' @param frac_reassign_feed_survive the fraction of mosquitoes that fall in the 'feed and die' category that are reassigned to 'feed and survive'
  #' @param indoor_net_protection fraction of indoor bites occurring during net-protected hours
  #' @param create_itn_kill_block_plots flag indicating whether the series of plots showing the blocking and killing rates for the DTK and feeding-attempt outcomes

  sim_itn_start_folder = paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage_mort')
  sim_itn_output_folder = paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage_block_kill')
  input_files = list.files(sim_itn_start_folder, pattern = '.csv')
  ifelse(!dir.exists(sim_itn_output_folder), dir.create(sim_itn_output_folder), FALSE)
  
  # iterate through all ITN files that have bioassay mortality and calculate the kill and block rates
  for(ii in 1:length(input_files)){
    itn_df = read.csv(paste0(sim_itn_start_folder, '/', input_files[ii]))
    if(nrow(itn_df)>0){
      itn_df = add_kill_block_to_df(itn_df, frac_reassign_feed_survive=frac_reassign_feed_survive, filename_cur=input_files[ii])
      itn_df$indoor_net_protection = indoor_net_protection
    } else{
      itn_df$block_initial = c()
      itn_df$kill_initial = c()
      itn_df$indoor_net_protection = c()
    }
    # save output
    write.csv(itn_df, paste0(sim_itn_output_folder, '/', input_files[ii]), row.names=FALSE)
  }
}



plot_kill_block_function_of_bioassay_mort = function(hbhi_dir, frac_reassign_feed_survive=1){
  ###################################################################################
  # plot hut & DTK mortality and BF as a function of permethrin mortality (bioassay)
  ###################################################################################
  fit_version_names = c('mortLLNash', 'mortLNash', 'mort_ave_L_LLNash')
  fit_version_BF_name='BFNash'  # BFExtract, BFNash
  for(ff in 1:length(fit_version_names)){
    fit_version_name = fit_version_names[ff]

    if(fit_version_name == 'mortLLNash'){
      fit_version = 'loglogistic_Nash'
    } else if (fit_version_name == 'mortLNash'){
      fit_version = 'logistic_Nash' 
    } else if(fit_version_name == 'mort_ave_L_LLNash' ){
      fit_version = 'ave_loglog_log_Nash' 
    } else warning('Did not recognize that mortality fit name')
    
    if(fit_version_BF_name == 'BFExtract'){
      fit_version_BF = 'ITN_extraction'
    } else if (fit_version_BF_name == 'BFNash'){
      fit_version_BF = 'Nash_2021'
    } else warning('Did not recognize that BF fit name')
    
    ifelse(!dir.exists(paste0(hbhi_dir, '/ento/blood feeding and mortality')), dir.create(paste0(hbhi_dir, '/ento/blood feeding and mortality')), FALSE)
    pdf(paste0(hbhi_dir, "/ento/blood feeding and mortality/from_bioassay_to_params_X",round(frac_reassign_feed_survive*100),"_",fit_version_name, "_", fit_version_BF_name,".pdf"), width=6, height=4.5)
    bioassay_mortality_values = seq(0,1,0.01)
    
    # non-PBO values
    hut_mortality = get_hut_mortality_from_bioassay_mortality(bioassay_mortality_values, fit_version=fit_version)  # loglogistic_Nash
    hut_BF = get_hut_BF_from_hut_mortality(k1_fit,k2_fit,hut_mortality, fit_version_BF=fit_version_BF)  # ITN_extraction
    params_and_fractions = get_dtk_block_kill_from_hut_mort_BF(hut_mortality, hut_BF, frac_reassign_feed_survive)
    
    # PBO values
    PBO_equivalent_bioassay = get_PBO_bioassay_from_bioassay(bioassay_mortality_values)
    PBO_hut_mortality = get_hut_mortality_from_bioassay_mortality(PBO_equivalent_bioassay, fit_version=fit_version)
    PBO_hut_BF = get_hut_BF_from_hut_mortality(k1_fit,k2_fit,PBO_hut_mortality, fit_version_BF=fit_version_BF)
    PBO_params_and_fractions = get_dtk_block_kill_from_hut_mort_BF(PBO_hut_mortality, PBO_hut_BF, frac_reassign_feed_survive)
    
    # IG2 values
    if(use_new_ig2){
      # new (4/19/2023) way of doing IG2 based on expected hut BF and mortality rather than dtk block and kill parameters
      IG2_hut_mortality = (PBO_hut_mortality+3*ig2_kill_rate)/4 #rep(ig2_kill_rate, length(bioassay_mortality_values))
      IG2_hut_BF = (hut_BF + 2*PBO_hut_BF)/3
      IG2_params_and_fractions = get_dtk_block_kill_from_hut_mort_BF(IG2_hut_mortality, IG2_hut_BF, frac_reassign_feed_survive=frac_reassign_feed_survive)
      IG2_params = IG2_params_and_fractions[[1]]
    } else{
      #(killing is set to a constant value and blocking is based on the values for other net types)
      IG2_params_and_fractions = params_and_fractions
      IG2_params_and_fractions[[1]][[2]] = rep(ig2_kill_rate, length(IG2_params_and_fractions[[1]][[2]]))
      IG2_params = IG2_params_and_fractions[[1]]
    }
    
    
    
    
    # ===== 
    # =====  plot showing hut-trial equivalent mortality and blood-feeding (both PBO and traditional LLINs) as functions of bioassay permethrin mortality
    # ===== 
    plot(NA, xlim=c(0,1), ylim=c(0,1), xlab='permethrin mortality (bioassay)', ylab='hut-trial equivalent mortality & blood-feeding', bty='L')
    lines(bioassay_mortality_values, hut_mortality, lwd=3, col=rgb(0.7,0.5,1))
    lines(bioassay_mortality_values, hut_BF, lwd=3, col=rgb(0.1,0.5,1))
    lines(bioassay_mortality_values, PBO_hut_mortality, lty=2, lwd=3, col=rgb(0.7,0.5,1, 0.5))
    lines(bioassay_mortality_values, PBO_hut_BF, lty=2, lwd=3, col=rgb(0.1,0.5,1, 0.5))
    legend('topleft',c('mortality','mortality (PBO)', 'BF', 'BF (PBO)'), lty=c(1,2,1,2), lwd=3, bty='n', col=c(rgb(0.7,0.5,1),rgb(0.7,0.5,1,0.5), rgb(0.1,0.5,1), rgb(0.1,0.5,1,0.5)))
    
    # ===== 
    # =====  plot showing fraction of hut-trial mosquitos that have each of four outcomes as functions of bioassay permethrin mortality
    # ===== 
    # non-PBO values
    hut_outcomes = params_and_fractions[[3]]
    # plot(NA, xlim=c(0,1), ylim=c(0,1), xlab='permethrin mortality (bioassay)', ylab='fraction of mosquitoes', main = 'Expected hut-trial outcomes', bty='L')
    # lines(bioassay_mortality, hut_outcomes[[1]], lwd=3, col=rgb(0.7,0.5,1))  # hut_frac_BF_survive
    # lines(bioassay_mortality, hut_outcomes[[2]], lwd=3, col=rgb(0.1,0.5,1))  # hut_frac_BF_die
    # lines(bioassay_mortality, hut_outcomes[[3]], lwd=3, col=rgb(0.9,0.8,1))   # hut_frac_noBF_survive
    # lines(bioassay_mortality, hut_outcomes[[4]], lwd=3, col=rgb(0.6,0.8,1))  # hut_frac_noBF_die
    # legend('right',c('BF & survive','BF & die', 'no BF & survive', 'no BF & die'), lty=1, lwd=3, bty='n', col=rep(c(rgb(0.7,0.5,1), rgb(0.1,0.5,1)), each=2))
    # stacked
    plot(NA, xlim=c(0,1), ylim=c(0,1), xlab='permethrin mortality (bioassay)', ylab='fraction of mosquitoes', main = 'Expected hut-trial outcomes', bty='L')
    polygon(c(0,bioassay_mortality_values,1), c(0,hut_outcomes[[1]]+hut_outcomes[[2]]+hut_outcomes[[3]]+hut_outcomes[[4]],0), col=rgb(0.3,0.2,0.2),border=NA)  # hut_frac_noBF_die
    polygon(c(0,bioassay_mortality_values,1), c(0,hut_outcomes[[1]]+hut_outcomes[[2]]+hut_outcomes[[3]],0), col=rgb(0.4,0,0.5),border=NA)  # hut_frac_noBF_survive
    polygon(c(0,bioassay_mortality_values,1), c(0,hut_outcomes[[1]]+hut_outcomes[[2]],0), col=rgb(0.45,0.6,0.9),border=NA)  # hut_frac_BF_die
    polygon(c(0,bioassay_mortality_values,1), c(0,hut_outcomes[[1]],0), col=rgb(0.5,1,0.8),border=NA)  # hut_frac_BF_survive
    text(x=0.6,y=0.8, 'die without BF', col=rgb(1,0.9,0.9))
    text(x=0.5,y=0.5, 'survive without BF', col=rgb(0.9,0.8,1))
    text(x=0.38,y=0.27, 'die with BF', col=rgb(0,0,0.4))
    text(x=0.2,y=0.1, 'survive with BF', col=rgb(0.05,0.3,0.2))
    
    # ===== 
    # =====  plot showing fraction of dtk mosquitos that have each of three outcomes as functions of bioassay permethrin mortality
    # ===== 
    dtk_outcomes = params_and_fractions[[2]]
    # stacked
    plot(NA, xlim=c(0,1), ylim=c(0,1), xlab='permethrin mortality (bioassay)', ylab='fraction of mosquitoes', main = 'DTK indoor-feeding outcomes', bty='L')
    polygon(c(0,bioassay_mortality_values,1), c(0,dtk_outcomes[[1]]+dtk_outcomes[[2]]+dtk_outcomes[[3]],0), col=rgb(0.3,0.2,0.2),border=NA)  # frac_noBF_die
    polygon(c(0,bioassay_mortality_values,1), c(0,dtk_outcomes[[1]]+dtk_outcomes[[2]],0), col=rgb(0.4,0,0.5),border=NA)  # frac_noBF_survive
    polygon(c(0,bioassay_mortality_values,1), c(0,dtk_outcomes[[1]],0), col=rgb(0.5,1,0.8),border=NA)  # frac_BF_survive
    text(x=0.6,y=0.8, 'die without BF', col=rgb(1,0.9,0.9))
    text(x=0.4+0.1*frac_reassign_feed_survive,y=0.45, 'survive without BF', col=rgb(0.9,0.8,1))
    text(x=0.2,y=0.1, 'survive with BF', col=rgb(0.05,0.3,0.2))
    
    # ===== 
    # =====  plot showing dtk blocking and killing parameters as a function of bioassay permethrin mortality
    # ===== 
    dtk_block_kill = params_and_fractions[[1]]
    PBO_dtk_block_kill = PBO_params_and_fractions[[1]]
    IG2_dtk_block_kill = IG2_params_and_fractions[[1]]
    plot(NA, xlim=c(0,1), ylim=c(0,1), xlab='permethrin mortality (bioassay)', ylab='blocking or killing rate', main = 'DTK parameters', bty='L')
    lines(bioassay_mortality_values, dtk_block_kill[[1]], lwd=3, col=rgb(0,0.52,0.44))
    lines(bioassay_mortality_values, dtk_block_kill[[2]], lwd=3, col=rgb(0.65,0.38,0.1))
    lines(bioassay_mortality_values, PBO_dtk_block_kill[[1]], lwd=3, col=rgb(0,0.52,0.44,0.5), lty=2)
    lines(bioassay_mortality_values, PBO_dtk_block_kill[[2]], lwd=3, col=rgb(0.65,0.38,0.1,0.5), lty=2)
    lines(bioassay_mortality_values, IG2_dtk_block_kill[[1]], lwd=3, col=rgb(0,0.52,0.44,0.5), lty=3)
    lines(bioassay_mortality_values, IG2_dtk_block_kill[[2]], lwd=3, col=rgb(0.65,0.38,0.1,0.5), lty=3)
    legend('bottomright',c('blocking rate', 'blocking rate (PBO)', 'blocking rate (IG2)','killing rate','killing rate (PBO)','killing rate (IG2)'), lty=c(1,2,3,1,2,3), lwd=3, bty='n', col=c(rgb(0,0.52,0.44), rgb(0,0.52,0.44,0.5), rgb(0,0.52,0.44,0.5), rgb(0.65,0.38,0.1), rgb(0.65,0.38,0.1), rgb(0.65,0.38,0.1,0.5)))
    
    dev.off()
  }
}
 


create_feeding_outcome_plot_for_GR = function(hbhi_dir, frac_reassign_feed_survive=1, indoor_net_protection=0.75){
  fit_version_mort_names = c('mortLNash','mortLLNash', 'mort_ave_L_LLNash')
  fit_version_BF_names = c('BFNash', 'BFExtract')

  num_net_types=3
  ifelse(!dir.exists(paste0(hbhi_dir, '/ento/blood feeding and mortality')), dir.create(paste0(hbhi_dir, '/ento/blood feeding and mortality')), FALSE)
  # png(paste0(hbhi_dir, "/ento/blood feeding and mortality/from_bioassay_survival_to_outcome_X",round(frac_reassign_feed_survive*100),"_noLabels.png"), width=4*0.79*num_net_types, height=0.82*3*length(fit_version_mort_names) *length(fit_version_BF_names), unit='in', res=1800)
  png(paste0(hbhi_dir, "/ento/blood feeding and mortality/from_bioassay_survival_to_outcome_X",round(frac_reassign_feed_survive*100),".png"), width=4*0.79*num_net_types, height=0.82*3*length(fit_version_mort_names) *length(fit_version_BF_names), unit='in', res=1800)
  par(mfrow=c(length(fit_version_mort_names) *length(fit_version_BF_names),num_net_types), mgp=c(2.3,1,0))
  
  
  for(mm in 1:length(fit_version_mort_names)){
    for(bb in 1:length(fit_version_BF_names)){
      fit_version_name = fit_version_mort_names[mm]
      fit_version_BF_name = fit_version_BF_names[bb]
      
      if(fit_version_name == 'mortLLNash'){
        fit_version = 'loglogistic_Nash' 
      } else if (fit_version_name == 'mortLNash'){
        fit_version = 'logistic_Nash' 
      } else if(fit_version_name == 'mort_ave_L_LLNash' ){
        fit_version = 'ave_loglog_log_Nash' 
      } else warning('Did not recognize that mortality fit name')
      
      if(fit_version_BF_name == 'BFExtract'){
        fit_version_BF = 'ITN_extraction' 
      } else if (fit_version_BF_name == 'BFNash'){
        fit_version_BF = 'Nash_2021' 
      } else warning('Did not recognize that BF fit name')
      

      
      bioassay_mortality_values = seq(0,1,0.01)
      
      # non-PBO values
      hut_mortality = get_hut_mortality_from_bioassay_mortality(bioassay_mortality_values, fit_version=fit_version)  # loglogistic_Nash
      hut_BF = get_hut_BF_from_hut_mortality(k1_fit,k2_fit,hut_mortality, fit_version_BF=fit_version_BF)  # ITN_extraction
      params_and_fractions = get_dtk_block_kill_from_hut_mort_BF(hut_mortality, hut_BF, frac_reassign_feed_survive)
      pyr_params = params_and_fractions[[1]]

      # PBO values
      PBO_equivalent_bioassay = get_PBO_bioassay_from_bioassay(bioassay_mortality_values)
      PBO_hut_mortality = get_hut_mortality_from_bioassay_mortality(PBO_equivalent_bioassay, fit_version=fit_version)
      PBO_hut_BF = get_hut_BF_from_hut_mortality(k1_fit,k2_fit,PBO_hut_mortality, fit_version_BF=fit_version_BF)
      PBO_params_and_fractions = get_dtk_block_kill_from_hut_mort_BF(PBO_hut_mortality, PBO_hut_BF, frac_reassign_feed_survive)
      PBO_params = PBO_params_and_fractions[[1]]
      
      # IG2 values
      if(use_new_ig2){
        # new (4/19/2023) way of doing IG2 based on expected hut BF and mortality rather than dtk block and kill parameters
        IG2_hut_mortality = (PBO_hut_mortality+3*ig2_kill_rate)/4 #rep(ig2_kill_rate, length(bioassay_mortality_values))
        IG2_hut_BF = (hut_BF + 2*PBO_hut_BF)/3
        IG2_params_and_fractions = get_dtk_block_kill_from_hut_mort_BF(IG2_hut_mortality, IG2_hut_BF, frac_reassign_feed_survive=frac_reassign_feed_survive)
        IG2_params = IG2_params_and_fractions[[1]]
      } else{
        #(killing is set to a constant value and blocking is based on the values for other net types)
        IG2_params_and_fractions = params_and_fractions
        IG2_params_and_fractions[[1]][[2]] = rep(ig2_kill_rate, length(IG2_params_and_fractions[[1]][[2]]))
        IG2_params = IG2_params_and_fractions[[1]]
      }

      
      
      all_net_type_params = list(pyr_params, PBO_params, IG2_params)
      net_type_names = c('Pyr', 'PBO', 'IG2')
      for(ii in 1:length(all_net_type_params)){
        cur_net_params = all_net_type_params[[ii]]
        dtk_outcomes = get_dtk_outcome_from_dtk_block_kill(dtk_killing_rate=cur_net_params[[2]], dtk_blocking_rate=cur_net_params[[1]])
        # =====  plot showing fraction of dtk mosquitos that have each of three outcomes as functions of bioassay permethrin mortality
        # stacked
        plot(NA, xlim=c(0,100), ylim=c(0,1), xlab='bioassay survival (%)', ylab='fraction of mosquitoes', main = paste0(net_type_names[ii], ' - ', fit_version_name,'_',fit_version_BF_name), bty='L')
        polygon((1-c(0,bioassay_mortality_values,1))*100, c(0,dtk_outcomes[[1]]+dtk_outcomes[[2]]+dtk_outcomes[[3]],0), col=rgb(0.3,0.2,0.2),border=NA)  # frac_noBF_die
        polygon((1-c(0,bioassay_mortality_values,1))*100, c(0,dtk_outcomes[[1]]+dtk_outcomes[[2]],0), col=rgb(0.2,0.3,0.6),border=NA)  # frac_noBF_survive
        polygon((1-c(0,bioassay_mortality_values,1))*100, c(0,dtk_outcomes[[1]],0), col=rgb(0.6,0.9,1),border=NA)  # frac_BF_survive
        text(x=25,y=0.8, 'die without feeding', col=rgb(1,0.9,0.9))
        text(x=45,y=0.44, 'survive without feeding', col=rgb(0.7,0.8,1))
        text(x=60,y=0.1, 'survive with feeding', col=rgb(0.05,0.2,0.3))
        
      }
      
      # # =====  plot showing fraction of dtk mosquitos that have each of three outcomes as functions of bioassay permethrin mortality
      # dtk_outcomes = params_and_fractions[[2]]
      # # stacked
      # plot(NA, xlim=c(0,100), ylim=c(0,1), xlab='bioassay survival (%)', ylab='fraction of mosquitoes', main = paste0(fit_version_name,'_',fit_version_BF_name), bty='L')
      # polygon((1-c(0,bioassay_mortality_values,1))*100, c(0,dtk_outcomes[[1]]+dtk_outcomes[[2]]+dtk_outcomes[[3]],0), col=rgb(0.3,0.2,0.2),border=NA)  # frac_noBF_die
      # polygon((1-c(0,bioassay_mortality_values,1))*100, c(0,dtk_outcomes[[1]]+dtk_outcomes[[2]],0), col=rgb(0.2,0.3,0.6),border=NA)  # frac_noBF_survive
      # polygon((1-c(0,bioassay_mortality_values,1))*100, c(0,dtk_outcomes[[1]],0), col=rgb(0.6,0.9,1),border=NA)  # frac_BF_survive
      # text(x=25,y=0.8, 'die without feeding', col=rgb(1,0.9,0.9))
      # text(x=45,y=0.44, 'survive without feeding', col=rgb(0.7,0.8,1))
      # text(x=60,y=0.1, 'survive with feeding', col=rgb(0.05,0.2,0.3))
    }
  }
  dev.off()
}






library(ggplot2)
library(reshape2)
library(dplyr)

data_dir = 'C:/Users/moniqueam/OneDrive - Bill & Melinda Gates Foundation/projects/HBHI_general/ITN_parameters'
script_dir = 'C:/Users/moniqueam/Documents/malaria-nga-snt22'
source(paste0(script_dir,'/data_processing/setup_inputs/add_hut_trial_mort_block.R'))

ht_results = read.csv(paste0(data_dir, '/Nash et al_2021/Nash_supp3_hut_outcomes.csv'))



####################################################################################
############## plot comparison of different net types from Nash et al. 2021 #################
####################################################################################
ht_results$p_24h_mort = ht_results$N_24h_mort / ht_results$N_total
ht_results$p_72h_mort = ht_results$N_72h_mort / ht_results$N_total


# compare 72 hour and 24 hour mortality when both are available
# hour_mort = reshape(data=ht_results[,c('Trial.Code', 'Insecticide_class', 'N_total', 'N_24h_mort', 'N_72h_mort')], direction='long', varying=c('N_24h_mort', 'N_72h_mort'), v.names='h_mort')# ,idvar=c('Trial.Code', 'Insecticide_class', 'N_total'))
ht_results_bothTimes = ht_results[!is.na(ht_results$N_24h_mort) & !is.na(ht_results$N_72h_mort),]
ggplot(ht_results_bothTimes, aes(x=p_24h_mort, y=p_72h_mort, color=Insecticide_class))+
  geom_point() + 
  geom_abline(slope=1, intercept=0)


# look at mortality for control huts (untreated net)
ggplot(ht_results[ht_results$Intervention=='Untreated',], aes(x=p_24h_mort))+
  geom_histogram(bins=20, fill='grey', color='black')
  # geom_density() 



# look at mortality for different types of nets

# get mean mortality for each insecticide class, weighted by number of mosquitoes in huts
ht_results_mean_mort = ht_results %>% group_by(Insecticide_class) %>%
  summarise(N_total = sum(N_total, na.rm=TRUE),
            N_24h_mort = sum(N_24h_mort, na.rm=TRUE)) %>% 
  ungroup()
ht_results_mean_mort$mort_rate = ht_results_mean_mort$N_24h_mort / ht_results_mean_mort$N_total

# create density plot with means
insecticide_classes_included = c('pyrethroid', 'pyrethroid + PBO',  'pyrethroid + pyrrole', 'untreated')
ggplot()+
  geom_density(data = ht_results[ht_results$Insecticide_class %in% insecticide_classes_included,], aes(x=p_24h_mort, color=Insecticide_class, fill=Insecticide_class), alpha=0.2)+
  geom_point(data=ht_results_mean_mort[ht_results_mean_mort$Insecticide_class %in% insecticide_classes_included,], aes(x=mort_rate, y=0.2, color=Insecticide_class), shape='|', size=10)
# create violin plot with dots
ggplot(data = ht_results[ht_results$Insecticide_class %in% insecticide_classes_included,], aes(y=p_24h_mort, x=Insecticide_class, color=Insecticide_class, fill=Insecticide_class))+
  geom_violin(alpha=0.2)+
  geom_jitter(aes(size=N_total/100), position=position_jitter(0.2), alpha=0.2)
  


# subset to sites where IG2 nets were tested
ig2_trial_subset = unique(ht_results$Trial.Code[which(ht_results$Insecticide_class == 'pyrethroid + pyrrole')])
ht_results_mean_mort_ig2 = ht_results[(ht_results$Insecticide_class %in% insecticide_classes_included) & (ht_results$Trial.Code %in% ig2_trial_subset),] %>% group_by(Insecticide_class) %>%
  summarise(N_total = sum(N_total, na.rm=TRUE),
            N_24h_mort = sum(N_24h_mort, na.rm=TRUE)) %>% 
  ungroup()
ht_results_mean_mort_ig2$mort_rate = ht_results_mean_mort_ig2$N_24h_mort / ht_results_mean_mort_ig2$N_total

# create violin plot with dots for trials with IG2 nets
ggplot(data = ht_results[(ht_results$Insecticide_class %in% insecticide_classes_included) & (ht_results$Trial.Code %in% ig2_trial_subset),], aes(y=p_24h_mort, x=Insecticide_class, color=Insecticide_class, fill=Insecticide_class))+
  geom_violin(alpha=0.2)+
  geom_jitter(aes(size=N_total/100), position=position_jitter(0.2), alpha=0.2)




# subset to sites where PBO nets were tested
pbo_trial_subset = unique(ht_results$Trial.Code[which(ht_results$Insecticide_class == 'pyrethroid + PBO')])
ht_results_mean_mort_pbo = ht_results[(ht_results$Insecticide_class %in% insecticide_classes_included) & (ht_results$Trial.Code %in% pbo_trial_subset),] %>% group_by(Insecticide_class) %>%
  summarise(N_total = sum(N_total, na.rm=TRUE),
            N_24h_mort = sum(N_24h_mort, na.rm=TRUE)) %>% 
  ungroup()
ht_results_mean_mort_pbo$mort_rate = ht_results_mean_mort_pbo$N_24h_mort / ht_results_mean_mort_pbo$N_total

# create violin plot with dots for trials with IG2 nets
ggplot(data = ht_results[(ht_results$Insecticide_class %in% insecticide_classes_included) & (ht_results$Trial.Code %in% pbo_trial_subset),], aes(y=p_24h_mort, x=Insecticide_class, color=Insecticide_class, fill=Insecticide_class))+
  geom_violin(alpha=0.2)+
  geom_jitter(aes(size=N_total/100), position=position_jitter(0.2), alpha=0.2)










############################################################################################################################
############## compare different net types used in same trials for hut  mortality and hut bloodfeeding  #################
############################################################################################################################

####===============================================================================###
####========                  reference datasets                     ==============###
####===============================================================================###
###========== add studies from Nash et al. 2021 ==========###

# get studies that have ig2 or pbo net types
trial_codes_cur = unique(ht_results$Trial.Code[grepl('PBO',ht_results$Insecticide_class) | grepl('pyrrole',ht_results$Insecticide_class)])
ht_results_comp = ht_results[ht_results$Trial.Code %in% trial_codes_cur,]
# subset to unwashed nets
ht_results_comp = ht_results_comp[ht_results_comp$is_washed==0,]

results_df = data.frame()
# get matched values for pyrethroid+ and pyrethroid-only versions of net
for(tt in 1:length(unique(ht_results_comp$Trial.Code))){
  cur_trial = ht_results_comp[ht_results_comp$Trial.Code==unique(ht_results_comp$Trial.Code)[tt],]
  # which pyrethroid+ type(s)
  plus_rows = cur_trial[grepl('pyrethroid +', cur_trial$Insecticide_class),]
  if(nrow(plus_rows)>0){
    
    # get the matching pyrethroid-only value for each row
    for(rr in 1:nrow(plus_rows)){
      cur_plus_row = plus_rows[rr,]
      
      # get the control mortality for this study (without treated net), if available
      control_row = cur_trial[(cur_trial$Insecticide=='untreated') & (cur_trial$Mosquito_species == cur_plus_row$Mosquito_species),]
      if(nrow(control_row)==1){
        control_mort = control_row$N_24h_mort / control_row$N_total
      } else {
        warning(paste0('PROBLEM: control study not found for trial ', unique(ht_results_comp$Trial.Code)[tt]))
        control_mort=0
      }
      
      # insecticide without the plus
      cur_insecticide = gsub(' +.*', '',cur_plus_row$Insecticide)
      match_pyr_row = cur_trial[(cur_trial$Insecticide == cur_insecticide) & (cur_trial$Mosquito_species == cur_plus_row$Mosquito_species),]
      if(nrow(match_pyr_row)!=1){
        # sometimes there are multiple net brands (often sharing the first few letters of name across net types)
        cur_brand = substr(x=cur_plus_row$Intervention, start=1,stop=4)
        # subset to the pyr rows matching the brand of the plus net
        match_pyr_row = match_pyr_row[grepl(cur_brand, match_pyr_row$Intervention),]
      }      
      if(nrow(match_pyr_row)!=1){
        warning(paste0('PROBLEM: pyrethroid-only row with same insecticide not found for ', unique(ht_results_comp$Trial.Code)[tt]))
      } else{
        if(grepl('PBO', cur_plus_row$Insecticide_class)){
          plus_type = 'PBO'
        } else if(grepl('pyrrole', cur_plus_row$Insecticide_class)){
          plus_type = 'IG2'
        } else{
          warning('PROBLEM: insecticide type not recognized')
          plus_type=NA
        }
        df_cur = data.frame('study' = unique(ht_results_comp$Trial.Code)[tt],
                            'pyr_mort' = match_pyr_row$N_24h_mort / match_pyr_row$N_total - control_mort,
                            'plus_mort' =  cur_plus_row$N_24h_mort / cur_plus_row$N_total - control_mort,
                            'pyr_fed' = match_pyr_row$pc_fed/100,
                            'plus_fed' =  cur_plus_row$pc_fed/100,
                            'plus_type' = plus_type)
        if(nrow(results_df)>0){
          results_df = rbind(results_df, df_cur)
        }
        else{
          results_df = df_cur
        }
      }
    }
  }
}


###========== add studies (mostly) not included in Nash et al. 2021 ==========##

###### PBO versus pyr ########

## Akoton et al. 2018
# datasets for each net type 
pyr_mort = c(.28,.11,.77,.49)
pbo_mort = c(.74,.73,1,.89)
pyr_block = 1-c(.18,.24,.22,.29)
pbo_block = 1-c(.1,.18,.04,0)
study = rep('Akoton_2018', length(pyr_mort))
df_cur = data.frame('study' = study,
                    'pyr_mort' = pyr_mort,
                    'plus_mort' =  pbo_mort,
                    'pyr_fed' = 1-pyr_block,
                    'plus_fed' =  1-pbo_block,
                    'plus_type' = 'PBO')
results_df = rbind(results_df, df_cur)


# ## Corbel et al. 2010 - already in Nash et al
# # datasets for each net type 
# control_net_mort1 = 0.04
# control_net_mort2 = 0.13
# control_net_mort3 = 0.05
# pyr_mort = c(.89-control_net_mort1, 0.83-control_net_mort2,.44-control_net_mort3)
# pbo_mort = c(0.97-control_net_mort1,0.94-control_net_mort2,0.78-control_net_mort3)
# pyr_block = 1-c(.04,.15,.35)
# pbo_block = 1-c(0,.28,21)
# pyr_mort_all = c(pyr_mort_all, pyr_mort)
# pbo_mort_all = c(pbo_mort_all, pbo_mort)
# pyr_block_all = c(pyr_block_all, pyr_block)
# pbo_block_all = c(pbo_block_all, pbo_block)
# study = c(study, rep('Corbel_2010', length(pyr_mort)))



## Bayili et al. 2019
# datasets for each net type 
control_net_mort = 0.11
pyr_mort = c(.15,.15)-control_net_mort
pbo_mort = c(.22,.43)-control_net_mort
pyr_block = 1-c(.6,.6)
pbo_block = 1-c(.43,.34)
study = rep('Bayili_2019', length(pyr_mort))
df_cur = data.frame('study' = study,
                    'pyr_mort' = pyr_mort,
                    'plus_mort' =  pbo_mort,
                    'pyr_fed' = 1-pyr_block,
                    'plus_fed' =  1-pbo_block,
                    'plus_type' = 'PBO')
results_df = rbind(results_df, df_cur)

## Menze et al. 2020
# datasets for each net type 
control_net_mort1 = 0.05
control_net_mort2 = 0.37
pyr_mort = c(c(.12,.1)-control_net_mort1, c(.51,.4)-control_net_mort2)
pbo_mort = c(c(.3,.25)-control_net_mort1, c(.66,.68)-control_net_mort2)
pyr_block = 1-c(.2,.15,.16,.16)
pbo_block = 1-c(.16,.17,.09,.2)
study = rep('Menze_2020', length(pyr_mort))
df_cur = data.frame('study' = study,
                    'pyr_mort' = pyr_mort,
                    'plus_mort' =  pbo_mort,
                    'pyr_fed' = 1-pyr_block,
                    'plus_fed' =  1-pbo_block,
                    'plus_type' = 'PBO')
results_df = rbind(results_df, df_cur)

## Syme et al. 2023
# datasets for each net type 
control_net_mort = 0.02
pyr_mort = c(.23)-control_net_mort
pbo_mort = c(.56)-control_net_mort
pyr_block = 1-c(.39)
pbo_block = 1-c(.13)
study = rep('Syme_2023', length(pyr_mort))
df_cur = data.frame('study' = study,
                    'pyr_mort' = pyr_mort,
                    'plus_mort' =  pbo_mort,
                    'pyr_fed' = 1-pyr_block,
                    'plus_fed' =  1-pbo_block,
                    'plus_type' = 'PBO')
results_df = rbind(results_df, df_cur)




###### IG2 vs pyr ########


## Bayili et al. 2017
# datasets for each net type 
control_net_mort = 0.13
pyr_mort = c(.27)-control_net_mort
ig2_mort = c(.81)-control_net_mort
pyr_block = 1-c(.52)
ig2_block = 1-c(.38)
study = rep('Bayili_2017', length(pyr_mort))
df_cur = data.frame('study' = study,
                    'pyr_mort' = pyr_mort,
                    'plus_mort' =  ig2_mort,
                    'pyr_fed' = 1-pyr_block,
                    'plus_fed' =  1-ig2_block,
                    'plus_type' = 'IG2')
results_df = rbind(results_df, df_cur)

## Syme et al. 2023
# datasets for each net type 
control_net_mort = 0.02
pyr_mort = c(.23)-control_net_mort
ig2_mort = c(.84)-control_net_mort
pyr_block = 1-c(.39)
ig2_block = 1-c(.2)
study = rep('Syme_2023', length(pyr_mort))
df_cur = data.frame('study' = study,
                    'pyr_mort' = pyr_mort,
                    'plus_mort' =  ig2_mort,
                    'pyr_fed' = 1-pyr_block,
                    'plus_fed' =  1-ig2_block,
                    'plus_type' = 'IG2')
results_df = rbind(results_df, df_cur)

# ## Camara et al. 2018 - already in Nash et al
# # datasets for each net type 
# control_net_mort = 0
# pyr_mort = c(.1)-control_net_mort
# ig2_mort = c(.87)-control_net_mort
# pyr_block = 1-c(.54)
# ig2_block = 1-c(.36)
# pyr_mort_all = c(pyr_mort_all, pyr_mort)
# ig2_mort_all = c(ig2_mort_all, ig2_mort)
# pyr_block_all = c(pyr_block_all, pyr_block)
# ig2_block_all = c(ig2_block_all, ig2_block)
# study = c(study, rep('Camara_2018', length(pyr_mort)))




###### IG2 vs PBO ########

pbo_mort_all = c()
pbo_block_all = c()
ig2_mort_all = c()
ig2_block_all = c()
study=c()
## Syme et al. 2023
# datasets for each net type
control_net_mort = 0.02
pbo_mort = c(.56)-control_net_mort
ig2_mort = c(.84)-control_net_mort
pbo_block = 1-c(.13)
ig2_block = 1-c(.2)
pbo_mort_all = c(pbo_mort_all, pbo_mort)
ig2_mort_all = c(ig2_mort_all, ig2_mort)
pbo_block_all = c(pbo_block_all, pbo_block)
ig2_block_all = c(ig2_block_all, ig2_block)
study = c(study, rep('Syme_2023', length(pbo_mort)))





####===============================================================================###
####========                  modeled relationships                  ==============###
####===============================================================================###

# fixed parameters for calculations - obtained from fit_blood_feeding_params.R
k1_fit = 0.27
k2_fit = 1.3
frac_reassign_feed_survive = 1
fit_version_mort_names = c('mortLNash','mortLLNash', 'mort_ave_L_LLNash')
fit_version_BF_names = c('BFNash', 'BFExtract')
ig2_version='stronger'  # 'stronger' , 'original'
if(ig2_version=='stronger'){
  ig2_kill_rate = 0.75 # 0.9, 0.75
} else{
  ig2_kill_rate = 0.75
}
use_new_ig2 = TRUE
num_net_types=3
fit_version_name = 'mortLNash'
fit_version_BF_name =  'BFNash'

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

# pyrethroid-only values
hut_mortality = get_hut_mortality_from_bioassay_mortality(bioassay_mortality_values, fit_version=fit_version)  # loglogistic_Nash
hut_BF = get_hut_BF_from_hut_mortality(k1_fit,k2_fit,hut_mortality, fit_version_BF=fit_version_BF)  # ITN_extraction
params_and_fractions1 = get_dtk_block_kill_from_hut_mort_BF(hut_mortality, hut_BF, frac_reassign_feed_survive=1)
params_and_fractions0 = get_dtk_block_kill_from_hut_mort_BF(hut_mortality, hut_BF, frac_reassign_feed_survive=0)
pyr_params1 = params_and_fractions1[[1]]
pyr_params0 = params_and_fractions0[[1]]

# PBO values
PBO_equivalent_bioassay = get_PBO_bioassay_from_bioassay(bioassay_mortality_values)
PBO_hut_mortality = get_hut_mortality_from_bioassay_mortality(PBO_equivalent_bioassay, fit_version=fit_version)
PBO_hut_BF = get_hut_BF_from_hut_mortality(k1_fit,k2_fit,PBO_hut_mortality, fit_version_BF=fit_version_BF)
PBO_params_and_fractions1 = get_dtk_block_kill_from_hut_mort_BF(PBO_hut_mortality, PBO_hut_BF, frac_reassign_feed_survive=1)
PBO_params_and_fractions0 = get_dtk_block_kill_from_hut_mort_BF(PBO_hut_mortality, PBO_hut_BF, frac_reassign_feed_survive=0)
PBO_params1 = PBO_params_and_fractions1[[1]]
PBO_params0 = PBO_params_and_fractions0[[1]]

# IG2 values
if(use_new_ig2){
  # new (4/19/2023) way of doing IG2 based on expected hut BF and mortality rather than dtk block and kill parameters
  IG2_hut_mortality = (PBO_hut_mortality+3*ig2_kill_rate)/4 #rep(ig2_kill_rate, length(bioassay_mortality_values))
  IG2_hut_BF = (hut_BF + 2*PBO_hut_BF)/3
  IG2_params_and_fractions1 = get_dtk_block_kill_from_hut_mort_BF(IG2_hut_mortality, IG2_hut_BF, frac_reassign_feed_survive=1)
  IG2_params_and_fractions0 = get_dtk_block_kill_from_hut_mort_BF(IG2_hut_mortality, IG2_hut_BF, frac_reassign_feed_survive=0)
  IG2_params1 = IG2_params_and_fractions1[[1]]
  IG2_params0 = IG2_params_and_fractions0[[1]]
  
} else{
# version 1 of specifying IG2 values, based on the dtk block and kill instead of hut equivalents 
  # IG2 values (blocking is set to the same as pyrethroid-only and the kill rate constant regardless of resistance level)
  IG2_params_and_fractions1 = params_and_fractions1
  IG2_params_and_fractions0 = params_and_fractions0
  if(ig2_version=='stronger'){
    IG2_params_and_fractions1[[1]][[1]] = (params_and_fractions1[[1]][[1]] + 2*PBO_params_and_fractions1[[1]][[1]]) / 3
    IG2_params_and_fractions0[[1]][[1]] = (params_and_fractions0[[1]][[1]] + 2*PBO_params_and_fractions0[[1]][[1]]) / 3
  } 
  IG2_params_and_fractions1[[1]][[2]] = rep(ig2_kill_rate, length(IG2_params_and_fractions1[[1]][[2]]))
  IG2_params_and_fractions0[[1]][[2]] = rep(ig2_kill_rate, length(IG2_params_and_fractions0[[1]][[2]]))
  IG2_params1 = IG2_params_and_fractions1[[1]]
  IG2_params0 = IG2_params_and_fractions0[[1]]
  # backtrack to hut_mortality and hut_BF
  IG2_hut_mortality = rep(NA,length(bioassay_mortality_values))
  IG2_hut_BF = rep(NA,length(bioassay_mortality_values))
  for(ii in 1:length(bioassay_mortality_values)){
    hut_results = get_reverse_hut_mort_BF_from_dtk_block_kill(dtk_blocking_rate=IG2_params_and_fractions1[[1]][[1]][ii], dtk_killing_rate=IG2_params_and_fractions1[[1]][[2]][ii], frac_reassign_feed_survive=frac_reassign_feed_survive)
    IG2_hut_mortality[ii] = hut_results[1]
    IG2_hut_BF[ii] = hut_results[2]
  }
  # get_reverse_hut_mort_BF_from_dtk_block_kill(dtk_blocking_rate=IG2_params_and_fractions0[[1]][[1]][ii], dtk_killing_rate=IG2_params_and_fractions0[[1]][[2]][ii], frac_reassign_feed_survive=0)
  # get_reverse_hut_mort_BF_from_dtk_block_kill(dtk_blocking_rate=IG2_params_and_fractions1[[1]][[1]][ii], dtk_killing_rate=IG2_params_and_fractions1[[1]][[2]][ii], frac_reassign_feed_survive=1)
}



##### get hut-trial equivalent values assumed in the model for each net type across bioassay permethrin mortality to allow net-to-net comparisons
dtk_pyr_hut = list(hut_mortality, hut_BF)
dtk_pbo_hut = list(PBO_hut_mortality, PBO_hut_BF)
dtk_ig2_hut = list(IG2_hut_mortality, IG2_hut_BF)


# ##### get outcomes from attempted feed with different net types (for same bioassay permethrin mortality to allow net-to-net comparisons)
# dtk_pyr_outcomes1 = get_dtk_outcome_from_dtk_block_kill(dtk_killing_rate=pyr_params1[[2]], dtk_blocking_rate=pyr_params1[[1]])
# dtk_pyr_outcomes0 = get_dtk_outcome_from_dtk_block_kill(dtk_killing_rate=pyr_params0[[2]], dtk_blocking_rate=pyr_params0[[1]])
# dtk_pbo_outcomes1 = get_dtk_outcome_from_dtk_block_kill(dtk_killing_rate=PBO_params1[[2]], dtk_blocking_rate=PBO_params1[[1]])
# dtk_pbo_outcomes0 = get_dtk_outcome_from_dtk_block_kill(dtk_killing_rate=PBO_params0[[2]], dtk_blocking_rate=PBO_params0[[1]])
# dtk_ig2_outcomes1 = get_dtk_outcome_from_dtk_block_kill(dtk_killing_rate=IG2_params1[[2]], dtk_blocking_rate=IG2_params1[[1]])
# dtk_ig2_outcomes0 = get_dtk_outcome_from_dtk_block_kill(dtk_killing_rate=IG2_params0[[2]], dtk_blocking_rate=IG2_params0[[1]])






####===============================================================================###
###========== create plots ==========###
####===============================================================================###

#### plot filled rectangle of outcomes
png(filename=paste0(data_dir, '/outcome_of_attempted_feed_', fit_version_name, '_', fit_version_BF_name, '_', ig2_version,'IG2v2.png'), width=4.9*0.85, height=11*0.85, units='in', res=900)
par(mfrow=c(3,1))
all_net_type_params = list(pyr_params0, PBO_params0, IG2_params0)
# all_net_type_params = list(pyr_params1, PBO_params1, IG2_params1)
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
dev.off()




### plot expected hut outcomes
par(mfrow=c(1,2))
# hut mortality by bioassay
plot(NA,ylim=c(0,1), xlim=c(0,1), xlab='bioassay mortality', ylab='expected hut mortality')
lines(bioassay_mortality_values, hut_mortality, col='black')
lines(bioassay_mortality_values, PBO_hut_mortality, col='green')
lines(bioassay_mortality_values, IG2_hut_mortality, col='purple')
# hut BF by bioassay
plot(NA,ylim=c(0,1), xlim=c(0,1), xlab='bioassay mortality', ylab='expected hut BF')
lines(bioassay_mortality_values, hut_BF, col='black')
lines(bioassay_mortality_values, PBO_hut_BF, col='green')
lines(bioassay_mortality_values, IG2_hut_BF, col='purple')
par(mfrow=c(1,1))



### plot dtk parameters
par(mfrow=c(2,2))
# hut mortality by bioassay - X=0
plot(NA,ylim=c(0,1), xlim=c(0,1), xlab='bioassay mortality', ylab='block param')
lines(bioassay_mortality_values, pyr_params0[[1]], col='black')
lines(bioassay_mortality_values, PBO_params0[[1]], col='green')
lines(bioassay_mortality_values, IG2_params0[[1]], col='purple')
# hut BF by bioassay - X=0
plot(NA,ylim=c(0,1), xlim=c(0,1), xlab='bioassay mortality', ylab='kill param')
lines(bioassay_mortality_values, pyr_params0[[2]], col='black')
lines(bioassay_mortality_values, PBO_params0[[2]], col='green')
lines(bioassay_mortality_values, IG2_params0[[2]], col='purple')
# hut mortality by bioassay - X=1
plot(NA,ylim=c(0,1), xlim=c(0,1), xlab='bioassay mortality', ylab='block param')
lines(bioassay_mortality_values, pyr_params1[[1]], col='black')
lines(bioassay_mortality_values, PBO_params1[[1]], col='green')
lines(bioassay_mortality_values, IG2_params1[[1]], col='purple')
# hut BF by bioassay - X=1
plot(NA,ylim=c(0,1), xlim=c(0,1), xlab='bioassay mortality', ylab='kill param')
lines(bioassay_mortality_values, pyr_params1[[2]], col='black')
lines(bioassay_mortality_values, PBO_params1[[2]], col='green')
lines(bioassay_mortality_values, IG2_params1[[2]], col='purple')
par(mfrow=c(1,1))





###### plot dtk expected hut equivalent against hut trial observation
png(filename=paste0(data_dir, '/compare_net_hut_results_', fit_version_name, '_', fit_version_BF_name, '_', ig2_version,'IG2v2.png'), width=9, height=7, units='in', res=900)

par(mfcol=c(2,3), mgp=c(2.3,1,0))
## combined pyr versus PBO across studies
# mortality
plot(dtk_pyr_hut[[1]], dtk_pbo_hut[[1]], type='l', xlim=c(0,1), ylim=c(0,1), ylab='pbo', xlab='pyr', main="mortality")
abline(a=0,b=1, col=rgb(0.9,0.9,0.9,0.9), lwd=2)
points(results_df$pyr_mort[results_df$plus_type=='PBO'], results_df$plus_mort[results_df$plus_type=='PBO'], col=as.factor(results_df$study[results_df$plus_type=='PBO']), cex=1.5, pch=20)
# bloodfed
plot(dtk_pyr_hut[[2]], dtk_pbo_hut[[2]], type='l', xlim=c(0,1), ylim=c(0,1), ylab='pbo', xlab='pyr', main="bloodfed")
abline(a=0,b=1, col=rgb(0.9,0.9,0.9,0.9), lwd=2)
points(results_df$pyr_fed[results_df$plus_type=='PBO'], results_df$plus_fed[results_df$plus_type=='PBO'], col=as.factor(results_df$study[results_df$plus_type=='PBO']), cex=1.5, pch=20)

## combined pyr versus IG2 across studies
# mortality
plot(dtk_pyr_hut[[1]], dtk_ig2_hut[[1]], type='l', xlim=c(0,1), ylim=c(0,1), ylab='ig2', xlab='pyr', main="mortality")
abline(a=0,b=1, col=rgb(0.9,0.9,0.9,0.9), lwd=2)
points(results_df$pyr_mort[results_df$plus_type=='IG2'], results_df$plus_mort[results_df$plus_type=='IG2'], col=as.factor(results_df$study[results_df$plus_type=='IG2']), cex=1.5, pch=20)
# bloodfed
plot(dtk_pyr_hut[[2]], dtk_ig2_hut[[2]], type='l', xlim=c(0,1), ylim=c(0,1), ylab='ig2', xlab='pyr', main="bloodfed")
abline(a=0,b=1, col=rgb(0.9,0.9,0.9,0.9), lwd=2)
points(results_df$pyr_fed[results_df$plus_type=='IG2'], results_df$plus_fed[results_df$plus_type=='IG2'], col=as.factor(results_df$study[results_df$plus_type=='IG2']), cex=1.5, pch=20)

## combined PBO versus IG2 across studies
# mortality
plot(dtk_pbo_hut[[1]], dtk_ig2_hut[[1]], type='l', xlim=c(0,1), ylim=c(0,1), ylab='ig2', xlab='pbo', main="mortality")
abline(a=0,b=1, col=rgb(0.9,0.9,0.9,0.9), lwd=2)
points(pbo_mort_all, ig2_mort_all, col=as.factor(study), cex=1.5, pch=20)
# bloodfed
plot(dtk_pbo_hut[[2]], dtk_ig2_hut[[2]], type='l', xlim=c(0,1), ylim=c(0,1), ylab='ig2', xlab='pbo', main="bloodfed")
abline(a=0,b=1, col=rgb(0.9,0.9,0.9,0.9), lwd=2)
points(1-pbo_block_all, 1-ig2_block_all, col=as.factor(study), cex=1.5, pch=20)

dev.off()





# ###### plot dtk expected outcome against hut trial observation
# arbitrary_mortality_ref_to_emod_rescaler = 1  # Given the EMOD constraints, it may be appropriate for EMOD mortality outcomes to be lower than those in the reference data, since in EMOD the killed mosquitoes didn't get to feed (and potentially transmit any current infections) before dying in the model, while they may be transmitting in the real world 
# png(filename=paste0(data_dir, '/compare_net_outcomes_', fit_version_name, '_', ig2_version,'IG2v2_refMortRescale', round(100*arbitrary_mortality_ref_to_emod_rescaler),'.png'), width=9, height=9, units='in', res=900)
# 
# par(mfcol=c(3,3), mgp=c(2.3,1,0))
# ## combined pyr versus PBO across studies
# # feed
# plot(dtk_pyr_outcomes1[[1]], dtk_pbo_outcomes1[[1]], type='l', xlim=c(0,1), ylim=c(0,1), ylab='pbo', xlab='pyr', main="fraction feed (and survive in model)")
# lines(dtk_pyr_outcomes0[[1]], dtk_pbo_outcomes0[[1]], type='l', col='grey')
# abline(a=0,b=1, col=rgb(0.9,0.9,0.9,0.9), lwd=2)
# points(results_df$pyr_fed[results_df$plus_type=='PBO'], results_df$plus_fed[results_df$plus_type=='PBO'], col=as.factor(results_df$study[results_df$plus_type=='PBO']), cex=1.5, pch=20)
# # don't feed
# plot(1-dtk_pyr_outcomes1[[1]], 1-dtk_pbo_outcomes1[[1]], type='l', xlim=c(0,1), ylim=c(0,1), ylab='pbo', xlab='pyr', main="fraction don't feed (survive or die)")
# lines(1-dtk_pyr_outcomes0[[1]], 1-dtk_pbo_outcomes0[[1]], type='l', col='grey')
# abline(a=0,b=1, col=rgb(0.9,0.9,0.9,0.9), lwd=2)
# points(1-results_df$pyr_fed[results_df$plus_type=='PBO'], 1-results_df$plus_fed[results_df$plus_type=='PBO'], col=as.factor(results_df$study[results_df$plus_type=='PBO']), cex=1.5, pch=20)
# # die
# plot(dtk_pyr_outcomes1[[3]], dtk_pbo_outcomes1[[3]], type='l', xlim=c(0,1), ylim=c(0,1), ylab='pbo', xlab='pyr', main='fraction die (without feeding in model)')
# lines(dtk_pyr_outcomes0[[3]], dtk_pbo_outcomes0[[3]], type='l', col='grey')
# abline(a=0,b=1, col=rgb(0.9,0.9,0.9,0.9), lwd=2)
# points(results_df$pyr_mort[results_df$plus_type=='PBO']*arbitrary_mortality_ref_to_emod_rescaler, results_df$plus_mort[results_df$plus_type=='PBO']*arbitrary_mortality_ref_to_emod_rescaler, col=as.factor(results_df$study[results_df$plus_type=='PBO']), cex=1.5, pch=20)
# 
# ## combined pyr versus IG2 across studies
# # feed
# plot(dtk_pyr_outcomes1[[1]], dtk_ig2_outcomes1[[1]], type='l', xlim=c(0,1), ylim=c(0,1), ylab='ig2', xlab='pyr', main="fraction feed (and survive in model)")
# lines(dtk_pyr_outcomes0[[1]], dtk_ig2_outcomes0[[1]], type='l', col='grey')
# abline(a=0,b=1, col=rgb(0.9,0.9,0.9,0.9), lwd=2)
# points(results_df$pyr_fed[results_df$plus_type=='IG2'], results_df$plus_fed[results_df$plus_type=='IG2'], col=as.factor(results_df$study[results_df$plus_type=='IG2']), cex=1.5, pch=20)
# # don't feed
# plot(1-dtk_pyr_outcomes1[[1]], 1-dtk_ig2_outcomes1[[1]], type='l', xlim=c(0,1), ylim=c(0,1), ylab='ig2', xlab='pyr', main="fraction don't feed (survive or die)")
# lines(1-dtk_pyr_outcomes0[[1]], 1-dtk_ig2_outcomes0[[1]], type='l', col='grey')
# abline(a=0,b=1, col=rgb(0.9,0.9,0.9,0.9), lwd=2)
# points(1-results_df$pyr_fed[results_df$plus_type=='IG2'], 1-results_df$plus_fed[results_df$plus_type=='IG2'], col=as.factor(results_df$study[results_df$plus_type=='IG2']), cex=1.5, pch=20)
# # die
# plot(dtk_pyr_outcomes1[[3]], dtk_ig2_outcomes1[[3]], type='l', xlim=c(0,1), ylim=c(0,1), ylab='ig2', xlab='pyr', main='fraction die (without feeding in model)')
# lines(dtk_pyr_outcomes0[[3]], dtk_ig2_outcomes0[[3]], type='l', col='grey')
# abline(a=0,b=1, col=rgb(0.9,0.9,0.9,0.9), lwd=2)
# points(results_df$pyr_mort[results_df$plus_type=='IG2']*arbitrary_mortality_ref_to_emod_rescaler, results_df$plus_mort[results_df$plus_type=='IG2']*arbitrary_mortality_ref_to_emod_rescaler, col=as.factor(results_df$study[results_df$plus_type=='IG2']), cex=1.5, pch=20)
# 
# ## combined PBO versus IG2 across studies
# # feed
# plot(dtk_pbo_outcomes1[[1]], dtk_ig2_outcomes1[[1]], type='l', xlim=c(0,1), ylim=c(0,1), ylab='ig2', xlab='pbo', main="fraction feed (and survive in model)")
# lines(dtk_pbo_outcomes0[[1]], dtk_ig2_outcomes0[[1]], type='l', col='grey')
# abline(a=0,b=1, col=rgb(0.9,0.9,0.9,0.9), lwd=2)
# points(1-pbo_block_all, 1-ig2_block_all, col=as.factor(study), cex=1.5, pch=20)
# # don't feed
# plot(1-dtk_pbo_outcomes1[[1]], 1-dtk_ig2_outcomes1[[1]], type='l', xlim=c(0,1), ylim=c(0,1), ylab='ig2', xlab='pbo', main="fraction don't feed (survive or die)")
# lines(1-dtk_pbo_outcomes0[[1]], 1-dtk_ig2_outcomes0[[1]], type='l', col='grey')
# abline(a=0,b=1, col=rgb(0.9,0.9,0.9,0.9), lwd=2)
# points(pbo_block_all, ig2_block_all, col=as.factor(study), cex=1.5, pch=20)
# # die
# plot(dtk_pbo_outcomes1[[3]], dtk_ig2_outcomes1[[3]], type='l', xlim=c(0,1), ylim=c(0,1), ylab='ig2', xlab='pbo', main='fraction die (without feeding in model)')
# lines(dtk_pbo_outcomes0[[3]], dtk_ig2_outcomes0[[3]], type='l', col='grey')
# points(pbo_mort_all*arbitrary_mortality_ref_to_emod_rescaler, ig2_mort_all*arbitrary_mortality_ref_to_emod_rescaler, col=as.factor(study), cex=1.5, pch=20)
# abline(a=0,b=1, col=rgb(0.9,0.9,0.9,0.9), lwd=2)
# 
# dev.off()

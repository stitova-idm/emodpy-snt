############################################################################################################################################
# compare_smc_coverage_mncp_mc.R
# May 2020
# Monique Ambrose

# plot comparison of SMC coverages reported in Malaria Consortium document for sampled regions and the coverage reported using NMCP reports
#   - for the number of doses given reported by NMCP, adjust to account for the fraction of children that vomit/spit out their first dose and are given a second dose
#   - for the population size, use the PNLP 2016 population sizes, but increased to account for a 2.86% population growth rate
############################################################################################################################################

library(data.table)
library(ggplot2)
library(dplyr)
library(rootSolve)

# --------- directories for input and output files --------- #
user = Sys.getenv("USERNAME")
user_path = file.path("C:/Users",user)
box_filepath = paste(user_path, '/Dropbox (IDM)/NU_collaboration', sep='')
ds_pop_filepath = paste0(box_filepath, '/hbhi_nigeria/snt_2022/admin_pop_archetype.csv')
smc_mc_filepath = paste0(box_filepath, '/hbhi_nigeria/SMC/MC_extracted_SMC_coverages.csv')
smc_nmcp_filepath_pre2019 = paste0(box_filepath, '/nigeria_who/Interventions/SMC/smc_distribution_by_LGA.csv')
# smc_nmcp_filepath_post2019 = paste0(box_filepath, '/nigeria_who/NGA_2022_SNT/_Submitted_data/Routine data/SMC_2020_2021_cleaned.csv')  # earlier version - only at state level
smc_nmcp_filepath_post2019 = paste0(box_filepath, '/nigeria_who/NGA_2022_SNT/NGA_Intervention_data_2021_31102022_SMC.csv')
smc_old_sims_filepath = paste0(box_filepath, '/hbhi_nigeria/simulation_inputs/projection_csvs/2010_2020_LGA_intervention_files/SMC/smc_df2_fin_by_round_v3.csv')
save_plot_filepath = paste(box_filepath, '/hbhi_nigeria/snt_2022/interventions/SMC/plots', sep='')
new_sim_input_filepath = paste0(box_filepath, '/hbhi_nigeria/snt_2022/simulation_inputs/interventions_2010_toPresent/smc_2010_toPresent.csv')
# new_sim_coverages_filepath = paste0(box_filepath, '/hbhi_burkina/simulation_inputs/_scenarios_ms/SMC_2010_2020.csv')
# new_sim_coverages_2020_filepath = paste0(box_filepath, '/hbhi_burkina/simulation_inputs/_scenarios_ms/SMC_projection_allAge5_allRounds4.csv')
script_dir = paste0(user_path, '/Documents/malaria-nga-snt22')
source(paste0(script_dir,'/standardize_admin_names.R'))

# start day for simulations
sim_startday = as.Date('01-01-2010', format='%m-%d-%Y')
smc_adherence = 0.8

# --------- read in files --------- #
ds_pop = read.csv(ds_pop_filepath)
smc_mc = read.csv(smc_mc_filepath)
smc_nmcp_pre2019_0 = fread(smc_nmcp_filepath_pre2019)
smc_nmcp_post2019_0 = fread(smc_nmcp_filepath_post2019)
smc_old_sims = fread(smc_old_sims_filepath)
old_high_access_frac = 0.5
# smc_new_sims = fread(new_sim_coverages_filepath)
# smc_new_sims_2020 = fread(new_sim_coverages_2020_filepath)

plot_year = 2020
pop_file_year = 2018
pop_increase_rescale=1.2

# plot parameters
png_height = 4 *3/5
png_width = 5*3/5

# --------- set parameters --------- #
# population growth rate assumed in population projections
annual_pop_growth = 1.026
# fraction of population U5
frac_U5 = 0.166 # BFA: 0.1677
# fraction of population 5-10 years
frac_5_10 = 0.146 # BFA: 0.1487
# fraction of children who spit out or vomit dose within first 30 minutes after administration and are given second dose
#   if data is in terms of children treated instead of total doses distributed, put frac_spit = 0. Otherwise, obtain value from MC. E.g., 0.14
frac_spit = 0

# set maximum U5 SMC coverage
max_smc_per_round_cov = 0.8
# maximum coverage in the high-access group
max_high_access_coverage = 0.98
# maximum coverage in the low-access group
max_low_access_coverage = 0.8

# what fraction of individuals ages 5-10 were treated (used when adjusting for doses given off-target)
frac_treated_5_10 = 0.3
# fraction of individuals 5-10 assumed to be in the high-access group
high_access_5_10 = 0.2
# _effective_ coverage in 5-10 year olds relative to U5 coverage (includes that underdosing may not be as effective)
coverage_multiplier_5_10 = 0.1


# # compare different population sizes estimates
# library(reshape2)
# ds_pop_long = melt(ds_pop[,c('LGA','population','geopode.pop','geopode.pop.0.4')], id.vars=c('LGA'))
# ggplot(ds_pop_long, aes(x=as.factor(variable), y=value))+
#   geom_violin(alpha=0.5)+
#   # geom_dotplot(aes(fill=LGA), dotsize=0.1)+
#   theme(legend.position='none')


# --------- reformat/prepare tables --------- #
standardize_column_names = function(df){
  colnames(df)[which(colnames(df) %in% c('state', 'adm1', 'Administrative.levels.0.and.1'))] = 'State'
  colnames(df)[which(colnames(df) %in% c('LGA', 'DS', 'adm2', 'admin2'))] = 'admin_name'
  colnames(df)[which(colnames(df) %in% c('Year'))] = 'year'
  colnames(df)[which(colnames(df) %in% c('Number.of.children.treated.in.1st.cycle', 'smc1_num'))] = 'smc_treat.1'
  colnames(df)[which(colnames(df) %in% c('Number.of.children.treated.in.2nd.cycle', 'smc2_num'))] = 'smc_treat.2'
  colnames(df)[which(colnames(df) %in% c('Number.of.children.treated.in.3rd.cycle', 'smc3_num'))] = 'smc_treat.3'
  colnames(df)[which(colnames(df) %in% c('Number.of.children.treated.in.4th.cycle', 'smc4_num'))] = 'smc_treat.4'
  colnames(df)[which(colnames(df) %in% c('smc1_targ'))] = 'smc_targ.1'
  colnames(df)[which(colnames(df) %in% c('smc2_targ'))] = 'smc_targ.2'
  colnames(df)[which(colnames(df) %in% c('smc3_targ'))] = 'smc_targ.3'
  colnames(df)[which(colnames(df) %in% c('smc4_targ'))] = 'smc_targ.4'
  colnames(df)[which(colnames(df) %in% c('cov1', 'smc1_cov'))] = 'smc_cov.1'
  colnames(df)[which(colnames(df) %in% c('cov2', 'smc2_cov'))] = 'smc_cov.2'
  colnames(df)[which(colnames(df) %in% c('cov3', 'smc3_cov'))] = 'smc_cov.3'
  colnames(df)[which(colnames(df) %in% c('cov4', 'smc4_cov'))] = 'smc_cov.4'
  
  return(df)
}


# standardize column names
smc_mc = standardize_column_names(smc_mc)
smc_old_sims = standardize_column_names(smc_old_sims)
smc_nmcp_pre2019_0 = standardize_column_names(smc_nmcp_pre2019_0)
smc_nmcp_post2019_0 = standardize_column_names(smc_nmcp_post2019_0)
smc_nmcp_0 = merge(smc_nmcp_pre2019_0, smc_nmcp_post2019_0, all=TRUE)
smc_nmcp_0$admin_name[which(smc_nmcp_0$admin_name %in% c('', ' ', 'Mmc', 'Birnin'))] = NA
smc_nmcp_0$admin_name[intersect(which(smc_nmcp_0$admin_name == 'Bassa'), which(smc_nmcp_0$State=='Kogi'))] = 'Bassa1'
smc_nmcp_0$admin_name[intersect(which(smc_nmcp_0$admin_name == 'Bassa'), which(smc_nmcp_0$State=='Plateau'))] = 'Bassa2'
smc_nmcp_0$admin_name[intersect(which(smc_nmcp_0$admin_name %in% c('Nassarawa', 'Nasarawa')), which(smc_nmcp_0$State=='Kano'))] = 'Nasarawa1'
smc_nmcp_0$admin_name[intersect(which(smc_nmcp_0$admin_name %in% c('Nassarawa', 'Nasarawa')), which(smc_nmcp_0$State=='Nasarawa'))] = 'Nasarawa2'
smc_nmcp_0 = as.data.frame(smc_nmcp_0)
smc_old_sims = as.data.frame(smc_old_sims)




######################################################################
# process NMEP data and turn into simulation inputs
######################################################################
# divide into dataframes with state-level and admin-level data
smc_nmcp_admin_0 = smc_nmcp_0[!is.na(smc_nmcp_0$admin_name),]
smc_nmcp_state_0 = smc_nmcp_0[is.na(smc_nmcp_0$admin_name),]
if(nrow(smc_nmcp_admin_0) + nrow(smc_nmcp_state_0) != nrow(smc_nmcp_0)) warning('There was an issue with partitioning the SMC data between entries with admin-level and only state-level information')

# remove state rows when there are admin rows for that state in that year already
duplicate_rows = c()
distinct_state_years_from_admin_level = distinct(smc_nmcp_admin_0[,which(colnames(smc_nmcp_admin_0) %in% c('State', 'year'))])
distinct_state_years_from_state_level = smc_nmcp_state_0[,which(colnames(smc_nmcp_state_0) %in% c('State', 'year'))]
for(rr in 1:nrow(smc_nmcp_state_0)){
  if(nrow(merge(distinct_state_years_from_state_level[rr,], distinct_state_years_from_admin_level))>0){
    duplicate_rows = c(duplicate_rows, rr)
  }
}
if(length(duplicate_rows)>0) smc_nmcp_state_0 = smc_nmcp_state_0[-duplicate_rows,]

# standardize admin name strings
smc_nmcp_admin_0 = standardize_admin_names_in_df(target_names_df=ds_pop, origin_names_df=smc_nmcp_admin_0, target_names_col='admin_name', origin_names_col='admin_name')
smc_old_sims = standardize_admin_names_in_df(target_names_df=ds_pop, origin_names_df=smc_old_sims, target_names_col='admin_name', origin_names_col='admin_name')

# check whether a particular admin-state are duplicated
duplicates = smc_nmcp_admin_0[,c('admin_name', 'year')][duplicated(smc_nmcp_admin_0[,c('admin_name', 'year')]),]
if(nrow(duplicates)>0){
  duplicate_df = merge(duplicates, smc_nmcp_admin_0, all.x=TRUE, all.y=FALSE)
  warning('some LGA-years are recorded across multiple rows')
}

# convert dataframes from wide to long format, with each round in its own column
smc_nmcp_admin = reshape(smc_nmcp_admin_0, direction='long', 
                        varying=matrix(c("smc_treat.1", "smc_treat.2", "smc_treat.3", "smc_treat.4", 
                                  "smc_targ.1", "smc_targ.2", "smc_targ.3", "smc_targ.4",
                                  "smc_cov.1", "smc_cov.2", "smc_cov.3", "smc_cov.4"), nrow=3, byrow=TRUE), 
                        timevar='round',
                        sep='.',
                        times=c('1', '2', '3', '4'),
                        v.names=c('smc_treat', 'smc_targ', 'smc_cov'),
                        idvar=c('admin_name', 'year'))
smc_nmcp_state = reshape(smc_nmcp_state_0, direction='long', 
                        varying=matrix(c("smc_treat.1", "smc_treat.2", "smc_treat.3", "smc_treat.4", 
                                         "smc_targ.1", "smc_targ.2", "smc_targ.3", "smc_targ.4",
                                         "smc_cov.1", "smc_cov.2", "smc_cov.3", "smc_cov.4"), nrow=3, byrow=TRUE), 
                        timevar='round',
                        sep='.',
                        times=c('1', '2', '3', '4'),
                        v.names=c('smc_treat', 'smc_targ', 'smc_cov'),
                        idvar=c('State', 'year'))


# create population size projections for all LGAs
ds_pop$pop_size_rescaled = ds_pop$pop_size * pop_increase_rescale
ds_pop = ds_pop %>% dplyr::select(admin_name, State, pop_size, pop_size_rescaled)
min_year = min(smc_nmcp_0$year, na.rm=TRUE)
max_year = max(smc_nmcp_0$year, na.rm=TRUE)
pop_df = data.frame('admin_name'=c(), 'State'=c(), 'year'=c(), 'population'=c())
for(yy in min_year:max_year){
  projection_multiplier = (annual_pop_growth)^(yy - pop_file_year)
  projected_pop_df = ds_pop
  projected_pop_df$population = projected_pop_df$pop_size_rescaled * projection_multiplier
  projected_pop_df$year = yy
  pop_df = rbind(pop_df, data.frame('admin_name'=projected_pop_df$admin_name, 'State'=projected_pop_df$State, 'year'=projected_pop_df$year, 'population'=projected_pop_df$population))
}
# ggplot(pop_df, aes(x=year, y=population, color=admin_name))+
#   geom_line()+
#   theme(legend.position='none')

# add column with total population size, U5 population size, and 5-10 population size for each DS
if(length(smc_nmcp_admin$admin_name[which(!(smc_nmcp_admin$admin_name %in% ds_pop$admin_name))])>0) warning('Some admin names are not matched between files')
smc_nmcp_admin = merge(smc_nmcp_admin, pop_df, by=c("admin_name", "year", 'State'))
state_pop = pop_df %>% group_by(State, year) %>%
  summarise(population = sum(population))
smc_nmcp_state = merge(smc_nmcp_state, state_pop, by=c("State", "year"))
smc_nmcp_admin$pop_U5 = smc_nmcp_admin$population * frac_U5
smc_nmcp_admin$pop_5_10 = smc_nmcp_admin$population * frac_5_10
smc_nmcp_state$pop_U5 = smc_nmcp_state$population * frac_U5
smc_nmcp_state$pop_5_10 = smc_nmcp_state$population * frac_5_10


# combine LGA-level data into state-level coverage to make comparisons with MC and population estimates
# expand the State-level coverages to LGAs
smc_nmcp_admin_to_state = smc_nmcp_admin[,c('State','year','round','population','pop_U5', 'pop_5_10','smc_treat','smc_targ')] %>% group_by(State, year, round) %>%
  summarise_all(sum)
smc_nmcp_admin_to_state$ data_level = 'lga'
smc_nmcp_state$data_level = 'state'
smc_nmcp_state = merge(smc_nmcp_admin_to_state, smc_nmcp_state, all=TRUE)
smc_nmcp_state$smc_cov = smc_nmcp_state$smc_treat/smc_nmcp_state$smc_targ

# NMCP: get number of doses / U5 population size (not accounting for spitting)
smc_nmcp_admin$doses_per_U5 = smc_nmcp_admin$smc_treat/smc_nmcp_admin$pop_U5
smc_nmcp_state$doses_per_U5 = smc_nmcp_state$smc_treat/smc_nmcp_state$pop_U5

# NMCP: get number of successful doses
#   given that some fraction of children vomit or spit out a dose within the first 30 minutes after administration and need to be given a second dose, adjust dose numbers to 
#   not include the wasted doses: total_doses_given = number_children_who_get_smc * ((1-frac_spit) * 1 + frac_spit * 2) --> number_children_who_get_smc = total_doses_given/(1+frac_spit)
# estimate number of treated children, accounting for wasted doses
smc_nmcp_admin$smc_treat_adj = smc_nmcp_admin$smc_treat/(1+frac_spit)
smc_nmcp_state$smc_treat_adj = smc_nmcp_state$smc_treat/(1+frac_spit)

# NMCP: get number of successful treatment doses / U5 population size (accounting for spitting)
smc_nmcp_admin$doses_per_U5_adj = smc_nmcp_admin$smc_treat_adj/smc_nmcp_admin$pop_U5
smc_nmcp_state$doses_per_U5_adj = smc_nmcp_state$smc_treat_adj/smc_nmcp_state$pop_U5

# NMCP: if we assume that the 5-10 coverage is frac_treated_5_10 in all regions with SMC and adjust the NMCP doses accordingly, what fraction of U5 would we estimate receive SMC?
# calculate number of doses in O5 and subtract them from the total number of doses. use new number to calculate fraction of U5 covered
smc_nmcp_admin$adjusted_nmcp_frac_U5_received = sapply((smc_nmcp_admin$smc_treat_adj - frac_treated_5_10 * smc_nmcp_admin$pop_5_10) / smc_nmcp_admin$pop_U5, max,0)
smc_nmcp_state$adjusted_nmcp_frac_U5_received = sapply((smc_nmcp_state$smc_treat_adj - frac_treated_5_10 * smc_nmcp_state$pop_5_10) / smc_nmcp_state$pop_U5, max,0)

# NMCP: cap per-round coverage at maximum
smc_nmcp_admin$cov_u5_estimate = sapply(smc_nmcp_admin$adjusted_nmcp_frac_U5_received, min, max_smc_per_round_cov)
smc_nmcp_state$cov_u5_estimate = sapply(smc_nmcp_state$adjusted_nmcp_frac_U5_received, min, max_smc_per_round_cov)




#### merge and insert placeholder state-level data
# expand the State-level coverages to LGAs
smc_nmcp_state_level_data = smc_nmcp_state[smc_nmcp_state$data_level == 'state',]
smc_nmcp_state_to_admin = merge(smc_nmcp_state_level_data[,-which(colnames(smc_nmcp_state_level_data) == 'admin_name')], ds_pop[,c('State', 'admin_name')], all=TRUE)
smc_nmcp_state_to_admin = smc_nmcp_state_to_admin[!is.na(smc_nmcp_state_to_admin$cov_u5_estimate),]
smc_nmcp_state_to_admin$data_level = 'state'
smc_nmcp_admin$data_level = 'lga'
smc_nmcp_all = merge(smc_nmcp_state_to_admin, smc_nmcp_admin, all=TRUE)

# since we do not have data from 2019, assume the same coverages and admins as in 2018
if(!(2019 %in% smc_nmcp_all$year)){
  smc_nmcp_all_2019 = smc_nmcp_all[smc_nmcp_all$year==2018,]
  smc_nmcp_all_2019$year = 2019
  smc_nmcp_all = rbind(smc_nmcp_all, smc_nmcp_all_2019)
}


################################################
# plots of coverage adjustment and timeseries
################################################

# children treated divided by estimated U5 population
ggplot(data = smc_nmcp_all[smc_nmcp_all$year==2020,], aes(x = round, y = doses_per_U5_adj, fill = as.factor(round)))+
  scale_fill_viridis_d( option = "D")+
  geom_violin(alpha=0.4, color='grey')+
  # geom_violin(alpha=0.4, position = position_dodge(width = .75),size=1,color="black") +
  # geom_boxplot(notch = TRUE,  outlier.size = -1, color="black",lwd=1.2, alpha = 0.7)+
  geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=1)+
  ylab(c(expression("number treated / \n estimated number U5")))  +
  xlab(c("SMC cycle"))+
  coord_cartesian(ylim=c(0,2))+
  theme_minimal()+ 
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(), legend.position = "none", axis.title.y = element_text(margin = margin(t = 0, r = 4, b = 0, l = 10)))
ggsave(paste0(save_plot_filepath, '/0_admin_smc_treat_over_est_pop_NMEP.png'), width=png_width, height=png_height, units='in', dpi=900)



# if we assume that the 5-10 coverage is frac_treated_5_10 in all regions with SMC and adjust the NMCP doses accordingly, what fraction of U5 would we estimate receive SMC?
# calculate number of doses in O5 and subtract them from the total number of doses. use new number to calculate fraction of U5 covered
ggplot(data = smc_nmcp_all[smc_nmcp_all$year==2020,], aes(x = round, y = adjusted_nmcp_frac_U5_received, fill = as.factor(round)))+
  scale_fill_viridis_d( option = "D")+
  geom_violin(alpha=0.4, color='grey')+
  # geom_violin(alpha=0.4, position = position_dodge(width = .75),size=1,color="black") +
  # geom_boxplot(notch = TRUE,  outlier.size = -1, color="black",lwd=1.2, alpha = 0.7)+
  geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=1)+
  ylab(c(expression("estimated number U5 treated / \n number of children U5")))  +
  xlab(c("SMC cycle"))+
  coord_cartesian(ylim=c(0,2))+
  theme_minimal()+ 
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(), legend.position = "none", axis.title.y = element_text(margin = margin(t = 0, r = 4, b = 0, l = 10)))
ggsave(paste0(save_plot_filepath, '/1_admin_smc_adj_offtarget_NMEP.png'), width=png_width, height=png_height, units='in', dpi=900)


# now cap maximum per-round coverage
ggplot(data = smc_nmcp_all[smc_nmcp_all$year==2020,], aes(x = round, y = cov_u5_estimate, fill = as.factor(round)))+
  scale_fill_viridis_d( option = "D")+
  geom_violin(alpha=0.4, color='grey')+
  # geom_violin(alpha=0.4, position = position_dodge(width = .75),size=1,color="black") +
  # geom_boxplot(notch = TRUE,  outlier.size = -1, color="black",lwd=1.2, alpha = 0.7)+
  geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=1)+
  ylab(c(expression("estimated number U5 treated / \n number of children U5")))  +
  xlab(c("SMC cycle"))+
  coord_cartesian(ylim=c(0,2))+
  theme_minimal()+ 
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(), legend.position = "none", axis.title.y = element_text(margin = margin(t = 0, r = 4, b = 0, l = 10)))
ggsave(paste0(save_plot_filepath, '/2_admin_smc_adj_offtarget_NMEP_wMaxCap.png'), width=png_width, height=png_height, units='in', dpi=900)


# plot average coverage in each LGA-year across rounds
smc_nmcp_year_ave = smc_nmcp_all %>% group_by(admin_name, year) %>%
  summarise(average_u5_coverage = mean(cov_u5_estimate))
# put zero coverage in years without data for admins that are ever present in the dataset
admins_with_smc = unique(smc_nmcp_year_ave$admin_name)
smc_nmcp_year_ave = merge(smc_nmcp_year_ave, pop_df[which(pop_df$admin_name %in% admins_with_smc), c('admin_name','year', 'State')], all=TRUE)
smc_nmcp_year_ave$average_u5_coverage[is.na(smc_nmcp_year_ave$average_u5_coverage)] = 0
# add a tiny bit of noise for jitter
smc_nmcp_year_ave$average_u5_coverage_jittered = smc_nmcp_year_ave$average_u5_coverage+runif(nrow(smc_nmcp_year_ave),0,0.05)
ggplot(smc_nmcp_year_ave, aes(x=year, y=average_u5_coverage_jittered, color=admin_name))+
  geom_point()+
  geom_line()+
  ylab('average per-cycle coverage in U5') +
  theme_bw()+
  theme(legend.position='none') +
  facet_wrap('State', nrow=5)
ggsave(paste0(save_plot_filepath, '/smc_coverage_timeseries_estimates.png'), width=png_width*3, height=png_height*3, units='in', dpi=900)



######################################################################
# get high/low access fractions
######################################################################


# use the number that receive 0, 1, 2, 3, or 4 doses to calculate high/low-access fractions
num_doses = c(0,1,2,3,4)
frac_received_num_doses = c(0.06, 0.052, 0.133, 0.153, 0.602)

# give the fraction of individuals in high-low access groups and the per-round coverage for each group, calculate expected number of individuals receiving 0, 1, 2, 3, or 4 total doses
calculate_expected_dist_num_doses = function(frac_high_access, frac_low_access, per_round_cov_high, per_round_cov_low){
  expected_frac_num_doses_high = c((1-per_round_cov_high)^4, 4*(1-per_round_cov_high)^3*per_round_cov_high, 6*(1-per_round_cov_high)^2*per_round_cov_high^2, 4*(1-per_round_cov_high)*per_round_cov_high^3, per_round_cov_high^4)
  expected_frac_num_doses_low = c((1-per_round_cov_low)^4, 4*(1-per_round_cov_low)^3*per_round_cov_low, 6*(1-per_round_cov_low)^2*per_round_cov_low^2, 4*(1-per_round_cov_low)*per_round_cov_low^3, per_round_cov_low^4)
  expected_frac_num_doses_total = frac_high_access * expected_frac_num_doses_high + frac_low_access * expected_frac_num_doses_low
  return(expected_frac_num_doses_total)
}

# goal will be to minimize the difference between the observed and expected dose distribution
calculate_difference = function(frac_received_num_doses, frac_high_access, frac_low_access, per_round_cov_high, per_round_cov_low){
  expected_frac_num_doses_total = calculate_expected_dist_num_doses(frac_high_access, frac_low_access, per_round_cov_high, per_round_cov_low)
  # sum_square_diff = sum((expected_frac_num_doses_total - frac_received_num_doses)^2)
  abs_diff = sum(abs(expected_frac_num_doses_total - frac_received_num_doses))
  # abs_rel_diff = sum(abs((expected_frac_num_doses_total - frac_received_num_doses)/frac_received_num_doses))
  return(abs_diff)
}


# # expected fraction of individuals receiving each number of doses when coverage is set at per-round max of 80%
# per_round_cov_est = get_high_low_access_coverage(high_access_fraction=0.62, total_coverage=0.8, max_high_access_coverage=0.98)
# calculate_expected_dist_num_doses(frac_high_access=0.62, frac_low_access=(1-0.62), per_round_cov_high=per_round_cov_est[1], per_round_cov_low=per_round_cov_est[2])
# # fraction receiving at least one dose
# 1-calculate_expected_dist_num_doses(frac_high_access=0.62, frac_low_access=(1-0.62), per_round_cov_high=per_round_cov_est[1], per_round_cov_low=per_round_cov_est[2])[1]


# # do a full sweep exploring the fraction of individuals in low/high access and the coverage for each to match observations
# frac_high_access_v = seq(0,1,0.01)
# per_round_cov_high_v = seq(0,1,0.01)
# per_round_cov_low_v = seq(0,1,0.01)
# diff_array = array(data=NA, dim=c(length(frac_high_access_v), length(per_round_cov_high_v), length(per_round_cov_low_v)))
# for(i1 in 1:length(frac_high_access_v)){
#   for(i2 in 1:length(per_round_cov_high_v)){
#     for(i3 in 1:length(per_round_cov_low_v)){
#       frac_high_access = frac_high_access_v[i1]
#       frac_low_access = 1 - frac_high_access
#       per_round_cov_high = per_round_cov_high_v[i2]
#       per_round_cov_low = per_round_cov_low_v[i3]
#       
#       diff_array[i1,i2,i3] = calculate_difference(frac_received_num_doses, frac_high_access, frac_low_access, per_round_cov_high, per_round_cov_low)
#     }
#   }
# }
# best_match = which(diff_array == min(diff_array), arr.ind = TRUE)[1,]
# print('Version where high-access coverage is flexible:')
# print(calculate_expected_dist_num_doses(frac_high_access=frac_high_access_v[best_match[1]], frac_low_access=(1-frac_high_access_v[best_match[1]]), per_round_cov_high=per_round_cov_high_v[best_match[2]], per_round_cov_low=per_round_cov_low_v[best_match[3]]))
# print(frac_received_num_doses)
# print(c(frac_high_access_v=max(frac_high_access_v[best_match[1]], 1-frac_high_access_v[best_match[1]]), per_round_cov_high=max(per_round_cov_high_v[best_match[2]], per_round_cov_low_v[best_match[3]]), per_round_cov_low=min(per_round_cov_high_v[best_match[2]], per_round_cov_low_v[best_match[3]])))

# if we assume max_high_access_coverage coverage in the high-access group, what are the best matches and how well do they compare to the unconstrained version?
frac_high_access_v = seq(0,1,0.01)
per_round_cov_high_v = max_high_access_coverage
per_round_cov_low_v = seq(0,1,0.01)
diff_array = array(data=NA, dim=c(length(frac_high_access_v), length(per_round_cov_high_v), length(per_round_cov_low_v)))
for(i1 in 1:length(frac_high_access_v)){
  for(i2 in 1:length(per_round_cov_high_v)){
    for(i3 in 1:length(per_round_cov_low_v)){
      frac_high_access = frac_high_access_v[i1]
      frac_low_access = 1 - frac_high_access
      per_round_cov_high = per_round_cov_high_v[i2]
      per_round_cov_low = per_round_cov_low_v[i3]
      
      diff_array[i1,i2,i3] = calculate_difference(frac_received_num_doses, frac_high_access, frac_low_access, per_round_cov_high, per_round_cov_low)
    }
  }
}
best_match = which(diff_array == min(diff_array), arr.ind = TRUE)[1,]
print('Version where high-access coverage is fixed:')
print(calculate_expected_dist_num_doses(frac_high_access=frac_high_access_v[best_match[1]], frac_low_access=(1-frac_high_access_v[best_match[1]]), per_round_cov_high=per_round_cov_high_v[best_match[2]], per_round_cov_low=per_round_cov_low_v[best_match[3]]))
print(frac_received_num_doses)
print(c(frac_high_access_v=max(frac_high_access_v[best_match[1]], 1-frac_high_access_v[best_match[1]]), per_round_cov_high=max(per_round_cov_high_v[best_match[2]], per_round_cov_low_v[best_match[3]]), per_round_cov_low=min(per_round_cov_high_v[best_match[2]], per_round_cov_low_v[best_match[3]])))
# the fixed version looks good, so we'll use that for allocating coverage for the simulation input assumptions

estimated_high_access_fraction = max(frac_high_access_v[best_match[1]], 1-frac_high_access_v[best_match[1]])

######################################################################
# format files for simulation input
######################################################################
# remove rows with NA or 0 coverage
smc_nmcp_all = smc_nmcp_all[!is.na(smc_nmcp_all$cov_u5_estimate),]
smc_nmcp_all = smc_nmcp_all[(smc_nmcp_all$cov_u5_estimate>0),]

# add in effective coverage among off-target children aged 5-10
smc_nmcp_all$cov_5_10_estimate = smc_nmcp_all$cov_u5_estimate * coverage_multiplier_5_10

# add in high-access and low-access fractions
smc_nmcp_all$high_access_U5 = estimated_high_access_fraction
smc_nmcp_all$high_access_5_10 = high_access_5_10


# function to calculate the coverage for the high and low access groups, given a total coverage and the fraction of individuals in each group
get_high_low_access_coverage = function(high_access_fraction, total_coverage, max_high_access_coverage, max_low_access_coverage){
  # assign coverage first to the high-access group, up to max_high_access_coverage, then assign the remaining coverage to the low-access group, up to the maximum
  if(total_coverage > high_access_fraction*max_high_access_coverage){
    high_access_coverage = max_high_access_coverage
    low_access_coverage = min(((total_coverage - high_access_fraction*max_high_access_coverage) / (1-high_access_fraction)), max_low_access_coverage)
  } else{
    high_access_coverage = total_coverage / high_access_fraction
    low_access_coverage = 0
  }
  return(c(high_access_coverage, low_access_coverage))
}
smc_nmcp_all$coverage_high_access_U5 = NA
smc_nmcp_all$coverage_low_access_U5 = NA
smc_nmcp_all$coverage_high_access_5_10 = NA
smc_nmcp_all$coverage_low_access_5_10 = NA
for(ii in 1:nrow(smc_nmcp_all)){
  # U5
  high_low_coverages_U5 = get_high_low_access_coverage(high_access_fraction=smc_nmcp_all$high_access_U5[ii], total_coverage=smc_nmcp_all$cov_u5_estimate[ii], max_high_access_coverage=max_high_access_coverage, max_low_access_coverage=max_low_access_coverage)
  smc_nmcp_all$coverage_high_access_U5[ii] = high_low_coverages_U5[1]
  smc_nmcp_all$coverage_low_access_U5[ii] = high_low_coverages_U5[2]
  # 5-10
  high_low_coverages_5_10 = get_high_low_access_coverage(high_access_fraction=smc_nmcp_all$high_access_5_10[ii], total_coverage=smc_nmcp_all$cov_5_10_estimate[ii], max_high_access_coverage=max_high_access_coverage, max_low_access_coverage=max_low_access_coverage)
  smc_nmcp_all$coverage_high_access_5_10[ii] = high_low_coverages_5_10[1]
  smc_nmcp_all$coverage_low_access_5_10[ii] = high_low_coverages_5_10[2]
}


# add in the fifth round for relevant LGA-years when relevant. there should be 4 rounds unless a month is present for smc5_month
# we weren't provided with coverage for the fifth round, so assume it is the same as the fourth
fifth_round_subset = smc_nmcp_all[intersect(which(smc_nmcp_all$round == 4), which(!is.na(smc_nmcp_all$smc5_month) & !(smc_nmcp_all$smc5_month %in% c('', ' ')))),]
fifth_round_subset$round = 5
smc_nmcp_all = rbind(smc_nmcp_all, fifth_round_subset)



# add in distribution dates
# function to extract distribution month from a wide range of possible date formats
get_date = function(date_str, year, placeholder_month=7){
  month_names = c('Jan', 'Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
  date_string_endings = c('st','nd','rd','th')
  blank_date_strings = c(' ', '', NA)
  # default values that will be returned unless a proper date can be extracted
  month_numeric = NA
  day_numeric = 15
  
  if(any(sapply(month_names, grepl, date_str))){
    month_numeric = as.numeric(which(sapply(month_names, grepl, date_str)))
    if(length(month_numeric) > 1) month_numeric = month_numeric[!is.na(month_numeric)][1]
    if(any(sapply(date_string_endings, grepl, date_str))){
      cur_ending = date_string_endings[as.numeric(which(sapply(date_string_endings, grepl, date_str)))]
      day_numeric = as.numeric(str_extract(date_str, paste0("\\d+(?=",cur_ending, ")")))
      if(length(day_numeric) > 1) day_numeric = day_numeric[!is.na(day_numeric)][1]
      if((day_numeric<1) | (day_numeric>31) | is.na(day_numeric)) day_numeric = 15
    }
  }else if (grepl('/', date_str)){
    try_date = tryCatch( as.Date(date_str, tryFormats = c('%m/%d/%Y', '%m/%d/%y')), error=function(cond){return(NA)})
    if(!is.na(try_date)){
      month_numeric = lubridate::month(try_date)
      day_numeric = lubridate::day(try_date)
    } else{
      if(!(date_str %in% blank_date_strings)) warning(paste0('Did not recognize date format for: ', date_str, '. Assuming SMC begins in month ', placeholder_month))
    }
  } else{
    if(!(date_str %in% blank_date_strings)) warning(paste0('Did not recognize date format for: ', date_str, '. Assuming SMC begins in month ', placeholder_month))
  }
  if(!is.na(month_numeric)){
    date = as.Date(paste0(month_numeric, '/', day_numeric, '/', year), tryFormats = c('%m/%d/%Y', '%m/%d/%y'))
  } else date = as.Date(paste0(placeholder_month, '/', day_numeric, '/', year), tryFormats = c('%m/%d/%Y', '%m/%d/%y'))
  return(date)
}


smc_nmcp_all$date_first_round = as.Date(NA)
for(rr in 1:nrow(smc_nmcp_all)){
  smc_nmcp_all$date_first_round[rr] = get_date(date_str=smc_nmcp_all$smc1_month[rr], year=smc_nmcp_all$year[rr])
}
# sort(unique(smc_nmcp_all$date_first_round))


# get the day of the simulated distribution based on the interval between the distribution and the start day of the sim
smc_nmcp_all$simday_first_round = as.numeric(smc_nmcp_all$date_first_round - sim_startday)
# sort(unique(smc_nmcp_all$simday_first_round))

# insert the day for distribution based on the day of the first round and the current round number
smc_nmcp_all$simday = smc_nmcp_all$simday_first_round + (as.numeric(smc_nmcp_all$round) - 1) * 30
# sort(unique(smc_nmcp_all$simday))

# insert adherence to later admodiaquine doses
smc_nmcp_all$adherence = smc_adherence
smc_sim_input = smc_nmcp_all[, c('admin_name', 'State', 'year', 'round', 'simday', 'date_first_round', 
                                 "coverage_high_access_U5", "coverage_low_access_U5", "coverage_high_access_5_10", "coverage_low_access_5_10", 
                                 "high_access_U5", "high_access_5_10", "adherence", "data_level")]

# check that each row is distinct (no duplicate campaigns on the same day)
if(nrow(smc_sim_input) == nrow(distinct(smc_sim_input))) {
  write.csv(smc_sim_input, file=new_sim_input_filepath, row.names=FALSE)
} else{
  warning('PROBLEM DETECTED: there are duplicate campaigns in the SMC input file.')
}










####################################################################################################
# compare average across rounds from MC with average across rounds from NMEP where both are available
####################################################################################################
smc_mc1 = smc_mc[grep('EoC', smc_mc$cycle),]
smc_mc_ave = smc_mc1[,c('State', 'coverage', 'cycle', 'year')] %>% group_by(State, year) %>%
  summarise(mean_mc_coverage = mean(coverage))
smc_nmcp_state_ave = smc_nmcp_state[c('State', 'round', 'adjusted_nmcp_frac_U5_received', 'year')] %>% group_by(State, year) %>%
  summarise(mean_nmcp_coverage = mean(adjusted_nmcp_frac_U5_received))

smc_mc_nmcp_ave = merge(smc_nmcp_state_ave, smc_mc_ave, by=c('State', 'year'), all=TRUE)
xy_max = 1.3 # min(5, max(smc_mc_nmcp_ave[,-c(1,2)], na.rm=TRUE))
plot(smc_mc_nmcp_ave$mean_nmcp_coverage, smc_mc_nmcp_ave$mean_mc_coverage, type='p', bty='L', ylim=c(0, xy_max), xlim=c(0, xy_max), ylab='Coverage estimate from MC survey', xlab='Coverage estimate from NMEP data', pch=20, col='darkblue', cex=1.1)
lines(c(0,xy_max), c(0, xy_max))


# --------- plots of population -------- #

# generally, the target population is larger than the estimate of the U5 population size from the other source
#  --> as a result, we rescaled the base population with pop_increase_rescale, now they should match more closely
max_xy = max(c(smc_nmcp_state$pop_U5, smc_nmcp_state$smc_targ), na.rm=TRUE)
plot(smc_nmcp_state$pop_U5, smc_nmcp_state$smc_targ, type='p', bty='L', ylim=c(0, max_xy), xlim=c(0, max_xy), ylab='NMEP target population', xlab='U5 population estimate')
lines(c(0,max_xy), c(0, max_xy))


# --------- plots of NMCP state data --------- #


# filter to desired year
smc_nmcp = smc_nmcp_state[smc_nmcp_state$year == plot_year,]
# smc_old_sims = smc_old_sims[smc_old_sims$year == plot_year,]
# smc_mc = smc_mc[smc_mc$year == plot_year,]

# raw numbers of doses divided by reported target population size
ggplot(data = smc_nmcp, aes(x = round, y = smc_cov, fill = as.factor(round)))+
  scale_fill_viridis_d( option = "D")+
  geom_violin(alpha=0.4, color='grey')+
  # geom_violin(alpha=0.4, position = position_dodge(width = .75),size=1,color="black") +
  # geom_boxplot(notch = TRUE,  outlier.size = -1, color="black",lwd=1.2, alpha = 0.7)+
  geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=1)+
  ylab(c(expression("reported fraction of \n target population covered")))  +
  xlab(c("SMC cycle"))+
  coord_cartesian(ylim=c(0,2))+
  theme_minimal()+ 
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(), legend.position = "none", axis.title.y = element_text(margin = margin(t = 0, r = 4, b = 0, l = 10)))
ggsave(paste0(save_plot_filepath, '/0_smc_NMEP_raw_reported_coverage_', plot_year, '.png'), width=png_width, height=png_height, units='in', dpi=900)



# raw numbers of doses divided by U5 population (not adjusted for spitting/comiting)
ggplot(data = smc_nmcp, aes(x = round, y = doses_per_U5, fill = as.factor(round)))+
  scale_fill_viridis_d( option = "D")+
  geom_violin(alpha=0.4, color='grey')+
  # geom_violin(alpha=0.4, position = position_dodge(width = .75),size=1,color="black") +
  # geom_boxplot(notch = TRUE,  outlier.size = -1, color="black",lwd=1.2, alpha = 0.7)+
  geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=1)+
  ylab(c(expression("number treated / \n estimated number U5")))  +
  xlab(c("SMC cycle"))+
  coord_cartesian(ylim=c(0,2))+
  theme_minimal()+ 
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(), legend.position = "none", axis.title.y = element_text(margin = margin(t = 0, r = 4, b = 0, l = 10)))
ggsave(paste0(save_plot_filepath, '/0_smc_treat_over_est_pop_NMEP.png'), width=png_width, height=png_height, units='in', dpi=900)


# number of doses, adjusted for spitting / vomiting first dose, divided by U5 population
ggplot(data = smc_nmcp, aes(x = round, y = doses_per_U5_adj, fill = as.factor(round)))+
  scale_fill_viridis_d( option = "D")+
  geom_violin(alpha=0.4, color='grey')+
  # geom_violin(alpha=0.4, position = position_dodge(width = .75),size=1,color="black") +
  # geom_boxplot(notch = TRUE,  outlier.size = -1, color="black",lwd=1.2, alpha = 0.7)+
  geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=1)+
  ylab(c(expression("number doses distributed (adjusted for spitting) / \n number of children under 5")))  +
  xlab(c("SMC cycle"))+
  coord_cartesian(ylim=c(0,2))+
  theme_minimal()+ 
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(), legend.position = "none", axis.title.y = element_text(margin = margin(t = 0, r = 4, b = 0, l = 10)))
ggsave(paste0(save_plot_filepath, '/0_smc_treat_over_est_pop_NMEP_wSpit.png'), width=png_width, height=png_height, units='in', dpi=900)



# if we assume that the 5-10 coverage is 20% in all regions with SMC and adjust the NMCP doses accordingly, what fraction of U5 would we estimate receive SMC?
# calculate number of doses in O5 and subtract them from the total number of doses. use new number to calculate fraction of U5 covered
ggplot(data = smc_nmcp, aes(x = round, y = adjusted_nmcp_frac_U5_received, fill = as.factor(round)))+
  scale_fill_viridis_d( option = "D")+
  geom_violin(alpha=0.4, color='grey')+
  # geom_violin(alpha=0.4, position = position_dodge(width = .75),size=1,color="black") +
  # geom_boxplot(notch = TRUE,  outlier.size = -1, color="black",lwd=1.2, alpha = 0.7)+
  geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=1)+
  ylab(c(expression("estimated number U5 treated / \n number of children U5")))  +
  xlab(c("SMC cycle"))+
  coord_cartesian(ylim=c(0,2))+
  theme_minimal()+ 
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(), legend.position = "none", axis.title.y = element_text(margin = margin(t = 0, r = 4, b = 0, l = 10)))
ggsave(paste0(save_plot_filepath, '/1_smc_adj_offtarget_NMEP.png'), width=png_width, height=png_height, units='in', dpi=900)


# now cap maximum per-round coverage
ggplot(data = smc_nmcp, aes(x = round, y = cov_u5_estimate, fill = as.factor(round)))+
  scale_fill_viridis_d( option = "D")+
  geom_violin(alpha=0.4, color='grey')+
  # geom_violin(alpha=0.4, position = position_dodge(width = .75),size=1,color="black") +
  # geom_boxplot(notch = TRUE,  outlier.size = -1, color="black",lwd=1.2, alpha = 0.7)+
  geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=1)+
  ylab(c(expression("estimated number U5 treated / \n number of children U5")))  +
  xlab(c("SMC cycle"))+
  coord_cartesian(ylim=c(0,2))+
  theme_minimal()+ 
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(), legend.position = "none", axis.title.y = element_text(margin = margin(t = 0, r = 4, b = 0, l = 10)))
ggsave(paste0(save_plot_filepath, '/1_smc_adj_offtarget_NMEP_wMaxCap.png'), width=png_width, height=png_height, units='in', dpi=900)





# compare against old version of simulation input
smc_old_sims$cov_u5_estimate = (smc_old_sims$coverage_high_access + smc_old_sims$coverage_low_access)/2
# plot average coverage in each LGA-year across rounds
smc_old_year_ave = smc_old_sims %>% group_by(admin_name, year) %>%
  summarise(average_u5_coverage = mean(cov_u5_estimate))
# put zero coverage in years without data for admins that are ever present in the dataset
admins_with_smc_old = unique(smc_old_year_ave$admin_name)
smc_old_year_ave = merge(smc_old_year_ave, pop_df[which(pop_df$admin_name %in% admins_with_smc_old), c('admin_name','year', 'State')], all=TRUE)
smc_old_year_ave = smc_old_year_ave[smc_old_year_ave$year<2020,]
smc_old_year_ave$average_u5_coverage[is.na(smc_old_year_ave$average_u5_coverage)] = 0
# add a tiny bit of noise for jitter
smc_old_year_ave$average_u5_coverage_jittered = smc_old_year_ave$average_u5_coverage+runif(nrow(smc_old_year_ave),0,0.05)
ggplot(smc_old_year_ave, aes(x=year, y=average_u5_coverage_jittered, color=admin_name))+
  geom_point()+
  geom_line()+
  ylab('average per-cycle coverage in U5') +
  theme(legend.position='none') +
  facet_wrap('State', nrow=5)
ggsave(paste0(save_plot_filepath, '/smc_coverage_timeseries_OLD2019_estimates.png'), width=png_width*3, height=png_height*3, units='in', dpi=900)



# # aggregate to admin 0 level
# 
# # NMCP: get sums of doses and populations within each admin 0
# smc_nmcp_admin0 = smc_nmcp %>% group_by(State) %>% 
#   summarise(sum(smc1_treat),
#             sum(smc2_treat), 
#             sum(smc3_treat), 
#             sum(smc4_treat),
#             sum(population),
#             sum(pop_U5), 
#             sum(pop_5_10)) %>% ungroup()
# colnames(smc_nmcp_admin0) = gsub('sum\\(', '', colnames(smc_nmcp_admin0))
# colnames(smc_nmcp_admin0) = gsub(')', '', colnames(smc_nmcp_admin0))
# 
# # NMCP: get number of successful doses
# #   given that some fraction of children vomit or spit out a dose within the first 30 minutes after administration and need to be given a second dose, adjust dose numbers to 
# #   not include the wasted doses: total_doses_given = number_children_who_get_smc * ((1-frac_spit) * 1 + frac_spit * 2) --> number_children_who_get_smc = total_doses_given/(1+frac_spit)
# smc_nmcp_admin0$smc1_treat = smc_nmcp_admin0$smc1_treat/(1+frac_spit)
# smc_nmcp_admin0$smc2_treat = smc_nmcp_admin0$smc2_treat/(1+frac_spit)
# smc_nmcp_admin0$smc3_treat = smc_nmcp_admin0$smc3_treat/(1+frac_spit)
# smc_nmcp_admin0$smc4_treat = smc_nmcp_admin0$smc4_treat/(1+frac_spit)
# 
# # NMCP: get number of successful doses / U5 population size
# smc_nmcp_admin0$doses_per_U5_1 = smc_nmcp_admin0$smc1_treat/smc_nmcp_admin0$pop_U5
# smc_nmcp_admin0$doses_per_U5_2 = smc_nmcp_admin0$smc2_treat/smc_nmcp_admin0$pop_U5
# smc_nmcp_admin0$doses_per_U5_3 = smc_nmcp_admin0$smc3_treat/smc_nmcp_admin0$pop_U5
# smc_nmcp_admin0$doses_per_U5_4 = smc_nmcp_admin0$smc4_treat/smc_nmcp_admin0$pop_U5
# 
# # NMCP: if we assume that the 5-10 coverage is 20% in all regions with SMC and adjust the NMCP doses accordingly, what fraction of U5 would we estimate receive SMC?
# # calculate number of doses in O5 and subtract them from the total number of doses. use new number to calculate fraction of U5 covered
# smc_nmcp_admin0$adjusted_nmcp_frac_U5_received_1 = sapply((smc_nmcp_admin0$smc1_treat - frac_treated_5_10 * smc_nmcp_admin0$pop_5_10) / smc_nmcp_admin0$pop_U5, max,0)
# smc_nmcp_admin0$adjusted_nmcp_frac_U5_received_2 = sapply((smc_nmcp_admin0$smc2_treat - frac_treated_5_10 * smc_nmcp_admin0$pop_5_10) / smc_nmcp_admin0$pop_U5, max,0)
# smc_nmcp_admin0$adjusted_nmcp_frac_U5_received_3 = sapply((smc_nmcp_admin0$smc3_treat - frac_treated_5_10 * smc_nmcp_admin0$pop_5_10) / smc_nmcp_admin0$pop_U5, max,0)
# smc_nmcp_admin0$adjusted_nmcp_frac_U5_received_4 = sapply((smc_nmcp_admin0$smc4_treat - frac_treated_5_10 * smc_nmcp_admin0$pop_5_10) / smc_nmcp_admin0$pop_U5, max,0)



# # merge with the malaria consortium coverage values
# smc_nmcp_mc = merge(smc_nmcp_admin0, smc_mc, by='State')
#
# # MC: what coverage would have been calculated for U5 if ignore that some doses go to O5? (using coverages from MC)
# #   (coverage_u5*num_U5 + coverage_o5*num_o5) / num_u5
# smc_nmcp_mc$adjusted_mc_doses_per_U5_1 = (smc_nmcp_mc$frac_U5_received_1 * smc_nmcp_mc$pop_U5 + smc_nmcp_mc$frac_O5_received_1 * smc_nmcp_mc$pop_5_10) / smc_nmcp_mc$pop_U5
# smc_nmcp_mc$adjusted_mc_doses_per_U5_2 = (smc_nmcp_mc$frac_U5_received_2 * smc_nmcp_mc$pop_U5 + smc_nmcp_mc$frac_O5_received_2 * smc_nmcp_mc$pop_5_10) / smc_nmcp_mc$pop_U5
# smc_nmcp_mc$adjusted_mc_doses_per_U5_3 = (smc_nmcp_mc$frac_U5_received_3 * smc_nmcp_mc$pop_U5 + smc_nmcp_mc$frac_O5_received_3 * smc_nmcp_mc$pop_5_10) / smc_nmcp_mc$pop_U5
# smc_nmcp_mc$adjusted_mc_doses_per_U5_4 = (smc_nmcp_mc$frac_U5_received_4 * smc_nmcp_mc$pop_U5 + smc_nmcp_mc$frac_O5_received_4 * smc_nmcp_mc$pop_5_10) / smc_nmcp_mc$pop_U5
# 
# 
# # NMCP: if we assume that the 5-10 coverage is 20% in all regions with SMC and adjust the NMCP doses accordingly, what fraction of U5 would we estimate receive SMC?
# # calculate number of doses in O5 and subtract them from the total number of doses. use new number to calculate fraction of U5 covered
# smc_nmcp_mc$adjusted_nmcp_frac_U5_received_1 = sapply((smc_nmcp_mc$smc1_treat - frac_treated_5_10 * smc_nmcp_mc$pop_5_10) / smc_nmcp_mc$pop_U5, max,0)
# smc_nmcp_mc$adjusted_nmcp_frac_U5_received_2 = sapply((smc_nmcp_mc$smc2_treat - frac_treated_5_10 * smc_nmcp_mc$pop_5_10) / smc_nmcp_mc$pop_U5, max,0)
# smc_nmcp_mc$adjusted_nmcp_frac_U5_received_3 = sapply((smc_nmcp_mc$smc3_treat - frac_treated_5_10 * smc_nmcp_mc$pop_5_10) / smc_nmcp_mc$pop_U5, max,0)
# smc_nmcp_mc$adjusted_nmcp_frac_U5_received_4 = sapply((smc_nmcp_mc$smc4_treat - frac_treated_5_10 * smc_nmcp_mc$pop_5_10) / smc_nmcp_mc$pop_U5, max,0)
# 



# ggplot(data = smc_nmcp_melt, aes(x = coverage_only_u5, y = adjusted_nmcp_frac_U5_received, fill = as.factor(round)))+
#   scale_fill_viridis_d( option = "D")+
#   # geom_violin(alpha=0.4, position = position_dodge(width = .75),size=1,color="black") +
#   # geom_boxplot(notch = TRUE,  outlier.size = -1, color="black",lwd=1.2, alpha = 0.7)+
#   geom_point( shape = 21,size=2, color="black",alpha=1)+
#   ylab(c("estimated number children U5 treated / number of children under 5"))  +
#   xlab(c("total number children treated / number of children under 5"))+ 
#   geom_abline(intercept = 0, slope=1) +
#   geom_abline(intercept = 1, slope=0, color=rgb(1,0,0,0.3), size=2, lty=2) +
#   theme_minimal()+
#   theme(legend.position = "none")
# ggsave(paste0(save_plot_filepath, '/2_smc_raw_adj_PNLP.png'), width=png_width, height=png_height, units='in', dpi=900)


# 
# # plot of old coverages used in sims versus adjusted coveraegs
# smc_old_sims = smc_old_sims %>% dplyr::select(admin_name, coverage_high_access, coverage_low_access, round, year)
# # get total coverage across high and low access children from input file used in old simulations
# smc_old_sims$coverage_old_sim = (smc_old_sims$coverage_high_access + smc_old_sims$coverage_low_access)/2
# smc_old_new = merge(smc_nmcp_melt, smc_old_sims, by=c('admin_name', 'round'))
# 
# ggplot(data = smc_old_new, aes(x = coverage_old_sim, y = adjusted_nmcp_frac_U5_received, fill = as.factor(admin_name)))+
#   geom_point( shape = 21,size=2, color="black",alpha=1)+
#   ylab(c("estimated number children U5 treated / number of children under 5"))  +
#   xlab(c("coverage estimate used in previous simulations"))+ 
#   geom_abline(intercept = 0, slope=1) +
#   geom_abline(intercept = 1, slope=0, color=rgb(1,0,0,0.3), size=2, lty=2)+
#   theme_minimal() +
#   theme(legend.position = "none")
# ggsave(paste0(save_plot_filepath, '/3_smc_old_adj_PNLP.png'), width=png_width, height=png_height, units='in', dpi=900)



# # --------- plots of MC data --------- #
# 
# # compare MC and NMEP values
# 
# 
# smc_mc_melt_U5 = melt(smc_mc, measure.vars = c("frac_U5_received_1", "frac_U5_received_2", "frac_U5_received_3", "frac_U5_received_4"),
#                       variable.name = "round", value.name = "U5_coverage")
# smc_mc_melt_U5$round = gsub('frac_U5_received_','',smc_mc_melt_U5$round)
# smc_mc_melt_O5 = melt(smc_mc, measure.vars = c("frac_O5_received_1", "frac_O5_received_2", "frac_O5_received_3", "frac_O5_received_4"),
#                       variable.name = "round", value.name = "O5_coverage")
# smc_mc_melt_O5$round = gsub('frac_O5_received_','',smc_mc_melt_O5$round)
# 
# ggplot(data = smc_mc_melt_U5, aes(x = round, y = U5_coverage, fill = as.factor(round)))+
#   scale_fill_viridis_d( option = "D")+
#   geom_violin(alpha=0.4, color='grey')+
#   # geom_violin(alpha=0.4, position = position_dodge(width = .75),size=1,color="black") +
#   # geom_boxplot(notch = TRUE,  outlier.size = -1, color="black",lwd=1.2, alpha = 0.7)+
#   geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=1)+
#   ylab(c("fraction of surveyed children \n U5 received SMC"))  +
#   xlab(c("SMC round"))+ 
#   ylim(0,1) +
#   theme_minimal()+ 
#   theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(), legend.position = "none")
# ggsave(paste0(save_plot_filepath, '/0_smc_raw_MC_U5.png'), width=png_width, height=png_height, units='in', dpi=900)
# 
# ggplot(data = smc_mc_melt_O5, aes(x = round, y = O5_coverage, fill = as.factor(round)))+
#   scale_fill_viridis_d( option = "D")+
#   geom_violin(alpha=0.4, color='grey')+
#   # geom_violin(alpha=0.4, position = position_dodge(width = .75),size=1,color="black") +
#   # geom_boxplot(notch = TRUE,  outlier.size = -1, color="black",lwd=1.2, alpha = 0.7)+
#   geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=1)+
#   ylab(c("fraction of surveyed children \n 5-10 received SMC"))  +
#   xlab(c("SMC round"))+ 
#   ylim(0,1) +
#   theme_minimal()+ 
#   theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(), legend.position = "none")
# ggsave(paste0(save_plot_filepath, '/0_smc_raw_MC_O5.png'), width=png_width, height=png_height, units='in', dpi=900)
# 
# 
# # --------- plots comparing MC and NMCP data --------- #
# # one plot per round in panel, region on x-axis, coverage on y-axis
# par(mfrow=c(2,2))
# for(rr in 1:4){
#   plot(NA, ylim=c(0,max(smc_nmcp_mc$doses_per_U5_4)), xlim=c(0,length(smc_nmcp_mc$doses_per_U5_4)), ylab=c('recorded doses per U5', 'or frac U5 received SMC'), xlab='Admin0', main=paste('round', rr), axes=FALSE)
#   axis(side=1, at=1:length(smc_nmcp_mc$State), labels = smc_nmcp_mc$State)
#   axis(side=2, at=sapply(seq(0,max(smc_nmcp_mc$doses_per_U5_4), length.out=4), round,2))
#   abline(h=1, lty=3, lwd=0.5, col='darkgrey')
#   # points for NMCP values
#   points(smc_nmcp_mc[,paste0('doses_per_U5_',rr)], col='red',pch=16)
#   # points for NMCP values, adjusted for 20% coverage in O5
#   points(smc_nmcp_mc[,paste0('adjusted_nmcp_frac_U5_received_',rr)], col='red',pch=22)
#   # assumed fraction of O5 with SMC for NMCP adjustments
#   points(rep(frac_treated_5_10, dim(smc_nmcp_mc)[1]), col='red',pch=24)
#   # points for U5 from MC
#   points(smc_nmcp_mc[,paste0('frac_U5_received_',rr)], col='blue',pch=15)
#   # points for O5 from MC
#   points(smc_nmcp_mc[,paste0('frac_O5_received_',rr)], col='blue',pch=17)
#   # points for the number of doses per U5 from MC
#   points(smc_nmcp_mc[,paste0('adjusted_mc_doses_per_U5_',rr)], col='blue',pch=1)
# }
# 
# 
# 
# # slide build sequence
# 
# # separate plots to show various comparisons
# png(paste0(save_plot_filepath, '/0_compare_smc_raw_PNLP_MC_r1.png'), width=png_width, height=png_height, units='in', res=900)
# par(mfrow=c(1,1))
# rr=1
# plot(NA, ylim=c(0,max(smc_nmcp_mc$doses_per_U5_4)), xlim=c(0,length(smc_nmcp_mc$doses_per_U5_4)), ylab=c('coverage'), xlab='Admin0', main=paste('round', rr), axes=FALSE)
# axis(side=1, at=1:length(smc_nmcp_mc$State), labels = smc_nmcp_mc$State)
# axis(side=2, at=sapply(seq(0,max(smc_nmcp_mc$doses_per_U5_4), length.out=4), round,2))
# abline(h=1, lty=3, lwd=2, col=rgb(0.4,0.4,0.4))
# # points for NMCP values
# points(smc_nmcp_mc[,paste0('doses_per_U5_',rr)], col='red',pch=16, cex=1.5)
# # points for U5 from MC
# points(smc_nmcp_mc[,paste0('frac_U5_received_',rr)], col='blue',pch=15, cex=1.5)
# # legend('bottomright', c('PNLP: children receiving SMC / pop U5', 'MC: fraction of surveyed U5 received SMC'), col=c('red', 'blue'), bty='n', pch=c(16,15), pt.cex=1.5)
# dev.off()
# 
# png(paste0(save_plot_filepath, '/1_compare_smc_raw_PNLP_MC_r1.png'), width=png_width, height=png_height, units='in', res=900)
# par(mfrow=c(1,1))
# rr=1
# plot(NA, ylim=c(0,max(smc_nmcp_mc$doses_per_U5_4)), xlim=c(0,length(smc_nmcp_mc$doses_per_U5_4)), ylab=c('coverage'), xlab='Admin0', main=paste('round', rr), axes=FALSE)
# axis(side=1, at=1:length(smc_nmcp_mc$State), labels = smc_nmcp_mc$State)
# axis(side=2, at=sapply(seq(0,max(smc_nmcp_mc$doses_per_U5_4), length.out=4), round,2))
# abline(h=1, lty=3, lwd=2, col=rgb(0.4,0.4,0.4))
# # points for NMCP values
# points(smc_nmcp_mc[,paste0('doses_per_U5_',rr)], col='red',pch=16, cex=1.5)
# # points for U5 from MC
# points(smc_nmcp_mc[,paste0('frac_U5_received_',rr)], col='blue',pch=15, cex=1.5)
# # points for O5 from MC
# points(smc_nmcp_mc[,paste0('frac_O5_received_',rr)], col='blue',pch=17)
# # legend('bottomright', c('PNLP: children receiving SMC / pop U5', 'MC: fraction of surveyed U5 received SMC', 'MC: fraction of surveyed O5 received SMC'), col=c('red', 'blue', 'blue'), bty='n', pch=c(16,15, 17), pt.cex=1.5)
# dev.off()
# 
# png(paste0(save_plot_filepath, '/2_compare_smc_raw_PNLP_adj_MC_r1.png'), width=png_width, height=png_height, units='in', res=900)
# par(mfrow=c(1,1))
# rr=1
# plot(NA, ylim=c(0,max(smc_nmcp_mc$doses_per_U5_4)), xlim=c(0,length(smc_nmcp_mc$doses_per_U5_4)), ylab=c('coverage'), xlab='Admin0', main=paste('round', rr), axes=FALSE)
# axis(side=1, at=1:length(smc_nmcp_mc$State), labels = smc_nmcp_mc$State)
# axis(side=2, at=sapply(seq(0,max(smc_nmcp_mc$doses_per_U5_4), length.out=4), round,2))
# abline(h=1, lty=3, lwd=2, col=rgb(0.4,0.4,0.4))
# # points for NMCP values
# points(smc_nmcp_mc[,paste0('doses_per_U5_',rr)], col='red',pch=16, cex=1.5)
# # points for the number of doses per U5 from MC
# points(smc_nmcp_mc[,paste0('adjusted_mc_doses_per_U5_',rr)], col='blue',pch=1, cex=1.4)
# # legend('bottomright', c('PNLP: children receiving SMC / pop U5', 'MC: children receiving SMC / pop U5'), col=c('red', 'blue'), bty='n', pch=c(16,1), pt.cex=1.5)
# dev.off()
# 
# 
# png(paste0(save_plot_filepath, '/3_compare_smc_adj_PNLP_raw_MC_r1.png'), width=png_width, height=png_height, units='in', res=900)
# par(mfrow=c(1,1))
# rr=1
# plot(NA, ylim=c(0,max(smc_nmcp_mc$doses_per_U5_4)), xlim=c(0,length(smc_nmcp_mc$doses_per_U5_4)), ylab=c('coverage'), xlab='Admin0', main=paste('round', rr), axes=FALSE)
# axis(side=1, at=1:length(smc_nmcp_mc$State), labels = smc_nmcp_mc$State)
# axis(side=2, at=sapply(seq(0,max(smc_nmcp_mc$doses_per_U5_4), length.out=4), round,2))
# abline(h=1, lty=3, lwd=2, col=rgb(0.4,0.4,0.4))
# # points for NMCP values, adjusted for 20% coverage in O5
# points(smc_nmcp_mc[,paste0('adjusted_nmcp_frac_U5_received_',rr)], col='red',pch=22, cex=1.5)
# # assumed fraction of O5 with SMC for NMCP adjustments
# points(rep(frac_treated_5_10, dim(smc_nmcp_mc)[1]), col='red',pch=24)
# # points for U5 from MC
# points(smc_nmcp_mc[,paste0('frac_U5_received_',rr)], col='blue',pch=15, cex=1.5)
# # points for O5 from MC
# points(smc_nmcp_mc[,paste0('frac_O5_received_',rr)], col='blue',pch=17)
# # legend('bottomleft', c('PNLP: est. fraction of U5 received SMC', 'PNLP: est. fraction of O5 received SMC', 'MC: fraction of surveyed U5 received SMC', 'MC: fraction of surveyed O5 received SMC'), col=c('red', 'red', 'blue', 'blue'), bty='n', pch=c(22,24,15,17), pt.cex=1.5)
# dev.off()
# 
# 
# png(paste0(save_plot_filepath, '/4_compare_smc_PNLP_MC.png'), width=png_width, height=png_height, units='in', res=900)
# par(mfrow=c(2,2))
# for(rr in 1:4){
#   plot(NA, ylim=c(0,max(smc_nmcp_mc$doses_per_U5_4)), xlim=c(0,length(smc_nmcp_mc$doses_per_U5_4)), ylab=c('recorded doses per U5', 'or frac U5 received SMC'), xlab='Admin0', main=paste('round', rr), axes=FALSE)
#   axis(side=1, at=1:length(smc_nmcp_mc$State), labels = smc_nmcp_mc$State) 
#   axis(side=2, at=sapply(seq(0,max(smc_nmcp_mc$doses_per_U5_4), length.out=4), round,2))
#   abline(h=1, lty=3, lwd=0.5, col='darkgrey')
#   # points for NMCP values
#   points(smc_nmcp_mc[,paste0('doses_per_U5_',rr)], col='red',pch=16)
#   # points for NMCP values, adjusted for 20% coverage in O5
#   points(smc_nmcp_mc[,paste0('adjusted_nmcp_frac_U5_received_',rr)], col='red',pch=22)
#   # assumed fraction of O5 with SMC for NMCP adjustments
#   points(rep(frac_treated_5_10, dim(smc_nmcp_mc)[1]), col='red',pch=24)
#   # points for U5 from MC
#   points(smc_nmcp_mc[,paste0('frac_U5_received_',rr)], col='blue',pch=15)
#   # points for O5 from MC
#   points(smc_nmcp_mc[,paste0('frac_O5_received_',rr)], col='blue',pch=17)
#   # points for the number of doses per U5 from MC
#   points(smc_nmcp_mc[,paste0('adjusted_mc_doses_per_U5_',rr)], col='blue',pch=1)
# }
# dev.off()
# 
# 
# 
# 
# 
# ###############################################
# # plot the post-adjusted SMC coverages used in simulations
# 
# png_width = 2
# png_height =2
# 
# # coverage for 2018; for high & low access as well as for all children in age group
# smc_new_sims_2018 = smc_new_sims[smc_new_sims[['year']]==2017,]
# # smc_new_sims_2018 = smc_new_sims_2020[smc_new_sims_2020[['year']]==2020,]  <- careful here! this is a lazy shortcut to look at 2020 coverage
# smc_new_sims_2018$total_U5_coverage = smc_new_sims_2018$coverage_high_access_U5 * smc_new_sims_2018$high_access_U5 +  smc_new_sims_2018$coverage_low_access_U5 * (1-smc_new_sims_2018$high_access_U5)
# smc_new_sims_2018$total_5_10_coverage = smc_new_sims_2018$coverage_high_access_5_10 * smc_new_sims_2018$high_access_5_10 +  smc_new_sims_2018$coverage_low_access_5_10 * (1-smc_new_sims_2018$high_access_5_10)
# 
# # children U5
# point_size = 1.5
# ggplot(data = smc_new_sims_2018, aes(x = round, y = coverage_high_access_U5, fill = as.factor(round)))+
#   scale_fill_viridis_d( option = "D")+
#   geom_violin(alpha=0.4, color='grey')+
#   # geom_violin(alpha=0.4, position = position_dodge(width = .75),size=1,color="black") +
#   # geom_boxplot(notch = TRUE,  outlier.size = -1, color="black",lwd=1.2, alpha = 0.7)+
#   geom_point( shape = 21,size=point_size, position = position_jitterdodge(), color="black",alpha=1)+
#   ggtitle('high-access U5')+
#   ylab(c("fraction received SMC"))  +
#   xlab(c("SMC round"))+ 
#   ylim(0,1)+
#   theme_minimal()+ 
#   theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(), legend.position = "none", plot.title = element_text(hjust = 0.5, size=12))
# ggsave(paste0(save_plot_filepath, '/smc_cov_new_sim_high_U5.png'), width=png_width, height=png_height, units='in', dpi=900)
# 
# ggplot(data = smc_new_sims_2018, aes(x = round, y = coverage_low_access_U5, fill = as.factor(round)))+
#   scale_fill_viridis_d( option = "D")+
#   geom_violin(alpha=0.4, color='grey')+
#   # geom_violin(alpha=0.4, position = position_dodge(width = .75),size=1,color="black") +
#   # geom_boxplot(notch = TRUE,  outlier.size = -1, color="black",lwd=1.2, alpha = 0.7)+
#   geom_point( shape = 21,size=point_size, position = position_jitterdodge(), color="black",alpha=1)+
#   ggtitle('low-access U5')+
#   ylab(c("fraction received SMC"))  +
#   xlab(c("SMC round"))+ 
#   ylim(0,1)+
#   theme_minimal()+ 
#   theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(), legend.position = "none", plot.title = element_text(hjust = 0.5, size=12))
# ggsave(paste0(save_plot_filepath, '/smc_cov_new_sim_low_U5.png'), width=png_width, height=png_height, units='in', dpi=900)
# 
# ggplot(data = smc_new_sims_2018, aes(x = round, y = total_U5_coverage, fill = as.factor(round)))+
#   scale_fill_viridis_d( option = "D")+
#   geom_violin(alpha=0.4, color='grey')+
#   # geom_violin(alpha=0.4, position = position_dodge(width = .75),size=1,color="black") +
#   # geom_boxplot(notch = TRUE,  outlier.size = -1, color="black",lwd=1.2, alpha = 0.7)+
#   geom_point( shape = 21,size=point_size, position = position_jitterdodge(), color="black",alpha=1)+
#   ggtitle('all children U5')+
#   ylab(c("fraction received SMC"))  +
#   xlab(c("SMC round"))+ 
#   ylim(0,1)+
#   theme_minimal()+ 
#   theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(), legend.position = "none", plot.title = element_text(hjust = 0.5, size=12))
# ggsave(paste0(save_plot_filepath, '/smc_cov_new_sim_all_U5.png'), width=png_width, height=png_height, units='in', dpi=900)
# 
# 
# 
# # children 5-10
# 
# ggplot(data = smc_new_sims_2018, aes(x = round, y = coverage_high_access_5_10, fill = as.factor(round)))+
#   scale_fill_viridis_d( option = "D")+
#   geom_violin(alpha=0.4, color='grey')+
#   # geom_violin(alpha=0.4, position = position_dodge(width = .75),size=1,color="black") +
#   # geom_boxplot(notch = TRUE,  outlier.size = -1, color="black",lwd=1.2, alpha = 0.7)+
#   geom_point( shape = 21,size=point_size, position = position_jitterdodge(), color="black",alpha=1)+
#   ggtitle('high-access 5-10')+
#   ylab(c("fraction received SMC"))  +
#   xlab(c("SMC round"))+ 
#   ylim(0,1)+
#   theme_minimal()+ 
#   theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(), legend.position = "none", plot.title = element_text(hjust = 0.5, size=12))
# ggsave(paste0(save_plot_filepath, '/smc_cov_new_sim_high_5_10.png'), width=png_width, height=png_height, units='in', dpi=900)
# 
# ggplot(data = smc_new_sims_2018, aes(x = round, y = coverage_low_access_5_10, fill = as.factor(round)))+
#   scale_fill_viridis_d( option = "D")+
#   geom_violin(alpha=0.4, color='grey')+
#   # geom_violin(alpha=0.4, position = position_dodge(width = .75),size=1,color="black") +
#   # geom_boxplot(notch = TRUE,  outlier.size = -1, color="black",lwd=1.2, alpha = 0.7)+
#   geom_point( shape = 21,size=point_size, position = position_jitterdodge(), color="black",alpha=1)+
#   ggtitle('low-access 5-10')+
#   ylab(c("fraction received SMC"))  +
#   xlab(c("SMC round"))+ 
#   ylim(0,1)+
#   theme_minimal()+ 
#   theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(), legend.position = "none", plot.title = element_text(hjust = 0.5, size=12))
# ggsave(paste0(save_plot_filepath, '/smc_cov_new_sim_low_5_10.png'), width=png_width, height=png_height, units='in', dpi=900)
# 
# ggplot(data = smc_new_sims_2018, aes(x = round, y = total_5_10_coverage, fill = as.factor(round)))+
#   scale_fill_viridis_d( option = "D")+
#   geom_violin(alpha=0.4, color='grey')+
#   # geom_violin(alpha=0.4, position = position_dodge(width = .75),size=1,color="black") +
#   # geom_boxplot(notch = TRUE,  outlier.size = -1, color="black",lwd=1.2, alpha = 0.7)+
#   geom_point( shape = 21,size=point_size, position = position_jitterdodge(), color="black",alpha=1)+
#   ggtitle('all children 5-10')+
#   ylab(c("fraction received SMC"))  +
#   xlab(c("SMC round"))+ 
#   ylim(0,1)+
#   theme_minimal()+ 
#   theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(), legend.position = "none", plot.title = element_text(hjust = 0.5, size=12))
# ggsave(paste0(save_plot_filepath, '/smc_cov_new_sim_all_5_10.png'), width=png_width, height=png_height, units='in', dpi=900)
# 
# 
# mean(smc_new_sims_2018$total_U5_coverage)
# mean(smc_new_sims_2018$total_5_10_coverage)
# 
# mean(smc_new_sims_2018$total_U5_coverage[smc_new_sims_2018$round==1])
# mean(smc_new_sims_2018$total_U5_coverage[smc_new_sims_2018$round==2])
# mean(smc_new_sims_2018$total_U5_coverage[smc_new_sims_2018$round==3])
# mean(smc_new_sims_2018$total_U5_coverage[smc_new_sims_2018$round==4])
# 
# mean(smc_new_sims_2018$total_5_10_coverage[smc_new_sims_2018$round==1])
# mean(smc_new_sims_2018$total_5_10_coverage[smc_new_sims_2018$round==2])
# mean(smc_new_sims_2018$total_5_10_coverage[smc_new_sims_2018$round==3])
# mean(smc_new_sims_2018$total_5_10_coverage[smc_new_sims_2018$round==4])





#' ########################################################
#' # high versus low access groups - old version
#' 
#' ff = function(fl, coverages, p_at_least_one){
#'   #' @description helper function used in solving for the fraction of the population in the low-access group (fl). Returns zero
#'   #'              when the correct value is reached.
#'   #' @param fl Fraction of the population in the lower-access group
#'   #' @param coverages Vector of length 4 with the coverages from each of the four rounds of SMC
#'   #' @param p_at_least_one The probability that an individual received at least one round of SMC
#'   print(fl)
#'   # return -999 if invalid value
#'   if((((coverages[1]-1+fl)/fl)<0) | (((coverages[2]-1+fl)/fl)<0) | (((coverages[3]-1+fl)/fl)<0) | (((coverages[4]-1+fl)/fl)<0)){
#'     zz=-999
#'   }else{
#'     # root for probability of receiving at least one SMC dose
#'     # equation equals zero for appropriate fl
#'     zz = 1 - p_at_least_one - ((1-(coverages[1]-1+fl)/fl) * 
#'                                  (1-(coverages[2]-1+fl)/fl) * 
#'                                  (1-(coverages[3]-1+fl)/fl) * 
#'                                  (1-(coverages[4]-1+fl)/fl)) * fl
#'   }
#'   return(zz)
#' }
#' 
#' coverage_U5 = c(0.892, 0.911,0.927, 0.861)
#' p_at_least_one_U5 = 0.939
#' coverage_O5 = c(0.1879, 0.1913, 0.1930, 0.1896)
#' p_at_least_one_O5 = 0.1963
#' xx_U5=uniroot(f=ff, lower = 0, upper = 1, coverage = coverage_U5, p_at_least_one=p_at_least_one_U5)$root
#' 1-xx_U5
#' xx_O5=uniroot(f=ff, lower = 0, upper = 1, coverage = coverage_O5, p_at_least_one=p_at_least_one_O5)$root
#' 1-xx_O5
#' 
#' coverage_test = c(0.5,0.5,0.5,0.5)
#' uniroot(f=ff, lower = 0, upper = 1, coverage = coverage_test, p_at_least_one=0.937)
#' 
#' # different regions
#' coverage_U5 = c(.88,94,96,91)
#' 1-uniroot(f=ff, lower = 0, upper = 1, coverage = coverage_U5, p_at_least_one=.978)$root  # high-access: 88%
#' coverage_U5 = c(.96,.96,.97,.98)
#' 1-uniroot(f=ff, lower = 0, upper = 1, coverage = coverage_U5, p_at_least_one=.99)$root  # high-access: 95%
#' coverage_U5 = c(.83,82,83,87)
#' 1-uniroot(f=ff, lower = 0, upper = 1, coverage = coverage_U5, p_at_least_one=.908)$root  # high-access: 83%
#' coverage_U5 = c(.95,.95,.97,.98)
#' 1-uniroot(f=ff, lower = 0, upper = 1, coverage = coverage_U5, p_at_least_one=.9908)$root  # high-access: 95%
#' coverage_U5 = c(.87,.90,.90,.91)
#' 1-uniroot(f=ff, lower = 0, upper = 1, coverage = coverage_U5, p_at_least_one=.96)$root  # high-access: 86%
#' 
#' 
#' 
#' # 2019
#' coverage_U5 = c(0.93, 0.94, 0.947, 0.893)
#' p_at_least_one_U5 = 0.969
#' xx_U5_2019=uniroot(f=ff, lower = 0, upper = 1, coverage = coverage_U5, p_at_least_one=p_at_least_one_U5)$root
#' xx_U5_2019
#' 
#' # probability of receiving all four doses
#' #  should be 82.8% for U5 and 18.2% for O5
#' fl=xx_U5
#' coverage = coverage_U5
#' (1-fl) * (1) + (fl) * (prod((coverage-1+fl)/fl))
#' fl=xx_O5
#' coverage = coverage_O5
#' (1-fl) * (1) + (fl) * (prod((coverage-1+fl)/fl))
#' 
#' 
#' 


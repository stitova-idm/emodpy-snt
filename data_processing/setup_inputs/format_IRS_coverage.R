# format_IRS_coverage.R

library(ggplot2)
library(dplyr)


# parameters for testing
# data_dir = 'C:/Users/mambrose/Dropbox (IDM)/Malaria Team Folder/data'
# hbhi_dir = 'C:/Users/mambrose/Dropbox (IDM)/Malaria Team Folder/projects/burundi_hbhi'
# irs_who_filename = paste0(data_dir, '/Burundi/WHO/Interventions/IRS_2016-2020_v2.csv')
# ds_pop_df_filename = paste(hbhi_dir, '/admin_pop_archetype.csv', sep='')
# mean_household_size=4.8
# irs_first_round_month=3.5
# irs_second_round_month=10.5
# sim_start_year=2010
# sim_end_year=2020
# burn_start_year=1961
# calib_start_year=2011
# max_effective_coverage_irs=0.75
# create_plots = FALSE



format_irs_inputs = function(irs_who_filename, ds_pop_df_filename,mean_household_size=4.8, irs_first_round_month=3.5,irs_second_round_month=10.5,  
                             sim_start_year=2010, burn_start_year=1961, calib_start_year=2011, 
                             effective_irs_coverage_multiplier=0.78, max_effective_coverage_irs=0.75, num_samples=50, irs_kill_seed_multipliers=1, create_plots = FALSE){
  # read in csv with households sprayed, number of rounds, and insecticide used in each admin
  irs_who = read.csv(irs_who_filename)
  irs_who$individuals_covered = irs_who$irs_hhs_num * mean_household_size
  irs_who$individuals_covered_round2 = irs_who$irs_hhs_num2 * mean_household_size
  admin_pop = read.csv(ds_pop_df_filename)
  
  irs_pop_df = merge(irs_who, admin_pop[,which(colnames(admin_pop) %in% c('admin_name', 'pop_size','cluster_id', 'seasonality_archetype'))], by.x='adm2', by.y='admin_name', all.x=TRUE)
  irs_pop_df$coverage_r1 = irs_pop_df$individuals_covered / irs_pop_df$pop_size
  irs_pop_df$coverage_r2 = irs_pop_df$individuals_covered_round2 / irs_pop_df$pop_size
  irs_pop_df$max_r1_r2 = apply(irs_pop_df[,which(colnames(irs_pop_df) %in% c('coverage_r1', 'coverage_r2'))], 1, max, na.rm=TRUE)
  irs_pop_df$max_r1_r2[irs_pop_df$max_r1_r2<0] = NA
  
  # get rescale value so that coverage in IRS districts in 2017 is 0.78
  irs_2017 = irs_pop_df[irs_pop_df$year==2017,]
  reported_coverage = sum(irs_2017$pop_size* irs_2017$max_r1_r2)/sum(irs_2017$pop_size)
  coverage_rescalar = 0.78/reported_coverage
  
  # get rescaled coverages in each round, capped at the maximum effective coverage
  irs_pop_df$effect_coverage_r1 = sapply((irs_pop_df$individuals_covered / irs_pop_df$pop_size * coverage_rescalar), min, max_effective_coverage_irs)
  irs_pop_df$effect_coverage_r2 = sapply((irs_pop_df$individuals_covered_round2 / irs_pop_df$pop_size * coverage_rescalar), min, max_effective_coverage_irs)
  irs_pop_df$max_r1_r2_rescaled = irs_pop_df$max_r1_r2 * coverage_rescalar
  irs_pop_df$max_r1_r2_rescaled_effective = sapply(irs_pop_df$max_r1_r2_rescaled, min, max_effective_coverage_irs)
  
  if(create_plots){
    # raw coverages
    ggplot(irs_pop_df, aes(x=year, y=max_r1_r2, col=adm2)) + 
    geom_line() +
    ylab('coverage') + 
    geom_point() + 
    theme_classic()
    # after 2017 rescaling
    ggplot(irs_pop_df, aes(x=year, y=max_r1_r2_rescaled, col=adm2)) + 
      geom_line() +
      ylab('coverage, rescaled') + 
      ylim(0,1.5) +
      geom_point() + 
      theme_classic()
    # after setting maximum effective coverage
    ggplot(irs_pop_df, aes(x=year, y=max_r1_r2_rescaled_effective, col=adm2)) + 
      geom_line() +
      ylab('effective coverage') + 
      ylim(0,1.5) +
      geom_point() + 
      theme_classic()
  }
  
  # split round 2 into separate rows
  irs_pop_df_r1 = irs_pop_df
  irs_pop_df_r1$round = 1
  irs_pop_df_r1$month = irs_first_round_month
  irs_pop_df_r1$effective_coverage = irs_pop_df_r1$effect_coverage_r1
  irs_pop_df_r2 = irs_pop_df[irs_pop_df$irs_rounds == 2,]
  irs_pop_df_r2$round = 2
  irs_pop_df_r2$month = irs_second_round_month
  irs_pop_df_r2$effective_coverage = irs_pop_df_r2$effect_coverage_r2
  irs_pop_df2 = rbind(irs_pop_df_r1, irs_pop_df_r2)
  irs_pop_df2 = irs_pop_df2[,which(colnames(irs_pop_df2) %in% c('adm2', 'year', 'month', 'round','insecticide','effective_coverage','pop_size','cluster_id','seasonality_archetype'))]
  # replace NAs with the average of all non-na coverages in that year in the dataset (if none in that year, use previous year, otherwise use average across all years in dataset)
  for(yy in unique(irs_pop_df2$year)){
    mean_yy = mean(irs_pop_df2$effective_coverage[which(irs_pop_df2$year == yy)], na.rm=TRUE)
    if(!is.na(mean_yy)){
      replacement_value_yy = mean_yy
    } else if(any(irs_pop_df2$year == (yy-1))){
        mean_yy_minus1 = mean(irs_pop_df2$effective_coverage[which(irs_pop_df2$year == (yy-1))], na.rm=TRUE)
        if(!is.na(mean_yy_minus1)){
          replacement_value_yy = mean_yy_minus1
        } else{
          # if there are no coverage values for this year, take average across all IRS years
          replacement_value_yy = mean(irs_pop_df2$effective_coverage, na.rm=TRUE)
        } 
    } else{
      # if there are no coverage values for this year, take average across all IRS years
      replacement_value_yy = mean(irs_pop_df2$effective_coverage, na.rm=TRUE)
    } 
    irs_pop_df2$effective_coverage[intersect(which(irs_pop_df2$year == yy), which(is.na(irs_pop_df2$effective_coverage)))] = replacement_value_yy
  }
  
  # add insecticide initial kill, block duration, and decay
  insecticide_names = c('bendiocarb','bendiocarb - fludora fusion', 'malathion', 'bendiocarb & malathion', 'actellic', 'fludora fusion', 'sumithion', '')  # if none listed, assume bendiocarb
  insecticide_kills = c(0.6, 0.6, 0.85, 0.75, 0.85, 0.6, 0.75, 0.6)
  insecticide_halflife = c(70, 100, 255, 160, 255, 100, 160, 70)
  insecticide_mean = insecticide_halflife / log(2)
  irs_pop_df2$initial_kill = NA
  irs_pop_df2$mean_duration = NA
  # if insecticide name not recognized, replace with ''
  irs_pop_df2$insecticide[!(irs_pop_df2$insecticide %in% insecticide_names)] = ''
  for(ii in 1:length(insecticide_names)){
    irs_pop_df2$initial_kill[irs_pop_df2$insecticide == insecticide_names[ii]] = insecticide_kills[ii]
    irs_pop_df2$mean_duration[irs_pop_df2$insecticide == insecticide_names[ii]] = insecticide_mean[ii]
  }
  
  
  # calibrations - separate entry for each archetype
  irs_pop_df2$num_covered = irs_pop_df2$pop_size * irs_pop_df2$effective_coverage
  irs_input_calib = irs_pop_df2 %>% group_by(year, cluster_id, seasonality_archetype, round, month) %>% 
    summarise(num_covered_arch = sum(num_covered),
              initial_kill = weighted.mean(initial_kill, pop_size),
              mean_duration = weighted.mean(mean_duration, pop_size))
  # get population sizes in each archetype
  ach_pop = admin_pop %>% group_by(seasonality_archetype, cluster_id) %>%
    summarise(pop_size = sum(pop_size))
  irs_input_calib2 = merge(irs_input_calib, ach_pop, by=c('seasonality_archetype', 'cluster_id'), all.x=TRUE)
  irs_input_calib2$effective_coverage = irs_input_calib2$num_covered_arch / irs_input_calib2$pop_size
  irs_input_calib2$simday = (irs_input_calib2$year - calib_start_year) * 365 + round(irs_input_calib2$month * 30)
  colnames(irs_input_calib2)[which(colnames(irs_input_calib2)=='seasonality_archetype')] = 'admin_name'
  irs_input_calib2 = irs_input_calib2[,which(colnames(irs_input_calib2) %in% c('admin_name', 'year','round','initial_kill','mean_duration', 'effective_coverage','simday'))]
  # rescale kill rates based on first value in irs_kill_seed_multipliers 
  irs_input_calib2$initial_kill = irs_input_calib2$initial_kill * irs_kill_seed_multipliers[1]
  
  
  # main simulations - separate entry for each admin that receives IRS
  # rescale kill rates based on irs_kill_seed_multipliers 
  if(length(irs_kill_seed_multipliers)>1){  # determine whether there will be different values for different seeds
    num_repeats = ceiling((num_samples+1) / length(irs_kill_seed_multipliers))
    multiplier_each_seed = rep(irs_kill_seed_multipliers, num_repeats)[1:(num_samples+1)]
    seeds_df = data.frame('seed'=1:(num_samples+1), 'kill_multiplier'=multiplier_each_seed)
    irs_pop_df3 = merge(irs_pop_df2, seeds_df, all=TRUE)
    irs_pop_df3$initial_kill = irs_pop_df3$initial_kill * irs_pop_df3$kill_multiplier
  } else if(length(irs_kill_seed_multipliers) == 1){
    irs_pop_df3=irs_pop_df2
    irs_pop_df3$initial_kill = irs_pop_df3$initial_kill * irs_kill_seed_multipliers
  } else{
    warning('PROBLEM: need to specify a valid parameter value for irs_kill_seed_multipliers.')
  }
  
  # format for simulation
  irs_input_main = irs_pop_df3
  irs_input_main$simday = (irs_input_main$year - sim_start_year) * 365 + round(irs_input_main$month * 30)
  colnames(irs_input_main)[colnames(irs_input_main) == 'adm2'] = 'admin_name'
  irs_input_main = irs_input_main[,which(colnames(irs_input_main) %in% c('admin_name', 'year','round','initial_kill','mean_duration', 'effective_coverage','simday', 'insecticide', 'seed'))]
  # rescale coverage to get 'effective coverage' accounting for building/modifying dwellings, cleaning walls, etc.
  irs_input_main$effective_coverage = irs_input_main$effective_coverage * effective_irs_coverage_multiplier
  irs_input_calib2$effective_coverage = irs_input_calib2$effective_coverage * effective_irs_coverage_multiplier
  # save csv files for simulation input
  write.csv(irs_input_calib2, paste0(hbhi_dir, '/simulation_inputs/interventions_calib/IRS_season_calib.csv'), row.names = FALSE)
  write.csv(irs_input_main, paste0(hbhi_dir, '/simulation_inputs/interventions_2010_toPresent/IRS_2010_toPresent.csv'), row.names = FALSE)
}



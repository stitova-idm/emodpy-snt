# plot impact of IPTp on population-level PfPR


library(data.table)
library(tidyverse)
library(lubridate)
library(reshape2)
library(plyr)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
#                     set filepaths, read in data, load needed functions
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###

##### User-specified information on scenarios to use: ############################################################################################### 
user = Sys.getenv("USERNAME")
user_path = file.path("C:/Users",user)

# set name of country to run
country_name = 'Nigeria'
if(country_name =='Burkina'){
  box_hbhi_filepath = paste(user_path, '/Box/hbhi_burkina', sep='')
  # filepath for simulation output
  # sim_output_filepath = paste(box_hbhi_filepath, '/simulation_output/2010_to_2020/_v18', sep='')
  sim_output_filepath = paste(box_hbhi_filepath, '/simulation_output/2020_to_2025/BF ms projection v2', sep='')
  sim_output_base_filepath_2010_allInter = paste(box_hbhi_filepath, '/simulation_output/2010_to_2020/_v18/BF 2010_2020 allInterventions', sep='')
  simname_stem = 'BF NSP projection '
} else if (country_name=='Nigeria'){
  box_hbhi_filepath = paste(user_path, '/Box/hbhi_nigeria', sep='')
  # filepath for simulation output
  sim_output_filepath = paste(box_hbhi_filepath, '/simulation_output/2010_to_2020_v10', sep='')
  # sim_output_filepath = paste(box_hbhi_filepath, '/simulation_output/2020_to_2030_v2', sep='')
  sim_output_base_filepath_2010_allInter = paste(box_hbhi_filepath, '/simulation_output/2010_to_2020_v10/NGA 2010-20 burnin_hs+itn+smc', sep='')
} else{
  warning('Country name not supported')
}

# set timespan of simulations
first_year = 2010
last_year = 2019
# first_year = 2020
# last_year = 2025

# filepath for scripts
script_base_filepath = paste(user_path, '/Documents/hbhi-nigeria/simulation/IPTp_mortality_postprocessing', sep='')
# load needed functions
source(paste0(script_base_filepath, '/malariaInPregnancy_functions.R'))
source(paste0(script_base_filepath, '/get_MiP_protection_coverage.R'))


# filename for estimated IPTp coverage
iptp_estimates_filename = paste(box_hbhi_filepath,'/IPTp/estimated_2010_2020_IPTp_each_DS.csv', sep='')
iptp_estimates_ci_l_filename = paste(box_hbhi_filepath,'/IPTp/estimated_ci_l_2010_2020_IPTp_each_DS.csv', sep='')
iptp_estimates_ci_u_filename = paste(box_hbhi_filepath,'/IPTp/estimated_ci_u_2010_2020_IPTp_each_DS.csv', sep='')
iptp_project_trajectory_filename = paste(box_hbhi_filepath,'/IPTp/projected_trajectory_2020_2030_IPTp_each_DS.csv', sep='')

# filename for estimated number of IPTp doses taken by IPTp recipients in each year
iptp_dose_number_filename = paste(box_hbhi_filepath,'/IPTp/estimated_2010_2020_num_doses.csv', sep='')
# read in estimated IPTp coverage
iptp_coverage_df = read.csv(iptp_estimates_filename, as.is=TRUE)[,-1]
ds_names = iptp_coverage_df$DS
# read in the number of IPTp doses taken by people who receive IPTp
iptp_dose_number = read.csv(iptp_dose_number_filename, row.names=1)


# get filenames and information for all of the scenarios currently being processed
scenario_adjusmtnet_info_filepath = paste(sim_output_filepath,  '/scenario_adjustment_info.csv', sep='')
if(file.exists(scenario_adjusmtnet_info_filepath)){
  scenario_adjustment_info = read.csv(scenario_adjusmtnet_info_filepath, as.is=TRUE)
  # check whether scenario name column already exists. if not, create it
  if(!('ScenarioName' %in% colnames(scenario_adjustment_info))){
    scenario_adjustment_info$ScenarioName = paste(simname_stem, scenario_adjustment_info$Scenario_no)
  }
} else{
  warning('csv file containing information about set of simulations not found')
}



for(i_scen in 1:length(scenario_adjustment_info$ScenarioName)){
  
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  #               read in scenario information and simulation output files
  ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ###
  
  # get filename and coverage information for this scenario
  scenarioName_cur = scenario_adjustment_info$ScenarioName[i_scen]
  scenarioName_IPTi = scenario_adjustment_info$IPTi_scenario_name[i_scen]
  sim_output_base_filepath = paste(sim_output_filepath, '/', scenarioName_cur, sep='')
  future_projection_flag = scenario_adjustment_info$future_projection_flag[i_scen]
  coverage_string = scenario_adjustment_info$IPTp_string[i_scen]
  name_for_pm_still_iptp = 'still100'
  
  # filename for simulation output on number of new infections in each month by age group
  sim_output_filename = paste(sim_output_base_filepath, "/newInfections_PfPR_cases_monthly_byAgeGroup.csv", sep='')
  sim_output_allAgeMonthly_filename_cur = paste(sim_output_base_filepath, "/All_Age_Monthly_Cases.csv", sep='')
  sim_output_allAgeMonthly_filename_MiP = gsub( '.csv', paste('_withMiP', coverage_string, '_', name_for_pm_still_iptp,'.csv', sep=''), sim_output_allAgeMonthly_filename_cur)
  
  
  # check whether output for this scenario exists
  if(file.exists(sim_output_filename)){
    print(paste('Starting scenario: ', scenarioName_cur, ' (', i_scen, ' out of ', length(scenario_adjustment_info$ScenarioName),')', sep=''))

    # get the IPTp coverages for this scenario (will be added to plot)
    IPTp_coverages = get_IPTp_coverages(iptp_estimates_filename=iptp_estimates_filename, iptp_dose_number_filename=iptp_dose_number_filename, 
                                        future_projection_flag=future_projection_flag, 
                                        first_year=first_year, last_year=last_year, coverage_string=coverage_string,
                                        iptp_estimates_ci_l_filename=iptp_estimates_ci_l_filename, iptp_estimates_ci_u_filename=iptp_estimates_ci_u_filename,
                                        iptp_project_trajectory_filename=iptp_project_trajectory_filename)
    iptp_coverage_df = IPTp_coverages[[1]]
    iptp_dose_number = IPTp_coverages[[2]]
    ds_names = IPTp_coverages[[3]]
    
    # among people who receive IPTp in a given year, what fraction gets each number of doses
    frac_iptp_years = first_year:last_year
    # read in from saved file (based on DHS/MIS data)
    frac_iptp_1 = as.numeric(iptp_dose_number[1,])
    frac_iptp_2 = as.numeric(iptp_dose_number[2,])
    frac_iptp_3 = as.numeric(iptp_dose_number[3,])
    
    # read in simulation file  post-adjustment (but different format from pre-adjustmant file)
    sim_output_MiP_adjusted = fread(sim_output_allAgeMonthly_filename_MiP)
    if('LGA' %in% colnames(sim_output_MiP_adjusted)){
      sim_ds_colname = 'LGA'
    } else{
      sim_ds_colname = 'DS_Name'
    }

    # iterate through runs, adding entries for new columns in the relevant run & month & ds row with the unadjusted PfPR and with IPTp coverage
    all_runs = unique(sim_output_MiP_adjusted$Run_Number)
    for(rr in 1:length(all_runs)){
      print(paste('currently adding new MiP column entries for run index', rr))
      
      cur_run = all_runs[rr]

      # calculate overall PfPR for population without IPTp
      # read in pregnancy-related values for this simulation
      popUnder15_size = fread(paste(sim_output_base_filepath, '/mortality_LBW/popUnder15_size_run',cur_run, '.csv', sep=''))
      pop1530_size = fread(paste(sim_output_base_filepath, '/mortality_LBW/pop1530_size_run',cur_run, '.csv', sep=''))
      pop3050_size = fread(paste(sim_output_base_filepath, '/mortality_LBW/pop3050_size_run',cur_run, '.csv', sep=''))
      popOver50_size = fread(paste(sim_output_base_filepath, '/mortality_LBW/popOver50_size_run',cur_run, '.csv', sep=''))
      pfpr_under15 = fread(paste(sim_output_base_filepath, '/mortality_LBW/pfpr_under15_run',cur_run, '.csv', sep=''))
      pfpr_1530 = fread(paste(sim_output_base_filepath, '/mortality_LBW/pfpr_1530_run',cur_run, '.csv', sep=''))
      pfpr_3050 = fread(paste(sim_output_base_filepath, '/mortality_LBW/pfpr_3050_run',cur_run, '.csv', sep=''))
      pfpr_over50 = fread(paste(sim_output_base_filepath, '/mortality_LBW/pfpr_over50_run',cur_run, '.csv', sep=''))
      
      pfpr_allAges_noIPTp = (pfpr_under15 * popUnder15_size + 
                             pfpr_1530 * pop1530_size + 
                             pfpr_3050 * pop3050_size + 
                             pfpr_over50 * popOver50_size) / 
                            (popUnder15_size + pop1530_size + pop3050_size + popOver50_size)
      
      
      # convert from wide to long format
      pfpr_allAges_noIPTp_long = melt(pfpr_allAges_noIPTp, id.vars=c('year', 'month'))
      colnames(pfpr_allAges_noIPTp_long)[colnames(pfpr_allAges_noIPTp_long)=='variable'] =  sim_ds_colname
      colnames(pfpr_allAges_noIPTp_long)[colnames(pfpr_allAges_noIPTp_long)=='value'] =  'unadjusted_mic_PfPR'
      pfpr_allAges_noIPTp_long$Run_Number = cur_run
      
      # add the values from this run to the cumulative data table
      if(rr == 1){
        pfpr_allAges_noIPTp_all = pfpr_allAges_noIPTp_long
      } else{
        pfpr_allAges_noIPTp_all = rbind(pfpr_allAges_noIPTp_all, pfpr_allAges_noIPTp_long)
      }
    } # finish iterating through runs
    
    # change admin names to all caps for both dataframes to allow merging
    sim_output_MiP_adjusted[[sim_ds_colname]] = toupper(sim_output_MiP_adjusted[[sim_ds_colname]])
    pfpr_allAges_noIPTp_all[[sim_ds_colname]] = toupper(pfpr_allAges_noIPTp_all[[sim_ds_colname]])
    
    # merge unadjusted values into adjusted dataframe
    compare_PfPR = merge(sim_output_MiP_adjusted, pfpr_allAges_noIPTp_all, by=c('year','month','Run_Number',sim_ds_colname))
    
    # add column for IPTp coverage in that admin/year
    iptp_coverage_long = melt(iptp_coverage_df, id.vars=c('DS'))
    colnames(iptp_coverage_long)[colnames(iptp_coverage_long)=='variable'] =  'year'
    colnames(iptp_coverage_long)[colnames(iptp_coverage_long)=='value'] =  'IPTp_coverage'
    colnames(iptp_coverage_long)[colnames(iptp_coverage_long)=='DS'] =  sim_ds_colname
    iptp_coverage_long[[sim_ds_colname]] = toupper(iptp_coverage_long[[sim_ds_colname]])
    iptp_coverage_long$year = gsub('X','', iptp_coverage_long$year)
    iptp_coverage_long$year = as.integer(iptp_coverage_long$year)
    
    compare_PfPR_2 = merge(compare_PfPR, iptp_coverage_long, by=c('year',sim_ds_colname),  sort=FALSE, all=TRUE)
    # compare_PfPR_2 = join(compare_PfPR, iptp_coverage_long, by=c('year',sim_ds_colname), type='left')
    
    fwrite(compare_PfPR, paste(sim_output_base_filepath, "/All_Age_Monthly_Cases_with_without_IPTp_adjustd_PfR.csv", sep=''))
  }
}

plot(compare_PfPR$unadjusted_mic_PfPR[5100:5200], compare_PfPR$PfPR_MiP_adjusted[5100:5200])
      
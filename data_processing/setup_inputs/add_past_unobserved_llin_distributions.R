# add_past_unobserved_llin_distributions.R

# sometimes we know there were mass distribution during certain years, but we have no information on coverage obtained, so we use coverage from previous  mass distribution in simulations
# given information about the coverage and timing of past distributions, this function adds campaigns to the input file. It assumes the same timing for all admins, so is not appropriate for
#    countries that have mass campaigns in different years at a subnational level

add_past_unobserved_llin_distributions = function(past_itn_mass_coverage_filepath,
                                                  sim_start_year = 2010,
                                                  additional_itn_years = c(2017, 2020), 
                                                  additional_itn_months = c(10, 7),
                                                  copied_itn_coverage_year = 2014  # assume coverage for unobserved mass LLIN distribution years is the same as coverage for this LLIN distribution
                                                  ){

  # read in ITN coverage csv for observed years
  itn_cov_df = read.csv(past_itn_mass_coverage_filepath)
  
  # copy rows matching the selected duplicated year
  for(yy in 1:length(additional_itn_years)){
    itn_cov_df_yy = itn_cov_df[itn_cov_df$dist_year == copied_itn_coverage_year,]
    itn_cov_df_yy$dist_year = additional_itn_years[yy]
    itn_cov_df_yy$month = additional_itn_months[yy]
    itn_cov_df_yy$simday = (itn_cov_df_yy$dist_year - sim_start_year) * 365 + round(itn_cov_df_yy$month*30)
    itn_cov_df = rbind(itn_cov_df, itn_cov_df_yy)
  }
  write.csv(itn_cov_df,paste0(hbhi_dir, '/simulation_inputs/intermediate_files/ITN_coverage/itn_mass_coverages_2010_2020.csv'), row.names=FALSE)  # intermediate files that only have coverage information
}






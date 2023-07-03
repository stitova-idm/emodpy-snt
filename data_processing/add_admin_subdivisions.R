# add_admin_subdivisions.R

# sometimes, the admins in a country are divided further to create additional admins
#   for example, in Burundi, three new DS have been created between the 2020 and 2023 analysis periods
#   to account for this in simulations, we re-run everything using the new set of admins
#   but data is often reported with the original admins. Therefore, this function can be used to 
#   add additional data rows for the new admins using the values from the original admin
#   for the years in which they were combined

add_rows_for_new_admins = function(filename_new_admins, intervention_filename, intervention_file_admin_name='admin_name'){
  new_admin_df = read.csv(filename_new_admins)
  intervention_df = read.csv(intervention_filename)
  original_df = intervention_df
  changes_made = FALSE
  
  # for every new admin, check whether it is already in the intervention file. If not, duplicate all rows from the admin it is splitting from (but replace with the new admin name) 
  for(aa in 1:nrow(new_admin_df)){
    cur_admin = new_admin_df$new_admin[aa]
    origin_admin = new_admin_df$original_admin[aa]
    if(!(cur_admin %in% intervention_df[[intervention_file_admin_name]])){
      if(origin_admin %in% intervention_df[[intervention_file_admin_name]]){
        # create copy of original admin's rows to assign to the new/split admin
        cur_copy = intervention_df[which(intervention_df[[intervention_file_admin_name]] == origin_admin),]
        cur_copy[[intervention_file_admin_name]] = cur_admin
        intervention_df = rbind(intervention_df, cur_copy)
        changes_made = TRUE
      } else{
        warning(paste0('For new admin ', cur_admin, ', the original / source admin ', origin_admin, ' was not found.'))
      }
    }
  }
  if(changes_made){
    # save the original version of the file
    write.csv(original_df, file=gsub('.csv', '_original.csv', intervention_filename), row.names=FALSE)
    # replace the intervention file with the one updated to have the newly added admins
    write.csv(intervention_df, file=intervention_filename, row.names=FALSE)
  }
}


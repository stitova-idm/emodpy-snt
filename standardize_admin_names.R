# standardize_admin_names.R
# November 2022
# contact: mambrose
#
# Goal: the names of LGAs are often different across different files and we want to make sure that they are used consistently
#    This function takes a 'target' naming system and a 'origin' set of names and matches the origin LGA names with the target names.
#    It then updates the origin names so that the same LGA names can be used consistently across all files.
library(stringr)

create_reference_name_match = function(lga_name){
  lga_name = str_replace_all(lga_name, pattern=' ', replacement='-')
  lga_name = str_replace_all(lga_name, pattern='/', replacement='-')
  lga_name = str_replace_all(lga_name, pattern='_', replacement='-')
  lga_name = str_replace_all(lga_name, pattern="'", replacement='')
  lga_name = toupper(lga_name)
  
  # first value (the name) is the one that is replaced by the second value
  replace_list = list(
                      # NGA 
                      'GANYE' = 'GANAYE',
                      'GIREI' = 'GIRERI',
                      'DAMBAM' = 'DAMBAN',
                      'MARKURDI' = 'MAKURDI',
                      'MUNICIPAL-AREA-COUNCIL' = 'ABUJA-MUNICIPAL',
                      'SHONGOM' = 'SHOMGOM',
                      # 'BIRNIN' = 'BIRINIWA', 'BIRNI-KUDU'?
                      'BIRNIN-KUDU' = 'BIRNI-KUDU',
                      'BIRNIWA' = 'BIRINIWA',
                      'KIRI-KASAMA' = 'KIRI-KASAMMA',
                      'MALAM-MADURI' = 'MALAM-MADORI',
                      'SULE-TANKAKAR' = 'SULE-TANKARKAR',
                      'MAKARFI' = 'MARKAFI',
                      'ZANGON-KATAF' = 'ZANGO-KATAF',
                      'DANBATTA' = 'DAMBATTA',
                      'TAKAI' = 'TAKALI',
                      'UNGONGO' = 'UNGOGO',
                      'DAN' = 'DAN-MUSA',
                      'DUTSIN' = 'DUTSIN-MA',
                      'MATAZU' = 'MATAZUU',
                      'BASSA' = 'BASSA1',
                      'IFELODUN' = 'IFELODUN1',
                      'IREPODUN' = 'IREPODUN1',
                      'PATIGI' = 'PATEGI',
                      'NASARAWA' = 'NASARAWA2',
                      'OBI' = 'OBI2',
                      'MUNYA' = 'MUYA',
                      'PAIKORO' = 'PAILORO',
                      'SURULERE' = 'SURULERE2',
                      'BARKIN-LADI' = 'BARIKIN-LADI',
                      'DANGE-SHUNI' = 'DANGE-SHNSI',
                      'GWADABAWA' = 'GAWABAWA',
                      'WAMAKKO' = 'WAMAKO',
                      'ARDO' = 'ARDO-KOLA',
                      'KARIM' = 'KURMI',
                      'BARDE' = 'BADE',
                      'BORSARI' = 'BURSARI',
                      'TARMUWA' = 'TARMUA',
                      # BDI
                      'BUJUMBURA-CENTRE' = 'ZONE-CENTRE',
                      'BUJUMBURA-NORD' = 'ZONE-NORD',
                      'BUJUMBURA-SUD' = 'ZONE-SUD'
                      # 'Bukinanyana' = 'Cibitoke',  # Bukinanyana is a commune in Cibitoke province - may now be a new DS
                      # 'Gisuru' = 'Ruyigi',  # Gisuru is a commune in Ruyigi province - may now be a new DS
                      # 'Rutovu' = 'Bururi'  # Rutovu is a commune in Bururi province - may now be a new DS

  )
  if(lga_name %in% names(replace_list)){
    lga_name = replace_list[lga_name][[1]]
  }
  return(lga_name)
}

standardize_admin_names_in_df = function(target_names_df, origin_names_df, target_names_col='admin_name', origin_names_col='admin_name'){

  target_names_df$matched_name = sapply(target_names_df[[target_names_col]], create_reference_name_match)
  origin_names_df$matched_name = sapply(origin_names_df[[origin_names_col]], create_reference_name_match)
  
  # check whether all names from the origin source are now matched to one of the target names
  if(!all(origin_names_df$matched_name %in% target_names_df$matched_name)){
    warning('Some of the source admin names could not be matched with a target admin name')
    View(distinct(origin_names_df[which(!(origin_names_df$matched_name %in% target_names_df$matched_name)), c('matched_name', 'State')]))
    View(target_names_df[,c('matched_name', 'State')])
  } 
  if('data.table' %in% class(origin_names_df)){  # changes how indexing works for columns in dataframe
    origin_names_df = as.data.frame(origin_names_df)
    was_data_table=TRUE
  } else was_data_table = FALSE
  if('data.table' %in% class(target_names_df)){  
    target_names_df = as.data.frame(target_names_df)
  }
  # remove the original admin name column from the origin dataframe, and also the target name column if applicable (for merging)
  origin_names_df = origin_names_df[,-(which(colnames(origin_names_df)==origin_names_col))]
  if(target_names_col %in% colnames(origin_names_df)){
    origin_names_df = origin_names_df[,-which(colnames(origin_names_df)==target_names_col)]
  }

  # add the updated admin names to the dataframe, under the original admin column name
  target_names_df =  target_names_df[,c('matched_name', target_names_col)]
  updated_names_df = merge(origin_names_df, target_names_df, all.x=TRUE, by='matched_name')
  colnames(updated_names_df)[colnames(updated_names_df)==target_names_col] = origin_names_col
  
  if(was_data_table) updated_names_df = as.data.table(updated_names_df)
  return(updated_names_df)
}



standardize_admin_names_in_vector = function(target_names, origin_names){
  
  target_matched_name = sapply(target_names, create_reference_name_match)
  origin_matched_name = sapply(origin_names, create_reference_name_match)
  
  # check whether all names from the origin source are now matched to one of the target names
  if(!all(origin_matched_name %in% target_matched_name)){
    warning('Some of the source admin names could not be matched with a target admin name')
  } 
  
  # create dataframes to merge
  target_df = data.frame('matched_name'=target_matched_name, 'target_name'=target_names)
  origin_df = data.frame('matched_name'=origin_matched_name, 'original_name'=origin_names)

  # add the updated admin names to the dataframe, under the original admin column name
  updated_names_df = merge(origin_df, target_df, all.x=TRUE, by='matched_name', sort=FALSE)
  if(all(updated_names_df$original_name == origin_names)){
    return(updated_names_df$target_name)
  } else{
    warning('merge did not maintain original order, need to fix standardize_admin_names_in_vector function.')
    return(NA)
  }
}

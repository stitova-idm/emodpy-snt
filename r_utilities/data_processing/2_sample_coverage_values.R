# using the DHS/MIS output, sample coverage parameters with probability proportional to likelihood

library(stats)

draw_timeseries = function(years, num_total_vec, num_true_vec, num_samples=100, plot_flag=FALSE){
  # using the number surveyed and the number observed positive in a set of years, draw sets of coverages proportional to the likelihoods for those years
  #   assuming the drawn coverages for the survey years, interpolate to get coverages in intermediate years
  sample_values = matrix(NA, ncol=length(min(years):max(years)),nrow=(num_samples+1))
  for(yy in 1:length(years)){
    sample_values[,years[yy]-min(years)+1] = c((num_true_vec[yy]/num_total_vec[yy]), sample(x=seq(0,1,0.001), size=num_samples, prob=dbinom(x=round(num_true_vec[yy]), size=round(num_total_vec[yy]), prob=seq(0,1,0.001)), replace=TRUE))
  }
  if(length(years)>1){
    for(ii in 1:(num_samples+1)){
      # fill in missing years with linear interpolated value
      sample_values[ii,]=approxfun(x=min(years):max(years), y=sample_values[ii,])(min(years):max(years))
    }
  }
  # if(plot_flag){  
  #   matplot(x=min(years):max(years), y=t(sample_values), type = "l", xlab='year', ylab='sampled coverage params', bty='L')
  # }
  colnames(sample_values) = min(years):max(years)
  return(sample_values) # first value is the mean value, all following are sampled values
}


# draw timeseries for baseline transmission calibration seeds (parameter values at each of 5 quantiles in the posterior: 0.1, 0.3, 0.5, 0.7, 0.9)
draw_quantile_timeseries = function(years, num_total_vec, num_true_vec, quantiles=c(0.1, 0.3, 0.5, 0.7, 0.9)){
  # using the number surveyed and the number observed positive in a set of years, get the coverage at each of five quantiles
  #   assuming the drawn coverages for the survey years, interpolate to get coverages in intermediate years
  sample_values = matrix(NA, ncol=length(min(years):max(years)),nrow=(length(quantiles)))
  for(yy in 1:length(years)){
    # sample_values[,years[yy]-min(years)+1] = c((num_true_vec[yy]/num_total_vec[yy]), sample(x=seq(0,1,0.001), size=num_samples, prob=dbinom(x=num_true_vec[yy], size=num_total_vec[yy], prob=seq(0,1,0.001)), replace=TRUE))
    sample_values[,years[yy]-min(years)+1] = qbinom(p=quantiles, size=num_total_vec[yy], prob=num_true_vec[yy]/num_total_vec[yy])/num_total_vec[yy]
  }
  if(length(years)>1){
    for(ii in 1:(length(quantiles))){
      # fill in missing years with linear interpolated value
      sample_values[ii,]=approxfun(x=min(years):max(years), y=sample_values[ii,])(min(years):max(years))
    }
  }
  # if(plot_flag){  
  #   matplot(x=min(years):max(years), y=t(sample_values), type = "l", xlab='year', ylab='sampled coverage params', bty='L')
  # }
  colnames(sample_values) = min(years):max(years)
  return(sample_values) # first value is the mean value, all following are sampled values
}





#################################################################################################################
# main function to sample parameter values from posterior
#################################################################################################################
sample_params_from_DHS_posterior = function(hbhi_dir, years, variables= c('mic','itn_u5','itn_5_10','itn_10_15','itn_15_20','itn_o20','iptp','cm','blood_test'),
                                            num_samples=50, quantiles=NA, sample_quantiles=FALSE, plot_flag=FALSE, sim_start_year=2010, last_sim_year=2020, min_num_total=30){
  if(!sample_quantiles){
    sample_type_name = ''
  } else{
    sample_type_name = 'quantiles_'
  }
  # iterate through variables and create csvs with sampled parameters to use in simulations
  for(i_var in 1:length(variables)){
    # iterate through DHS/MIS years, creating new dataframes with sampled coverages
    num_total_df = data.frame()
    num_true_df = data.frame()
    # dates_df = data.frame()
    cur_years = c()
    for(yy in 1:length(years)){
      admin_sums = read.csv(paste0(hbhi_dir, '/estimates_from_DHS/DHS_admin_minN',min_num_total,'_', years[yy], '.csv'))[,-1]
      if(paste0(variables[i_var], '_num_total') %in% colnames(admin_sums)){
        cur_years = c(cur_years, years[yy])
        num_total_df_cur = admin_sums[, c('NOMDEP',paste0(variables[i_var], '_num_total'))]
        num_true_df_cur = admin_sums[, c('NOMDEP',paste0(variables[i_var], '_num_true'))]
        # if('mean_date' %in% colnames(admin_sums)) dates_df_cur = admin_sums[, c('NOMDEP','mean_date')]
        colnames(num_total_df_cur)[2] = years[yy]
        colnames(num_true_df_cur)[2] = years[yy]
        # colnames(dates_df_cur)[2] = years[yy]
        if(nrow(num_total_df)>0){
          num_total_df = merge(num_total_df, num_total_df_cur, by='NOMDEP')
          num_true_df = merge(num_true_df, num_true_df_cur, by='NOMDEP')
          # if('mean_date' %in% colnames(admin_sums)) dates_df = merge(dates_df, dates_df_cur, by='NOMDEP')
        } else {
          num_total_df = num_total_df_cur
          num_true_df = num_true_df_cur
          # if('mean_date' %in% colnames(admin_sums)) dates_df = dates_df_cur
        }
      }
    }
    if(nrow(num_total_df)>0){
      # error message if admin names are in different orders
      if(all(num_total_df$NOMDEP == num_true_df$NOMDEP)) {
      
        # create dataframe where parameter value for each seed are in separate columns and each row is a different admin-year pair
        sample_df = data.frame()
        # itarate through admins, obtaining num_samples timeseries of parameter values (or obtaining a timeseries for each quantile)
        for(i_chief in 1:nrow(num_total_df)){
          if(!sample_quantiles){
            sample_df_cur = as.data.frame(t(draw_timeseries(years=cur_years, num_total_vec=as.numeric(num_total_df[i_chief,-1]), num_true_vec=as.numeric(num_true_df[i_chief,-1]), num_samples=num_samples)))
          } else{
            sample_df_cur = as.data.frame(t(draw_quantile_timeseries(years=cur_years, num_total_vec=as.numeric(num_total_df[i_chief,-1]), num_true_vec=as.numeric(num_true_df[i_chief,-1]), quantiles=quantiles)))
          }
          sample_df_cur$admin_name = num_total_df$NOMDEP[i_chief]
          sample_df_cur$year = min(cur_years):max(cur_years)
          sample_df = rbind(sample_df, sample_df_cur)
        }
        
        # ggplot(sample_df, aes(x=year, y=V1, color=admin_name)) + 
        #   geom_line()+
        #   theme(legend.position='none')
        
        
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
        # standardize years to go from sim_start_year to last_sim_year
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
        # if first datapoint is before sim_start_year, truncate at 2010, otherwise assume coverage is constant between 2010 and the first measure
        if(min(cur_years)<=sim_start_year){
          sample_df_2010 = sample_df[sample_df$year>=sim_start_year,]
        } else{
          sample_df_2010 = sample_df
          temp_df = sample_df[sample_df$year==min(cur_years),]
          for(yy in sim_start_year:(min(cur_years)-1)){
            temp_df$year = yy
            sample_df_2010 = rbind(sample_df_2010, temp_df)
          }
        }
        
        # if final datapoint is last_sim_year, truncate at last_sim_year, otherwise assume constant coverage between final datapoint and last_sim_year
        sample_df_2010_2020 = sample_df_2010
        if(max(cur_years)<last_sim_year){
          temp_df = sample_df_2010_2020[sample_df_2010_2020$year==max(cur_years),]
          for(yy in (max(cur_years)+1):last_sim_year){
            temp_df$year = yy
            sample_df_2010_2020 = rbind(sample_df_2010_2020, temp_df)
          }  
        } 
        
        # plot(NA, xlim=c(2010,2020), ylim=c(0,1))
        # for(i_chief in 1:nrow(num_total_df)){
        #   sample_df_cur = sample_df_2010_2020[sample_df_2010_2020$admin_name == num_total_df$NOMDEP[i_chief],]
        #   year_order = order(sample_df_cur$year)
        #   lines(sample_df_cur$year[year_order], sample_df_cur$V1[year_order])
        # }
        sample_df_2010_2020 = sample_df_2010_2020[,c(which(colnames(sample_df_2010_2020)=='admin_name'), which(colnames(sample_df_2010_2020)=='year'), grep('V',colnames(sample_df_2010_2020)))]
        colnames(sample_df_2010_2020)[grep('V',colnames(sample_df_2010_2020))] = gsub('V','sample_', colnames(sample_df_2010_2020)[grep('V',colnames(sample_df_2010_2020))])
        write.csv(sample_df_2010_2020, paste0(hbhi_dir, '/estimates_from_DHS/DHS_sampled_params_', sample_type_name, variables[i_var],'.csv'))
        
      } else {
        warning(paste0('for variable ', variables[i_var],', need to standardize admin order'))
      }
    }
  }
  
  
  if(plot_flag){
    ################################################################################################
    # plots for slides showing sampled timeseries for some/all admins
    ################################################################################################
    if(!dir.exists(paste0(hbhi_dir,'/simulation_inputs/plots'))) dir.create(paste0(hbhi_dir,'/simulation_inputs/plots'))
    for(i_var in 1:length(variables)){
      filename_cur = paste0(hbhi_dir, '/estimates_from_DHS/DHS_sampled_params_', sample_type_name, variables[i_var],'.csv')
      if(file.exists(filename_cur)){
        sample_df_2010_2020 = read.csv(filename_cur)[,-1]
        png(paste0(hbhi_dir,'/simulation_inputs/plots/', sample_type_name, variables[i_var], '.png'), width=8.5, height=5.5, units='in',res=300)
        par(mfrow=c(3,3))
        for(i_ds in seq(1,length(unique(sample_df_2010_2020$admin_name))/2, length.out=9)){
          # subset to a single admin
          single_admin_df = sample_df_2010_2020[sample_df_2010_2020$admin_name==unique(sample_df_2010_2020$admin_name)[i_ds],]
          # plot example
          year_order = order(single_admin_df$year)
          matplot(x=single_admin_df$year[year_order], y=single_admin_df[year_order,grep('sample_',colnames(single_admin_df))], type = "l", xlab='year', ylab='sampled coverage params', bty='L', ylim=c(0,1), main=paste0(variables[i_var], ' - admin: ', unique(sample_df_2010_2020$admin_name)[i_ds]), lty=1)
          if(!sample_quantiles){
            lines(x=single_admin_df$year[year_order], y=single_admin_df$sample_1[year_order], col='black', lwd=2)
          }
        }
        dev.off()
      }
    }

    # # single admin
    # sample_df_2010_2020 = read.csv(paste0(hbhi_dir, '/estimates_from_DHS/DHS_sampled_params_iptp.csv'))[,-1]
    # png(paste0(hbhi_dir,'/simulation_inputs/plots/iptp_Gisagara.png'), width=5, height=5, units='in',res=300)
    # par(mfrow=c(1,1))
    # # subset to a single admin
    # single_admin_df = sample_df_2010_2020[sample_df_2010_2020$admin_name=='Gisagara',]
    # # plot example
    # year_order = order(single_admin_df$year)
    # matplot(x=single_admin_df$year[year_order], y=single_admin_df[year_order,grep('sample_',colnames(single_admin_df))], type = "l", xlab='year', ylab='IPTp coverage used in simulation', bty='L', ylim=c(0,0.7), main='IPTp coverage for Gisagara simulations', lty=1, col=rgb(0.2,0.7,0.4))
    # lines(x=single_admin_df$year[year_order], y=single_admin_df$sample_1[year_order], col='black', lwd=2)
    # dev.off()


    # large grid with all admins
    for(i_var in 1:length(variables)){
      filename_cur = paste0(hbhi_dir, '/estimates_from_DHS/DHS_sampled_params_', sample_type_name, variables[i_var],'.csv')
      if(file.exists(filename_cur)){
        sample_df_2010_2020 = read.csv(paste0(hbhi_dir, '/estimates_from_DHS/DHS_sampled_params_', sample_type_name, variables[i_var],'.csv'))[,-1]
        nrow_plot = round(length(unique(sample_df_2010_2020$admin_name))^(1/2))
        ncol_plot = ceiling(length(unique(sample_df_2010_2020$admin_name))/nrow_plot)
        png(paste0(hbhi_dir,'/simulation_inputs/plots/', sample_type_name, variables[i_var],'_allAdmins.png'), width=3*ncol_plot, height=2.5*nrow_plot, units='in',res=300)
        par(mfrow=c(nrow_plot, ncol_plot))
        for(i_ds in seq(1,length(unique(sample_df_2010_2020$admin_name)))){
          # subset to a single admin
          single_admin_df = sample_df_2010_2020[sample_df_2010_2020$admin_name==unique(sample_df_2010_2020$admin_name)[i_ds],]
          # plot example
          year_order = order(single_admin_df$year)
          matplot(x=single_admin_df$year[year_order], y=single_admin_df[year_order,grep('sample_',colnames(single_admin_df))], type = "l", xlab='year', ylab='simulation coverage', bty='L', ylim=c(0,1), main=unique(sample_df_2010_2020$admin_name)[i_ds], lty=1, col=rgb(0.2,0.7,0.4))
          if(!sample_quantiles){
            lines(x=single_admin_df$year[year_order], y=single_admin_df$sample_1[year_order], col='black', lwd=2)        }
        }
        dev.off()
      }
    }
    par(mfrow=c(1,1))
  }
}







# 
# 
# 
# 
# 
# # IPTp figures for presentation explaining why and how we sample from these distributions
# par(mfrow=c(3,1), mar=c(4,4,1,1))
# NN=150
# ptrue=0.75
# plot(seq(1,NN), dbinom(seq(1,NN), size=NN, prob=ptrue), ylab='probability density', xlab='number reporting ITN use', bty='L', type='l')
# plot(seq(1,NN)/NN, dbinom(seq(1,NN), size=NN, prob=ptrue), ylab='probability density', xlab='estimated use rate', bty='L', type='l')
# abline(v=ptrue, col='red')
# rbinom(10,NN,ptrue)
# pobs=0.70
# plot(seq(1,NN)/NN, dbinom(seq(1,NN), size=NN, prob=pobs), ylab=c('fraction of simulations', ' run with this value'), xlab='use rate parameter', bty='L', type='l')
# abline(v=ptrue, col='red')
# par(mfrow=c(1,1))
# hist(rbinom(100, size=NN, prob=pobs)/NN, ylab=c('fraction of simulations', ' run with this value'), xlab='use rate parameter', bty='L', type='l', main='')
# 
# # IPTp Gisigara
# par(mfrow=c(1,3), mar=c(4,4,1,1))
# NN=71
# XX=0
# plot(seq(0,1,0.01), dbinom(XX, size=NN, prob=seq(0,1,0.01)), ylab=c('fraction of simulations', ' run with this value'), xlab='use rate parameter', bty='L', type='l', col=rgb(0.2,0.7,0.3), lwd=4)
# NN=38
# XX=0
# plot(seq(0,1,0.01), dbinom(XX, size=NN, prob=seq(0,1,0.01)), ylab=c('fraction of simulations', ' run with this value'), xlab='use rate parameter', bty='L', type='l', col=rgb(0.2,0.7,0.3), lwd=4)
# NN=101
# XX=28
# plot(seq(0,1,0.01), dbinom(XX, size=NN, prob=seq(0,1,0.01)), ylab=c('fraction of simulations', ' run with this value'), xlab='use rate parameter', bty='L', type='l', col=rgb(0.2,0.7,0.3), lwd=4)
# 
# # ITN U5 Gisigara
# par(mfrow=c(1,3), mar=c(4,4,1,1))
# NN=132
# XX=72
# plot(seq(0,1,0.01), dbinom(XX, size=NN, prob=seq(0,1,0.01)), ylab=c('fraction of simulations', ' run with this value'), xlab='use rate parameter', bty='L', type='l', col=rgb(0.2,0.7,0.3), lwd=4)
# NN=69
# XX=22
# plot(seq(0,1,0.01), dbinom(XX, size=NN, prob=seq(0,1,0.01)), ylab=c('fraction of simulations', ' run with this value'), xlab='use rate parameter', bty='L', type='l', col=rgb(0.2,0.7,0.3), lwd=4)
# NN=180
# XX=21
# plot(seq(0,1,0.01), dbinom(XX, size=NN, prob=seq(0,1,0.01)), ylab=c('fraction of simulations', ' run with this value'), xlab='use rate parameter', bty='L', type='l', col=rgb(0.2,0.7,0.3), lwd=4)
# 
# par(mfrow=c(1,3), mar=c(4,4,1,1))
# NN=132
# XX=72
# plot(seq(0,1,0.01), dbinom(XX, size=NN, prob=seq(0,1,0.01)), ylab=c('probability density'), xlab='use rate parameter', bty='L', type='l', col=rgb(0.2,0.7,0.3), lwd=4)
# plot(seq(0,1,0.01), dbinom(XX, size=NN, prob=seq(0,1,0.01)*0.662/0.855),ylab=c('probability density'), xlab='use rate parameter', bty='L', type='l', col=rgb(0.2,0.7,0.3), lwd=4)
# 
# 
# # blood-test  Murore
# par(mfrow=c(1,3), mar=c(4,4,1,1))
# NN=51
# XX=10
# plot(seq(0,1,0.01), dbinom(XX, size=NN, prob=seq(0,1,0.01)), ylab=c('fraction of simulations', ' run with this value'), xlab='CM parameter', bty='L', type='l', col=rgb(0.2,0.7,0.3), lwd=4)
# NN=43
# XX=19
# plot(seq(0,1,0.01), dbinom(XX, size=NN, prob=seq(0,1,0.01)), ylab=c('fraction of simulations', ' run with this value'), xlab='CM parameter', bty='L', type='l', col=rgb(0.2,0.7,0.3), lwd=4)
# NN=129
# XX=92
# plot(seq(0,1,0.01), dbinom(XX, size=NN, prob=seq(0,1,0.01)), ylab=c('fraction of simulations', ' run with this value'), xlab='CM parameter', bty='L', type='l', col=rgb(0.2,0.7,0.3), lwd=4)
# # burnin samples of blood-test, Murore
# quantiles = c(0.1, 0.3, 0.5, 0.7, 0.9)
# par(mfrow=c(1,3), mar=c(4,4,1,1))
# NN=51
# XX=10
# plot(seq(0,1,0.01), dbinom(XX, size=NN, prob=seq(0,1,0.01)), ylab='probability density', xlab='CM parameter', bty='L', type='l', col=rgb(0.2,0.7,0.3), lwd=4)
# for(qq in 1:length(quantiles)){
#   abline(v=qbinom(p=quantiles[qq], size=NN, prob=XX/NN)/NN, col=rgb(0.3,0.3,0.4), lwd=0.5)
# }
# NN=43
# XX=19
# plot(seq(0,1,0.01), dbinom(XX, size=NN, prob=seq(0,1,0.01)), ylab='probability density', xlab='CM parameter', bty='L', type='l', col=rgb(0.2,0.7,0.3), lwd=4)
# for(qq in 1:length(quantiles)){
#   abline(v=qbinom(p=quantiles[qq], size=NN, prob=XX/NN)/NN, col=rgb(0.3,0.3,0.4), lwd=0.5)
# }
# NN=129
# XX=92
# plot(seq(0,1,0.01), dbinom(XX, size=NN, prob=seq(0,1,0.01)), ylab='probability density', xlab='CM parameter', bty='L', type='l', col=rgb(0.2,0.7,0.3), lwd=4)
# for(qq in 1:length(quantiles)){
#   abline(v=qbinom(p=quantiles[qq], size=NN, prob=XX/NN)/NN, col=rgb(0.3,0.3,0.4), lwd=0.5)
# }
# 
# 
# 
# 
# 
# 
# 
# 
# # how many samples need to be taken to approximately converge on full distribution, depending on survey size?
# par(mfrow=c(3,3))
# nn=seq(10,90,10)
# xx = nn*0.9
# for(ii in 1:length(nn)){
#   hist(sample(x=seq(0,1,0.001), size=50, prob=dbinom(xx[ii], size=nn[ii], prob=seq(0,1,0.001)), replace=TRUE), freq=FALSE, main=paste0('N=',nn[ii]), breaks=seq(0,1,0.01), xlab='estimated coverage', ylab='prob density')
#   lines(seq(0,1,0.01), dbinom(xx[ii], size=nn[ii], prob=seq(0,1,0.01))/sum( dbinom(xx[ii], size=nn[ii], prob=seq(0,1,0.01)))*100, col='red')
# }
# par(mfrow=c(1,1))
# 
# 
# 
# 
# 
# 
# # # create array where dim 1 is year, dim2 is seed, dim3 is admin Dimensions: c('year','seed','admin')
# # sample_array = array(NA, dim=c((max(years)-min(years)+1), num_samples, nrow(num_total_df)), dimnames = list(min(years):max(years), 1:num_samples, num_total_df$NOMDEP))
# # sample_array[,,i_chief] = t(draw_timeseries(years=years, num_total_vec=as.numeric(num_total_df[i_chief,-1]), num_true_vec=as.numeric(num_true_df[i_chief,-1]), num_samples=num_samples, plot_flag=FALSE))

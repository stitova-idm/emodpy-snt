# explore_net_retention_choices.R







# examine different net-retention distributions
# imgine true median time is 1.7 years

# current distribution:
r1 = c(rexp(10000000*.9,1/1.7), rexp(10000000*.1,1/10))
mean(r1)
median(r1)
# if we switched to calculating mean from median and got rid of dual
r2 = (rexp(10000000,1/1.7*log(2)))
mean(r2)
median(r2)
# if we kept dual and adjusted first distribution mean downward to obtain desired median
r3 = c(rexp(10000000*.9,1/2.2), rexp(10000000*.1,1/10))
mean(r3)
median(r3)


# function Amelia used for smooth-compact loss function:
#   Loss(t,tau) = exp(18 - 18/(1-(t/tau)^2)), where t is time
# so if the halflife is 1.7, the tau should satisfy
#   0.5 = exp(18 - 18/(1-(1.7/tau)^2))
#  we can solve for tau
#    ln(0.5) = 18 - 18/(1-(1.7/tau)^2)
#    18/(1-(1.7/tau)^2) = 18 - ln(0.5)
#    (1-(1.7/tau)^2)  = 18 / (18 - ln(0.5))
#    (1.7/tau)^2 = 1 - 18 / (18 - ln(0.5))
#    1.7/tau = (1 - 18 / (18 - ln(0.5))) ^ (1/2)
#    1.7 * 1 - 18 / (18 - ln(0.5))) ^ (-1/2) = tau
1.7 * (1 - 18 / (18 - log(0.5)))^(-1/2)
# check that the algebra worked...
# exp(18 - 18 / (1 - (1.7/8.828)^2))
#
# Nigeria (used half-life of 1.7 years for comparison with WHO report, but note that the Bertozzi-villa estimate for Nigeria was 2.22 years)
halflife = 1.7#2.22
tau = halflife * (1 - 18 / (18 - log(0.5)))^(-1/2)
# which of the distributions looks most like this one?
tt = seq(0,10,0.1)  # seq(0,5,0.1)
plot(NA, xlim=c(min(tt), max(tt)), ylim=c(0,1), xlab='time since distribution (years)', ylab='fraction retaining net', bty='L')
lines(tt, exp(18 - 18 / (1 - (tt/tau)^2)), type='l', lwd=3)  # Amelia's function
# mean_est = halflife/log(2)
# lines(tt, exp(-(1/mean_est) * tt), col='blue')
# lines(tt, exp(-(1/1.7) * tt), col='green')
# lines(tt, 0.9*exp(-(1/mean_est) * tt) + 0.1*exp(-1/100*tt), col='dodgerblue')  # dual exponential
# lines(tt, 1-plnorm(tt, meanlog=log(halflife), sdlog=0.5), col='red')  # lognormal
lines(tt, 1-plnorm(tt, meanlog=log(halflife), sdlog=0.8), col='darkred', lwd=2)  # lognormal
# try the weibull distribution instead of lognormal
w_shape = 2.3
w_scale = halflife / log(2)^(1/w_shape)
lines(tt, 1-pweibull(tt, shape=w_shape, scale=w_scale), col='green', lwd=2)
legend('topright', c("Amelia's model",'lognormal, sd=0.8', 'weibull, shape=2.3'), lty=1, col=c('black','darkred', 'green'), lwd=c(3,2,2), bty='n')
# legend('topright', c("Amelia's model", 'dual exponential currently used', 'exponential with matching median', 'lognormal, sd=0.5', 'lognormal, sd=0.8', 'weibull, shape=2.3'), lty=1, col=c('black', 'dodgerblue','blue','red','darkred', 'green'), lwd=c(3,1,1,1,2), bty='n')
# legend('topright', c('Amelia model', 'exponential with matching median','exponential with mean as median'), lty=1, col=c('black','blue','green'), lwd=c(3,1,1), bty='n')
# lines(tt, 1-ppois(tt*365, lambda=2.5*365), col='purple')  # poisson

mean(x=sample(seq(0.0001,10,0.001), size=10000, prob=exp(18 - 18 / (1 - (seq(0.0001,10,0.001)/tau)^2)), replace=TRUE))

# lines(tt, 0.9*exp(-(1/mean_est*2) * tt) + 0.1*exp(-1/100*tt), col='dodgerblue')  # dual exponential

# plot mean and other sampled retention shapes
tt = seq(0,5,0.1)
sampled_retentions = rnorm(n=(num_samples), mean=median_llin_retention_est, sd=((median_llin_retention_est-median_llin_retention_lowerCI)/2))
plot(tt, 1-plnorm(tt, meanlog=log(median_llin_retention_est), sdlog=0.8), type='l', lwd=2, xlab='time since distribution (years)', ylab='fraction retaining net', bty='L', col=rgb(0.6,0.1,1))
for(ss in 1:num_samples){
  lines(tt, 1-plnorm(tt, meanlog=log(sampled_retentions[ss]), sdlog=0.8), type='l', lwd=1, col=rgb(0.8,0.7,1))  # col=rgb(0.6,0.1,1,0.5))
}
lines(tt, 1-plnorm(tt, meanlog=log(median_llin_retention_est), sdlog=0.8), type='l', lwd=2, xlab='time since distribution (years)', ylab='fraction retaining net', bty='L', col=rgb(0.6,0.1,1))



# 
# ### - - - - - - - - - - - - ###
# # # Burundi values
# ### - - - - - - - - - - - - ###
# halflife = 1.31
# tau = halflife * (1 - 18 / (18 - log(0.5)))^(-1/2)
# # plot comparison between Amelia's net retention function and our assumption
# tt = seq(0,6,0.1)  # seq(0,5,0.1)
# plot(tt, exp(18 - 18 / (1 - (tt/tau)^2)), type='l', lwd=3, xlab='time since distribution (years)', ylab='fraction retaining net', bty='L') # Amelia's function
# lines(tt, 1-plnorm(tt, meanlog=log(halflife), sdlog=0.8), col='darkred', lwd=2)  # lognormal
# legend('topright', c("function used in Bertozzi-Villa fit", 'function used here'), lty=1, col=c('black', 'darkred'), lwd=c(3,2), bty='n')
# # legend('topright', c('Amelia model', 'exponential with matching median','exponential with mean as median'), lty=1, col=c('black','blue','green'), lwd=c(3,1,1), bty='n')
# # lines(tt, 1-ppois(tt*365, lambda=2.5*365), col='purple')  # poisson
# 
# lines(tt, 1-plnorm(tt*365, meanlog=log(halflife*365), sdlog=0.8), col='green', lwd=2)  # lognormal... in days
# halflife*365
# mean(rlnorm(10000,meanlog=log(halflife*365),sd=0.8))
# median(rlnorm(10000,meanlog=log(halflife*365),sd=0.8))
# hist(sapply(rlnorm(10000,meanlog=log(halflife*365),sd=0.8),min,800), breaks=500)
# 
# 
# # match on means instead of on medians
#   # 1-cum_dist(tt) = exp(18 - 18 / (1 - (tt/tau)^2))
#   # cum_dist(tt) = 1 - exp(18 - 18 / (1 - (tt/tau)^2))
# halflife_days = halflife * 365
# tau_days = halflife_days * (1 - 18 / (18 - log(0.5)))^(-1/2)
# tt_days = tt*365
# plot(tt_days, exp(18 - 18 / (1 - (tt_days/tau_days)^2)), type='l', lwd=3, xlab='time since distribution (days)', ylab='fraction retaining net', bty='L') # Amelia's function
# 
# ss = seq(0,6*365,1)
# cd = 1 - exp(18 - 18 / (1 - (ss/tau_days)^2))
# dd = cumsum(cd[1:(length(cd)-1)]) - cumsum(cd[2:(length(cd))])
# dd = cd[2:(length(cd))] - cd[1:(length(cd)-1)]
# dens = dd/sum(dd)
# amelia_mean = sum(dens*ss[-1])  # mean of distribution
# print(paste0("mean is estimated at: ", round(amelia_mean)))
# ss[which.min(abs(cumsum(dens) - 0.5))+1]  # should be close to the halflife
# 
# 
# # mean of lognormal distribution is  exp(mu + sigma^2/2) while median is exp(mu)
# #   when we were basing it off the median, we calculated mu = log(median)
# #   to base the parameterization off the mean, mu = log(mean) - sigma^2/2
# lines(tt_days, 1-plnorm(tt_days, meanlog=(log(amelia_mean)-0.8^2/2), sdlog=0.8), col='darkred', lwd=2)  # lognormal
# legend('topright', c("function used in Bertozzi-Villa fit", 'lognormal function used here'), lty=1, col=c('black', 'darkred'), lwd=c(3,2), bty='n')
# mean(rlnorm(10000, meanlog=(log(amelia_mean)-0.8^2/2), sdlog=0.8))  # should be ethe same as the mean of Amelia's distribution
# 
# 
# # calculate the density function (taking the derivative of the CDF of the function in Amelia's paper) and use that to calculate the mean
# plot(ss[-1], dens, type='l', main='density functions describing net retention times', bty='L', ylim=c(0,0.002))
# # derivative of CDF:
# lines(tt_days, 18*exp(18 - 18 / (1 - (tt_days/tau_days)^2)) * 2*tt_days / tau_days^2 / (1 - (tt_days/tau_days)^2)^2, col='black', lwd=2)  # it should overlap the plot of the density calculated by subtracting the subsequent cumsums
# # calculate expected value (distribution mean)
# xx = seq(0,(6*365),0.1)
# exp_value = weighted.mean(xx,  18*exp(18 - 18 / (1 - (xx/tau_days)^2)) * 2*xx / tau_days^2 / (1 - (xx/tau_days)^2)^2)
# abline(v=exp_value, col='blue')  # expected value
# abline(v=halflife * 365, col='grey')  # median
# 
# 
# lines(tt_days, dlnorm(tt_days, meanlog=(log(exp_value)-0.8^2/2), sdlog=0.8), col='blue', lty=2)
# lines(tt_days, dlnorm(tt_days, meanlog=(log(halflife*365)), sdlog=0.8), col='grey', lty=2)
# legend('topright', c("function from Amelia's paper", "lognormal with matching mean", "lognormal with matching median"), col=c("black","blue","grey"), lwd=c(2,1,1), bty='n')
# 
# 
# 
# 
# # plot lognormal retention time distribution for quantiles
# tt = seq(0,5*365,1)
# sampled_retentions = c(5.805830331,5.858027862,5.892646209,5.926106118,5.972523465)
# plot(NA, xlim=c(0,max(tt)/365), ylim=c(0,1), type='l', lwd=2, xlab='time since distribution (years)', ylab='fraction retaining net', bty='L', col=rgb(0.6,0.1,1))
# for(ss in 1:length(sampled_retentions)){
#   lines(tt/365, 1-plnorm(tt, meanlog=sampled_retentions[ss], sdlog=0.8), type='l', lwd=1, col=rgb(0.6,0.1,1))  # col=rgb(0.6,0.1,1,0.5))
# }


# # plot lognormal retention time distribution for sampled values
# tt = seq(0,5*365,1)
# sampled_retentions = read.csv(paste0(hbhi_dir, '/simulation_inputs/itn_discard_decay_params.csv'))$net_life_lognormal_mu
# plot(NA, xlim=c(0,max(tt)/365), ylim=c(0,1), type='l', lwd=2, xlab='time since distribution (years)', ylab='fraction retaining net', bty='L', col=rgb(0.6,0.1,1))
# for(ss in 1:length(sampled_retentions)){
#   lines(tt/365, 1-plnorm(tt, meanlog=sampled_retentions[ss], sdlog=0.8), type='l', lwd=1, col=rgb(0.6,0.1,1))  # col=rgb(0.6,0.1,1,0.5))
# }

#MCMC_samps is the posterior samples from the fit in analysis_script.R
library(dlimIM)
library(ggplot2)

load("m_hdl_c48_mcmc_samps.rda") #from analysis_script.R

plot_bdlim(x = MCMC_samps,
           m_star = seq(0.1, 0.6, 0.1),
           burnin = burnin,
           thin = thin,
           type = "by_modifier",
           exp_times = -2:22) + xlab("Months since LMP")+
           ylab(paste("Change in", response_name, "per exposure unit"))


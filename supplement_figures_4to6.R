#MCMC_samps is the posterior samples from the fit in analysis_script.R
library(dlimIM)
library(ggplot2)

                  ### Supplemental Figure 4 ###

load("m_ldl_c48_mcmc_samps.rda") #from analysis_script.R

#MCMC info
niter <- nrow(MCMC_samps$chain1)
burnin <- attr(MCMC_samps, "burnin")
thin <- 10
inf_ind <- seq(burnin, niter, thin)

plot_bdlim(x = MCMC_samps,
           m_star = m_star,
           burnin = burnin,
           thin = thin,
           type = "cumulative") + 
xlab("Months since LMP") +
ylab(paste("Cumulative change in LDL \n per exposure unit"))

                    ### Supplemental Figure 5 ###

load("m_ldl_c48_mcmc_samps.rda") #from analysis_script.R

#MCMC info
niter <- nrow(MCMC_samps$chain1)
burnin <- attr(MCMC_samps, "burnin")
thin <- 10
inf_ind <- seq(burnin, niter, thin)

est <- colMeans(MCMC_samps$chain1[inf_idx,weight_idx])
est_m_star <- M%*%est
m_star_quarts <- round(quantile(est_m_star, probs=c(0.1, 0.25, 0.3, 0.4, 0.75, 0.9)),2)
plot_bdlim(x = MCMC_samps,
           m_star = m_star_quarts,#c(0.1, 0.3, 0.5),#seq(0.1, 0.6, 0.1),
           burnin = burnin,
           thin = thin,
           type = "by_modifier",
           exp_times = -2:22) + 
xlab("Months since LMP") +
ylab(paste("Change in LDL per exposure unit"))

                     ### Supplemental Figure 6 ###

load("m_hdl_c48_mcmc_samps.rda") #from analysis_script.R

#MCMC info
niter <- nrow(MCMC_samps$chain1)
burnin <- attr(MCMC_samps, "burnin")
thin <- 10
inf_ind <- seq(burnin, niter, thin)

plot_bdlim(x = MCMC_samps,
           m_star = m_star,
           burnin = burnin,
           thin = thin,
           type = "cumulative") + 
  xlab("Months since LMP") +
  ylab(paste("Cumulative change in HDL \n per exposure unit"))
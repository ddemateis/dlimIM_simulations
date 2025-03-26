#saves two .rda files with MCMC samples
#hdl_mcmc_samps.rda and ldl_mcmc_samps.rda
#use in maintext_figure3.R and supplement_figures_13to16.R

library(dlimIM)
library(ggplot2)
library(gridExtra)

#df
df <- 3
df_m <- df
df_l <- df

#hyperparameters
tau2 <- 100
xi2 <- 110
a <- 1
b <- 0.001

#MCMC specifications
niter <- 30000
burnin <- 10000
thin <- 1
inf_idx <- seq(burnin, niter, thin)

#data
modifiers <- c("epdsT", "nleT", "pssT", "sum_staiT")
covariates <- c("mother_age00","smk", "marital", "parity",
                "mother_pre_bmi2T", "SES_3cat2T", "lmpseason",
                "alcohol_preg48", "flagmed48") #static covariates

set.seed(101124)

responses <- c(HDL = "hdl", LDL = "ldl")
for(response in responses){
  #load the data set
  load(paste0(response, "_complete_data.rda"))
  
  #model specifications
  x <- as.matrix(dta_complete[,grep("X", colnames(dta_complete))])
  y <- dta_complete[,1]
  M <- as.matrix(dta_complete[,modifiers])
  z <- dta_complete[,covariates]
  n <- nrow(x)
  
  #scale modifiers
  M <- apply(M, 2, scale_01)
  
  #fit
  MCMC_samps <- MCMC_sampler_m(x = x,
                               y = y,
                               M = M,
                               z = z,
                               df_m = df_m,
                               df_l = df_l,
                               niter = niter,
                               burnin = burnin,
                               tau2 = tau2,
                               xi2 = xi2,
                               a = a,
                               b = b,
                               var_select = F,
                               weights_prior = rep(5, 4),
                               nchains = 5)
  
  save(MCMC_samps, file = paste0(response,"_mcmc_samps.rda"))
}

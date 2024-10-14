                      #######################
                      ### Stress Analysis ###
                      #######################
                      
#saves two .rda files with MCMC samples
#m_hdl_c48_mcmc_samps.rda and m_ldl_c48_mcmc_samps.rda
#use in maintext_figure3.R, supplement_figure4.R,
#supplement_figure5.R, and supplement_figure6.R
                      
library(dlimIM)
                      
responses <- c("m_hdl_c48", "m_ldl_c48")

#modifiers, covariates, exposure
modifiers <- c("epdsT", "nleT", "pssT", "sum_staiT")
covariates <- c("mother_age00","smk", "marital", "parity",
                "mother_pre_bmi2T", "SES_3cat2T", "lmpseason",
                "alcohol_preg48", "flagmed48") #static covariates
exposure <- colnames(data96_clean)[grep("X", colnames(data96_clean))]

#df
df <- 3
df_m <- df
df_l <- df

dir.create(paste0("df", df))

#hyperparameters
tau2 <- 100
xi2 <- 110
a <- 1
b <- 0.001

#MCMC specifications
niter <- 30000
burnin <- 10000
thin <- 10
inf_idx <- seq(burnin, niter, thin)

set.seed(101124)

for(response in responses){
  
  pdf(paste0("df", df,"/results_", response, "_df", df, ".pdf"))
  
  #subset for complete cases
  dta <- full_data[,c(response,
                      exposure, 
                      modifiers, 
                      covariates)]
  dta$SES_3cat2T <- as.factor(dta$SES_3cat2T)
  dta_complete <- dta[complete.cases(dta),]#272 for HDL48 LDL48
  
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
                               weights_prior = rep(5, 4))
  
  save(MCMC_samps, file = paste0(response,"_mcmc_samps.rda"))
}

                      #############################
                      ### EnviroScreen Analysis ###
                      #############################

#saves two .rda files with MCMC samples
#analysis_MCMC_samps_Dir1_df3_chain1.rda and 
#analysis_MCMC_samps_Dir1_df4_chain1.rda
#use in maintext_figure2.R, maintext_table2.R, supplement_figure2.R,
#supplement_figure3.R, and supplement_table2.R

### Data cleaning ###

#Ander: add cleaning code you used before giving me the data set.

#load data set, name it co_birth_index
#data already subset to the front range region and elevation under 6000ft
load("co_birth_index.rda")

#extract PM exposures
pm_pred <- (as.matrix(co_birth_index[,grep("pm25", colnames(co_birth_index), ignore.case=T)]))[,1:37]
pm_pred <- scale(pm_pred)

#remove PM columns
dat <- co_birth_index[,-grep("PM25", colnames(co_birth_index), ignore.case=T)]

#save
save(dat, file = "bw_es_dat_index.rda", compress = "xz")
save(pm_pred, file = "bw_es_pm_pred_index.rda", compress = "xz")


### Set up data for model fitting ###

library(dlimIM)
load("bw_es_dat_index.rda") #dat
load("bw_es_pm_pred_index.rda") #pm_pred

#percentiles making up each component score
sens_pop <- c(Asthma = "As___P", Cancer = "Cn__P", Diabetes = "Db__P", Heart_disease = "H____", Life_expectancy = "Lf__P", Mental_health = "M___P", Percent_64 = "P__64_", Percent_5 = "P__5_")
demographics <- c(Housing_cost = "H___P", Disability = "Pr__P", HS= "P______", Lingustic_iso = "Prcnt_ln__", Low_income = "Prcnt_lw__", POC = "Prc____P")
modifiers <- dat[,colnames(dat) %in% c(sens_pop, demographics)]
modifiers <- modifiers/100 #scale modifier percentiles to lie between 0 and 1

#covariates
z <- dat[, colnames(dat) %in% c("MOC" , "YOC", 
                                "MatAge", "MotherBMI", "MEduc", 
                                "MotherHeightIn",  "PriorWeight", "race", 
                                "hispanic", "Income", "Marital2", 
                                "elev_feet_tract", "fipscoor", 
                                "tmmx_tri1", "tmmx_tri2", "tmmx_tri3", 
                                "PrenatalCare", "Smk")]
z$MOC <- as.factor(z$MOC)
z$YOC <- as.factor(z$YOC)
z$MatAge2 <- z$MatAge^2

#create full data set 
full_dat <- data.frame(x = pm_pred, 
                       y = dat$bwgaz, 
                       mods = modifiers, 
                       Z = z)


### Fit the DLIM-IM ###

#hyperparameters
tau2 <- 100
xi2 <- 110
a <- 1
b <- 0.001
df <- 3
df_m <- df
df_l <- df
dir_param <- 1
niter <- 50000
burnin <- 30000

#set seed for replication
set.seed(040124+0)
rm(dat)
rm(pm_pred)

#specify arguments for MCMCM
x <- as.matrix(full_dat[,grep("pm25", colnames(full_dat))])
y <- full_dat$y
M <- as.matrix(full_dat[,grep("mods", colnames(full_dat))])
z <- full_dat[,grep("Z.", colnames(full_dat))]
n <- nrow(x)

#run MCMC sampler
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
                             model_type = "ns",
                             weights_prior = rep(dir_param,ncol(modifiers)),
                             var_select = T)
save(MCMC_samps, file = paste0("analysis_MCMC_samps_Dir",
                               dir_param,
                               "_df",
                               df,
                               "_chain",
                               1,
                               ".rda"))

### Fit the DLIM-IM with 4 degrees of freedom (for supplement)

df <- 4
df_m <- df
df_l <- df

#run MCMC sampler
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
                             model_type = "ns",
                             weights_prior = rep(dir_param,ncol(modifiers)),
                             var_select = T)
save(MCMC_samps, file = paste0("analysis_MCMC_samps_Dir",
                               dir_param,
                               "_df",
                               df,
                               "_chain",
                               1,
                               ".rda"))


### Fit the DLIM with linear modification and fixed modification index ###

library(dlim)
library(mgcv)

linear_fit <- dlim(y = y, 
                   x = x, 
                   modifier = dat$Hl__S_F_S,
                   z = z, 
                   df_m = 20,
                   df_l = 20, 
                   penalize = TRUE, 
                   model_type = "linear",
                   fit_fn = "bam") 

save(linear_fit, file = "linear_fit.rda")


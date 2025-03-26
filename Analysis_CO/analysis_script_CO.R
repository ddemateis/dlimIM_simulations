                      #############################
                      ### EnviroScreen Analysis ###
                      #############################

#saves 10 .rda files with MCMC samples
#analysis_MCMC_samps_Dir1_df3_chain1.rda,
#analysis_MCMC_samps_Dir1_df3_chain2.rda,
#analysis_MCMC_samps_Dir1_df3_chain3.rda,
#analysis_MCMC_samps_Dir1_df3_chain4.rda,
#analysis_MCMC_samps_Dir1_df3_chain5.rda,
#analysis_MCMC_samps_Dir1_df4_chain1.rda,
#analysis_MCMC_samps_Dir1_df4_chain2.rda,
#analysis_MCMC_samps_Dir1_df4_chain3.rda,
#analysis_MCMC_samps_Dir1_df4_chain4.rda, and
#analysis_MCMC_samps_Dir1_df4_chain5.rda
#also save .rda file with DLIM fit, linear_fit.rda
#use in maintext_figure2_table2.R, supplement_figures2_3_table_2.R
                      
library(dlimIM)
                      

### Data cleaning ###

#To be used on cleaned data. See read me for more details. 

#load data set, name it co_birth_danielle_index
#data already subset to the front range region and elevation under 6000ft
load("co_birth_danielle_index.rda")

#extract PM exposures
pm_pred <- (as.matrix(co_birth_danielle_index[,grep("pm25", colnames(co_birth_danielle_index), ignore.case=T)]))[,1:37]
pm_pred <- scale(pm_pred)

#remove PM columns
dat <- co_birth_danielle_index[,-grep("PM25", colnames(co_birth_danielle_index), ignore.case=T)]

### Set up data for model fitting ###

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


### Fit the DLIM-IM with 3 and 4 degrees of freedom ###

## Run the next section on HPC using the array 
#command #SBATCH --array=0-9 to create separate jobs for each 
#degrees of freedom (2) and chain (5) combination 


#read in array index
array_idx <- as.numeric(commandArgs(trailingOnly=TRUE)) #read in array values from the shell file, which is just the task number

#chain, prior, df
prior_grid <- 1 
df_grid <- c(3, 4)
chain_grid <- 1:5
n_prior <- length(prior_grid)
n_df <- length(df_grid)
n_chain <- length(chain_grid)
dir_param <- prior_grid[floor(array_idx/(n_df*n_chain))+1]
df <- df_grid[(floor(array_idx/n_chain))%%n_df+1]
chain <- chain_grid[floor(array_idx)%%n_chain+1]

#create directory to save files
dir.create(paste0("analysis_results"))
dir.create(paste0("analysis_results/Dir", dir_param))
dir.create(paste0("analysis_results/Dir", dir_param, "/df", df))

#hyperparameters
tau2 <- 100
xi2 <- 110
a <- 1
b <- 0.001
df_m <- df
df_l <- df
niter <- 50000
burnin <- 30000

#set seed for replication
set.seed(040124+array_idx)
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
save(MCMC_samps, file = paste0("analysis_results/Dir", dir_param, "/df", df,
                               "/analysis_MCMC_samps_Dir",
                               dir_param,
                               "_df",
                               df,
                               "_chain",
                               chain,
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


library(dlimIM)
library(dlim)
library(mgcv)

#create directories
dir.create("analysis_results")
dir.create("analysis_results/single_modifiers")

#read in array index
array_idx <- as.numeric(commandArgs(trailingOnly=TRUE)) #read in array values from the shell file, which is just the task number

#To be used on cleaned data. See read me for more details. 

#load data set, name it co_birth_danielle_index
#data already subset to the front range region and elevation under 6000ft
load("co_birth_danielle_index.rda")

#extract PM exposures
pm_pred <- (as.matrix(co_birth_danielle_index[,grep("pm25", colnames(co_birth_danielle_index), ignore.case=T)]))[,1:37]
pm_pred <- scale(pm_pred)

#remove PM columns
dat <- co_birth_danielle_index[,-grep("PM25", colnames(co_birth_danielle_index), ignore.case=T)]

#percentiles making up each component score
sens_pop <- c(Asthma = "As___P", Cancer = "Cn__P", Diabetes = "Db__P", Heart_disease = "H____", Life_expectancy = "Lf__P", Mental_health = "M___P", Percent_64 = "P__64_", Percent_5 = "P__5_")
demographics <- c(Housing_cost = "H___P", Disability = "Pr__P", HS= "P______", Lingustic_iso = "Prcnt_ln__", Low_income = "Prcnt_lw__", POC = "Prc____P")

#pick out the modifier variables
modifiers <- dat[,colnames(dat) %in% c(sens_pop, demographics)]
modifiers <- modifiers/100

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

#create full data set and subset to complete observations
full_dat <- data.frame(x = pm_pred, 
                       y = dat$bwgaz, 
                       mods = modifiers, 
                       Z = z)
full_dat <- full_dat[complete.cases(full_dat),]

#hyperparameters
tau2 <- 100
xi2 <- 110
a <- 1
b <- 0.001
df <- 3
df_m <- df
df_l <- df
niter <- 50000
burnin <- 30000
dir_param <- 1

#set seed for replication
set.seed(040124+array_idx)#original seed

#specify arguments for MCMCM
x <- as.matrix(full_dat[,grep("pm25", colnames(full_dat))])
y <- full_dat$y
M <- as.matrix(full_dat[,grep("mods", colnames(full_dat))])
m <- matrix(M[,array_idx + 1], ncol=1)
z <- full_dat[,grep("Z.", colnames(full_dat))]
n <- nrow(x)

#DLIM-IM
MCMC_samps <- MCMC_sampler_m(x = x,
                             y = y,
                             M = m,
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
                             weights_prior = dir_param)
save(MCMC_samps, file = paste0("analysis_results/single_modifiers/analysis_MCMC_samps_dlimIM_modifier_",
                               array_idx, ".rda"))


## Run on HPC using provided submit_all.sh script for using on HPC
## or remove the array_idx line below and loop array_idx over 0 to 971

## run modifier_correlation_structures.R first to create the correlation matrices
## creates a simulation results folder with .rda files to be used in
## maintext_figure1_table1.rda, supplemental_figure2_table1.rda

library(dlim)
library(MASS)
library(dlimIM)

#read in array index
array_idx <- as.numeric(commandArgs(trailingOnly=TRUE)) #read in array values from the shell file, which is just the task number

#create directory to save files
dir.create("simulation_results")
dir.create("simulation_traceplots")

#set seed
set.seed(080123)

#MCMC specifications
niter <- 10000 
burnin <- 8000 

#number of simulated data sets
total_sims <- 200
max_jobs <- 1000 #max number of jobs I can run on HPC

#set up grids for simulation scenarios
SNR_grid <- c(0.1, 0.5, 1) #length 3
model_list <- list(Bayes_select = "Bayes_select",
                   Bayes_Dir1 = "Bayes_Dir1",
                   Freq_ps = "Freq_ps")
true_weight_list <- list(sparse_50 = c(rep(1,3)/3, rep(0,47)),#for first 50 data sets
                         sparse_10 = c(c(0.5, 0.4, 0.1), rep(0,7)),
                         equal_3 = rep(1,3)/3,
                         diff_3 = c(0.5, 0.4, 0.1),
                         equal_10 = rep(1,10)/10,
                         diff_10 = c(0.3, 0.2, rep(0.1,2), rep(0.05,6)) #length 9
)

#set up simulated data sets grid
n_scenarios <- length(true_weight_list)*length(SNR_grid)*length(model_list) #number of simulation scenario groups, floor to not exceed max_jobs
max_parallel_runs <- floor(max_jobs/n_scenarios) #max number of parallel runs per simulation scenario
min_sims_per_run <- floor(total_sims/max_parallel_runs) #minimum simulations per job
runs_with_extra_sims <- total_sims %% max_parallel_runs #number of left over simulated data sets

if(min_sims_per_run == 0){
  max_parallel_runs <- total_sims
}
sim_distribution <- list()
num_distributed <- 0
for(i in 1:max_parallel_runs){
  if(i <= runs_with_extra_sims){
    num_in_run <- min_sims_per_run + 1
  }else{
    num_in_run <- min_sims_per_run
  }
  sim_distribution[[i]] <- seq(num_distributed+1, num_distributed + num_in_run)
  num_distributed <- num_distributed + num_in_run
}

#number of each scenario
n_SNR <- length(SNR_grid)
n_model <- length(model_list)
n_weight <- length(true_weight_list)
n_sim_groups <- length(sim_distribution)
#n_SNR * n_model * n_weight * n_sim_groups

#choose type, SNR, and weights
true_weights <- true_weight_list[[floor(array_idx/(n_sim_groups*n_model*n_SNR))+1]]
SNR <- SNR_grid[floor(array_idx/(n_sim_groups*n_model))%%n_SNR+1]
model <- model_list[[(floor(array_idx/n_sim_groups))%%n_model+1]]
sims <- sim_distribution[[array_idx %% max_parallel_runs + 1]]

#specify weight scenario name
weight_scenario <- names(true_weight_list)[[floor(array_idx/(n_sim_groups*n_model*n_SNR))+1]]

#model specifications
model_type <- "ns" #natural splines for bases
df_m <- 5 #set to 20 for penalized DLIM below
df_l <- 5 #set to 20 for penalized DLIM below

#data specifications
x <- exposure #exposure from library(dlim)
n <- nrow(x) #number of observations
type <- 4 #non-linear modification structure

#hyperparameters
tau2 <- 100
xi2 <- 110
a <- 1
b <- 0.001

#other parameters
true_gammas <- round(runif(length(true_weights), -1, 1),1)
n_mod <- length(true_weights)

#set up for storage
table_results <- c() #uncomment for HPC
weights_table <- c() #uncomment for HPC
all_betas <- vector(mode = "list", length = length(sims))#uncomment for HPC
WAIC <- c() #uncomment for HPC
chains <- list()

#start simulation loop
for(i in sims){
  
  #set seed specific to each data set
  set.seed(080123 + i)
  
  ### Fit ###
  
  #load modifier-specific covariance structure
  #constructed in modifier_correlation_structure.R
  load(paste0("modifier_Sigma_",n_mod,".rda"))
  
  #generate modifiers
  M <- mvrnorm(n=n, 
               mu=rep(0, n_mod), 
               Sigma = Sigma)
  M <- apply(M, 2, scale_01)
  
  #simulate data
  dta <- sim_data_m(x=x, 
                    M=M, 
                    w=true_weights,
                    SNR = SNR, 
                    type=type, 
                    gamma = true_gammas) 
  y <- dta$y
  z <- dta$Z
  pred_seq <- seq(0.2, 0.8, length = 100)
  
  #Fit model
  if(model == "Bayes_select" | model == "Bayes_lin_select"){
    if(grepl("lin", model)){
      model_type <- "linear"
    }
    fit <- MCMC_sampler_m(x = x,
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
                          model_type = model_type,
                          var_select = T,
                          WAIC = F)
    fit_type <- "Bayes"
    ests <- pred_m(posterior_list = fit, 
                   burnin = burnin,
                   m_star = pred_seq)
    weight_idx <- grep("^weight",colnames(fit[[1]]))
    weight_ests <- colMeans(fit$chain1[-c(1:burnin),weight_idx])
    weight_PIPs <- colMeans(fit$chain1[-c(1:burnin),weight_idx]!=0)
    weight_LB <- apply(fit$chain1[-c(1:burnin),weight_idx], 2, quantile, probs=0.025)
    weight_UB <- apply(fit$chain1[-c(1:burnin),weight_idx], 2, quantile, probs=0.975)
    sigma2_LB <- quantile(fit$chain1[-c(1:burnin),grep("^sigma",colnames(fit[[1]]))], probs=0.025)
    sigma2_UB <- quantile(fit$chain1[-c(1:burnin),grep("^sigma",colnames(fit[[1]]))], probs=0.975)
  }else if(model == "Bayes_Dir1" | model == "Bayes_lin_Dir1"){
    if(grepl("lin", model)){
      model_type <- "linear"
    }
    fit <- MCMC_sampler_m(x = x,
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
                          model_type = model_type,
                          weights_prior = rep(1,n_mod),
                          WAIC = T)
    fit_type <- "Bayes"
    ests <- pred_m(posterior_list = fit, 
                   burnin = burnin,
                   m_star = pred_seq)
    weight_idx <- grep("^weight",colnames(fit[[1]]))
    weight_ests <- colMeans(fit$chain1[-c(1:burnin),weight_idx])
    weight_PIPs <- colMeans(fit$chain1[-c(1:burnin),weight_idx]!=0)
    weight_LB <- apply(fit$chain1[-c(1:burnin),weight_idx], 2, quantile, probs=0.025)
    weight_UB <- apply(fit$chain1[-c(1:burnin),weight_idx], 2, quantile, probs=0.975)
    sigma2_LB <- quantile(fit$chain1[-c(1:burnin),grep("^sigma",colnames(fit[[1]]))], probs=0.025)
    sigma2_UB <- quantile(fit$chain1[-c(1:burnin),grep("^sigma",colnames(fit[[1]]))], probs=0.975)
  }else if(model == "Bayes_Dir0.1"| model == "Bayes_lin_Dir0.1"){
    if(grepl("lin", model)){
      model_type <- "linear"
    }
    fit <- MCMC_sampler_m(x = x,
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
                          model_type = model_type,
                          weights_prior = rep(0.1,n_mod),
                          WAIC = T)
    fit_type <- "Bayes"
    ests <- pred_m(posterior_list = fit, 
                   burnin = burnin,
                   m_star = pred_seq)
    weight_idx <- grep("^weight",colnames(fit[[1]]))
    weight_ests <- colMeans(fit$chain1[-c(1:burnin),weight_idx])
    weight_PIPs <- colMeans(fit$chain1[-c(1:burnin),weight_idx]!=0)
    weight_LB <- apply(fit$chain1[-c(1:burnin),weight_idx], 2, quantile, probs=0.025)
    weight_UB <- apply(fit$chain1[-c(1:burnin),weight_idx], 2, quantile, probs=0.975)
    sigma2_LB <- quantile(fit$chain1[-c(1:burnin),grep("^sigma",colnames(fit[[1]]))], probs=0.025)
    sigma2_UB <- quantile(fit$chain1[-c(1:burnin),grep("^sigma",colnames(fit[[1]]))], probs=0.975)
  }else if(model == "Bayes_DirInfoCorrect"){
    fit <- MCMC_sampler_m(x = x,
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
                          model_type = model_type,
                          weights_prior = c(rep(10,2), rep(1,n_mod-2)),
                          WAIC = F)
    fit_type <- "Bayes"
    ests <- pred_m(posterior_list = fit, 
                   burnin = burnin,
                   m_star = pred_seq)
    weight_idx <- grep("^weight",colnames(fit[[1]]))
    weight_ests <- colMeans(fit$chain1[-c(1:burnin),weight_idx])
    weight_PIPs <- colMeans(fit$chain1[-c(1:burnin),weight_idx]!=0)
    weight_LB <- apply(fit$chain1[-c(1:burnin),weight_idx], 2, quantile, probs=0.025)
    weight_UB <- apply(fit$chain1[-c(1:burnin),weight_idx], 2, quantile, probs=0.975)
    sigma2_LB <- quantile(fit$chain1[-c(1:burnin),grep("^sigma",colnames(fit[[1]]))], probs=0.025)
    sigma2_UB <- quantile(fit$chain1[-c(1:burnin),grep("^sigma",colnames(fit[[1]]))], probs=0.975)
  }else if(model == "Bayes_DirInfoMis"){
    fit <- MCMC_sampler_m(x = x,
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
                          model_type = model_type,
                          weights_prior = c(rep(1,n_mod-2),rep(10,2)),
                          WAIC = F)
    fit_type <- "Bayes"
    ests <- pred_m(posterior_list = fit, 
                   burnin = burnin,
                   m_star = pred_seq)
    weight_idx <- grep("^weight",colnames(fit[[1]]))
    weight_ests <- colMeans(fit$chain1[-c(1:burnin),weight_idx])
    weight_PIPs <- colMeans(fit$chain1[-c(1:burnin),weight_idx]!=0)
    weight_LB <- apply(fit$chain1[-c(1:burnin),weight_idx], 2, quantile, probs=0.025)
    weight_UB <- apply(fit$chain1[-c(1:burnin),weight_idx], 2, quantile, probs=0.975)
    sigma2_LB <- quantile(fit$chain1[-c(1:burnin),grep("^sigma",colnames(fit[[1]]))], probs=0.025)
    sigma2_UB <- quantile(fit$chain1[-c(1:burnin),grep("^sigma",colnames(fit[[1]]))], probs=0.975)
  }else if(model == "Bayes_Dir0.5" | model == "Bayes_lin_Dir0.5"){
    if(grepl("lin", model)){
      model_type <- "linear"
    }
    fit <- MCMC_sampler_m(x = x,
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
                          model_type = model_type,
                          weights_prior = rep(0.5,n_mod),
                          WAIC = T)
    fit_type <- "Bayes"
    ests <- pred_m(posterior_list = fit, 
                   burnin = burnin,
                   m_star = pred_seq)
    weight_idx <- grep("^weight",colnames(fit[[1]]))
    weight_ests <- colMeans(fit$chain1[-c(1:burnin),weight_idx])
    weight_PIPs <- colMeans(fit$chain1[-c(1:burnin),weight_idx]!=0)
    weight_LB <- apply(fit$chain1[-c(1:burnin),weight_idx], 2, quantile, probs=0.025)
    weight_UB <- apply(fit$chain1[-c(1:burnin),weight_idx], 2, quantile, probs=0.975)
    sigma2_LB <- quantile(fit$chain1[-c(1:burnin),grep("^sigma",colnames(fit[[1]]))], probs=0.025)
    sigma2_UB <- quantile(fit$chain1[-c(1:burnin),grep("^sigma",colnames(fit[[1]]))], probs=0.975)
  }else if(model == "Bayes_DLM"){
    
    fit <- Gibbs_sampler(x = x,
                         y = y,
                         z = z,
                         df_l = df_l,
                         niter = niter,
                         burnin = burnin,
                         tau2 = tau2,
                         xi2 = xi2,
                         a = a, 
                         b = b)
    
    fit_type <- "Bayes"
    ests <- NULL
    weight_ests <- NA
    weight_PIPs <- NA
    weight_LB <- NA
    weight_UB <- NA
    sigma2_LB <- NA
    sigma2_UB <- NA
    
  }else if(model == "Freq_ns"){
    m_star_p <- M%*%matrix(rep(1,ncol(M))/ncol(M), ncol=1)
    fit <- dlim(y = y,
                x = x,
                modifiers = m_star_p,
                z = z, 
                df_m = df_m,
                df_l = df_l,
                penalize = F,
                method = "REML")
    fit_type <- "Freq"
    pred <- predict(fit, pred_seq)
    ests <- pred$est_dlim
  }else if(model == "Freq_ps"){
    m_star_p <- M%*%matrix(rep(1,ncol(M))/ncol(M), ncol=1)
    fit <- dlim(y = y,
                x = x,
                modifiers = m_star_p,
                z = z, 
                df_m = 20,
                df_l = 20,
                penalize = T,
                method = "REML")
    fit_type <- "Freq"
    pred <- predict(fit, pred_seq)
    ests <- pred$est_dlim
  }
  
  ### Summarize ###
  
  if(fit_type == "Bayes"){
    if(is.null(attr(fit, "WAIC"))){
      WAIC <- NA
    }else{
      WAIC <- attr(fit, "WAIC")
    }
    if(i < 10 & model != "Bayes_DLM"){
      chains[[i]] <- fit$chain1
    }
  }else{
    WAIC <- NA
  }
  
  #weights
  if(fit_type=="Bayes"){
    weights_table <- rbind(weights_table, 
                           data.frame(Model = rep(model,n_mod),
                                      Weights = rep(weight_scenario,n_mod),
                                      SNR = rep(as.character(SNR),n_mod),
                                      Type = rep(as.character(type),n_mod),
                                      Weight_id = 1:n_mod,
                                      Estimate = weight_ests,
                                      Truth = true_weights,
                                      PIP = weight_PIPs))
    
  }
  
  
  #true betas
  true_betas <- t(sim_dlf(L = 36,
                          modifiers = pred_seq,
                          type = type))
  cumul_truth <- rowSums(true_betas)
  
  #CI Width
  width_pw <- mean(ests$UB - ests$LB)
  width_cumul <- mean(ests$cumul_UB - ests$cumul_LB)
  width_weight <- ifelse(fit_type=="Bayes", mean(weight_UB - weight_LB), NA)
  #make table
  table_results <- rbind(table_results, 
                         data.frame(Model = rep(model,3),
                                    Weights = rep(weight_scenario,3),
                                    SNR = rep(as.character(SNR),3),
                                    Type = rep(as.character(type),3),
                                    Parameter = c("PW", "Cumul", "Weights"),
                                    Metric = rep("Width",3),
                                    Value = c(width_pw, width_cumul, width_weight)))
  
  #coverage
  cover_pw <- mean(ests$UB > true_betas & ests$LB < true_betas)
  cover_cumul <- mean(ests$cumul_UB > cumul_truth & ests$cumul_LB < cumul_truth)
  cover_weight <- ifelse(fit_type=="Bayes", mean(weight_UB > true_weights & weight_LB < true_weights), NA)
  cover_sigma2 <- ifelse(fit_type=="Bayes", sigma2_LB < dta$noise2 & dta$noise2 < sigma2_UB, NA)
  #make table
  table_results <- rbind(table_results, 
                         data.frame(Model = rep(model,4),
                                    Weights = rep(weight_scenario,4),
                                    SNR = rep(as.character(SNR),4),
                                    Type = rep(as.character(type),4),
                                    Parameter = c("PW", "Cumul", "Weights", "sigma2"),
                                    Metric = rep("Coverage",4),
                                    Value = c(cover_pw, cover_cumul, cover_weight, cover_sigma2)))
  
  #RMSE
  RMSE_pw <- sqrt(mean((ests$betas - true_betas)^2))
  RMSE_cumul <- sqrt(mean((ests$betas_cumul - cumul_truth)^2))
  RMSE_weight <- ifelse(fit_type=="Bayes", sqrt(mean((weight_ests - true_weights)^2)), NA)
  #make table
  table_results <- rbind(table_results, 
                         data.frame(Model = rep(model,3),
                                    Weights = rep(weight_scenario,3),
                                    SNR = rep(as.character(SNR),3),
                                    Type = rep(as.character(type),3),
                                    Parameter = c("PW", "Cumul", "Weights"),
                                    Metric = rep("RMSE",3),
                                    Value = c(RMSE_pw, RMSE_cumul, RMSE_weight)))
  
  #WAIC
  table_results <- rbind(table_results, 
                         data.frame(Model = model,
                                    Weights = weight_scenario,
                                    SNR = SNR,
                                    Type = type,
                                    Parameter = "None",
                                    Metric = "WAIC",
                                    Value = WAIC))
  
  #save point-wise estimates for spaghetti plots
  all_betas[[i]] <- ests$betas
  
}



# } #comment out for HPC

#save results in list
sim_results <- list(table_results = table_results, 
                    weights_table = weights_table,
                    all_betas = all_betas)
attr(sim_results, "Type") <- type
attr(sim_results, "SNR") <- SNR
attr(sim_results, "true_weights") <- true_weights
attr(sim_results, "niter") <- niter
attr(sim_results, "burnin") <- burnin
attr(sim_results, "df_l") <- df_l
attr(sim_results, "df_m") <- df_m
attr(sim_results, "model_type") <- model_type

save(sim_results, file=paste0("simulation_results/sim_results_", array_idx, ".rda"))


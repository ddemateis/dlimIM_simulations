library(dlim)
library(MASS)
library(dlimIM)

## Run on HPC using the array command #SBATCH --array=0-80
## or remove the array_idx line below and loop array_idx over 0 to 80

#read in array index
array_idx <- as.numeric(commandArgs(trailingOnly=TRUE)) #read in array values from the shell file, which is just the task number

#MCMC specifications
niter <- 10000 
burnin <- 8000 
sims <- 1:200  

#create directory to save files
dir.create("simulation_results")
dir.create("simulation_traceplots")

#set seed
set.seed(080123)

#set up grids
type_grid <- 4
SNR_grid <- c(0.1, 0.5, 1) 
model_list <- list(Bayes_select = "Bayes_select",
                   Bayes_Dir1 = "Bayes_Dir1",
                   Freq_ps = "Freq_ps")
true_weight_list <- list(sparse_50 = c(rep(1,3)/3, rep(0,47)),#for first 50 data sets
                         sparse_50 = c(rep(1,3)/3, rep(0,47)),#for first 50 data sets
                         sparse_50 = c(rep(1,3)/3, rep(0,47)),#for first 50 data sets
                         sparse_50 = c(rep(1,3)/3, rep(0,47)),#for second 50 data sets
                         sparse_10 = c(c(0.5, 0.4, 0.1), rep(0,7)),
                         equal_3 = rep(1,3)/3,
                         diff_3 = c(0.5, 0.4, 0.1),
                         equal_10 = rep(1,10)/10,
                         diff_10 = c(0.3, 0.2, rep(0.1,2), rep(0.05,6)) #length 9
)
n_type <- length(type_grid)
n_SNR <- length(SNR_grid)
n_weight <- length(true_weight_list)
n_model <- length(model_list)
weight_scenario <- names(true_weight_list)

x <- exposure
n <- nrow(x)
model_type <- "ns"
df_m <- 5 #set to 20 for penalized DLIM below
df_l <- 5 #set to 20 for penalized DLIM below

#hyperparameters
tau2 <- 100
xi2 <- 110
a <- 1
b <- 0.001

#choose type, SNR, and weights
type <- type_grid[floor(array_idx/(n_SNR*n_weight*n_model))+1]
SNR <- SNR_grid[(floor(array_idx/n_weight/n_model))%%n_SNR+1]
model <- model_list[[floor(array_idx/n_weight)%%n_model+1]]
true_weights <- true_weight_list[[array_idx%%n_weight+1]]

weight_scenario <- names(true_weight_list)[[array_idx%%n_weight+1]]

#parameters
true_gammas <- round(runif(length(true_weights), -1, 1),1)
n_mod <- length(true_weights)


table_results <- c() 
weights_table <- c() 
all_betas <- vector(mode = "list", length = length(sims))

if(weight_scenario=="sparse_50"){ #can only run 50 sims for sparse case in a day on HPC
  if(array_idx%%4==0){
    sims <- 1:50
  }else if(array_idx%%4==1){
    sims <- 51:100
  }else if(array_idx%%4==2){
    sims <- 101:150
  }else if(array_idx%%4==3){
    sims <- 151:200
  }
}

chains <- list()

for(i in sims){
  
  #set seed specific to each data set
  set.seed(080123 + i)
  
  ### Fit ###
  
  #modifiers 
  M <- mvrnorm(n=n, 
               mu=rep(0, n_mod), 
               Sigma = diag(0.5, n_mod) + matrix(rep(0.5, n_mod^2), ncol=n_mod))
  M <- apply(M, 2, scale_01)
  
  #simulate data
  dta <- sim_data_m(x=x, 
                    M=M, 
                    w=true_weights,
                    SNR = SNR, 
                    type=type, 
                    ncovariates = 3, 
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
  }else if(model == "Freq_ps"){
    m_star_p <- M%*%matrix(rep(1,ncol(M))/ncol(M), ncol=1)
    fit <- dlim(y = y,
                x = x,
                modifiers = m_star_p,
                z = z, 
                df_m = 20,
                df_l = 20,
                penalize = T)
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

#save trace plots

if(fit_type == "Bayes"){
  pdf(paste0("simulation_traceplots/type", type, "_SNR", SNR, "_", model, "_", weight_scenario,".pdf"))
  for(p in 1:ncol(fit$chain1)){
    par(mfrow=c(3,3))
    for(i in 1:9){
      idx <- seq(burnin, niter, 1)
      plot(idx,chains[[i]][idx,p],type="l", ylab = colnames(fit$chain1)[p], xlab = "Iteration")
    }
  }
  dev.off()
}
#MCMC_samps is the posterior samples from the fit in analysis_script.R
library(dlimIM)
library(ggplot2)
library(gridExtra)

#load model fit for LDL from analysis_script_PROGRESS.R
load("ldl_mcmc_samps.rda") 
response_name <- "LDL"
response <- "ldl"

#prediction grid
m_star <- seq(0,1,0.01)
inf_idx <- seq(10000, 30000, 1)

                                ## Figure 13 ##

#cumulative plot
pdf(paste0(response,"_cumulative.pdf"), width=6, height=6)
print(plot_bdlim(x = MCMC_samps,
                 m_star = m_star,
                 sel = inf_idx,
                 type = "cumulative") + 
        xlab("Months since LMP") +
        ylab(paste("Cumulative change in", response_name, "\n per exposure unit")))
dev.off()

                                ## Figure 14 ##

#weight estimates
combined_post <- c()
for(ch in 1:length(MCMC_samps)){
  combined_post <- rbind(combined_post,
                         MCMC_samps[[ch]][inf_idx,])
}
weight_idx <- grep("^weight",colnames(combined_post))
est <- colMeans(combined_post[,weight_idx])

#modiifers
M <- attr(MCMC_samps, "M")

#prediction values, specific quantiles of estimate modifier index
est_m_star <- M%*%est
m_star_quarts <- round(quantile(est_m_star, probs=c(0.1, 0.25, 0.3, 0.4, 0.75, 0.9)),2)

#PW plots
pdf(paste0(response,"_pw.pdf"), width=6, height=6)
print(plot_bdlim(x = MCMC_samps,
                 m_star = m_star_quarts,
                 sel = inf_idx,
                 type = "by_modifier",
                 exp_times = -2:22) + 
        xlab("Months since LMP") +
        ylab(paste("Change in", response_name, "per exposure unit")))
dev.off()

                             ## Figure 15 ##

#load model fit for HDL from analysis_script_PROGRESS.R
load("hdl_mcmc_samps.rda") 
response_name <- "HDL"
response <- "hdl"

#cumulative plot
pdf(paste0(response,"_cumulative.pdf"), width=6, height=6)
print(plot_bdlim(x = MCMC_samps,
                 m_star = m_star,
                 sel = inf_idx,
                 type = "cumulative") + 
        xlab("Months since LMP") +
        ylab(paste("Cumulative change in", response_name, "\n per exposure unit")))
dev.off()

                            ## Figure 16 ##

#first run the single-modifier models, then plot

#load data
data96_clean <- read.csv("data96_clean2.csv") #change this to be the current version of data across all stages

#load stress modifiers
stress_data <- read.csv("MM_Sandra_Feb27_2023_stress.csv")

#merge data with stress modifiers
full_data <- merge(data96_clean, stress_data, by = "folio")

#combine trimester metrics into one stress measurement during pregnancy
full_data$epdsT <- rowSums(cbind(full_data$epds2T, full_data$epds3T), na.rm = T)
full_data$nleT <- rowSums(cbind(full_data$nle2T, full_data$nle3T), na.rm = T)
full_data$pssT <- rowSums(cbind(full_data$pss2T, full_data$pss3T), na.rm = T)
full_data$sum_staiT <- rowSums(cbind(full_data$sum_stai2T, full_data$sum_stai3T), na.rm = T)

#response
response <- c(HDL = "m_hdl_c48")
response_name <- "HDL"

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

dir.create("Stress")

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

set.seed(101124)

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

#load original multi-modifier model from analysis_script_PROGRESS.R
load("hdl_mcmc_samps.rda")

#predictions with multi-modifier model
m_star <- c( 0.25, 0.31, 0.45) #25, 40, 75th percentiles based on original index presented in main text
pred_full <- pred_m(MCMC_samps, sel = inf_idx, m_star = m_star)

#fit
for(i in 1:ncol(M)){
  modifier <- matrix(M[,i], ncol=1)
  MCMC_samps <- MCMC_sampler_m(x = x,
                               y = y,
                               M = modifier,
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
                               weights_prior = rep(5, 1))

  mod <- modifiers[i]

  #predict for single modifier model
  pred_sing <- pred_m(MCMC_samps, sel=inf_idx, m_star = m_star)
  
  #PW plots
  pw_df <- data.frame(Effect = c(t(pred_full$betas),t(pred_sing$betas)),
                      LB = c(t(pred_full$LB), t(pred_sing$LB)),
                      UB = c(t(pred_full$UB), t(pred_sing$UB)),
                      Week = rep(1:25, length(m_star)*2),
                      Index = rep(rep(m_star, each = 25),2),
                      Model = rep(c("Multi-modifier", "Single-modifier"), each=25*length(m_star)))
  
  pdf(paste0("sing_mod_stress_", mod,".pdf"), width = 6, height = 4)
  print(ggplot(pw_df, aes(x = Week, y = Effect, color = Model, fill = Model)) +
          geom_hline(yintercept = 0)+
          geom_ribbon(aes(ymin=LB, ymax=UB), alpha=0.5, color=FALSE)+
          geom_line(aes(linetype = Model))+
          facet_wrap(vars(Index), labeller = "label_both") +
          xlab("Months since LMP") +
          ylab("Change in HDL per exposure unit") +
          theme_classic() +
          theme(legend.position="bottom",
                text = element_text(size = 18)) +
          scale_fill_manual(values=c("#CDC1D1", "#C3DCD7"))+
          scale_color_manual(values=c("#421052", "#3BA68A")))
  dev.off()
}



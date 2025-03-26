#this file is to be used with the output from analysis_script_CO.R

library(dplyr)
library(xtable)
library(ggplot2)
library(dlimIM)

#load penalized DLIM-linear model
load("linear_fit.rda") #fit on newest data set

#view results from HPC
MCMC_samples <- c() #initialize posterior samples, post-burnin
files <- list.files("analysis_results/Dir1/df3")
for(file in files){
  load(paste0("analysis_results/Dir1/df3/",file))
  MCMC_samples <- rbind(MCMC_samples,
                        MCMC_samps$chain1[attr(MCMC_samps, "burnin"):nrow(MCMC_samps$chain1),])
  rm(MCMC_samps)
}

niter <- 50000
burnin <- 30000
inf_idx <- seq(burnin, niter, thin = 1)

                                ## Table 2 ##

#Estimates and PIPs
sens_pop <- c(Asthma = "As___P", 
              Cancer = "Cn__P", 
              Diabetes = "Db__P", 
              Heart_disease = "H____",
              Life_expectancy = "Lf__P", 
              Mental_health = "M___P", 
              Percent_64 = "P__64_", 
              Percent_5 = "P__5_")
demographics <- c(Housing_cost = "H___P", 
                  Disability = "Pr__P", 
                  HS= "P______", 
                  Lingustic_iso = "Prcnt_ln__", 
                  Low_income = "Prcnt_lw__", 
                  POC = "Prc____P")
weight_idx <- grep("^weight",colnames(MCMC_samples))
PIP <- colMeans(MCMC_samples[,weight_idx]!=0) 
est <- colMeans(MCMC_samples[,weight_idx])
sd <- apply(MCMC_samples[,weight_idx], 2, sd)
weight_df <- data.frame(est = est,
                        sd = sd,
                        PIP = PIP)
rownames(weight_df) <- c(names(sens_pop), names(demographics))
sorted_weight_df <- weight_df[order(weight_df$est, decreasing = TRUE),]
xtable(sorted_weight_df,3, digits=3)

                            ## Figure 2a ##

#load results from HPC
files <- list.files("analysis_results/Dir1/df3")
load(paste0("analysis_results/Dir1/df3/",files[1]))
MCMC_samples <- MCMC_samps
rm(MCMC_samps)
for(file in files[-1]){
  load(paste0("analysis_results/Dir1/df3/",file))
  MCMC_samples[[which(file == files)]] <- MCMC_samps$chain1
  rm(MCMC_samps)
}

#predict
m_star <- seq(15,80,1)
DLIM_pred <- predict(linear_fit, newdata = m_star)
pred <- pred_m(MCMC_samples, sel = inf_idx, m_star = m_star/100)

#cumulative plot
cumul_df <- data.frame(Effect = c(t(pred$betas_cumul), t(DLIM_pred$est_dlim$betas_cumul)),
                       LB = c(t(pred$cumul_LB), t(DLIM_pred$est_dlim$cumul_LB)),
                       UB = c(t(pred$cumul_UB), t(DLIM_pred$est_dlim$cumul_UB)),
                       HSF = rep(m_star, 2),
                       Model = rep(c("Proposed", "Fixed-index"), each=length(m_star)))

pdf("COES_cumulative.pdf", width=6, height=4)
ggplot(cumul_df, aes(x = HSF, y = Effect, color = Model, fill = Model))+
  geom_hline(yintercept = 0)+
  geom_ribbon(aes(ymin=LB, ymax=UB), alpha=0.5, color=FALSE)+
  geom_line(aes(linetype = Model))+
  xlab("HSF") +
  ylab("Cumulative change in BWGAZ \n per exposure unit") +
  theme_classic() +
  theme(legend.position="bottom") +
  scale_fill_manual(values=c("#CDC1D1", "#C3DCD7"))+
  scale_color_manual(values=c("#421052", "#3BA68A"))
dev.off()  

                           ## Figure 2b ##

mods <- c(32, 46, 59)
mod_idx <- which(m_star %in% mods)
pw_df <- data.frame(Effect = c(t(pred$betas[mod_idx,]), t(DLIM_pred$est_dlim$betas[mod_idx,])),
                    LB = c(t(pred$LB[mod_idx,]), t(DLIM_pred$est_dlim$LB[mod_idx,])),
                    UB = c(t(pred$UB[mod_idx,]), t(DLIM_pred$est_dlim$UB[mod_idx,])),
                    Week = rep(1:37, length(mod_idx)*2),
                    HSF = rep(rep(mods, each = 37),2),
                    Model = rep(c("Proposed", "Fixed-index"), each=37*length(mods)))
pdf("COES_PW.pdf", width=6, height=4)
ggplot(pw_df, aes(x = Week, y = Effect, color = Model, fill = Model)) +
  geom_hline(yintercept = 0)+
  geom_ribbon(aes(ymin=LB, ymax=UB), alpha=0.5, color=FALSE)+
  geom_line(aes(linetype = Model))+
  facet_wrap(vars(HSF), labeller = "label_both") +
  xlab("Week of gestation") +
  ylab("Change in BWGAZ per exposure unit") +
  theme_classic() +
  theme(legend.position="bottom") +
  scale_fill_manual(values=c("#CDC1D1", "#C3DCD7"))+
  scale_color_manual(values=c("#421052", "#3BA68A"))
dev.off()
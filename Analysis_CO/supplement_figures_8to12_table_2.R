#this file is to be used with the output from analysis_script_CO.R
#and analysis_script_CO_single_modifier_models.R


library(dplyr)
library(xtable)
library(ggplot2)
library(dlimIM)

niter <- 50000
burnin <- 30000
thin <- 1 
inf_idx <- seq(burnin, niter, thin)

#load 3 dfresults from HPC
files <- list.files("analysis_results/Dir1/df3")
load(paste0("analysis_results/Dir1/df3/",files[1]))
MCMC_samples <- MCMC_samps
rm(MCMC_samps)
for(file in files[-1]){
  load(paste0("analysis_results/Dir1/df3/",file))
  MCMC_samples[[which(file == files)]] <- MCMC_samps$chain1
  rm(MCMC_samps)
}

#load 4 df results from HPC
files <- list.files("analysis_results/Dir1/df4")
load(paste0("analysis_results/Dir1/df4/",files[1]))
MCMC_samples_4 <- MCMC_samps
rm(MCMC_samps)
for(file in files[-1]){
  load(paste0("analysis_results/Dir1/df4/",file))
  MCMC_samples_4[[which(file == files)]] <- MCMC_samps$chain1
  rm(MCMC_samps)
}

m_star <- seq(15,80,1)
pred <- pred_m(MCMC_samples, sel = inf_idx, m_star = m_star/100)
pred_df4 <- pred_m(MCMC_samples_4, sel = inf_idx, m_star = m_star/100)

                             ## Figure 8 ##

#cumulative plot
cumul_df <- data.frame(Effect = c(t(pred$betas_cumul), t(pred_df4$betas_cumul)),
                       LB = c(t(pred$cumul_LB), t(pred_df4$cumul_LB)),
                       UB = c(t(pred$cumul_UB), t(pred_df4$cumul_UB)),
                       HSF = rep(m_star, 2),
                       DF = rep(c("3", "4"), each=length(m_star)))

pdf("COES_cumulative_df4.pdf", width=6, height=4)
ggplot(cumul_df, aes(x = HSF, y = Effect, color = DF, fill = DF))+
  geom_hline(yintercept = 0)+
  geom_ribbon(aes(ymin=LB, ymax=UB), alpha=0.5, color=FALSE)+
  geom_line(aes(linetype = DF))+
  xlab("HSF") +
  ylab("Cumulative change in BWGAZ \n per exposure unit") +
  theme_classic() +
  theme(legend.position="bottom") +
  scale_fill_manual(values=c("darkgray", "lightgray"))+
  scale_color_manual(values=c("black", "black"))
dev.off() 

                             ## Figure 9 ##

mods <- c(32, 46, 59)
mod_idx <- which(m_star %in% mods)
pw_df <- data.frame(Effect = c(t(pred$betas[mod_idx,]), t(pred_df4$betas[mod_idx,])),
                    LB = c(t(pred$LB[mod_idx,]), t(pred_df4$LB[mod_idx,])),
                    UB = c(t(pred$UB[mod_idx,]), t(pred_df4$UB[mod_idx,])),
                    Week = rep(1:37, length(mod_idx)*2),
                    HSF = rep(rep(mods, each = 37),2),
                    DF = rep(c("3", "4"), each=37*length(mods)))
pdf("COES_PW_df4.pdf", width=6, height=4)
ggplot(pw_df, aes(x = Week, y = Effect, color = DF, fill = DF)) +
  geom_hline(yintercept = 0)+
  geom_ribbon(aes(ymin=LB, ymax=UB), alpha=0.5, color=FALSE)+
  geom_line(aes(linetype = DF))+
  facet_wrap(vars(HSF), labeller = "label_both") +
  xlab("Week of gestation") +
  ylab("Change in BWGAZ per exposure unit") +
  theme_classic() +
  theme(legend.position="bottom") +
  scale_fill_manual(values=c("darkgrey", "lightgrey"))+
  scale_color_manual(values=c("black", "black"))
dev.off()

                            ## Figures 10-12 ##

#multi-modifier model for comparison
load("analysis_results/Dir1/df3/analysis_MCMC_samps_Dir1_df3_chain1.rda") #new data
MCMC_samps_full <- MCMC_samps
m_star <- seq(15,80,1)
pred_full <- pred_m(MCMC_samps_full, sel = inf_idx, m_star = m_star/100)

#single-modifier models
for(i in 0:13){
  load(paste0("analysis_results/single_modifiers/analysis_MCMC_samps_dlimIM_modifier_", i, ".rda"))
  MCMC_samps_sing <- MCMC_samps
  rm(MCMC_samps)
  
  pred_sing <- pred_m(MCMC_samps_sing, sel = inf_idx, m_star = m_star/100)
  
  mods <- c(32, 46, 59)
  mod_idx <- which(m_star %in% mods)
  pw_df <- data.frame(Effect = c(t(pred_full$betas[mod_idx,]),t(pred_sing$betas[mod_idx,])),
                      LB = c(t(pred_full$LB[mod_idx,]), t(pred_sing$LB[mod_idx,])),
                      UB = c(t(pred_full$UB[mod_idx,]), t(pred_sing$UB[mod_idx,])),
                      Week = rep(1:37, length(mod_idx)*2),
                      HSF = rep(rep(mods, each = 37),2),
                      Model = rep(c("Multi-modifier", "Single-modifier"), each=37*length(mods)))
  
  pdf(paste0("single_modifier_COES_",i,".pdf"), width=6, height=4)
  print(ggplot(pw_df, aes(x = Week, y = Effect, color = Model, fill = Model)) +
          geom_hline(yintercept = 0)+
          geom_ribbon(aes(ymin=LB, ymax=UB), alpha=0.5, color=FALSE)+
          geom_line(aes(linetype = Model))+
          facet_wrap(vars(HSF), labeller = "label_both") +
          xlab("Week of gestation") +
          ylab("Change in BWGAZ \n per exposure unit") +
          theme_classic() +
          theme(legend.position="bottom",
                text = element_text(size = 20)) +
          scale_fill_manual(values=c("#CDC1D1", "#C3DCD7"))+
          scale_color_manual(values=c("#421052", "#3BA68A")))
  dev.off()
}

                               ## Table 2 ##

#3 df
#view results from HPC
MCMC_samples <- c() #initialize posterior samples, post-burnin
files <- list.files("analysis_results/Dir1/df3")
for(file in files){
  load(paste0("analysis_results/Dir1/df3/",file))
  MCMC_samples <- rbind(MCMC_samples,
                        MCMC_samps$chain1[attr(MCMC_samps, "burnin"):nrow(MCMC_samps$chain1),])
  rm(MCMC_samps)
}

#4 df
MCMC_samples_4 <- c() #initialize posterior samples, post-burnin
files <- list.files("analysis_results/Dir1/df4")
for(file in files){
  load(paste0("analysis_results/Dir1/df4/",file))
  MCMC_samples_4 <- rbind(MCMC_samples_4,
                          MCMC_samps$chain1[attr(MCMC_samps, "burnin"):nrow(MCMC_samps$chain1),])
  rm(MCMC_samps)
}

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

#3 df
weight_idx <- grep("^weight",colnames(MCMC_samples))
PIP <- colMeans(MCMC_samples[,weight_idx]!=0) 
est <- colMeans(MCMC_samples[,weight_idx])
sd <- apply(MCMC_samples[,weight_idx], 2, sd)
weight_df <- data.frame(est = est,
                        sd = sd,
                        PIP = PIP)

#4 df
weight_idx <- grep("^weight",colnames(MCMC_samples_4))
PIP <- colMeans(MCMC_samples_4[,weight_idx]!=0)
est <- colMeans(MCMC_samples_4[,weight_idx])
sd <- apply(MCMC_samples_4[,weight_idx], 2, sd)
df_df4 <- data.frame(est = est,
                     sd = sd,
                     PIP = PIP)
df <- data.frame(weight_df$est,
                 df_df4$est,
                 weight_df$sd,
                 df_df4$sd,
                 weight_df$PIP,
                 df_df4$PIP)
rownames(df) <- c(names(sens_pop), names(demographics))
xtable(round(df,3))
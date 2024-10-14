## Plot
library(ggplot2)
library(dlimIM)
library(xtable)

#MCMC_samps is the posterior samples for the DLIM-IM fit in 
#analysis_script.R

#load DLIM-IM fit
load("analysis_MCMC_samps_Dir1_df3_chain1.rda") #file saved from analysis_script.R
MCMC_samps_df3 <- MCMC_samps
load("analysis_MCMC_samps_Dir1_df4_chain1.rda") #file saved from analysis_script.R
MCMC_samps_df4 <- MCMC_samps #used for supplement weights table

#MCMC info
niter <- nrow(MCMC_samps$chain1)
burnin <- attr(MCMC_samps, "burnin")
thin <- 10
inf_ind <- seq(burnin, niter, thin)
m_star <- seq(15,80,1)

                      ### Supplemental Figure 2 ###
pred <- pred_m(MCMC_samps_df3, burnin = burnin, thin = thin, m_star = m_star/100)
pred_df4 <- pred_m(MCMC_samps_df4, burnin = burnin, thin = thin, m_star = m_star/100)

#cumulative plot
cumul_df <- data.frame(Effect = c(t(pred$betas_cumul), t(pred_df4$betas_cumul)),
                       LB = c(t(pred$cumul_LB), t(pred_df4$cumul_LB)),
                       UB = c(t(pred$cumul_UB), t(pred_df4$cumul_UB)),
                       HSF = rep(m_star, 2),
                       DF = rep(c("3", "4"), each=length(m_star)))

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

                      ### Supplemental Figure 3 ###

mods <- c(32, 46, 59)
mod_idx <- which(m_star %in% mods)
pw_df <- data.frame(Effect = c(t(pred$betas[mod_idx,]), t(pred_df4$betas[mod_idx,])),
                    LB = c(t(pred$LB[mod_idx,]), t(pred_df4$LB[mod_idx,])),
                    UB = c(t(pred$UB[mod_idx,]), t(pred_df4$UB[mod_idx,])),
                    Week = rep(1:37, length(mod_idx)*2),
                    HSF = rep(rep(mods, each = 37),2),
                    DF = rep(c("3", "4"), each=37*length(mods)))
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

                      ### Supplemental Table 2 ###

weight_idx <- grep("^weight",colnames(MCMC_samps_df3[[1]]))
PIP <- colMeans(MCMC_samps_df3$chain1[inf_ind,weight_idx]!=0)
est <- colMeans(MCMC_samps_df3$chain1[inf_ind,weight_idx])
sd <- apply(MCMC_samps_df3$chain1[inf_ind,weight_idx], 2, sd)
df_df3 <- data.frame(est = est,
                     sd = sd,
                     PIP = PIP)
weight_idx <- grep("^weight",colnames(MCMC_samps_df4[[1]]))
PIP <- colMeans(MCMC_samps_df4$chain1[inf_ind,weight_idx]!=0)
est <- colMeans(MCMC_samps_df4$chain1[inf_ind,weight_idx])
sd <- apply(MCMC_samps_df4$chain1[inf_ind,weight_idx], 2, sd)
df_df4 <- data.frame(est = est,
                     sd = sd,
                     PIP = PIP)
df <- data.frame(df_df3$est,
                 df_df4$est,
                 df_df3$sd,
                 df_df4$sd,
                 df_df3$PIP,
                 df_df4$PIP)
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
rownames(df) <- c(names(sens_pop), names(demographics))
xtable(round(df,3))

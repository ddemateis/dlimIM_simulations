library(ggplot2)
library(dlimIM)
library(xtable)
#MCMC_samps is the posterior samples for the DLIM-IM fit in 
#analysis_script.R

                                ### Figure 2 ###

#load DLIM-IM fit
load("analysis_MCMC_samps_Dir1_df3_chain1.rda") #file saved from analysis_script.R

#load DLIM fit
load("linear_fit.rda")#file saved from analysis_script.R

#MCMC info
niter <- nrow(MCMC_samps$chain1)
burnin <- attr(MCMC_samps, "burnin")
thin <- 10
inf_ind <- seq(burnin, niter, thin)

#predict
m_star <- seq(15,80,1)
DLIM_pred <- predict(linear_fit, newdata = m_star) #fit from DLIM-linear with fixed HSF
pred <- pred_m(MCMC_samps, burnin = burnin, thin = thin, m_star = m_star/100)

#cumulative plot
cumul_df <- data.frame(Effect = c(t(pred$betas_cumul), t(DLIM_pred$est_dlim$betas_cumul)),
                       LB = c(t(pred$cumul_LB), t(DLIM_pred$est_dlim$cumul_LB)),
                       UB = c(t(pred$cumul_UB), t(DLIM_pred$est_dlim$cumul_UB)),
                       HSF = rep(m_star, 2),
                       Model = rep(c("Proposed", "Fixed-index"), each=length(m_star)))

ggplot(cumul_df, aes(x = HSF, y = Effect, color = Model, fill = Model))+
  geom_hline(yintercept = 0)+
  geom_ribbon(aes(ymin=LB, ymax=UB), alpha=0.5, color=FALSE)+
  geom_line(aes(linetype = Model))+
  xlab("HSF") +
  ylab("Cumulative change in BWGAZ \n per exposure unit") +
  theme_classic() +
  theme(legend.position="bottom") +
  scale_fill_manual(values=c("grey", "lightblue"))+
  scale_color_manual(values=c("black", "blue"))

#PW plots
mods <- c(32, 46, 59)
mod_idx <- which(  m_star %in% mods)
pw_df <- data.frame(Effect = c(t(pred$betas[mod_idx,]), t(DLIM_pred$est_dlim$betas[mod_idx,])),
                    LB = c(t(pred$LB[mod_idx,]), t(DLIM_pred$est_dlim$LB[mod_idx,])),
                    UB = c(t(pred$UB[mod_idx,]), t(DLIM_pred$est_dlim$UB[mod_idx,])),
                    Week = rep(1:37, length(mod_idx)*2),
                    HSF = rep(rep(mods, each = 37),2),
                    Model = rep(c("Proposed", "Fixed-index"), each=37*length(mods)))
ggplot(pw_df, aes(x = Week, y = Effect, color = Model, fill = Model)) +
  geom_hline(yintercept = 0)+
  geom_ribbon(aes(ymin=LB, ymax=UB), alpha=0.5, color=FALSE)+
  geom_line(aes(linetype = Model))+
  facet_wrap(vars(HSF), labeller = "label_both") +
  xlab("Week of gestation") +
  ylab("Change in BWGAZ per exposure unit") +
  theme_classic() +
  theme(legend.position="bottom") +
  scale_fill_manual(values=c("grey", "lightblue"))+
  scale_color_manual(values=c("black", "blue"))

                                 ### Table 2 ###

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
weight_idx <- grep("^weight",colnames(MCMC_samps[[1]]))
PIP <- colMeans(MCMC_samps$chain1[inf_ind,weight_idx]!=0) 
est <- colMeans(MCMC_samps$chain1[inf_ind,weight_idx])
sd <- apply(MCMC_samps$chain1[inf_ind,weight_idx], 2, sd)
weight_df <- data.frame(est = est,
                        sd = sd,
                        PIP = PIP)
rownames(weight_df) <- c(names(sens_pop), names(demographics))
xtable(round(weight_df,3))

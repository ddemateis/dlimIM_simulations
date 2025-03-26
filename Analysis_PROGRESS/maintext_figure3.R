#MCMC_samps is the posterior samples from the fit in analysis_script.R
library(dlimIM)
library(ggplot2)

load("hdl_mcmc_samps.rda") #from analysis_script.R

                                ## Figure 3 ##

#MCMC samples for inference
inf_idx <- seq(10000, 30000, 1)

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

#plot
pdf(paste0("hdl_pw.pdf"), width=6, height=6)
print(plot_bdlim(x = MCMC_samps,
                 m_star = m_star_quarts,
                 sel = inf_idx,
                 type = "by_modifier",
                 exp_times = -2:22) + 
        xlab("Months since LMP") +
        ylab(paste("Change in HDL per exposure unit")))
dev.off()

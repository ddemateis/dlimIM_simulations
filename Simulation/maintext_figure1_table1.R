library(ggplot2)
library(dplyr)
library(xtable)
library(viridis)
library(MASS)
library(dlimIM)

#folder_name should be a folder containing saved simulation results
#from the simulation_script.R

                                ### Table 1 ###

#read in simulation data
weights_table <- c()
table_results <- c()
files <- list.files("simulation_results")
for(fl in files){
  load(paste0("simulation_results/", fl))
  
  table_results <- rbind(table_results, sim_results[[1]])
  weights_table <- rbind(weights_table, sim_results[[2]])
}

#post-processing

#m* post-processing
#post-processing
m_star_table <- c()
post_m_star_mean <- c()
type=4
for(weight in unique(weights_table$Weights)){
  for(snr in unique(weights_table$SNR)){
    for(model in c("Bayes_select", "Bayes_Dir1")){
      
      sub_data <- weights_table %>% filter(Weights == weight, 
                                           Model == model,
                                           SNR == snr,
                                           Type == type)
      
      n_mod <- length(unique(sub_data$Weight_id))
      true_weights <- sub_data$Truth[1:n_mod]
      
      for(i in 1:200){
        #set same seed to recover modifiers
        set.seed(080123 + i)
        
        #generate modifiers
        load(paste0("modifier_Sigma_",n_mod,".rda"))
        M <- mvrnorm(n=1000, 
                     mu=rep(0, n_mod), 
                     Sigma = Sigma)
        M <- apply(M, 2, scale_01)
        
        #construct true m_star
        true_m_star <- M %*% true_weights
        
        #obtain estimated weights
        est_weights <- sub_data$Estimate[(n_mod*(i-1)+1):(n_mod*(i-1)+1 + (n_mod-1))]
        freq_weights <- rep(1,n_mod)/n_mod
        est_m_star <- M%*%est_weights
        freq_m_star <- M%*%freq_weights
        
        #obtain mean/sd estimated m_star
        m_star_mean <- mean(est_m_star)
        m_star_sd <- sd(est_m_star)
        m_star_5 <- quantile(est_m_star, probs = 0.05)
        m_star_95 <- quantile(est_m_star, probs = 0.95)
        
        #bias and RMSE
        bias <- mean(abs(true_m_star - est_m_star))
        bias_f <- mean(abs(true_m_star - freq_m_star))
        rmse <- sqrt(mean((true_m_star - est_m_star)^2))
        rmse_f <- sqrt(mean((true_m_star - freq_m_star)^2))
        
        #add to table
        m_star_table <- rbind(m_star_table,
                              data.frame(Model = c(model, "Freq_ps"),
                                         Weights = rep(weight,2),
                                         SNR = rep(snr,2),
                                         Abs_bias = c(bias, bias_f),
                                         RMSE = c(rmse, rmse_f)))
        post_m_star_mean <- rbind(post_m_star_mean,
                                  data.frame(Model = c(model, "Freq_ps"),
                                             Weights = rep(weight,2),
                                             SNR = rep(snr,2),
                                             Metric = c("mean", "sd", "q5", "q95"),
                                             Value = c(m_star_mean, m_star_sd, m_star_5, m_star_95)))
        
      }#end for loop
    }#end model
  }#end SNR
}#end weights

#change weights to scenario
table_results$Scenario <- ifelse(table_results$Weights == "equal_3", "1", NA)
table_results$Scenario[table_results$Weights=="diff_3"] <- "2"
table_results$Scenario[table_results$Weights=="sparse_50"] <- "3"
table_results$Scenario <- factor(table_results$Scenario, levels = c("1", "2", "3", NA))

paper_table <- table_results %>% filter(Model == "Bayes_Dir1" | Model == "Bayes_select" | Model == "Freq_ps", 
                                        Weights == "equal_3" | Weights == "diff_3" | Weights == "sparse_50",
                                        Type == "4",
                                        Parameter == "Cumul" | Parameter == "PW",
                                        Metric == "RMSE" | Metric == "Coverage" | Metric == "Width")

#summary table for paper
cumul_rmse_table <- paper_table %>% filter(Parameter == "Cumul", Metric == "RMSE") %>% group_by(Weights, SNR, Model) %>%
  summarise(across(where(is.numeric), .fns = 
                     list(RMSE = mean))) 
cumul_rmse_table$Value_RMSE <- round(cumul_rmse_table$Value_RMSE,3)

cumul_cover_table <- paper_table %>% filter(Parameter == "Cumul",
                                            Metric == "Coverage") %>% group_by(Weights, SNR, Model) %>%
  summarise(across(where(is.numeric), .fns = 
                     list(Coverage = mean))) 
cumul_cover_table$Value_Coverage <- round(cumul_cover_table$Value_Coverage,2)

cumul_width_table <- paper_table %>% filter(Parameter == "Cumul",
                                            Metric == "Width") %>% group_by(Weights, SNR, Model) %>%
  summarise(across(where(is.numeric), .fns = 
                     list(Width = mean))) 
cumul_width_table$Value_Width <- round(cumul_width_table$Value_Width,3)


pw_rmse_table <- paper_table %>% filter(Parameter == "PW",
                                        Metric == "RMSE") %>% group_by(Weights, SNR, Model) %>%
  summarise(across(where(is.numeric), .fns = 
                     list(RMSE = mean))) 
pw_rmse_table$Value_RMSE <- round(pw_rmse_table$Value_RMSE,3)

pw_cover_table <- paper_table %>% filter(Parameter == "PW",
                                         Metric == "Coverage") %>% group_by(Weights, SNR, Model) %>%
  summarise(across(where(is.numeric), .fns = 
                     list(Coverage = mean))) 
pw_cover_table$Value_Coverage <- round(pw_cover_table$Value_Coverage,2)

pw_width_table <- paper_table %>% filter(Parameter == "PW",
                                         Metric == "Width") %>% group_by(Weights, SNR, Model) %>%
  summarise(across(where(is.numeric), .fns = 
                     list(Width = mean))) 
pw_width_table$Value_Width <- round(pw_width_table$Value_Width,3)

#construct table
sum_table <- data.frame(cumul_rmse_table,
                        Cumul_coverage = cumul_cover_table$Value_Coverage,
                        Cumul_width = cumul_width_table$Value_Width,
                        PW_RMSE = pw_rmse_table$Value_RMSE,
                        PW_coverage = pw_cover_table$Value_Coverage,
                        PW_width = pw_width_table$Value_Width)
#latex code for table
#xtable(sum_table)

#view the table
View(sum_table)

## additional first column of table for m* summaries ##

#summary table for m* (we include RMSE, this is the "Index" column in Table 1)
m_star_sum <- m_star_table %>% filter(Model == "Bayes_Dir1" | Model == "Bayes_select" | Model == "Freq_ps",
                                      Weights == "equal_3" | Weights == "diff_3" | Weights == "sparse_50") %>% group_by(Weights, SNR, Model) %>%
  summarise(across(where(is.numeric), .fns = 
                     list(Avg = mean))) 
View(m_star_sum)
#xtable(m_star_sum, digits=4)


                           ### Figure 1 ###

#read in simulation data
weights_table <- c()
table_results <- c()
files <- list.files("simulation_results")
for(fl in files){
  load(paste0("simulation_results/", fl))
  
  table_results <- rbind(table_results, sim_results[[1]])
  weights_table <- rbind(weights_table, sim_results[[2]])
}

#post-processing

#change weights to scenario
weights_table$Scenario <- ifelse(weights_table$Weights == "equal_3", "1", NA)
weights_table$Scenario[weights_table$Weights=="diff_3"] <- "2"
weights_table$Scenario[weights_table$Weights=="sparse_50"] <- "3"
weights_table$Scenario <- factor(weights_table$Scenario, levels = c("1", "2", "3", NA))

#process the weights in the sparse case, aggregate all zero weights
weights_table$Weight_id[weights_table$Weights=="sparse_50" & (weights_table$Weight_id != 1 & weights_table$Weight_id != 2 & weights_table$Weight_id != 3)] <- 4

#weights box plot
pdf("weights.pdf")
ggplot(weights_table %>% filter(Type == 4, 
                                Scenario == "1" | Scenario == "2" | Scenario == "3", 
                                Model == "Bayes_Dir1" |  Model == "Bayes_select"), aes(x = as.factor(Weight_id), y = Estimate, color = Model)) + geom_boxplot() + 
  geom_segment(aes(x=as.numeric(Weight_id) - 0.5,xend=as.numeric(Weight_id)+0.5,y=Truth,yend=Truth), color="black")+ 
  facet_grid(rows = vars(SNR), cols = vars(Scenario), labeller = "label_both", scale="free") + 
  ylab("Modifier Index Weight Estimate") + 
  xlab("Modifier Weight Number") + 
  theme_bw()+
  theme(legend.position="bottom") +
  scale_color_viridis(name = "", labels = c("Model without modifier selection", "Model with modifier selection"), discrete = TRUE, end = 0.5)
dev.off()

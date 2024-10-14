library(ggplot2)
library(dplyr)
library(xtable)

#folder_name should be a folder containing saved simulation results
#from the simulation_script.R

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
ggplot(weights_table %>% filter(Type == 4, 
                                Scenario == "1" | Scenario == "2" | Scenario == "3", 
                                Model == "Bayes_Dir1" |  Model == "Bayes_select"), aes(x = as.factor(Weight_id), y = Estimate, color = Model)) + geom_boxplot() + 
  geom_segment(aes(x=as.numeric(Weight_id) - 0.5,xend=as.numeric(Weight_id)+0.5,y=Truth,yend=Truth), color="black")+ 
  facet_grid(rows = vars(SNR), cols = vars(Scenario), labeller = "label_both", scale="free") + 
  ylab("Estimate") + 
  xlab("Weight Number") + 
  theme_bw()+
  theme(legend.position="bottom") +
  scale_color_discrete(name = "Selection", labels = c("No", "Yes"))


                                ### Table 1 ###

#post-processing

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

#summary table
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


sum_table <- data.frame(cumul_rmse_table,
                        Cumul_coverage = cumul_cover_table$Value_Coverage,
                        Cumul_width = cumul_width_table$Value_Width,
                        PW_RMSE = pw_rmse_table$Value_RMSE,
                        PW_coverage = pw_cover_table$Value_Coverage,
                        PW_width = pw_width_table$Value_Width)
xtable(sum_table)

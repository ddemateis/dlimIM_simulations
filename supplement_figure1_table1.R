#folder_name should be a folder containing saved simulation results
#from the simulation_script.R

                              ### Figure 1 ###

files <- list.files(paste0("HPC_results/", folder_name))
for(fl in files){
  load(paste0("HPC_results/", folder_name, "/", fl))
  weights_table <- rbind(weights_table, sim_results[[2]])
}

ggplot(weights_table %>% filter(Type == 4, 
                                Scenario == "1" | Scenario == "2" | Scenario == "3", 
                                Model == "Bayes_select"), aes(x = as.factor(Weight_id), y = PIP)) + 
geom_boxplot() + 
geom_segment(aes(x=as.numeric(Weight_id) - 0.5,xend=as.numeric(Weight_id)+0.5,y=ifelse(Truth!=0, 1, 0),yend=ifelse(Truth!=0, 1, 0)), color="red")+ 
geom_hline(yintercept = 0.5)+
facet_grid(rows = vars(SNR), cols = vars(Scenario), 
           labeller = "label_both", scale="free") +
ylab("PIP") + 
xlab("Weight Number") +
theme_bw()

                               ### Table 1 ###

#read in simulation data
weights_table <- c()
files <- list.files(paste0("HPC_results/", folder_name))
for(fl in files){
  load(paste0("HPC_results/", folder_name, "/", fl))
  
  table_results <- rbind(table_results, sim_results[[1]])
  weights_table <- rbind(weights_table, sim_results[[2]])
  WAIC_table <- rbind(WAIC_table, sim_results$WAIC)
}

#post-processing

#change weights to scenario
weights_table$Scenario <- ifelse(weights_table$Weights == "equal_3", "1", NA)
weights_table$Scenario[weights_table$Weights=="diff_3"] <- "2"
weights_table$Scenario[weights_table$Weights=="sparse_50"] <- "3"
weights_table$Scenario <- factor(weights_table$Scenario, levels = c("1", "2", "3", NA))

table_results$Scenario <- ifelse(table_results$Weights == "equal_3", "1", NA)
table_results$Scenario[table_results$Weights=="diff_3"] <- "2"
table_results$Scenario[table_results$Weights=="sparse_50"] <- "3"
table_results$Scenario <- factor(table_results$Scenario, levels = c("1", "2", "3", NA))


#draft summary table for supplement 
cumul_rmse_table <- table_results %>% filter(Model == "Bayes_Dir1" | Model == "Bayes_select" | Model == "Freq_ps",
                                             Type == "4",
                                             Weights == "equal_10" | Weights == "diff_10" | Weights == "sparse_10",
                                             Parameter == "Cumul", 
                                             Metric == "RMSE") %>% group_by(Weights, SNR, Model) %>%
  summarise(across(where(is.numeric), .fns = 
                     list(RMSE = mean))) 
cumul_rmse_table$Value_RMSE <- round(cumul_rmse_table$Value_RMSE,3)

cumul_cover_table <- table_results %>% filter(Model == "Bayes_Dir1" | Model == "Bayes_select" | Model == "Freq_ps",
                                              Type == "4",
                                              Weights == "equal_10" | Weights == "diff_10" | Weights == "sparse_10",
                                              Parameter == "Cumul", 
                                              Metric == "Coverage") %>% group_by(Weights, SNR, Model) %>%
  summarise(across(where(is.numeric), .fns = 
                     list(Coverage = mean))) 
cumul_cover_table$Value_Coverage <- round(cumul_cover_table$Value_Coverage,2)

cumul_width_table <- table_results %>% filter(Model == "Bayes_Dir1" | Model == "Bayes_select" | Model == "Freq_ps",
                                              Type == "4",
                                              Weights == "equal_10" | Weights == "diff_10" | Weights == "sparse_10",
                                              Parameter == "Cumul", 
                                              Metric == "Width") %>% group_by(Weights, SNR, Model) %>%
  summarise(across(where(is.numeric), .fns = 
                     list(Width = mean))) 
cumul_width_table$Value_Width <- round(cumul_width_table$Value_Width,3)


pw_rmse_table <- table_results %>% filter(Model == "Bayes_Dir1" | Model == "Bayes_select" | Model == "Freq_ps",
                                          Type == "4",
                                          Weights == "equal_10" | Weights == "diff_10" | Weights == "sparse_10",
                                          Parameter == "PW", 
                                          Metric == "RMSE") %>% group_by(Weights, SNR, Model) %>%
  summarise(across(where(is.numeric), .fns = 
                     list(RMSE = mean))) 
pw_rmse_table$Value_RMSE <- round(pw_rmse_table$Value_RMSE,3)

pw_cover_table <- table_results %>% filter(Model == "Bayes_Dir1" | Model == "Bayes_select" | Model == "Freq_ps",
                                           Type == "4",
                                           Weights == "equal_10" | Weights == "diff_10" | Weights == "sparse_10",
                                           Parameter == "PW", 
                                           Metric == "Coverage") %>% group_by(Weights, SNR, Model) %>%
  summarise(across(where(is.numeric), .fns = 
                     list(Coverage = mean))) 
pw_cover_table$Value_Coverage <- round(pw_cover_table$Value_Coverage,2)

pw_width_table <- table_results %>% filter(Model == "Bayes_Dir1" | Model == "Bayes_select" | Model == "Freq_ps",
                                           Type == "4",
                                           Weights == "equal_10" | Weights == "diff_10" | Weights == "sparse_10",
                                           Parameter == "PW", 
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
kable(sum_table)
View(sum_table)
xtable(sum_table)
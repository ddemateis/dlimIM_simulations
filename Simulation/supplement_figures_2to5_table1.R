library(dlimIM)
library(dplyr)
library(dlim)
library(ggplot2)

#folder_name should be a folder containing saved simulation results
#from the simulation_script.R

                              ### Figure 2 ###

folder_name <- "simulation_results"
weights_table <- c()

#read in all simulation files
files <- list.files(folder_name)
for(fl in files){
  load(paste0(folder_name, "/", fl))
  weights_table <- rbind(weights_table, sim_results[[2]])
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

#process the weights in the sparse case, aggregate all zero weights
weights_table$Weight_id[weights_table$Weights=="sparse_50" & (weights_table$Weight_id != 1 & weights_table$Weight_id != 2 & weights_table$Weight_id != 3)] <- 4

pip_table <- weights_table %>% filter(Type == 4, 
                                      Scenario == "1" | Scenario == "2" | Scenario == "3", 
                                      Model == "Bayes_select")
tmp_weight_id <- case_when(
  pip_table$Weight_id == "1" ~ "1",
  pip_table$Weight_id == "2" ~ "2",
  pip_table$Weight_id == "3" ~ "3",
  pip_table$Weight_id == "4-50" ~ "4")

#plot of posterior inclusion probabilities (PIPs)
pdf("PIPs.pdf")
ggplot(pip_table, aes(x = as.factor(Weight_id), y = PIP)) +
  geom_hline(yintercept = 0.5)+
  geom_boxplot() + 
  geom_segment(aes(x=as.numeric(tmp_weight_id) - 0.5,xend=as.numeric(tmp_weight_id)+0.5,y=ifelse(Truth!=0, 1, 0),yend=ifelse(Truth!=0, 1, 0)), color="red")+
  facet_grid(rows = vars(SNR), cols = vars(Scenario), 
             labeller = "label_both", scale="free") +
  ylab("PIP") + 
  xlab("Weight Number") +
  theme_bw()
dev.off()

                              ### Figures 3-5 ###

#data frame for saving plotting info
plot_df <- c() #Estimate, model, modifier, data set

#true DLF
truth <- sim_dlf(36, c(0.3, 0.5, 0.7), type=4)


#post-processing
for(fl in files){
  load(paste0(folder_name, "/", fl))
  
  #Modification Type, SNR, Model, Weight Scenario
  SNR <- attr(sim_results, "SNR")
  model <- unique(sim_results$table_results$Model)
  weights <- unique(sim_results$table_results$Weights)
  
  if(((model == "Bayes_Dir1" | model == "Freq_ps" | model == "Bayes_select") & SNR == 1 & (weights == "equal_3" | weights == "diff_3" | weights == "sparse_50")) & !(unique(paste(model, weights)) %in% unique(paste(plot_df$Model, plot_df$Scenario)))){
    
    #estimated DLF for 10 data sets, for modifiers 0.3, 0.5, 0.7
    for(i in 1:10){
      mod0_3 <- sim_results$all_betas[[i]][18,]
      mod0_5 <- sim_results$all_betas[[i]][51,]
      mod0_7 <- sim_results$all_betas[[i]][84,]
      
      if(!is.null(mod0_3) & !is.null(mod0_5) & !is.null(mod0_7)){
        #save to df
        plot_df <- rbind(plot_df,
                         data.frame(Estimate = c(mod0_3, 
                                                 mod0_5, 
                                                 mod0_7),
                                    Time = rep(1:length(mod0_3), 3),
                                    Model = model,
                                    Scenario = weights,
                                    Modifier = c(rep(0.3, length(mod0_3)),
                                                 rep(0.5, length(mod0_5)),
                                                 rep(0.7, length(mod0_7))),
                                    Dataset = as.character(i)))
      }
      
    }
  }
}

plot_df$Model <- case_when(
  plot_df$Model %in% "Bayes_Dir1" ~ "DLIM-IM",
  plot_df$Model %in% "Bayes_select" ~ "DLIM-IM sel",
  plot_df$Model %in% "Freq_ps" ~ "Fixed Index"
)

#create plots
for(scenario in unique(plot_df$Scenario)){
  pdf(paste0("spaghetti_,", scenario,".pdf"))
  print(ggplot(plot_df %>% filter(Scenario == scenario), aes(x=Time, y = Estimate, fill = Dataset)) + 
          geom_line(color="gray") + 
          geom_line(aes(x = Time, y = rep(truth,30)), color = "black")+
          facet_grid(Model ~ Modifier, labeller = "label_both")) 
  dev.off()
}


                               ### Table 1 ###

#read in simulation data
table_results <- c()
files <- list.files(folder_name)
for(fl in files){
  load(paste0(folder_name, "/", fl))
  table_results <- rbind(table_results, sim_results[[1]])
}

#post-processing

#change weights to scenario
table_results$Scenario <- ifelse(table_results$Weights == "equal_3", "1", NA)
table_results$Scenario[table_results$Weights=="diff_3"] <- "2"
table_results$Scenario[table_results$Weights=="sparse_50"] <- "3"
table_results$Scenario <- factor(table_results$Scenario, levels = c("1", "2", "3", NA))

#summary table for supplement 
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
View(sum_table)
#library(xtable)
#xtable(sum_table)
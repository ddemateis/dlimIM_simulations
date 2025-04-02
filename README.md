# dlimIM_simulations
The repository contains scripts to recreate the simulations and analyses in "Distributed Lag Interaction Model with Index Modification" submitted to _Biostatistics_. There are 3 folders, one containing simulation scripts and the others containing analysis scripts (one for each respective analysis).  

Simulation folder:

First, run the modifier_correlation_structures.R script. Then, run the simulation.R script on a HPC using the submit_all.sh script (in terminal, run $sbatch submit_all.sh). The rest of the files in this folder are used to create figures/tables for this simulation. 

Analysis_CO folder: 

First, run analysis_script_CO.R on a HPC using the submit_analysis.sh script (in terminal, run $sbatch submit_analysis.sh). Now, you can use maintext_figure2_table2.R to create figures/tables. For additional supplemental figures, run analysis_script_CO_single_modifier_models.R on a HPC using the submit_single_modifier.sh script (in terminal, fun $sbtach submit_single_modifier.sh). Now, the supplement_figures_8to12_table2.R will be able to create all supplemental figures/tables. 

Analysis_PROGRESS folder:

First, run analysis_script_PROGRESS.R. The rest of the files in this folder are used to create figures/tables for this simulation. 

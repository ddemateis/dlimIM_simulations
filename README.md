# dlimIM_simulations
The repository contains scripts to recreate the simulations and analyses in "Distributed Lag Interaction Model with Index Modification" submitted to Biostatistics. Files are named according to which figure/table they recreate. 

The analysis_script_PROGRESS.R, analysis_script_CO, and simulation_script.R contain code to replicate the two analyses (Colorado and PROGRESS cohorts) and simulation results presented in the main text and supplemental text. The other files contain code to recreate the figures and tables in the main text (file name starts with "maintext") or in the supplemental text (file starts with "supplement"). Run the analysis_script.R or simulation_script.R before running figure or table scripts. 

To run the simulation_script.R, use the submit_all.sh script on a HPC using the sbatch command. This is the main simulation script. 

To run the analysis_script_CO.R, use the submit_analysis.sh script on a HPC using the sbatch command. This is the main CO analysis script.

To run the analysis_script_CO_single_modifier_models.R, use the submit_single_modifier.sh script on a HPC using the sbatch command. This is the supplemental CO analysis script for single modifier models. 

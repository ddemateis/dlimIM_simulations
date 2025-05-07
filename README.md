# dlimIM_simulations
The repository contains scripts to recreate the simulations and analyses in "Distributed Lag Interaction Model with Index Modification" in _Biostatistics_. 

There are 3 folders, one containing simulation scripts and the others containing analysis scripts (one for each respective analysis).  
Information on data sharing and processing is available in the README files in the Analysis_CO/DataPrep and Analysis_PROGRESS folders.

## Simulation folder:

First, run the modifier_correlation_structures.R script. Then, run the simulation.R script on a HPC using the submit_all.sh script (in terminal, run $sbatch submit_all.sh). The rest of the files in this folder are used to create figures/tables for this simulation. 

## Analysis_CO folder: 

First, run analysis_script_CO.R on a HPC using the submit_analysis.sh script (in terminal, run $sbatch submit_analysis.sh). Now, you can use maintext_figure2_table2.R to create figures/tables. For additional supplemental figures, run analysis_script_CO_single_modifier_models.R on a HPC using the submit_single_modifier.sh script (in terminal, fun $sbtach submit_single_modifier.sh). Now, the supplement_figures_8to12_table2.R will be able to create all supplemental figures/tables. 

Information on obtaining data and creating the analysis file is in the folder Analysis_CO/DataPrep. The README file in that folder contains additional details.

## Analysis_PROGRESS folder:

First, run analysis_script_PROGRESS.R. The rest of the files in this folder are used to create figures/tables for this simulation. 

Information on data sharing for the PROGESS data is available in the README file in the folder Analysis_PROGRESS.



## Packages

This uses the following packages:
- dlimIM version 0.1.1
- dlim version 0.1.0
- mgcv version 1.9.1
- dplyr version 1.1.4
- xtable version 1.8.4
- ggplot2 version 3.5.2
- gridExtra version 2.3
- data.table version 1.17.0
- tidyverse version 2.0.0
- tidycensus version 1.7.1
- elevatr version 0.99.0
- tigris version 2.2.1
- sf version 1.0.20
- ncdf4 version 1.24

All packages are on CRAN with the exception of `dlimIM` which is available at [github.com/ddemateis/dlimIM](https://github.com/ddemateis/dlimIM) and [DOI:10.5281/zenodo.15354392](doi.org/10.5281/zenodo.15354392).



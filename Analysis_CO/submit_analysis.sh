#!/bin/bash

#SBATCH --partition=amilan 							# Partition (a.k.a.queue) to use
#SBATCH --qos=long								# Quality of service (either normal, long, or mem)
#SBATCH --nodes=1 									# Total number of nodes (a.k.a. servers) requested
#SBATCH --mem=15G 									# Memory allocation
#SBATCH --job-name=analysis							# Job name
#SBATCH --output=Logs/%x_%A_%a.out					# Name of stdout output file (%A jobId, %a array index, %x job name)
#SBATCH --error=Logs/%x_%A_%a.err					# Name of error file (%A jobId, %a array index, %x job name)
#SBATCH --time 6-23:59:00 							# Run time (days-hh:mm:ss)
#SBATCH --mail-type=ALL 							# Send emails on start, end and failure
#SBATCH --mail-user=Danielle.Demateis@colostate.edu # Address for sending emails
#SBATCH --array=0-9							 	# Array indices for true weights list 

echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

# load R
module load R/4.4.0

# Launch a serial job
echo "Starting @ "`date`
Rscript analysis_script_CO.R $SLURM_ARRAY_TASK_ID
echo "Completed @ "`date`
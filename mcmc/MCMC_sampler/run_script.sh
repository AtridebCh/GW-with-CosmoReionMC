#!/bin/bash
# run_script.sh - bash wrapper for Condor

# REQUIRED for your system
source /usr/share/Modules/init/bash
# load modules
module load anaconda3/2022.10 

# run your Python script
python /Users/users/achatterjee/Data/GW_with_CosmoReionMC/mcmc/MCMC_sampler/MCMC_run.py
#echo "Hello wold"

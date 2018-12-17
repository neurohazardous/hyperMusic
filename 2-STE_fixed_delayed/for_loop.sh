#!/bin/bash
# Simplest case - using for loop to submit a serial farm
# The input file table.dat contains individual cases - one case per line

export TABLE=table.dat

# Total number of cases (= number of jobs to submit):
N_cases=$(cat "$TABLE" | wc -l)

# Submitting one job per case using the for loop:
for ((i=1; i<=$N_cases; i++))
  do
  # Using environment variable I_FOR to communicate the case number to individual jobs:
  export I_FOR=$i
  sbatch job_script.sh
  done
  

#!/bin/bash
# Bundling up multiple short cases (simplest - static - workload distribution).
# The input file table.dat contains individual cases - one case per line
# Four environment variables are passed to the job script - TABLE, N_CASES, N_JOBS, I0.

export TABLE=table.dat

# Target bundle factor (number of cases per job; can be a bit bigger for individual jobs):
N_bundle=3

# Total number of cases (= number of jobs to submit):
export N_CASES=$(cat "$TABLE" | wc -l)

# Number of jobs (rounded to the smaller integer value - this ensures that the actual bundle factor is never smaller than $N_bundle):
export N_JOBS=$(($N_CASES / $N_bundle))

echo "Submitting $N_JOBS jobs..."

for ((i=1; i<=$N_JOBS; i++))
  do
  export I0=$i
  sbatch job_script.sh
  done
  

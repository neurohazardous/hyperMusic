#!/bin/bash
# Here you should provide the sbatch arguments to be used in all jobs in this serial farm
#SBATCH -t 0-04:30
#SBATCH --mem=6000
#SBATCH -A def-ibruce

# Simple way to achieve the best static workload distribution of $N_CASES cases over $N_JOBS jobs, 
# with the number of cases per job equal or larger that the target $N_bundle:
for ((i=$I0; i<=$N_CASES; i=i+$N_JOBS))
  do
  # Extracing the $i-th line from file $TABLE:
  LINE=`sed -n ${i}p "$TABLE"`

  # Echoing the command (optional), with the case number prepended:
  echo "$i; $LINE"

  # Executing the command:
  eval "$LINE"
  done
  

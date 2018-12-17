#!/bin/bash
# Here you should provide the sbatch arguments to be used in all jobs in this serial farm
#SBATCH -t 0-00:20
#SBATCH --mem=5000
#SBATCH -A def-ibruce

# Extracing the $I_FOR-th line from file $TABLE:
LINE=`sed -n ${I_FOR}p "$TABLE"`

# Echoing the command (optional), with the case number prepended:
echo "$I_FOR; $LINE"

# Executing the command:
eval "$LINE"

#!/usr/bin/env bash

#SBATCH --job-name=pyhmmer_array
#SBATCH --account=grp-org-sc-mds
#SBATCH --qos=jgi_normal
#SBATCH --time=16:00:00
#SBATCH --cpus-per-task=64
#SBATCH --array=0-39
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100G

readarray -t LIST < no_kofam.txt
FILE=${LIST[$SLURM_ARRAY_TASK_ID]}

srun python pyhmmer_array_search_all.py ${FILE} && echo ${FILE} >> RAN_PYHMMER

#!/usr/bin/env bash

#SBATCH --job-name=GENOMAD
#SBATCH --account=grp-org-sc-mds
#SBATCH --qos=jgi_normal
#SBATCH --time=30:00:00
#SBATCH --mem=256G
#SBATCH --nodes=12
#SBATCH --ntasks-per-node=1

srun --no-kill --ntasks=12 --wait=0 ./antonio_payload.sh genomad_list.txt 

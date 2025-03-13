#!/usr/bin/env bash

#SBATCH --job-name=pyhmmer
#SBATCH --account=grp-org-sc-mds
#SBATCH --qos=jgi_normal
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=64
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100G

srun jgi-trace python pyhmmer_array_search_all.py hmmsearch_outputs_5kb/new_isolates_and_rescues.faa

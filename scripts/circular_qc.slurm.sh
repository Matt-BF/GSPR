#!/usr/bin/env bash

#SBATCH --job-name=circular_qc
#SBATCH --account=grp-org-sc-mds
#SBATCH --mail-type=ALL
#SBATCH --mail-user mbfiamenghi@lbl.gov
#SBATCH --qos=jgi_normal
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G

srun --ntasks=1 seqkit seq assemblies_organized_5kb/*.gz | ./circular_sequence_qc.py  > qc.tsv

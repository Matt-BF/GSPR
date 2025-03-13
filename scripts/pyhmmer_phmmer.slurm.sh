#!/usr/bin/env bash

#SBATCH --job-name=phmmer
#SBATCH --account=grp-org-sc-mds
#SBATCH --mail-type=ALL
#SBATCH --mail-user mbfiamenghi@lbl.gov
#SBATCH --qos=jgi_normal
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50G

srun --ntasks=1 python pyhmmer_phmmer.py genomad_outputs_5kb/part_001/part_001_annotate/part_001_proteins.faa ../dbs/iceberg/all_ICEs.faa 


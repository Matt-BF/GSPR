#!/bin/bash

if [[ -z "${SLURM_NODEID}" ]]; then
    echo "need \$SLURM_NODEID set"
    exit
fi
if [[ -z "${SLURM_NNODES}" ]]; then
    echo "need \$SLURM_NNODES set"
    exit
fi
cat $1 |                                               \
awk -v NNODE="$SLURM_NNODES" -v NODEID="$SLURM_NODEID" \
'NR % NNODE == NODEID' |                               \
rush -j 1 'genomad end-to-end assemblies_organized/{} genomad_outputs/{^.fna.gz} ../dbs/genomad_db -s 4.8 --cleanup --enable-score-calibration --relaxed && echo {} >> DONE'

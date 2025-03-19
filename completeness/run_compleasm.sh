#!/bin/bash -l
#PBS -N casm
#PBS -l walltime=03:30:00
#PBS -l mem=35G
#PBS -l ncpus=8

cd $PBS_O_WORKDIR

#####

# Specify compleasm directory
COMPLEASMDIR=/home/stewarz2/various_programs/compleasm_kit

# Specify lineage to search
LINEAGE=vertebrata

# Specify computational parameters
CPUS=8

#####

COMPLEASMDB=${COMPLEASMDIR}/db

# Run clustering evaluation
FOLDER=vertebrates
RUNDIR=${FOLDER}/analysis/run_6d08041e9e70fa297f6d
BINGEFILE=${RUNDIR}/BINge_clustering_result.tsv;
RUNNAME=$(basename ${RUNDIR});
PROTEINFILE=${FOLDER}/representatives/${RUNNAME}/BINge_clustering_representatives.aa;
OUTDIR=${FOLDER}/busco;
mkdir -p ${OUTDIR};
python ${COMPLEASMDIR}/compleasm.py protein -p ${PROTEINFILE} \
    -L ${COMPLEASMDB} -l ${LINEAGE} -t ${CPUS} \
    -o ${OUTDIR};


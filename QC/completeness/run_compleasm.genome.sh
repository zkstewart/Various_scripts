#!/bin/bash -l
#PBS -N casm
#PBS -l walltime=06:00:00
#PBS -l mem=55G
#PBS -l ncpus=12

cd $PBS_O_WORKDIR

#####

# Specify compleasm directory
COMPLEASMDIR=/home/stewarz2/various_programs/compleasm_kit
COMPLEASMDB=${COMPLEASMDIR}/db

# Specify protein FASTA file location
FASTA=/location/to/genome.fasta

# Specify output location
OUTDIR=compleasm_genome_results

# Specify program parameters
LINEAGE=eudicotyledons
CPUS=12

#####

mkdir -p ${OUTDIR};
python ${COMPLEASMDIR}/compleasm.py run \
-a ${FASTA} \
-L ${COMPLEASMDB} \
-l ${LINEAGE} \
-t ${CPUS} \
-o ${OUTDIR};

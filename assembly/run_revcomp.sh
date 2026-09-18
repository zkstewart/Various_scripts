#!/bin/bash -l
#PBS -N revcomp
#PBS -l walltime=12:00:00
#PBS -l mem=60G
#PBS -l ncpus=8
#PBS -P <>

cd $PBS_O_WORKDIR

#################################

# Specify the location of the Various_scripts repository
VARSCRIPTDIR=/home/stewarz2/scripts/Various_scripts

# Specify the location of the reference genome
REFGENOME=/work/ePGL/sequencing/dna/nanopore/citrus/glauca_shearing_tests/diploid/Pete_Bush_28_10_2024/reference_assembly/hap1.final.fasta

# Specify how many CPUs to use
CPUS=8

# Specify behavioural params
OUTDIR=output
PRESET=asm20

#################################

# STEP 0: Specify the location of the input files
ASSEMBLYDIR=$(realpath ../scaffolding/output)
HAP1=${ASSEMBLYDIR}/hap1.fasta
HAP2=${ASSEMBLYDIR}/hap2.fasta

# STEP 1: Run the pipeline
python ${VARSCRIPTDIR}/assembly/revcomp_to_match_reference.py -i ${HAP1} ${HAP2} \
    -r ${REFGENOME} \
    -o ${OUTDIR} \
    --threads ${CPUS} --preset ${PRESET}

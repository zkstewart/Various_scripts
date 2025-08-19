#!/bin/bash -l
#PBS -N BINge
#PBS -l walltime=12:00:00
#PBS -l mem=50G
#PBS -l ncpus=12

cd $PBS_O_WORKDIR

####

# Specify the location of the BINge directory
BINGEDIR=/home/stewarz2/scripts/BINge

# Specify the locations of the genome and annotation files
GENOMEDIR=/home/stewarz2/plant_group/anuradha/upr_annotation/citrus_genomes
SEQDIR=/home/stewarz2/plant_group/anuradha/upr_annotation/liftoff

# Specify computational and behavioural parameters
THREADS=12
IDENTITY=0.85 # 85% identity for MMseqs unbinned clustering

# Specify the species to be analysed
PREFIX=binge_results

####

# STEP 1: Initialise directory
python ${BINGEDIR}/BINge.py initialise -d ${PREFIX} \
    -i ${GENOMEDIR}/species_genome.gff3,${GENOMEDIR}/species_genome.fasta \
    --ix ${SEQDIR}/species_sequences.mrna,${SEQDIR}/species_sequences.cds,${SEQDIR}/species_sequences.aa \
    --threads ${THREADS}

# STEP 2: Run BINge clustering [--gmapIdentity test]
python ${BINGEDIR}/BINge.py cluster -d ${PREFIX} \
    --threads ${THREADS} \
    --identity ${IDENTITY} \
    --debug;

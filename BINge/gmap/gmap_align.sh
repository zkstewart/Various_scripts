#!/bin/bash -l
#PBS -N gmapA6
#PBS -l walltime=48:00:00
#PBS -l mem=20G
#PBS -l ncpus=18

cd $PBS_O_WORKDIR

# Specify the location of the genome FASTA
## Note: Indexing is assumed to have happened, and the .gmap index is located at GENDIR
GENDIR=/home/stewarz2/banana_group/genomes
GENFILE=Musa_acuminata_pahang_v4.fasta

# Specify the location of the query FASTA file
QDIR=/home/stewarz2/banana_group/qcav_dge/transcriptome/results
QFILE=qcav_dge_okay-okalt.fasta

# Specify program behavioural parameters
F=2 # output format; 2 == GFF3
NUMPATHS=6 # number of places to align the sequence to
X=50 # chimera margin; amount unaligned sequence that triggers search for remainder
B=5 # batch mode; 5 == fast but memory intensive
INTRONMIDDLE=500000 # max allowed length for internal intron
INTRONENDS=50000 # max allowed length of first or last intron

# Specify computational parameters
CPUS=18

# Specify output prefix
PREFIX=qcav

#####

# > STEP 1: Automatically generate output prefix
PREFIX=${PREFIX}_to_A.${NUMPATHS}paths.gmap

# > STEP 2: Run GMAP aligner
gmap -D ${GENDIR} -d ${GENFILE}.gmap \
    -f ${F} -n ${NUMPATHS} -x ${X} \
    -B 5 -t ${CPUS} \
    --max-intronlength-middle=${INTRONMIDDLE} \
    --max-intronlength-ends=${INTRONENDS} \
    ${QDIR}/${QFILE} > ${PREFIX}.gff3

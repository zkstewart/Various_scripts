#!/bin/bash -l
#PBS -N gmapAlign
#PBS -l walltime=12:00:00
#PBS -l mem=20G
#PBS -l ncpus=12

cd $PBS_O_WORKDIR

#####

# Specify the location of the genome FASTA
## Note: Indexing is assumed to have happened, and the .gmap index is located at GENDIR
GENDIR=/home/stewarz2/telmatactis
GENFILE=telmatactis_HGAP.arr4.pil2.fasta

# Specify the location of the query FASTA file
QDIR=/home/stewarz2/telmatactis/gene_models/transcriptomes/evidentialgene/concatenated
QFILE=tel_hgap_okay-okalt.cds

# Specify program behavioural parameters
F=2 # output format; 2 == GFF3
NUMPATHS=12 # number of places to align the sequence to
X=50 # chimera margin; amount unaligned sequence that triggers search for remainder
B=5 # batch mode; 5 == fast but memory intensive
INTRONMIDDLE=500000 # max allowed length for internal intron
INTRONENDS=100000 # max allowed length of first or last intron

# Specify computational parameters
CPUS=12

# Specify the output prefix
PREFIX=telmatactis

#####


# > STEP 1: Automatically generate output prefix
PREFIX=${TXFILE}_n${NUMPATHS}_gmap

# > STEP 2: Run GMAP aligner
gmap -D ${GENDIR} -d ${GENFILE}.gmap \
    -f ${F} -n ${NUMPATHS} -x ${X} \
    -B 5 -t ${CPUS} \
    --max-intronlength-middle=${INTRONMIDDLE} \
    --max-intronlength-ends=${INTRONENDS} \
    ${QDIR}/${QFILE} > ${PREFIX}.gff3

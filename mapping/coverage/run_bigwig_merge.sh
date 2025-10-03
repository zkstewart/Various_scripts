#!/bin/bash -l
#PBS -N bigwigmerge
#PBS -l walltime=04:00:00
#PBS -l mem=40G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

####

# Specify the bedGraphToBigWig exe location
BGMERGE_EXE=/home/stewarz2/various_programs/ucsc-twobit-toolkit/bigWigMerge
BGTOBW_EXE=/home/stewarz2/various_programs/ucsc-twobit-toolkit/bedGraphToBigWig

# Specify the genome lengths file location
GENOMEDIR=/home/stewarz2/citrus/andrew_miles/powerpole/mapping/reference
GENOMELENGTH=Power_both_polished.fasta.lengths.tsv

# Specify the suffix to identify files by
BIGWIGSUFFIX=.bigwig

# Specify the prefix for the merged output
PREFIX=toxin1

####

# STEP 1: Merge bigwigs to bedgraph
${BGMERGE_EXE} *${BIGWIGSUFFIX} ${PREFIX}.cov

# STEP 2: Convert bedgraph back to bigwig
${BGTOBW_EXE} ${PREFIX}.cov ${GENOMEDIR}/${GENOMELENGTH} ${PREFIX}.bigwig;

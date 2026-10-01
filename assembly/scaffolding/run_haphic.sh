#!/bin/bash -l
#PBS -N haphic
#PBS -l walltime=12:00:00
#PBS -l mem=55G
#PBS -l ncpus=8

cd $PBS_O_WORKDIR

conda activate haphic

####

# Specify the HapHiC install location
HAPHICDIR=/home/stewarz2/various_programs/HapHiC

# Specify the hifiasm unitigs file
## should be a .gfa and .fasta following the prefix
UTGPREFIX=/work/ePGL/sequencing/dna/nanopore/citrus/NGS_743/minneola/assembly/no_hic/minneola.asm.bp.p_utg

# Specify the Hi-C mapping BAM file
BAMFILE=/scratch/stewarz2/assembly/minneola/mapping/NGS743_Minneola_S7.haphic.bam

# Specify how many scaffolds should be assembled
## This should be nchr * ploidy
NUMSCAFFOLDS=18 # 9 * 2

# Specify computational parameters
## I assume these run in addition, so THREAD+PROC==ncpus ?
THREAD=1
PROC=1
INFLATION=5.0

####

${HAPHICDIR}/haphic pipeline ${UTGPREFIX}.fasta ${BAMFILE} ${NUMSCAFFOLDS} --threads ${THREAD} --processes ${PROC} --gfa ${UTGPREFIX}.gfa --max_inflation ${INFLATION}

## Or, run it manually
# ${HAPHICDIR}/haphic cluster $UTGFASTA $BAMFILE 18 --min_inflation 1.1 --max_inflation 10.0
#!/bin/bash -l
#PBS -N gff2gtf
#PBS -l walltime=02:00:00
#PBS -l mem=30G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

## MANUAL SETUP BELOW

# >> SETUP: Specify conda env where Perl can be run from
CONDAENV=perl5

# >> SETUP: Specify AGAT bin location
AGATDIR=/home/stewarz2/various_programs/AGAT/bin

# >> SETUP: Specify gene annotation GFF3 location
GFFDIR=/home/stewarz2/plant_group/plant_rnaseq/annotations
GFFFILE=Min.fixed.gff3

# >> SETUP: Specify output file name
OUTFILE=Min.fixed.gtf

## MANUAL SETUP END

## RUN PROGRAM
# STEP 1: Activate conda env
conda activate ${CONDAENV}

# STEP 2: Run agat_convert_sp_gff2gtf.pl
${AGATDIR}/agat_convert_sp_gff2gtf.pl --gff ${GFFDIR}/${GFFFILE} -o ${OUTFILE}

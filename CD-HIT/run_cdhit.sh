#!/bin/bash -l
#PBS -N cdhit
#PBS -l walltime=78:00:00
#PBS -l mem=30G
#PBS -l ncpus=8

cd $PBS_O_WORKDIR

####

# SETUP: Specify the location of CD-HIT executable
CDHITDIR=/home/stewarz2/various_programs/cdhit-4.8.1
CDHITEXE=cd-hit
## Note: cd-hit-est for nucleotides, cd-hit for proteins

# SETUP: Specify the location of the input FASTA file
FASTADIR=/home/stewarz2/anemones/cassie/transcriptome/contaminant_removal/symbionts/cdhit
FASTA=concat_symbionts.fasta

# SETUP: Specify the output file prefix
PREFIX=symbionts

# SETUP: Specify CD-HIT parameters
IDENTITY=0.97
N=5
## Note: Make sure identity and N values are compatible
CPUS=8
MEM=25000
aS=0.9
aL=0.0

####


# STEP 1: Run CD-HIT
${CDHITDIR}/${CDHITEXE} -i ${FASTADIR}/${FASTA} -o ${PREFIX}_c${IDENTITY}_aS${aS}_aL${aL}.fasta -c ${IDENTITY} -n ${N} -d 0 -M ${MEM} -T ${CPUS} -aS ${aS} -aL ${aL}


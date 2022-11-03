#!/bin/bash -l
#PBS -N psy_cassie
#PBS -l walltime=100:00:00
#PBS -l mem=30G
#PBS -l ncpus=18

cd $PBS_O_WORKDIR

####

module load libsvm/3.22-intel-2016b

# SETUP: Specify psytrans dir
PSYTRANSDIR=/home/stewarz2/various_programs/psytrans

# SETUP: Specify host and symbiont files
HOSTDIR=/home/stewarz2/anemones/cassie/transcriptome/contaminant_removal/anemones/cdhit
HOSTFASTA=anemones_c0.95_aS0.9_aL0.0.fasta

SYMDIR=/home/stewarz2/anemones/cassie/transcriptome/contaminant_removal/symbionts/cdhit
SYMFASTA=symbionts_c0.95_aS0.9_aL0.0.fasta

# SETUP: Specify target transcriptome file
TARGETDIR=/home/stewarz2/anemones/cassie/transcriptome/transcriptomes/evidentialgene/concatenated
TARGETFASTA=cassie_okay-okalt.cds

# SETUP: Specify program parameters and output prefix
PREFIX=cassie
CPUS=18
TMP_DIR=main_results


####


mkdir -p ${TMP_DIR}
python2 ${PSYTRANSDIR}/psytrans.py ${TARGETDIR}/${TARGETFASTA} -A ${HOSTDIR}/${HOSTFASTA} -B ${SYMDIR}/${SYMFASTA} -p ${CPUS} -t ${TMP_DIR}

#!/bin/bash -l
#PBS -N norm
#PBS -l walltime=32:00:00
#PBS -l mem=800G
#PBS -l ncpus=16

cd $PBS_O_WORKDIT

module load Java/17.0.6
conda activate perl5

####

# Specify program location
TRINITYDIR=/home/stewarz2/various_programs/trinityrnaseq-v2.15.1

# Specify reads location
READSDIR=/home/stewarz2/citrus/andrew_miles/powerpole/reads/rna/toxin1
R1SUFFIX=.trimmed_1P.fq.gz
R2SUFFIX=.trimmed_2P.fq.gz

# Specify computational requirements
CPUS=16
MEM=700G

####

# STEP 1: Locate all reads files for argument formatting
declare -a R1FILES
i=0
for f in ${READSDIR}/*${R1SUFFIX}; do
    R1FILES[${i}]=$(echo "${f}");
    i=$((i+1));
done

declare -a R2FILES
i=0
for f in ${READSDIR}/*${R2SUFFIX}; do
    R2FILES[${i}]=$(echo "${f}");
    i=$((i+1));
done

# STEP 2: Run Trinity insilico normalisation
${TRINITYDIR}/util/insilico_read_normalization.pl \
    --seqType fq --max_cov 30 \
    --JM ${MEM} --CPU ${CPUS} \
    --left $(ls ${READSDIR}/*${R1SUFFIX} | paste -sd "," -) \
    --right $(ls ${READSDIR}/*${R2SUFFIX} | paste -sd "," -) \
    --pairs_together --PARALLEL_STATS \
    2>&1 >> Trinity.log

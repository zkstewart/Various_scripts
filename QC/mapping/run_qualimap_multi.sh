#!/bin/bash -l
#PBS -N multiqual
#PBS -l walltime=06:00:00
#PBS -l mem=15G
#PBS -l ncpus=1
#PBS -W depend=afterok:

cd $PBS_O_WORKDIR

conda activate qualimap

####

DIRSUFFIX=_bamqc
OUTDIR=multi_bamqc
MEM=10G

TSVFILE=qualimap_samples.tsv

####

# STEP 1: Find qualimap output directories
declare -a PREFIXES
declare -a QMDIRS
i=0
for DIR in *${DIRSUFFIX}; do
    PREFIXES[${i}]=$(echo "${DIR%%${DIRSUFFIX}}");
    QMDIRS[${i}]=$(echo "${DIR}");
    i=$((i+1));
done

# STEP 2: Generate a TSV file for input to qualimap
rm ${TSVFILE}
touch ${TSVFILE}
for x in ${!PREFIXES[@]}; do
    echo "${PREFIXES[$x]}	${QMDIRS[$x]}" >> ${TSVFILE};
done

# STEP 3: Run the multi-QC step
qualimap multi-bamqc -d ${TSVFILE} -outdir ${OUTDIR} --java-mem-size=${MEM}

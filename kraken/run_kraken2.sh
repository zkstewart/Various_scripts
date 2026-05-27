#!/bin/bash -l
#PBS -N k2classify
#PBS -l walltime=01:00:00
#PBS -l mem=160G
#PBS -l ncpus=8
#PBS -J 1-X

cd $PBS_O_WORKDIR

conda activate kraken2 # should have bioconda::krakentools and bioconda::krona in this environment

####

# Specify the database directory
DBDIR=/work/ePGL/databases/kraken2/citrus_plants_contaminants

# Specify reads dir
READSDIR=/scratch/stewarz2/devindee_mapping/rlsite/trimmed_reads
R1SUFFIX=.trimmed_1P.fq.gz
R2SUFFIX=.trimmed_2P.fq.gz

# Specify number of threads to use
CPUS=8

####

# STEP 1: Find file prefixes
declare -a RNAFILES
i=0
for f in ${READSDIR}/*${R1SUFFIX}; do
    RNAFILES[${i}]=$(echo "${f%%${R1SUFFIX}}");
    i=$((i+1));
done

# STEP 2: Get job details
ARRAY_INDEX=$((${PBS_ARRAY_INDEX}-1))
FILEPREFIX=${RNAFILES[${ARRAY_INDEX}]}
BASEPREFIX=$(basename ${FILEPREFIX})

# STEP 3: Classify the reads
mkdir -p ${BASEPREFIX}
k2 classify --threads ${CPUS} --db ${DBDIR} \
    --output ${BASEPREFIX}/${BASEPREFIX}_output.txt \
    --unclassified-out "${BASEPREFIX}/${BASEPREFIX}#_unclassified.fastq" \
    --classified-out "${BASEPREFIX}/${BASEPREFIX}#_classified.fastq" \
    --report ${BASEPREFIX}/${BASEPREFIX}_report.txt \
    --paired ${FILEPREFIX}${R1SUFFIX} ${FILEPREFIX}${R2SUFFIX}

# STEP 4: Generate a visual of the results
python $(which kreport2krona.py) \
    -r ${BASEPREFIX}/${BASEPREFIX}_report.txt \
    -o ${BASEPREFIX}/${BASEPREFIX}_report.krona

ktImportText ${BASEPREFIX}/${BASEPREFIX}_report.krona \
             -o ${BASEPREFIX}/${BASEPREFIX}_report.html

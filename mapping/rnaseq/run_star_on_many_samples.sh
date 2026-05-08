#!/bin/bash -l
#PBS -N star
#PBS -l walltime=04:00:00
#PBS -l mem=20G
#PBS -l ncpus=4
#PBS -J 1-X

cd $PBS_O_WORKDIR

module load STAR/2.7.11b

####

# Specify reference index directory
INDEXDIR=/work/ePGL/genomes/mango/indica/CATAS_Mindica_2.1/STAR_gtf_index

# Specify reads location & file suffix
READSDIR=/scratch/stewarz2/mapping_steph/trimmed_reads
SUFFIX=P.fq.gz

# Specify computational resources
CPUS=4

####

# STEP 1: Locate all read prefixes for mapping
declare -a RNAPREFIXES
i=0
for f in ${READSDIR}/*1${SUFFIX}; do
    RNAPREFIXES[${i}]=$(echo "${f%%1${SUFFIX}}");
    i=$((i+1));
done

# STEP 2: Handle batch submission variables
ARRAY_INDEX=$((${PBS_ARRAY_INDEX}-1))
PREFIX=${RNAPREFIXES[${ARRAY_INDEX}]}
BASENAME=$(basename ${PREFIX} _) # strip any _ trailing characters

# STEP 3: Setup directory for output files
mkdir -p ${BASENAME}
cd ${BASENAME}

# STEP 4: Run mapping procedure for each sample
## We set up cmds in an array so we can conditionally add something to it
cmd=(STAR --runThreadN ${CPUS} \
          --genomeDir ${INDEXDIR} \
          --readFilesIn ${PREFIX}1${SUFFIX} ${PREFIX}2${SUFFIX} \
          --outSAMtype BAM SortedByCoordinate \
          --outSAMunmapped Within \
          --quantMode TranscriptomeSAM GeneCounts)
if [[ ${SUFFIX} == *.gz ]] ; then
	cmd+=(--readFilesCommand zcat);
fi;
"${cmd[@]}"

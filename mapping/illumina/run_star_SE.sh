#!/bin/bash -l
#PBS -N star_SE
#PBS -l walltime=48:00:00
#PBS -l mem=50G
#PBS -l ncpus=1
#PBS -J 1-12

cd $PBS_O_WORKDIR

## MANUAL SETUP BELOW

# >> SETUP: Specify STAR executable location
STARDIR=/home/stewarz2/various_programs/STAR-2.7.10a/bin/Linux_x86_64_static

# >> SETUP: Specify genome location
GENDIR=/home/stewarz2/plant_group/leena/genome
GENFILE=NbLab360.genome.fasta

# >> SETUP: Specify gene annotation GTF location
## If you have a GFF3, convert it to GTF with agat_convert_sp_gff2gtf.pl
GTFDIR=/home/stewarz2/plant_group/leena/annotation
GTFFILE=NbLab360.v103.gtf

# >> SETUP: Specify reads location & file suffix
READSDIR=/home/stewarz2/plant_group/leena/trimmomatic
SUFFIX=.trimmed.fq

# >> SETUP: Specify computational resources
CPUS=1

# >> SETUP: Specify output prefix
OUTPREFIX=leena

## MANUAL SETUP END

## RUN PROGRAM
# STEP 1: Locate all reads for mapping
declare -a RNAFILES
i=0
for f in ${READSDIR}/*${SUFFIX}; do
    RNAFILES[${i}]=$(echo "${f}");
    i=$((i+1));
done

# STEP 2: Handle batch submission variables
ARRAY_INDEX=$((${PBS_ARRAY_INDEX}-1))
FILE=${RNAFILES[${ARRAY_INDEX}]}
BASENAME=$(basename ${FILE} ${SUFFIX})

# STEP 3: Setup directory for output files
mkdir -p ${BASENAME}
cd ${BASENAME}

# STEP 4: Run mapping procedure for each sample
${STARDIR}/STAR --runThreadN ${CPUS} \
	--genomeDir ${PBS_O_WORKDIR} \
	--readFilesIn ${FILE} \
	--outSAMtype BAM SortedByCoordinate \
	--outSAMunmapped Within \
	--quantMode TranscriptomeSAM GeneCounts

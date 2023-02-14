#!/bin/bash -l
#PBS -N fbayes_p
#PBS -l walltime=72:00:00
#PBS -l mem=55G
#PBS -l ncpus=18

cd $PBS_O_WORKDIR

module unload gcc/4.9.3-2.25
module unload gcccore/4.9.3
module unload binutils/2.25-gcccore-4.9.3
module unload zlib/1.2.8-foss-2016a
module load gcc/10.3.0

module unload libxml2/2.9.3-foss-2016a
module load libxml2/2.9.10-gcccore-10.3.0

#################################

# Specify the location of Freebayes and vcflib
FBDIR=/home/stewarz2/various_programs/freebayes_built/freebayes/build # should be executable called 'freebayes' in this dir
FBSCRIPT=/home/stewarz2/various_programs/freebayes_built/freebayes/scripts
VCFLIBDIR=/home/stewarz2/various_programs/freebayes_built/freebayes/vcflib/bin

# Specify the location of the genome FASTA
GENOMEDIR=/home/stewarz2/citrus/genome
GENOME=citrus.fasta

# Specify the location of the mapped BAM files
MAPDIR=/home/stewarz2/citrus/map

# Specify the suffix that identifies mapped BAM files
SUFFIX=.sorted.md.bam

# Specify computational parameters
CPUS=18
MAXDEPTH=500

# Specify output file prefix
PREFIX=citrus

#################################

# > STEP 1: Make sure freebayes and vcflib are in our path
export PATH="${VCFLIBDIR}:$PATH"
export PATH="${FBDIR}:$PATH"
export PATH="${FBSCRIPT}:$PATH"

# > STEP 2: Get our file list
declare -a BAMFILES
i=0
for f in ${MAPDIR}/*${SUFFIX}; do
    BAMFILES[${i}]=$(echo "${f}");
    i=$((i+1));
done

# > STEP 3: Get our input files argument
SEPARATOR=" "
INPUT_ARG="$( printf "${SEPARATOR}%s" "${BAMFILES[@]}" )"

# > STEP 4: Run parallelised freebayes
${FBSCRIPT}/freebayes-parallel <(${FBSCRIPT}/fasta_generate_regions.py ${GENOMEDIR}/${GENOME}.fai 100000) ${CPUS} \
    -f ${GENOMEDIR}/${GENOME} -g ${MAXDEPTH} ${INPUT_ARG} > ${PREFIX}.vcf

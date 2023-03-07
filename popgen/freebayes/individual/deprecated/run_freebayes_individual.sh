#!/bin/bash -l
#PBS -N fbayes_p
#PBS -l walltime=72:00:00
#PBS -l mem=55G
#PBS -l ncpus=1
#PBS -J 1-186

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
CPUS=1
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

# > STEP 3: Get our array index
declare -i index
index=${PBS_ARRAY_INDEX}-1

# > STEP 4: Get our file for analysis
INPUTFILE=${BAMFILES[${index}]}

# > STEP 5: Get our output file prefix
PREFIX=${INPUTFILE%%.sorted.bam}

# > STEP 4: Run freebayes
if [[ ! -f ${PREFIX}.vcf ]]; then
    ${FBDIR}/freebayes -f ${GENOMEDIR}/${GENOME} ${MAPDIR}/${INPUTFILE} > ${PREFIX}.vcf;
fi

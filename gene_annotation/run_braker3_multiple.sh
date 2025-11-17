#!/bin/bash -l
#PBS -N h1b3
#PBS -l walltime=48:00:00
#PBS -l mem=50G
#PBS -l ncpus=16
#PBS -J 1-26

cd $PBS_O_WORKDIR

SCRATCHDIR=/scratch/stewarz2

####

# Specify BRAKER3 singularity .sif location
B3SIF=/home/stewarz2/various_programs/BRAKER-3/braker3.sif

# Specify protein evidence file location
PROTSEQ=/work/ePGL/sequencing/dna/nanopore/citrus/glauca_shearing_tests/diploid/results/braker3/proteins/citrus.plus.viridiplantae.proteins.fasta

# Specify the directory of softmasked genome files
TARGETDIR=/work/ePGL/sequencing/dna/nanopore/citrus/glauca_shearing_tests/diploid/results/repeats/hap1
TARGETSUFFIX=_hap1.fasta.masked

# Specify RNAseq input details
RNADIR=/work/ePGL/genomes/external/citrus/australasica/trimmed_rnaseq
R1SUFFIX=.trimmed_1P.fq.gz

# Specify haplotype
HAP=hap1

# Specify parameters
CPUS=16
BUSCOLINEAGE=eudicots

####

# STEP 0: Get the RNAseq file prefixes
RNAIDS=""
DELIMITER=""
for f in ${RNADIR}/*${R1SUFFIX}; do
    RNAIDS="${RNAIDS}${DELIMITER}$(basename "${f}" "${R1SUFFIX}")" 
    DELIMITER=',' 
done

# STEP 1: Locate all softmasked genome files
declare -a TARGETFILES
i=0
for f in ${TARGETDIR}/*/*/*${TARGETSUFFIX}; do
    TARGETFILES[${i}]=$(echo "${f%%${TARGETSUFFIX}}");
    i=$((i+1));
done

# STEP 2: Get job details
ARRAY_INDEX=$((${PBS_ARRAY_INDEX}-1))
TARGETPREFIX=${TARGETFILES[${ARRAY_INDEX}]}
BASEPREFIX=$(basename ${TARGETPREFIX})

# STEP 3: Move to and setup scratch location
WORKDIR=${SCRATCHDIR}/${BASEPREFIX}_${HAP}_braker3
mkdir -p ${WORKDIR}
cd ${WORKDIR}

# STEP 4: Set up writable location for Augustus config files
mkdir -p augustus
singularity exec --bind /work/ePGL,${WORKDIR} ${B3SIF} cp -R /opt/Augustus/config .
mv config/ augustus/

# STEP 5: Run BRAKER3 pipeline
singularity exec --bind /work/ePGL,${SCRATCHDIR},${WORKDIR}/augustus/config:/augustus_config ${B3SIF} braker.pl --genome=${TARGETPREFIX}${TARGETSUFFIX} --species=${BASEPREFIX}_${HAP} \
    --rnaseq_sets_dir=${RNADIR} --rnaseq_sets_ids=${RNAIDS} \
    --prot_seq=${PROTSEQ} --busco_lineage=${BUSCOLINEAGE} \
    --threads=${CPUS} --gff3 --nocleanup \
    --AUGUSTUS_CONFIG_PATH=/augustus_config

# STEP 6: Move results back to origin
mv braker/braker.gff ${PBS_O_WORKDIR}/${BASEPREFIX}_${HAP}.braker.gff3
mv braker/braker.gtf ${PBS_O_WORKDIR}/${BASEPREFIX}_${HAP}.braker.gtf

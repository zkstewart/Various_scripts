#!/bin/bash -l
#PBS -N b3Aus
#PBS -l walltime=40:00:00
#PBS -l mem=80G
#PBS -l ncpus=24

cd $PBS_O_WORKDIR

SCRATCHDIR=/scratch/stewarz2

####

# Specify BRAKER3 singularity .sif location
B3SIF=/home/stewarz2/various_programs/BRAKER-3/braker3.sif

# Specify protein evidence file location
PROTSEQ=/work/ePGL/genomes/external/citrus/databases/citrus.plus.viridiplantae.proteins.fasta

# Specify the softmasked genome file
GENOMEDIR=/work/ePGL/genomes/external/citrus/australasica/repeats/Fingerlime_genome_softmask
GENOME=Fingerlime_genome.fna.masked

# Specify RNAseq input details
RNADIR=/work/ePGL/genomes/external/citrus/australasica/trimmed_rnaseq
R1SUFFIX=.trimmed_1P.fq.gz

# Specify parameters
CPUS=24
BUSCOLINEAGE=eudicots

# Specify output prefix
PREFIX=Fingerlime_genome

####

# STEP 0: Get the RNAseq file prefixes
RNAIDS=""
DELIMITER=""
for f in ${RNADIR}/*${R1SUFFIX}; do
    RNAIDS="${RNAIDS}${DELIMITER}$(basename "${f}" "${R1SUFFIX}")" 
    DELIMITER=',' 
done

# STEP 1: Move to and setup scratch location
WORKDIR=${SCRATCHDIR}/${PREFIX}_braker3
mkdir -p ${WORKDIR}
cd ${WORKDIR}

# STEP 4: Set up writable location for Augustus config files
mkdir -p augustus
singularity exec --bind /work/ePGL,${WORKDIR} ${B3SIF} cp -R /opt/Augustus/config .
mv config/ augustus/

# STEP 5: Run BRAKER3 pipeline
singularity exec --bind /work/ePGL,${SCRATCHDIR},${WORKDIR}/augustus/config:/augustus_config ${B3SIF} braker.pl --genome=${GENOMEDIR}/${GENOME} --species=${PREFIX} \
    --rnaseq_sets_dir=${RNADIR} --rnaseq_sets_ids=${RNAIDS} \
    --prot_seq=${PROTSEQ} --busco_lineage=${BUSCOLINEAGE} \
    --threads=${CPUS} --gff3 --nocleanup \
    --AUGUSTUS_CONFIG_PATH=/augustus_config

# STEP 6: Move results back to origin
mv braker/braker.gff ${PBS_O_WORKDIR}/${PREFIX}.braker.gff3
mv braker/braker.gtf ${PBS_O_WORKDIR}/${PREFIX}.braker.gtf

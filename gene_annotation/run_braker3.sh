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

# Specify input evidence file locations
RNASEQ=/work/ePGL/sequencing/dna/nanopore/citrus/glauca_shearing_tests/diploid/results/braker3/reads
PROTSEQ=/work/ePGL/sequencing/dna/nanopore/citrus/glauca_shearing_tests/diploid/results/braker3/proteins/citrus.plus.viridiplantae.proteins.fasta

# Specify the directory of softmasked genome files
TARGETDIR=/work/ePGL/sequencing/dna/nanopore/citrus/glauca_shearing_tests/diploid/results/repeats/hap1
TARGETSUFFIX=_hap1.fasta.masked

# Specify RNAseq IDs
RNAIDS="NGS-588-101_S94,NGS-588-102_S95,NGS-588-104_S97,NGS-588-37_S36,NGS-588-39_S37,NGS-588-40_S38,NGS-588-41_S39,NGS-588-42_S40,NGS-588-43_S41,NGS-588-44_S42,NGS-588-47_S43,NGS-588-48_S44,NGS-588-49_S45,NGS-588-50_S46,NGS-588-51_S47,NGS-588-52_S48,NGS-588-53_S49,NGS-588-54_S50,NGS-588-55_S51,NGS-588-56_S52,NGS-588-57_S53,NGS-588-58_S54,NGS-588-59_S55,NGS-588-60_S56,NGS-588-61_S57,NGS-588-62_S58,NGS-588-63_S59,NGS-588-64_S60,NGS-588-65_S61,NGS-588-66_S62,NGS-588-67_S63,NGS-588-68_S64,NGS-588-70_S65,NGS-588-71_S66,NGS-588-72_S67,NGS-588-73_S68,NGS-588-74_S69,NGS-588-75_S70,NGS-588-77_S71,NGS-588-78_S72,NGS-588-79_S73,NGS-588-80_S74,NGS-588-81_S75,NGS-588-82_S76,NGS-588-83_S77,NGS-588-84_S78,NGS-588-85_S79,NGS-588-86_S80,NGS-588-87_S81,NGS-588-88_S82,NGS-588-89_S83,NGS-588-90_S84,NGS-588-91_S85,NGS-588-92_S86,NGS-588-94_S87,NGS-588-95_S88,NGS-588-96_S89,NGS-588-97_S90,NGS-588-98_S91"

# Specify haplotype
HAP=hap1

# Specify parameters
CPUS=16
BUSCOLINEAGE=eudicots

####

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
    --rnaseq_sets_dir=${RNASEQ} --rnaseq_sets_ids=${RNAIDS} \
    --prot_seq=${PROTSEQ} --busco_lineage=${BUSCOLINEAGE} \
    --threads=${CPUS} --gff3 --nocleanup \
    --AUGUSTUS_CONFIG_PATH=/augustus_config

# STEP 6: Move results back to origin
mv braker/braker.gff ${PBS_O_WORKDIR}/${BASEPREFIX}_${HAP}.braker.gff3
mv braker/braker.gtf ${PBS_O_WORKDIR}/${BASEPREFIX}_${HAP}.braker.gtf

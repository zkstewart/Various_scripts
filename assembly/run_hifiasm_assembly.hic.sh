#!/bin/bash -l
#PBS -N asm
#PBS -l walltime=24:00:00
#PBS -l mem=90G
#PBS -l ncpus=14

cd $PBS_O_WORKDIR

####

# Specify the Genome_analysis_scripts location
GENSCRIPTDIR=/home/stewarz2/scripts/Genome_analysis_scripts

# Specify compleasm directory
COMPLEASMDIR=/home/stewarz2/various_programs/compleasm_kit
COMPLEASMDB=${COMPLEASMDIR}/db

# Specify computational and behavioural parameters
CPUS=14
PLOIDY=2 # if ploidy > 2, you need to manually add some changes below

# Specify lineage to search
LINEAGE=eudicotyledons

# Specify the Hi-C read prefix
HICPREFIX=BUSK_S2

####

# STEP 0: Get the file name and prefix
INDIR=$(dirname ${PBS_O_WORKDIR})/fastq_corrected/output
FILEIN=$(ls ${INDIR}/*.fasta)
PREFIX=$(basename ${FILEIN} _corrected_reads.fasta)
HICDIR=$(realpath ../hic)

# STEP 1: Run hifiasm
hifiasm --version > hifiasm.version # leave a record of the version used to create this assembly
hifiasm -o ${PREFIX}.asm \
    --n-hap ${PLOIDY} \
    -t ${CPUS} \
    --h1 ${HICDIR}/${HICPREFIX}.trimmed_1P.fq.gz --h2 ${HICDIR}/${HICPREFIX}.trimmed_2P.fq.gz \
    ${FILEIN}

# STEP 2: Generate the FASTA files
## Add hap3... if needed
awk '/^S/{print ">"$2;print $3}' ${PREFIX}.asm.bp.hap1.p_ctg.gfa > ${PREFIX}.asm.bp.hap1.p_ctg.fasta
awk '/^S/{print ">"$2;print $3}' ${PREFIX}.asm.bp.hap2.p_ctg.gfa > ${PREFIX}.asm.bp.hap2.p_ctg.fasta
awk '/^S/{print ">"$2;print $3}' ${PREFIX}.asm.bp.p_ctg.gfa > ${PREFIX}.asm.bp.p_ctg.fasta

# STEP 3: Get the assembly statistics
python ${GENSCRIPTDIR}/genome_stats.py -i ${PREFIX}.asm.bp.hap1.p_ctg.fasta -o ${PREFIX}.asm.bp.hap1.p_ctg.stats
python ${GENSCRIPTDIR}/genome_stats.py -i ${PREFIX}.asm.bp.hap2.p_ctg.fasta -o ${PREFIX}.asm.bp.hap2.p_ctg.stats
python ${GENSCRIPTDIR}/genome_stats.py -i ${PREFIX}.asm.bp.p_ctg.fasta -o ${PREFIX}.asm.bp.p_ctg.stats

# STEP 4: Run completeness scoring of assemblies
## Add hap3... if needed
for assembly in ${PREFIX}.asm.bp.hap1.p_ctg.fasta ${PREFIX}.asm.bp.hap2.p_ctg.fasta ${PREFIX}.asm.bp.p_ctg.fasta; do
    mkdir -p ${assembly}.compleasm
    python ${COMPLEASMDIR}/compleasm.py run -a ${assembly} \
        -L ${COMPLEASMDB} -l ${LINEAGE} -t ${CPUS} \
        -o ${assembly}.compleasm;
done;

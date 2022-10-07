#!/bin/bash -l
#PBS -N clair3
#PBS -l walltime=01:00:00
#PBS -l mem=25G
#PBS -l ncpus=4
#PBS -J 1-15

cd $PBS_O_WORKDIR

#################################

# Specify the name of the conda environment where clair3 has been installed
CLAIR3ENV=clair3

# Specify the location of the genome FASTA
GENOMEDIR=/home/stewarz2/plant_group/ryan_analysis/genome
GENOME=JZ.v1.0.genome.fa ## Make sure this is samtools faidx indexed first!

# Specify the location of the mapped BAM files
MAPDIR=/home/stewarz2/plant_group/ryan_analysis/mapping

# Specify the suffix that identifies mapped BAM files
SUFFIX=.sorted.bam

# Specify the location of the model to be used
MODELDIR=/home/stewarz2/various_programs/clair3/models/ont

# Specify the contig(s) to get variants from
CONTIGS=scaffold514_cov87

# Specify the computational resources to use
CPUS=4

# Specify the location to store the main outputs
MAINOUTDIR=main_results

#################################

# > STEP 1: Activate the conda environment
conda activate ${CLAIR3ENV}

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
PREFIX=$(basename ${INPUTFILE} ${SUFFIX})

# > STEP 6: Run Clair3
run_clair3.sh \
	--bam_fn=${INPUTFILE} \
	--ref_fn=${GENOMEDIR}/${GENOME} \
	--model_path=${MODELDIR} \
	--ctg_name=${CONTIGS} \
	--threads=${CPUS} \
	--platform="ont" \
	--output=${PREFIX} \
	--enable_phasing

# > STEP 7: Gunzip file for human reading
gunzip -c ${PREFIX}/phased_merge_output.vcf.gz > ${PREFIX}/phased_merge_output.vcf

# > STEP 8: Link to file in main results folder
mkdir -p ${MAINOUTDIR}
cd ${MAINOUTDIR}
ln -s ../${PREFIX}/phased_merge_output.vcf ${PREFIX}_clair3_output.vcf

# > STEP 9: bgzip, index, and reheader to enable later merging
bgzip -i ${PREFIX}_clair3_output.vcf
bcftools index ${PREFIX}_clair3_output.vcf.gz

echo ${PREFIX} > tmp_${PREFIX}_reheader.txt
bcftools reheader --samples tmp_${PREFIX}_reheader.txt -o ${PREFIX}_clair3_output_fix.vcf.gz ${PREFIX}_clair3_output.vcf.gz
rm tmp_${PREFIX}_reheader.txt
bcftools index ${PREFIX}_clair3_output_fix.vcf.gz


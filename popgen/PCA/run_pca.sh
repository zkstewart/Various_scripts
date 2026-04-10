#!/bin/bash -l
#PBS -N pca
#PBS -l walltime=04:00:00
#PBS -l mem=60G
#PBS -l ncpus=4

cd $PBS_O_WORKDIR

conda activate plink2

####

# Specify the location of the Various_scripts github repository
VARSCRIPTDIR=/home/stewarz2/scripts/Various_scripts

# Specify the VCF file location
VCF=/home/stewarz2/flgenomics/commercial/variant_calls/merged.filtered.vcf.gz

# Specify the FAM file location
## TSV file with five columns: FAMID SAMPLEID PARENT1ID PARENT2ID SEX PHENO
FAM=/home/stewarz2/flgenomics/commercial/metadata/flgenomics_commercial.fam

# Specify the prefix for output files
PREFIX=flgenomics_commercial

# Specify behavioural parameters
MAXSNPMISS=0.1 # 10% missing data allowed
MAF=0.01 # 1% minor allele frequency minimum
WINDOW=50 # Each genomic window is 50 SNPs in width
STEP=5 # Window steps 5 SNPs to the right then re-checks the window
LDTHRESHOLD=0.5 # Linkage disequilibrium (r-squared) allowed to be up to 0.5 before a linked SNP gets removed

####

# STEP 1: Run PLINK2 without setting --new-id-max-allele-len
PLINKERROR=${PREFIX}.tmpplink.error

if [[ ! -f ${PREFIX}.step1.ok ]]; then
    plink2 --vcf ${VCF} \
           --fam ${FAM} \
           --sort-vars --set-all-var-ids "@:#\$r,\$a" \
           --geno ${MAXSNPMISS} --maf ${MAF} \
           --make-pgen --out tmp_${PREFIX} 2> ${PLINKERROR} || touch ${PREFIX}.step1.ok;
fi;

# STEP 2: Identify the appropriate value to set for --new-id-max-allele-len
REGEX="has length ([0-9]+)\."
ERRORTEXT=$(cat ${PLINKERROR})

if [[ ${ERRORTEXT} =~ ${REGEX} ]]; then
    MAXALEN=$(echo "--new-id-max-allele-len ${BASH_REMATCH[1]} ");
else
    MAXALEN="";
fi;

# STEP 3: Generate PLINK2 representation of the data
if [[ ! -f ${PREFIX}.step3.ok ]]; then
    plink2 --vcf ${VCF} \
           --fam ${FAM} \
           --sort-vars --set-all-var-ids "@:#\$r,\$a" --rm-dup \
           ${MAXALEN} \
           --geno ${MAXSNPMISS} --maf ${MAF} \
           --make-pgen --out ${PREFIX} && touch ${PREFIX}.step3.ok;
fi;

# STEP 4: Perform LD pruning of variants
if [[ ! -f ${PREFIX}.step4.ok ]]; then
    plink2 --pfile ${PREFIX} \
           --indep-pairwise ${WINDOW} ${STEP} ${LDTHRESHOLD} \
           --make-founders \
           --out ${PREFIX} && touch ${PREFIX}.step4.ok;
fi;

# STEP 5: Extract pruned variants
if [[ ! -f ${PREFIX}.step5.ok ]]; then
    plink2 --pfile ${PREFIX} \
           --extract ${PREFIX}.prune.in \
           --make-pgen \
           --out ${PREFIX}_pruned && touch ${PREFIX}.step5.ok;
fi;

# STEP 6: Identify related samples using king-cutoff [0.177 filters first-degree relatives or closer)
if [[ ! -f ${PREFIX}.step6.ok ]]; then
    plink2 --pfile ${PREFIX}_pruned \
           --king-cutoff 0.177 \
           --out ${PREFIX}_kingFilt && touch ${PREFIX}.step6.ok;
fi;

# STEP 7: Obtain allele frequencies for the entire dataset
if [[ ! -f ${PREFIX}.step7.ok ]]; then
    plink2 --pfile ${PREFIX}_pruned \
           --nonfounders --freq counts \
           --out ${PREFIX}_counts && touch ${PREFIX}.step7.ok;
fi;

# STEP 8: Run PCA after removal of related samples
if [[ ! -f ${PREFIX}.step8.ok ]]; then
    plink2 --pfile ${PREFIX}_pruned \
           --keep ${PREFIX}_kingFilt.king.cutoff.in.id \
           --nonfounders --read-freq ${PREFIX}_counts.acount \
           --pca allele-wts 10 \
           --out ${PREFIX}_counts && touch ${PREFIX}.step8.ok;
fi;

# STEP 9: Derive the column indices needed for the next step
SCOREHEADER=($(head -n 1 ${PREFIX}_counts.eigenvec.allele))

for i in "${!SCOREHEADER[@]}"; do
    if [[ "${SCOREHEADER[$i]}" = "ID" ]]; then
        IDCOL=$((i+1));
    fi;
    if [[ "${SCOREHEADER[$i]}" = "A1" ]]; then
        EFFECTCOL=$((i+1));
    fi;
    if [[ "${SCOREHEADER[$i]}" = "PC1" ]]; then
        FIRSTPCCOL=$((i+1));
    fi;
done

LASTPCCOL=$((i+1));

# STEP 10: Project all samples into the PCA ordination space
if [[ ! -f ${PREFIX}.step10.ok ]]; then
    plink2 --pfile ${PREFIX}_pruned \
           --read-freq ${PREFIX}_counts.acount \
           --score ${PREFIX}_counts.eigenvec.allele ${IDCOL} ${EFFECTCOL} header-read no-mean-imputation variance-standardize \
           --score-col-nums ${FIRSTPCCOL}-${LASTPCCOL} --out ${PREFIX}_projected && touch ${PREFIX}.step10.ok;
fi;

# STEP 11: Plot the PCA
conda activate base
if [[ ! -f ${PREFIX}.step11.ok ]]; then
    python ${VARSCRIPTDIR}/GWAS/plot_PCA.py -s ${PREFIX}_projected.sscore \
                                            -e ${PREFIX}_counts.eigenval \
                                            -o ${PREFIX}_pca.html && touch ${PREFIX}.step11.ok;
fi;

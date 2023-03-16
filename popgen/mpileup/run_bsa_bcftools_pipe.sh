#!/bin/bash -l
#PBS -N bsaB_PaSe
#PBS -l walltime=80:00:00
#PBS -l mem=30G
#PBS -l ncpus=2

cd $PBS_O_WORKDIR

# Specify the location of the genome FASTA
GENOMEDIR=/home/stewarz2/flies/genome_based_2022/original/genome
GENOME=btrys06_freeze2.rename.fasta

# Specify the location of the GATK exe
GATKDIR=/home/stewarz2/various_programs/gatk-4.3.0.0

# Specify the prefix for output files
PREFIX=bris_orig_bsa_PaSe

# Specify computational parameters
CPUS=2

####

# > STEP 1.1: Merge bulk1 BAM files together
if [[ ! -f bulk1.unsorted.bam  ]]; then
    samtools merge bulk1.unsorted.bam bulk1/*.bam
fi

# > STEP 1.2: Merge bulk2 BAM files together
if [[ ! -f bulk2.unsorted.bam  ]]; then
    samtools merge bulk2.unsorted.bam bulk2/*.bam
fi

# > STEP 2.1: Sort bulk1 BAM file
if [[ ! -f bulk1.sorted.bam  ]]; then
    samtools sort \
	    -m 5G \
	    -@ ${CPUS} \
	    -o bulk1.sorted.bam \
	    bulk1.unsorted.bam;
    samtools index bulk1.sorted.bam;
fi

# > STEP 2.2: Sort bulk2 BAM file
if [[ ! -f bulk2.sorted.bam  ]]; then
    samtools sort \
            -m 5G \
            -@ ${CPUS} \
            -o bulk2.sorted.bam \
            bulk2.unsorted.bam;
    samtools index bulk2.sorted.bam;
fi

# > STEP 3: Run bcftools mpileup
if [[ ! -f ${PREFIX}.mpileup  ]]; then
    bcftools mpileup \
        -f ${GENOMEDIR}/${GENOME} \
        -q 10 -Q 20 \
	-a AD \
        --threads ${CPUS} \
	--ignore-RG \
        -o ${PREFIX}.mpileup \
        bulk1.sorted.bam bulk2.sorted.bam
fi

# > STEP 4: Run bcftools call & index file
if [[ ! -f ${PREFIX}.vcf.gz  ]]; then
    bcftools call -m -v -Oz \
        -o ${PREFIX}.vcf.gz \
        ${PREFIX}.mpileup;
    tabix ${PREFIX}.vcf.gz;
    tabix -C ${PREFIX}.vcf.gz;
fi

# > STEP 5: Filter vcf
MISSING=0.5 ## this means >50% of individuals need to have the site
MINQ=30 ## minimum SNP quality of 30
MAC=1 ## minor allele count must be >= 1
MAF=0.05 ## minor allele frequency greater than or equal to 5%
if [[ ! -f ${PREFIX}.filtered.vcf  ]]; then
    vcftools --gzvcf ${PREFIX}.vcf.gz --max-missing ${MISSING} --mac ${MAC} --minQ ${MINQ} --remove-filtered-all --recode --recode-INFO-all --maf ${MAF} --out ${PREFIX}.filtered.vcf;
    mv ${PREFIX}.filtered.vcf.recode.vcf ${PREFIX}.filtered.vcf;
fi

# > STEP 6: Remove indels and keep only biallelic sites
if [[ ! -f ${PREFIX}.filtered.noindels.vcf  ]]; then
    bcftools view --max-alleles 2 --exclude-types indels -Ov -o ${PREFIX}.filtered.noindels.vcf ${PREFIX}.filtered.vcf
fi

# > STEP 7: Skip any further filtering, just symbolic link
if [[ ! -f ${PREFIX}.final.vcf  ]]; then
    ln -s ${PREFIX}.filtered.noindels.vcf ${PREFIX}.final.vcf
fi

# > STEP 8: Convert to TSV for BSA
if [[ ! -f ${PREFIX}.final.tsv  ]]; then
    ${GATKDIR}/gatk VariantsToTable -V ${PREFIX}.final.vcf -O ${PREFIX}.final.tsv 
fi

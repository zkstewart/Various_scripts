#!/bin/bash

#######################################################

#### USER SPECIFIED VARIABLES ####

# Specify which conda environment we should be using for each analysis step
CUTADAPT_ENV=cutadapt
DEMULTIPLEX_ENV=base
MERGE_ENV=base
COMPLEMENT_ENV=base
MSA_ENV=base
CHIMERA_ENV=base
CRISPRESSO_ENV=crispresso2_env # make sure this env has had 'pip install biopython' run inside it

# Specify the location where you are running this analysis
WORKING_DIR=/mnt/f/banana_group/crispr_coding/inv_run

# Specify the location of the Various_scripts GitHub directory
## Note: 'git clone https://github.com/zkstewart/Various_scripts.git' somewhere to your computer
## Or, if you haven't done it in a while, go into that directory and do 'git pull' in case I've made updates
VARSCRIPTDIR=/mnt/c/git/Various_scripts

# Specify the location of the MUSCLE and MAFFT executable files
## Note: a version of MUSCLE 3 should be used here
MUSCLE=/usr/bin/muscle
MAFFT=/usr/local/bin/mafft

# Specify the location of the allele reference file(s)
## Note: one or more files can exist in the indicated directory
ALLELES_DIR=/mnt/f/banana_group/novogene_demultiplex/alleles # this is where you'll find the reference alleles FASTA files
ALLELES_SUFFIX=.fasta # this is the file extension of each reference FASTA file

# Specify the location of the adapters file
## Note: this should be a single file containing ALL forward AND reverse adapters with sequence IDs equal to the sequence itself e.g.,
## >ATGACGGAAGGGAGGGATC
## ATGACGGAAGGGAGGGATC
## >...
ADAPTERS_FILE=/mnt/f/banana_group/novogene_demultiplex/new_raw/01.RawData/adapters/both.fasta

# Specify the location of the metadata file
## Note: this should be a tab-separated values (TSV) file consisting of three columns and one header line.
## The header row contents do not matter as it will be ignored. The file itself should look a little like:
## header_1   | header_2  | header_3
## sample_id1 | ATCGTATAT | GATAGATCAGT
## sample_id2 | ATCCTTGAT | CCCAGATCAGT
## ... [note that the nucleotide sequences above are examples only]
METADATA_FILE=/mnt/f/banana_group/novogene_demultiplex/new_raw/01.RawData/demultiplexing_metadata_s1.txt

# Specify the location of the reads files
## Note: this expects a single forward and reverse file; you may need to concatenate lanes together beforehand
FORWARD=/mnt/f/banana_group/novogene_demultiplex/new_raw/01.RawData/S1_INV_TREAT/S1_INV_TREAT_1.fq
REVERSE=/mnt/f/banana_group/novogene_demultiplex/new_raw/01.RawData/S1_INV_TREAT/S1_INV_TREAT_2.fq

# Specify the location you want to contain results files
## Note: this directory should be empty and specific to this analysis run only
OUTDIR=S1_INV_TREAT_results

# Specify how many CPU cores are available for multithreading
CPUS=8

# Specify what Levenshtein ratio to use for reads filtration
## Note: a value of 0 (for no filtering) or a value in the range of 0.75-0.9 might be appropriate here
## you may want to set it at 0.85 and see how many reads get filtered at this value
LEV=0.85


#######################################################
## Note: nothing below this point should require your modification


#### STEP 0 - AUTOMATIC ANALYSIS SETUP ####

# Set up directories for intermediate files
mkdir -p ${WORKING_DIR}/1_cutadapt
mkdir -p ${WORKING_DIR}/2_demultiplex
mkdir -p ${WORKING_DIR}/3_merge
mkdir -p ${WORKING_DIR}/4_complement
mkdir -p ${WORKING_DIR}/5_alleles
mkdir -p ${WORKING_DIR}/6_chimera_filter
mkdir -p ${WORKING_DIR}/7_crispresso

# Get prefixes of alleles files
declare -a PREFIXES
i=0
for f in ${ALLELES_DIR}/*${ALLELES_SUFFIX}; do
    PREFIXES[${i}]=$(basename ${f} ${ALLELES_SUFFIX});
    i=$((i+1));
done


#### STEP 1 - CUTADAPT ####
source activate ${CUTADAPT_ENV}
cutadapt -e 1 --cores ${CPUS} \
    -g file:${ADAPTERS_FILE} \
    -G file:${ADAPTERS_FILE} \
    -o ${WORKING_DIR}/1_cutadapt/{name1}-{name2}.1.fastq -p ${WORKING_DIR}/1_cutadapt/{name1}-{name2}.2.fastq \
    ${FORWARD} ${REVERSE}
echo "# STEP 1 (CUTADAPT) DONE"

if [[ ! ${CUTADAPT_ENV} == "base" ]]; then
    conda deactivate;
fi;


#### STEP 2 - DEMULTIPLEXING ####
source activate ${DEMULTIPLEX_ENV}
python ${VARSCRIPTDIR}/QC/demultiplexing/demultiplex_cutadapt_results.py \
    -i ${WORKING_DIR}/1_cutadapt \
    -m ${METADATA_FILE} \
    -o ${WORKING_DIR}/2_demultiplex
echo "# STEP 2 (DEMULTIPLEX) DONE"

if [[ ! ${DEMULTIPLEX_ENV} == "base" ]]; then
    conda deactivate;
fi;


#### STEP 3 - MERGE ####
source activate ${MERGE_ENV}
for R1FILE in ${WORKING_DIR}/2_demultiplex/*.1.fastq; do
    SAMPLEID=$(basename ${R1FILE} .1.fastq);
    usearch \
        -fastq_mergepairs ${WORKING_DIR}/2_demultiplex/${SAMPLEID}.1.fastq \
        -reverse ${WORKING_DIR}/2_demultiplex/${SAMPLEID}.2.fastq \
        -fastqout ${WORKING_DIR}/3_merge/${SAMPLEID}.fastq;
done;
echo "# STEP 3 (MERGE) DONE"

if [[ ! ${MERGE_ENV} == "base" ]]; then
    conda deactivate;
fi;


#### STEP 4 - COMPLEMENT ####
source activate ${COMPLEMENT_ENV}
for SAMPLEID in ${PREFIXES[@]}; do
    # Skip this allele if it didn't appear in the demultiplexed output
    if [[ ! -f ${WORKING_DIR}/3_merge/${SAMPLEID}.fastq ]]; then
        echo "STEP 4: Skipping ${SAMPLEID} allele as it is not in your merged files";
        continue;
	fi;
    
    # Obtain a representative sequence for this allele
    REPALLELE=$(python ${VARSCRIPTDIR}/fasta_handling_master_code.py -f echoindex \
        -n 1 -i ${ALLELES_DIR}/${SAMPLEID}${ALLELES_SUFFIX});
    
    # Reverse complement sequences as appropriate
    python ${VARSCRIPTDIR}/QC/reads/complement_amplicons.py -a ${REPALLELE} \
        -i ${WORKING_DIR}/3_merge/${SAMPLEID}.fastq \
        -o ${WORKING_DIR}/4_complement/${SAMPLEID}.fastq \
        --minimum_lev ${LEV};
done;
echo "# STEP 4 (COMPLEMENT) DONE"

if [[ ! ${COMPLEMENT_ENV} == "base" ]]; then
    conda deactivate;
fi;


#### STEP 5 - MSA ####
source activate ${MSA_ENV}
python ${VARSCRIPTDIR}/Fasta_related/Alignment/Auto_aligners/mafftAlign.py \
    -i ${ALLELES_DIR} -o ${WORKING_DIR}/5_alleles \
    --mafft ${MAFFT} --alg einsi --maxiterate 10 \
    --threads ${CPUS} --fastaExtensions ${ALLELES_SUFFIX};
echo "# STEP 5 (MSA) DONE"

if [[ ! ${MSA_ENV} == "base" ]]; then
    conda deactivate;
fi;


#### STEP 6 - CHIMERA FILTER ####
source activate ${CHIMERA_ENV}
for FQFILE in ${WORKING_DIR}/4_complement/*.fastq; do
    SAMPLEID=$(basename ${FQFILE} .fastq);
    python ${VARSCRIPTDIR}/QC/reads/predict_chimeric_amplicons.py -fq ${FQFILE} \
        -fa ${WORKING_DIR}/5_alleles/${SAMPLEID}${ALLELES_SUFFIX} \
        -m ${MUSCLE} \
        -o ${WORKING_DIR}/6_chimera_filter \
        --noAmbiguity --turnOffSmoothing;
    mv ${WORKING_DIR}/6_chimera_filter/chimeras.fastq ${WORKING_DIR}/6_chimera_filter/${SAMPLEID}_chimeras.fastq;
    mv ${WORKING_DIR}/6_chimera_filter/non_chimeras.fastq ${WORKING_DIR}/6_chimera_filter/${SAMPLEID}_non_chimeras.fastq;
done;
echo "# STEP 6 (CHIMERA_FILTER) DONE"

if [[ ! ${CHIMERA_ENV} == "base" ]]; then
    conda deactivate;
fi;


#### STEP 7 - CRISPRESSO ####
source activate ${CRISPRESSO_ENV}
for FQFILE in ${WORKING_DIR}/4_complement/*.fastq; do
    SAMPLEID=$(basename ${FQFILE} .fastq);
    
    # Obtain the allele sequences and names as comma-delimited strings
    ALLELES_STR=$(python ${VARSCRIPTDIR}/fasta_handling_master_code.py -f crispressoalleles \
        -i ${ALLELES_DIR}/${SAMPLEID}${ALLELES_SUFFIX});
    NAMES_STR=$(python ${VARSCRIPTDIR}/fasta_handling_master_code.py -f crispressonames \
        -i ${ALLELES_DIR}/${SAMPLEID}${ALLELES_SUFFIX});
    
    # Run CRISPResso2
    CRISPResso -r1 ${WORKING_DIR}/6_chimera_filter/${SAMPLEID}_non_chimeras.fastq \
        -a ${ALLELES_STR} -an ${NAMES_STR} \
        -o ${WORKING_DIR}/7_crispresso \
        -w 2 -wc -3 --amplicon_min_alignment_score 0 --plot_window_size 70 \
        --write_cleaned_report --place_report_in_output_folder \
        --fastq_output --allele_plot_pcts_only_for_assigned_reference;
done;
echo "# STEP 7 (CRISPRESSO) DONE"

if [[ ! ${CRISPRESSO_ENV} == "base" ]]; then
    conda deactivate;
fi;

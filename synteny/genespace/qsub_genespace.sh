#!/bin/bash -l
#PBS -N allGENESPACE
#PBS -l walltime=36:00:00
#PBS -l mem=120G
#PBS -l ncpus=12

cd $PBS_O_WORKDIR

conda activate synteny

####

# Setup: Manual specification of program resources
CPUS=12

####

# STEP 1: Run OrthoFinder [with files in 'peptide']
orthofinder -t ${CPUS} -a ${CPUS} -f peptide -S mmseqs -X -o orthofinder

# STEP 2: Produce results folder for GENESPACE
mkdir -p results
for versusdir in orthofinder/Results_*/Orthologues/Orthologues_*; do cp ${versusdir}/* results; done
for blastfile in orthofinder/Results_*/WorkingDirectory/Blast*.txt; do gzip ${blastfile}; done
for gzfile in orthofinder/Results_*/WorkingDirectory/Blast*.txt.gz; do cp ${gzfile} results; done
cp orthofinder/Results_*/Phylogenetic_Hierarchical_Orthogroups/N0.tsv results
cp orthofinder/Results_*/Orthogroups/Orthogroups.tsv results
cp orthofinder/Results_*/WorkingDirectory/SequenceIDs.txt results
cp orthofinder/Results_*/WorkingDirectory/SpeciesIDs.txt results
cp orthofinder/Results_*/Species_Tree/SpeciesTree_rooted.txt results

# STEP 3: Hope that GENESPACE resumes correctly
Rscript run_genespace.R --no-save

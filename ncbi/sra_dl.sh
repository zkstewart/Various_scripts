#!/bin/bash -l
#PBS -N sra_852
#PBS -l ncpus=1
#PBS -l walltime=48:00:00
#PBS -l mem=30G

cd $PBS_O_WORKDIR

TOOLKITDIR=/home/stewarz2/various_programs/sratoolkit.3.0.0-ubuntu64/bin

$TOOLKITDIR/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRR24350314
$TOOLKITDIR/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRR24350315
$TOOLKITDIR/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRR24350316
$TOOLKITDIR/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRR24350317

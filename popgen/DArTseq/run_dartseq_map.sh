
#!/bin/bash -l
#PBS -N dart_map
#PBS -l walltime=27:00:00
#PBS -l mem=20G
#PBS -l ncpus=1
#PBS -J 1-128

cd $PBS_O_WORKDIR

eval $(cat cmd_dartseq_map.txt | head -n ${PBS_ARRAY_INDEX} | tail -n 1)

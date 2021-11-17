#! python3
# Script to handle DArTseq read mapping via BWA-MEM.
# It follows the process of https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4243306/ 

import os, argparse, subprocess
from pathlib import Path

# Define functions
def validate_args(args):
    # Validate input file location
    if not os.path.isfile(args.fastaFile):
        print('I am unable to locate the genome FASTA file (' + args.fastaFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.metadataCsv):
        print('I am unable to locate the metadata CSV file (' + args.metadataCsv + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isdir(args.demultiplexedDirectory):
        print('I am unable to locate the demultiplexed directory (' + args.demultiplexedDirectory + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()

def parse_metadata_csv(metadataCsv, speciesIdCol, genotypeCols):
    header = None
    speciesIds = []
    genotypeIndices = []
    genotypes = []
    
    with open(metadataCsv, "r") as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n ").split(",")
            # Handle header line
            if header == None:
                header = sl
                # Raise errors if columns don't exist
                if speciesIdCol not in header:
                    raise Exception("Species ID column does not exist")
                for gCol in genotypeCols:
                    if gCol not in header:
                        raise Exception("Genotype column(s) do not exist")
                    else:
                        # Get column index if it exists
                        genotypeIndices.append(header.index(gCol))
                # Get column indices if they exist
                speciesIdIndex = header.index(speciesIdCol)
            # Handle content lines
            else:
                genotype = []
                for gIndex in genotypeIndices:
                    genotype.append(sl[gIndex])
                    
                genotypes.append("_".join(genotype))
                speciesIds.append(sl[speciesIdIndex])
                
    return speciesIds, genotypes

def process_readgroups(speciesIds, genotypes, platform, library, unit):
    readgroups=[]
    for i in range(len(genotypes)):
        readgroups.append("@RG\tID:{0}\tSM:{1}\tPL:{2}\tLB:{3}\tPU:{4}".format(speciesIds[i], genotypes[i], platform, library, unit))
    return readgroups

def process_fqfiles(speciesIds, fqDir):
    fqFiles=[]
    for i in range(len(speciesIds)):
        fqFiles.append(os.path.join(fqDir, "{0}.fq.gz".format(speciesIds[i]))) # This is what dartseq_process will give us
    return fqFiles

def create_cmd_file(speciesIds, fqDir, readgroups, genomeFile, outputFileName="cmd_dartseq_map.txt"):
    fqFiles = process_fqfiles(speciesIds, fqDir)
    
    # Validations
    if os.path.isfile(outputFileName):
        raise FileExistsError("File name <{0}> already exists".format(outputFileName))
    assert len(fqFiles) == len(readgroups)
    
    # Write to file
    with open(outputFileName, "w") as fileOut:
        for i in range(len(fqFiles)):
            sid = speciesIds[i]
            fq = fqFiles[i]
            rg = readgroups[i]
            
            fileOut.write("bwa mem -R '{0}' -p {1} {2} > {3}.sam\n".format(rg, genomeFile, fq, sid))

def create_shell_script(cmdFile, numJobs, outputFileName="run_dartseq_map.sh"):
    # Specify hard-coded script features
    jobname = "dart_map"
    walltime = "12:00:00"
    mem = "10G"
    
    # Setup the script's contents
    formatStr = """
#!/bin/bash -l
#PBS -N {0}
#PBS -l walltime={1}
#PBS -l mem={2}
#PBS -l ncpus=1
#PBS -J 1-{3}

cd $PBS_O_WORKDIR

$(cat {4} | head -n ${{PBS_ARRAY_INDEX}} | tail -n 1)
"""

    # Write to file
    with open(outputFileName, "w") as fileOut:
        fileOut.write(formatStr.format(jobname, walltime, mem, numJobs, cmdFile))

def main():
    # User input
    usage = """%(prog)s enables easier running of BWA-MEM with readgroup
    tags for downstream GATK analysis. It will output a shell script amenable
    to batched job submission via PBS.
    
    It allows for fastq.gz file format as input.
    It assumes the bwa executable can be found in environment variables; otherwise
    the shell script will need to be modified.
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-d", dest="demultiplexedDirectory",
                   required=True,
                   help="Input multiplexed fastq file")
    p.add_argument("-f", dest="fastaFile",
                   required=True,
                   help="Input genome fasta file")
    p.add_argument("-csv", dest="metadataCsv",
                   required=True,
                   help="Input Metadata CSV file")
    p.add_argument("-id", dest="speciesIdCol",
                   required=True,
                   help="Column name where species ID is located")
    p.add_argument("-g", dest="genotypeCols", nargs="+",
                   required=True,
                   help="One or more columns to concatenate as a genotype name")
    p.add_argument("-p", dest="platform",
                   required=True,
                   choices=["illumina", "pacbio"], #incomplete list
                   default="illumina",
                   help="String to use for readgroup platform e.g., 'illumina' by default")
    p.add_argument("-l", dest="library",
                   required=True,
                   default="lib1",
                   help="String to use for library e.g., 'lib1' by default")
    p.add_argument("-u", dest="unit",
                   required=True,
                   default="unit1",
                   help="String to use for unit e.g., 'unit1' by default")
    args = p.parse_args()
    validate_args(args)

    # Parse CSV file columns
    speciesIds, genotypes = parse_metadata_csv(args.metadataCsv, args.speciesIdCol, args.genotypeCols)

    # Process for readgroups and fqfiles
    readgroups = process_readgroups(speciesIds, genotypes, args.platform, args.library, args.unit)
    # fqFiles = process_fqfiles(speciesIds, args.demultiplexedDirectory) # this is done internally to cmd generation
    
    # Create cmd file
    cmdFileName = "cmd_dartseq_map.txt"
    create_cmd_file(speciesIds, args.demultiplexedDirectory, readgroups, args.fastaFile, cmdFileName)
    
    # Create shell script
    create_shell_script(cmdFileName, len(readgroups))
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

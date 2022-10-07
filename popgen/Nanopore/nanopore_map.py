#! python3
# Script to handle Nanopore single-end read mapping via minimap2
# It generally follows the process of https://www.mdpi.com/1422-0067/21/23/9177/htm

import os, argparse, distutils.spawn

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
    if not os.path.isdir(args.fastqDirectory):
        print('I am unable to locate the FASTQ files directory (' + args.fastqDirectory + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate executable locations
    if not os.path.isfile(args.minimap2):
        args.minimap2 = distutils.spawn.find_executable("minimap2")
    if args.minimap2 == None or not os.path.isfile(args.minimap2):
        print('I am unable to locate the minimap2 executable file (' + args.minimap2 + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.nanofilt):
        args.nanofilt = distutils.spawn.find_executable("NanoFilt")
    if args.nanofilt == None or not os.path.isfile(args.nanofilt):
        print('I am unable to locate the NanoFilt executable file (' + args.nanofilt + ')')
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
                genotype = "_".join(genotype)
                
                assert genotype not in genotypes, "Genotype '{0}' isn't unique!".format(genotype)
                genotypes.append(genotype)
                speciesIds.append(sl[speciesIdIndex])
    
    return speciesIds, genotypes

def process_readgroups(speciesIds, genotypes, library, unit):
    readgroups=[]
    for i in range(len(genotypes)):
        readgroups.append("@RG\\tID:{0}\\tSM:{1}\\tPL:{2}\\tLB:{3}\\tPU:{4}".format(speciesIds[i], genotypes[i], "nanopore", library, unit))
    return readgroups

def process_fqfiles(speciesIds, fqDir, suffix=".fq.gz"):
    fqFiles=[]
    for i in range(len(speciesIds)):
        fqFiles.append(os.path.join(fqDir, "{0}{1}".format(speciesIds[i], suffix)))
    return fqFiles

def create_cmd_file(speciesIds, fqDir, readgroups, genomeFile, minimap2Exe, nanofiltExe, minLength, maxLength, outputFileName):
    fqFiles = process_fqfiles(speciesIds, fqDir)
    
    # Validations
    if os.path.isfile(outputFileName):
        raise FileExistsError("File name <{0}> already exists".format(outputFileName))
    assert len(fqFiles) == len(readgroups)
    
    # Generate the cmd file
    with open(outputFileName, "w") as fileOut:
        for i in range(len(fqFiles)):
            # Get relevant details for this sample
            sid = speciesIds[i]
            fq = fqFiles[i]
            rg = readgroups[i]
            
            # Format a command and write to file
            cmd = "gunzip -c {0} | {1} -l {2} --maxlength {3} | ".format(fq, nanofiltExe, minLength, maxLength)
            cmd += "{0} -ax map-ont {1} - > {2}.sam\n".format(minimap2Exe, genomeFile, sid)
            fileOut.write(cmd)

def create_shell_script(cmdFile, numJobs, outputFileName="run_nanopore_map.sh"):
    # Specify hard-coded script features
    jobname = "nanopore_map"
    walltime = "12:00:00"
    mem = "20G"
    
    # Setup the script's contents
    formatStr = """#!/bin/bash -l
#PBS -N {0}
#PBS -l walltime={1}
#PBS -l mem={2}
#PBS -l ncpus=1
#PBS -J 1-{3}

cd $PBS_O_WORKDIR

eval $(cat {4} | head -n ${{PBS_ARRAY_INDEX}} | tail -n 1)
"""

    # Write to file
    with open(outputFileName, "w") as fileOut:
        fileOut.write(formatStr.format(jobname, walltime, mem, numJobs, cmdFile))

def main():
    # User input
    usage = """%(prog)s enables easier running of minimap2 with readgroup
    tags for downstream popgen analysis of Nanopore reads. It will output a shell script
    amenable to batched job submission via PBS.
    
    FASTQ files must have .fq.gz file format as input. This is non-negiotable unfortunately.
    This script will also streamline any read filtration using Nanofilt which can be useful
    when dealing with amplicon-based variant calling.
    
    You must provide the minimap2 executable location as an input. If it's in your system's
    PATH, then providing "-m minimap2" will suffice. The same applies to NanoFilt.
    
    The species ID column should contain unique values for each sample. It should
    also be the prefix to each .fq.gz file! The genotype columns should, when concatenated,
    render a unique value for each sample as well.
    """
    p = argparse.ArgumentParser(description=usage)
    ## Required
    p.add_argument("-d", dest="fastqDirectory",
                   required=True,
                   help="Input directory where fastq files reside")
    p.add_argument("-f", dest="fastaFile",
                   required=True,
                   help="Input genome fasta file")
    p.add_argument("-csv", dest="metadataCsv",
                   required=True,
                   help="Input Metadata CSV file")
    p.add_argument("-m", dest="minimap2",
                   required=True,
                   help="""Input the full path to the minimap2 executable; if it 
                   is in your system path, just specify minimap2 here""")
    p.add_argument("-n", dest="nanofilt",
                   required=True,
                   help="""Input the full path to the NanoFilt executable; if it
                   is in your system path, just specify nanofilt here""")
    p.add_argument("-id", dest="speciesIdCol",
                   required=False,
                   help="Column name where species ID is located")
    p.add_argument("-g", dest="genotypeCols", nargs="+",
                   required=True,
                   help="One or more columns to concatenate as a genotype name")
    p.add_argument("--minLength", dest="minLength",
                   type=int,
                   required=True,
                   help="Specify a minimum read length for NanoFilt")
    p.add_argument("--maxLength", dest="maxLength",
                   type=int,
                   required=True,
                   help="Specify a maximum read length for NanoFilt")
    ## Optional
    p.add_argument("-l", dest="library",
                   required=False,
                   help="String to use for library e.g., 'lib1' by default",
                   default="lib1")
    p.add_argument("-u", dest="unit",
                   required=False,
                   help="String to use for unit e.g., 'unit1' by default",
                   default="unit1")
    args = p.parse_args()
    validate_args(args)
    
    # Parse CSV file columns
    speciesIds, genotypes = parse_metadata_csv(args.metadataCsv, args.speciesIdCol, args.genotypeCols)
    
    # Process for readgroups and fqfiles
    readgroups = process_readgroups(speciesIds, genotypes, args.library, args.unit)
    
    # Create cmd file
    cmdFileName = "cmd_nanopore_map.txt"
    create_cmd_file(speciesIds, args.fastqDirectory, readgroups, args.fastaFile, args.minimap2, args.nanofilt, args.minLength, args.maxLength, cmdFileName)
    
    # Create shell script
    create_shell_script(cmdFileName, len(readgroups))
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

#! python3
# Script to handle Nanopore single-end read mapping via minimap2
# It generally follows the process of https://www.mdpi.com/1422-0067/21/23/9177/htm

import os, argparse, distutils.spawn, gzip

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
        if args.minimap2 == None:
            print(f'I am unable to locate the minimap2 executable file')
        else:
            print(f'I am unable to locate the minimap2 executable file ({args.minimap2})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    
    if not os.path.isfile(args.qcProgram):
        args.qcProgram = distutils.spawn.find_executable(args.qcProgram)
    if args.qcProgram == None or not os.path.isfile(args.qcProgram):
        if args.minimap2 == None:
            print('I am unable to locate the NanoFilt / Chopper executable file')
        else:
            print('I am unable to locate the NanoFilt / Chopper executable file (' + args.qcProgram + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()

def is_gz(file):
    with gzip.open(file, 'r') as fh:
        try:
            fh.read(1)
            return True
        except OSError:
            return False

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
                    raise Exception(f"Species ID column '{speciesIdCol}' does not exist in header '{header}'")
                for gCol in genotypeCols:
                    if gCol not in header:
                        raise Exception(f"Genotype column '{gCol}' do not exist in header '{header}'")
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
        fqFileName = os.path.join(fqDir, "{0}{1}".format(speciesIds[i], suffix))
        
        if not os.path.isfile(fqFileName):
            print(f"I expected to find the file {fqFileName}, but failed to do so!")
            print("Either fix your metadata, or make sure the indicated file exists.")
            quit()
        
        fqFiles.append(fqFileName)
    return fqFiles

def create_cmd_file(speciesIds, fqDir, readgroups, genomeFile, minimap2Exe,
                    qcExe, minLength, maxLength, outputFileName, fqSuffix=".fq.gz"):
    fqFiles = process_fqfiles(speciesIds, fqDir, fqSuffix)
    
    # Validations
    if os.path.isfile(outputFileName):
        raise FileExistsError(f"File name '{outputFileName}' already exists!")
    
    if len(fqFiles) != len(readgroups):
        print(f"The number of samples indicated in your metadata {len(readgroups)} " +
              f"differs to the number of FQ files found {len(fqFiles)}")
        print("Make sure these values corroborate each other, then try again.")
        quit()
    
    # Generate the cmd file
    with open(outputFileName, "w") as fileOut:
        for i in range(len(fqFiles)):
            # Get relevant details for this sample
            sid = speciesIds[i]
            fq = fqFiles[i]
            rg = readgroups[i]
            isGzip = is_gz(fq)
            
            # Format a command and write to file
            if isGzip:
                cmd = f"gunzip -c {fq} | {qcExe} -l {minLength} --maxlength {maxLength} 2> {sid}.qc | "
            else:
                cmd = f"cat {fq} | {qcExe} -l {minLength} --maxlength {maxLength} 2> {sid}.qc | "
            cmd += f"{minimap2Exe} -R '{rg}' -ax map-ont {genomeFile} - > {sid}.sam\n"
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
    
    This script will also streamline any read filtration using Nanofilt or Chopper which can
    be useful when dealing with amplicon-based variant calling.
    
    You must provide the minimap2 executable location as an input. If it's in your system's
    PATH, then providing "-m minimap2" will suffice. The same applies to NanoFilt or
    Chopper.
    
    The species ID column should contain unique values for each sample. It should
    also be the prefix to each FASTQ file! The genotype columns should, when concatenated,
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
    p.add_argument("-qc", dest="qcProgram",
                   required=True,
                   help="""Input the full path to the executable being used for QC;
                   if it is in your system path, just specify the name here""")
    p.add_argument("-id", dest="speciesIdCol",
                   required=False,
                   help="Column name where species ID is located")
    p.add_argument("-g", dest="genotypeCols", nargs="+",
                   required=True,
                   help="One or more columns to concatenate as a genotype name")
    p.add_argument("--minLength", dest="minLength",
                   type=int,
                   required=True,
                   help="Specify a minimum read length")
    p.add_argument("--maxLength", dest="maxLength",
                   type=int,
                   required=True,
                   help="Specify a maximum read length")
    ## Optional
    p.add_argument("-l", dest="library",
                   required=False,
                   help="String to use for library e.g., 'lib1' by default",
                   default="lib1")
    p.add_argument("-u", dest="unit",
                   required=False,
                   help="String to use for unit e.g., 'unit1' by default",
                   default="unit1")
    p.add_argument("--fqSuffix", dest="fqSuffix",
                   required=False,
                   help="""Specify the file suffix which ALL your FASTQ files have
                   e.g., '.fq.gz' by default""",
                   default=".fq.gz")
    args = p.parse_args()
    validate_args(args)
    
    # Parse CSV file columns
    speciesIds, genotypes = parse_metadata_csv(args.metadataCsv, args.speciesIdCol, args.genotypeCols)
    
    # Process for readgroups and fqfiles
    readgroups = process_readgroups(speciesIds, genotypes, args.library, args.unit)
    
    # Create cmd file
    cmdFileName = "cmd_nanopore_map.txt"
    create_cmd_file(speciesIds, args.fastqDirectory, readgroups,
                    args.fastaFile, args.minimap2, args.qcProgram,
                    args.minLength, args.maxLength, cmdFileName, args.fqSuffix)
    
    # Create shell script
    create_shell_script(cmdFileName, len(readgroups))
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

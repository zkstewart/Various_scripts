#!/usr/bin/env python3
# contig_variant_caller.py
# Script to call variants using the bcftools mpileup -> call pipeline on
# a specific contig. That contig will be chunked in a variant calling-safe
# way to enable faster calling.

import argparse
import math
import os
import re
import subprocess
import sys

from Bio import SeqIO

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))) # 3 dirs up is where we find GFF3IO
from Function_packages import ZS_SeqIO, ZS_Utility

# Specify file names being produced by this script
BAM_LIST = "bamlist.txt"
CHUNK_LIST = "{0}_chunks.txt"
CALLING_SCRIPT = "{0}_call.sh"
NORMALISE_SCRIPT = "{0}_normalise.sh"
CONCAT_SCRIPT = "{0}_concatenation.sh"

####

# Define functions
def validate_args(args):
    global CHUNK_LIST
    CHUNK_LIST = CHUNK_LIST.format(args.contigName)
    
    global CALLING_SCRIPT
    CALLING_SCRIPT = CALLING_SCRIPT.format(args.contigName)
    
    global NORMALISE_SCRIPT
    NORMALISE_SCRIPT = NORMALISE_SCRIPT.format(args.contigName)
    
    global CONCAT_SCRIPT
    CONCAT_SCRIPT = CONCAT_SCRIPT.format(args.contigName)
    
    def _not_discoverable_error(program):
        raise FileNotFoundError(f"{program} not discoverable in your system PATH.")
    
    # Validate input file locations
    if not os.path.isfile(args.fastaFile):
        raise FileNotFoundError(f"I am unable to locate the genome FASTA file '{args.fastaFile}'")
    args.fastaFile = os.path.abspath(args.fastaFile)
    
    if not os.path.isdir(args.bamDirectory):
        raise FileNotFoundError(f"I am unable to locate the BAM directory '{args.bamDirectory}'")
    
    # Validate BAM suffix
    foundABAM = False
    for file in os.listdir(args.bamDirectory):
        if file.endswith(args.bamSuffix):
            foundABAM = True
            break
    if not foundABAM:
        raise FileNotFoundError(f"No BAM files with suffix '{args.bamSuffix}' found in directory '{args.bamDirectory}'")
    
    # Validate program discoverability
    if ZS_Utility.wsl_which("samtools") is None:
        _not_discoverable_error("samtools")
    
    if ZS_Utility.wsl_which("bcftools") is None:
        _not_discoverable_error("bcftools")
    
    if ZS_Utility.wsl_which("tabix") is None:
        _not_discoverable_error("tabix")
    
    if ZS_Utility.wsl_which("vt") is None:
        _not_discoverable_error("vt")
    
    # Validate behavioural inputs
    if not re.match(r"^\d+:\d{2}:\d{2}$", args.walltime):
        raise ValueError("--walltime should be in the format of 'HH:MM:SS'")
    if not re.match(r"^\d+G$", args.mem):
        raise ValueError("--mem should be in the format of '*G' where * is a number")
    if args.jobPrefix != "":
        if not args.jobPrefix.endswith("_"):
            args.jobPrefix += "_" # Add underscore if not present
    if args.afterok != None:
        if not re.match(r"^\d+(\[\])?\.\w+$", args.afterok):
            raise ValueError("--afterok should be a PBS job ID")
    
    # Validate numeric values
    if args.threads < 2:
        raise ValueError("-t must be given a value > 1, or there is nothing to parallelise")
    if args.windowSize < 0:
        raise ValueError("--window must be given a value >= 0")
    
    # Validate output file names
    if os.path.exists(CALLING_SCRIPT):
        raise FileExistsError(f"'{CALLING_SCRIPT}' already exists.")
    if os.path.exists(NORMALISE_SCRIPT):
        raise FileExistsError(f"'{NORMALISE_SCRIPT}' already exists.")

def get_contig_length(fastaFile, contigID):
    '''
    Parameters:
        fastaFile -- a string indicating the location of a FASTA file
        contigID -- a string indicating the contig to obtain a length for
    Returns:
        contigLength -- an integer of the contig sequence length
    '''
    with open(fastaFile, "r") as fastaIn:
        records = SeqIO.parse(fastaIn, "fasta")
        for record in records:
            if record.id == contigID:
                return len(record)
    raise ValueError(f"'{contigID}' not found within '{fastaFile}'!")

def get_chunking_points(numberToChunk, chunks, isNumOfChunks=True):
    '''
    This is a general purpose function to take in a number of "things"
    that you want to chunk, and find out how to chunk them evenly.
    
    The resulting list should be interpreted as the 0-based indices where
    a new chunk should form. You should check for this index at the start
    of a loop, and form a new file if your index == the value in this list.
    
    Also, this uses "allocated chunking" such that it will try to keep
    the number of things per chunk approximately equal. Even if you specify
    X number of things per chunk, it might be more optimal to have X-1 in each
    chunk so as to make sure the last chunk doesn't contain a single thing.
    This might not be what you want, but usually, allocated chunking leads to
    more optimal code (e.g., a major use of this function could be for
    parallel processing of the chunks).
    
    Params:
        numberToChunk -- an integer value, possibly derived from a list length as example.
        chunks -- an integer value for the desired number of chunks OR the number of
                  sequences to contain within each chunk, determined by
        isNumOfChunks -- a boolean indicating whether you want the number to be the number
                         of chunks (True), or the number of sequences within each chunk (False)
    '''
    assert isinstance(numberToChunk, int)
    assert isinstance(chunks, int)
    if numberToChunk < chunks:
        raise Exception(f"Chunking only valid if chunkSize <= chunks i.e., {chunks} <= {numberToChunk}")
    
    # Derive how many chunks we want to split the file into
    if isNumOfChunks:
        numChunks = chunks
    else:
        numChunks = math.ceil(numberToChunk / chunks)
    
    rawNum = numberToChunk / numChunks # This line is more relevant in the multithreading code I took this from, but it's okay to just leave it.
    numRoundedUp = round((rawNum % 1) * numChunks, 0) # By taking the decimal place and multiplying it by the num of chunks, we can figure out how many chunks need to be rounded up
    
    # Store positions at which to start a new chunk
    chunkPoints = []
    ongoingCount = 0
    for i in range(numChunks):
        
        # Determine where chunks begin in 0-based indexing
        if i < numRoundedUp: # decide if the number of sequences in this chunk should be rounded up
            point = math.ceil(rawNum) + ongoingCount # Round up the rawNum, and also add our ongoingCount which corresponds to the number of things already put into a chunk
            
            # Prevent chunking beyond the last index where a chunk should start
            if point >= numberToChunk: # Without this check, if we have more chunks than things to chunk, we can end up with "extra" numbers in the list (e.g., [1, 2, 3, 4, 5, 6, 6, 6, 6, 6]).
                break  # This doesn't actually affect program function, but for aesthetic reasons and for clarity of how this function works, I prevent this from occurring.
            
            chunkPoints.append(point)
            ongoingCount += math.ceil(rawNum)
        else:
            point = math.floor(rawNum) + ongoingCount # Round down the rawNum since we've already accounted for any extra uneven numbers
            
            if point >= numberToChunk:
                break
            
            chunkPoints.append(point)
            ongoingCount += math.floor(rawNum)
    
    return chunkPoints

def qsub(scriptFileName):
    '''
    Parameters:
        scriptFileName -- a string indicating the location of a script file to submit
                          to the HPC cluster
    Returns:
        jobID -- a string indicating the job ID of the submitted job
    '''
    qsubProcess = subprocess.Popen(f"qsub {scriptFileName}",
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE,
                                   shell=True)
    jobID, stderr = qsubProcess.communicate()
    jobID, stderr = jobID.decode(), stderr.decode()
    if stderr == "":
        return jobID.strip(" \r\n")
    else:
        raise Exception(f"qsub died with stderr == {stderr}")

def make_calling_script(argsContainer, PREFIX="", WALLTIME="72:00:00", MEM="40G"):
    '''
    Function to format a script for qsub submission to the HPC cluster to perform
    variant calling using the bcftools mpileup->call pipeline in parallel.
    
    Parameters:
        argsContainer -- an object containing the following keys:
            numJobs -- the number of jobs to submit to the HPC cluster;
                       should be equal to the number of contigs in the genome
            workingDir -- the directory to run the script in; should be the CWD
            indelsCns -- a boolean indicating whether to use the --indels-cns option
            ontSup -- a boolean indicating whether to use the -X ont-sup option
            gvcf -- a boolean indicating whether to output GVCF or variants-only VCF
            genomeDir -- the directory containing the genome FASTA file
            genome -- the name of the genome FASTA file
            outputFileName -- the name of the script file to write
        PREFIX -- OPTIONAL; a string to prepend to the job name; default == ""
        WALLTIME -- OPTIONAL; a string indicating how much walltime to give
                    the job; default == "72:00:00"
        MEM -- OPTIONAL; a string indicating how much memory to give the job;
               default == "40G"
    '''
    # Set quality parameters depending on user input
    if argsContainer.ontSup: # use defaults for -X ont-sup with bcftools version 1.20
        qualityLine = ("-B -Q1 --max-BQ 35 --delta-BQ 99 -F0.2 -o15 -e1 -h110 --del-bias 0.4 " + 
                       "--indel-bias 0.7 --poly-mqual --seqq-offset 130 --indel-size 80")
    else: # use long-standing ZKS defaults
        qualityLine = "-q 10 -Q 20"
    
    # Add on optional arguments which control mpileup behaviour
    if argsContainer.indelsCns:
        qualityLine += " --indels-cns"
    if argsContainer.gvcf:
        qualityLine += " --gvcf 0"
    
    # Set call parameters depending on user input
    if argsContainer.gvcf:
        variantLine = "--gvcf 0"
    else:
        variantLine = "-v"
    
    # Generate the script text
    scriptText = \
"""#!/bin/bash -l
#PBS -N 1_{contig}
#PBS -l walltime={WALLTIME}
#PBS -l mem={MEM}
#PBS -l ncpus=1
#PBS -J 1-{numJobs}
{afterokLine}

cd {workingDir}

####

# Specify the location of the genome FASTA
GENOMEDIR={genomeDir}
GENOME={genome}

# Specify the name of the chunks and bamlist file
CONTIG={contig}
CHUNK_LIST={CHUNK_LIST}
BAM_LIST={BAM_LIST}

####

# STEP 1: Get the chunk region to work on
CHUNK=$(cat ${{CHUNK_LIST}} | head -n ${{PBS_ARRAY_INDEX}} | tail -n 1)

# STEP 2: Run bcftools mpileup->call pipeline
bcftools mpileup -Ou \\
    -f ${{GENOMEDIR}}/${{GENOME}} \\
    -r ${{CHUNK}} \\
    --bam-list ${{BAM_LIST}} \\
    {qualityLine} \\
    -a AD,DP | bcftools call -m {variantLine} -Oz -o ${{CONTIG}}_${{PBS_ARRAY_INDEX}}.vcf.gz

# STEP 3: Index the VCF file
tabix ${{CONTIG}}_${{PBS_ARRAY_INDEX}}.vcf.gz
tabix -C ${{CONTIG}}_${{PBS_ARRAY_INDEX}}.vcf.gz

""".format(
    PREFIX=PREFIX,
    WALLTIME=WALLTIME,
    MEM=MEM,
    numJobs=argsContainer.numJobs,
    workingDir=argsContainer.workingDir,
    genomeDir=argsContainer.genomeDir,
    genome=argsContainer.genome,
    contig=argsContainer.contig,
    CONTIG_LIST=CONTIG_LIST, # global variable
    BAM_LIST=BAM_LIST, # global variable
    qualityLine=qualityLine,
    variantLine=variantLine,
    afterokLine = "#PBS -W depend=afterok:{0}".format(":".join(argsContainer.runningJobIDs)) if argsContainer.runningJobIDs != [] else ""
)

    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_normalise_script(argsContainer, PREFIX="", WALLTIME="08:00:00", MEM="15G"):
    '''
    Function to format a script for qsub submission to the HPC cluster to perform
    variant normalisation.
    
    Parameters:
        argsContainer -- an object containing the following keys:
            numJobs -- the number of jobs to submit to the HPC cluster;
                       should be equal to the number of contigs in the genome
            workingDir -- the directory to run the script in; should be the CWD
            genomeDir -- the directory containing the genome FASTA file
            genome -- the name of the genome FASTA file
            outputFileName -- the name of the script file to write
            runningJobIDs -- a list of strings indicating the job IDs of jobs that
                             must complete before this script can run
        PREFIX -- OPTIONAL; a string to prepend to the job name; default == ""
        WALLTIME -- OPTIONAL; a string indicating how much walltime to give
                    the job; default == "08:00:00"
        MEM -- OPTIONAL; a string indicating how much memory to give the job;
               default == "15G"
    '''
    
    scriptText = \
"""#!/bin/bash -l
#PBS -N 2_{contig}
#PBS -l walltime={WALLTIME}
#PBS -l mem={MEM}
#PBS -l ncpus=1
#PBS -J 1-{numJobs}
{afterokLine}

cd {workingDir}

####

# Specify the location of the genome FASTA
GENOMEDIR={genomeDir}
GENOME={genome}

# Specify the contig being worked on
CONTIG={contig}

####

# STEP 1: Get the file prefix to work on
PREFIX=${{CONTIG}}_${{PBS_ARRAY_INDEX}}

# STEP 2: Split multiallelic records to biallelic
bcftools norm -m- -Oz -o ${{PREFIX}}.split.vcf.gz -N ${{PREFIX}}.vcf.gz

# STEP 3: Rejoin biallic sites into multiallelic sites
bcftools norm -m+ -Oz -o ${{PREFIX}}.rejoin.vcf.gz -N ${{PREFIX}}.split.vcf.gz

# STEP 4: Left-align and normalise everything
bcftools norm -f ${{GENOMEDIR}}/${{GENOME}} -Ov -o ${{PREFIX}}.normalised.vcf ${{PREFIX}}.rejoin.vcf.gz

# STEP 5: vt decompose SNPs
vt decompose_blocksub ${{PREFIX}}.normalised.vcf > ${{PREFIX}}.decomposed.vcf

# STEP 6: Make file ready for next steps
bgzip -c ${{PREFIX}}.decomposed.vcf > ${{PREFIX}}.decomposed.vcf.gz
bcftools index ${{PREFIX}}.decomposed.vcf.gz

""".format(
    PREFIX=PREFIX,
    WALLTIME=WALLTIME,
    MEM=MEM,
    numJobs=argsContainer.numJobs,
    workingDir=argsContainer.workingDir,
    genomeDir=argsContainer.genomeDir,
    genome=argsContainer.genome,
    contig=argsContainer.contig,
    CONTIG_LIST=CONTIG_LIST, # global variable
    afterokLine = "#PBS -W depend=afterok:{0}".format(":".join(argsContainer.runningJobIDs)) if argsContainer.runningJobIDs != [] else ""
)

    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_concatenation_script(argsContainer, PREFIX="", WALLTIME="08:00:00", MEM="15G"):
    '''
    Function to format a script for qsub submission to the HPC cluster to perform
    concatenation of multiple VCF files produced in parallel for the same genome.
    
    Parameters:
        argsContainer -- an object containing the following keys:
            numJobs -- the number of jobs to submit to the HPC cluster;
                       should be equal to the number of contigs in the genome
            workingDir -- the directory to run the script in; should be the CWD
            outputFileName -- the name of the script file to write
            runningJobIDs -- a list of strings indicating the job IDs of jobs that
                             must complete before this script can run
        PREFIX -- OPTIONAL; a string to prepend to the job name; default == ""
        WALLTIME -- OPTIONAL; a string indicating how much walltime to give
                    the job; default == "08:00:00"
        MEM -- OPTIONAL; a string indicating how much memory to give the job;
               default == "15G"
    '''
    
    scriptText = \
"""#!/bin/bash -l
#PBS -N 3_{contig}
#PBS -l walltime={WALLTIME}
#PBS -l mem={MEM}
#PBS -l ncpus=1
{afterokLine}

cd {workingDir}

####

# Specify the contig being worked on
CONTIG={contig}

####

# STEP 1: Get our file list
declare -a VCFFILES
i=0
for f in ${{CONTIG}}_*.decomposed.vcf.gz; do
    VCFFILES[${{i}}]=$(echo "${{f}}");
    i=$((i+1));
done

# STEP 2: Get our input files argument
SEPARATOR=" "
VCFTOOLS_ARG="$( printf "${{SEPARATOR}}%s" "${{VCFFILES[@]}}" )"

# STEP 3: Merge individual VCFs
bcftools concat --allow-overlaps --rm-dups -Oz -o ${{CONTIG}}.vcf.gz ${{VCFTOOLS_ARG}}

# STEP 4: Index VCF
tabix ${{CONTIG}}.vcf.gz;
tabix -C ${{CONTIG}}.vcf.gz;

""".format(
    PREFIX=PREFIX,
    WALLTIME=WALLTIME,
    MEM=MEM,
    numJobs=argsContainer.numJobs,
    workingDir=argsContainer.workingDir,
    CONTIG_LIST=CONTIG_LIST, # global variable
    OUTPUT_PREFIX=OUTPUT_PREFIX, # global variable
    afterokLine = "#PBS -W depend=afterok:{0}".format(":".join(argsContainer.runningJobIDs)) if argsContainer.runningJobIDs != [] else ""
)

    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

class Container:
    def __init__(self, paramsDict):
        for key, value in paramsDict.items():
            self.__dict__[key] = value

## Main
def main():
    # User input
    usage = """%(prog)s accepts a reference genome and a directory containing one or more BAM
    files, one for each sample mapped to that reference genome. Unlike variant_calling_pipeline.py,
    this script will focus on only one contig, and accelerate the calling process by chunking
    that contig into variant calling-safe sections.
    
    It will produce scripts suitable for submission to a HPC cluster to perform variant calling
    using the bcftools mpileup->call pipeline for the contig in parallel. This is followed by a
    normalisation process and concatenation of the resulting VCF files. The final VCF file will
    be named <contig>.vcf.gz"""
    
    p = argparse.ArgumentParser(description=usage)
    # Reqs
    p.add_argument("-f", dest="fastaFile",
                   required=True,
                   help="Input genome FASTA file")
    p.add_argument("-c", dest="contigName",
                   required=True,
                   help="Specify the contig to call variants for")
    p.add_argument("-b", dest="bamDirectory",
                   required=True,
                   help="Input directory containing BAM file(s)")
    p.add_argument("-t", dest="threads",
                   required=True,
                   type=int,
                   help="Specify how many threads should be used")
    # Opts (file location)
    p.add_argument("--bamSuffix", dest="bamSuffix",
                   required=False,
                   help="""Indicate the suffix of the BAM files to look for;
                   default == '.sorted.bam'""",
                   default=".sorted.bam")
    # Opts (behavioural)
    p.add_argument("--window", dest="windowSize",
                   required=False,
                   type=int,
                   help="""Optionally specify how much of a buffer window (in bp) should
                   be allowed around the chunk points, in order to ensure that adjacent
                   regions have sufficient overlap to ensure variant calls are made
                   accurately at the junctions; default == 1000""",
                   default=1000)
    p.add_argument("--indels-cns", dest="indelsCns",
                   required=False,
                   action="store_true",
                   help="""Provide this flag to make use of the --indels-cns option as is 
                   currently recommended by bcftools for diploid genomes""",
                   default=False)
    p.add_argument("--gvcf", dest="gvcf",
                   required=False,
                   action="store_true",
                   help="""Provide this flag to call variants with gvcf output
                   which enables merging of samples from separate analyses at the
                   cost of larger files containing invariant sites""",
                   default=False)
    p.add_argument("--ont-sup", dest="ontSup",
                   required=False,
                   action="store_true",
                   help="""Provide this flag to make use of the --ont-sup option for
                   accurate ONT long-read variant calling""",
                   default=False)
    # Opts (PBS)
    p.add_argument("--afterok", dest="afterok",
                   required=False,
                   help="""Optionally, indicate a job ID to wait for before starting this job;
                   default is to start immediately""",
                   default=None)
    p.add_argument("--walltime", dest="walltime",
                   required=False,
                   help="""Optionally, specify how much walltime you want the variant
                   prediction job to have on the HPC cluster; default == '48:00:00'""",
                   default="48:00:00")
    p.add_argument("--mem", dest="mem",
                   required=False,
                   help="""Optionally, specify how much memory you want the variant
                   prediction job to have on the HPC cluster; default == '40G'""",
                   default="40G")
    p.add_argument("--jobPrefix", dest="jobPrefix",
                   required=False,
                   help="""Optionally, specify a prefix to append to each job;
                   default == '' (no prefix); suggest keeping this short""",
                   default="")
    
    args = p.parse_args()
    validate_args(args)
    
    # Index FASTA file (if necessary)
    if not os.path.exists(f"{args.fastaFile}.fai"):
        ZS_SeqIO.StandardProgramRunners.samtools_faidx(args.fastaFile, "samtools") # must be in PATH to pass arg validation
    
    # Get the length of the indicated FASTA contig
    contigLength = get_contig_length(args.fastaFile, args.contigName)
    
    # Get the chunking points for this contig, written as a list for qsub interpretation
    chunkPoints = get_chunking_points(contigLength, args.threads, isNumOfChunks=True)
    if os.path.exists(CHUNK_LIST):
        if not os.path.isfile(CHUNK_LIST):
            raise FileExistsError(f"ERROR: '{CHUNK_LIST}' exists and is not a file.")
        else:
            print(f"chunks file '{CHUNK_LIST}' already exists; skipping.")
    else:
        lastPosition = 1
        with open(CHUNK_LIST, "w") as fileOut:
            for nextPosition in chunkPoints:
                fileOut.write(f"{args.contigName}:{lastPosition}-{nextPosition}\n")
                lastPosition = nextPosition - args.windowSize
            # Write out the last chunk
            fileOut.write(f"{args.contigName}:{lastPosition}-{contigLength}\n")
    
    # Locate BAM files
    bamFiles = [
        os.path.join(args.bamDirectory, file)
        for file in os.listdir(args.bamDirectory)
        if file.endswith(args.bamSuffix)
    ]
    
    # Write out bamlist file
    if os.path.exists(BAM_LIST):
        if not os.path.isfile(BAM_LIST):
            raise FileExistsError(f"ERROR: '{BAM_LIST}' exists and is not a file.")
        else:
            print(f"bamlist file '{BAM_LIST}' already exists; skipping.")
    else:
        with open(BAM_LIST, "w") as fileOut:
            for bamFile in bamFiles:
                fileOut.write(f"{bamFile}\n")
    
    # Write and qsub mpileup->call pipeline script
    runningJobs = [args.afterok] if args.afterok != None else []
    make_calling_script(Container({
            "numJobs": len(chunkPoints),
            "workingDir": os.getcwd(),
            "genomeDir": os.path.dirname(args.fastaFile),
            "genome": os.path.basename(args.fastaFile),
            "contig": args.contigName,
            "indelsCns": args.indelsCns,
            "ontSup": args.ontSup,
            "gvcf": args.gvcf,
            "outputFileName": CALLING_SCRIPT,
            "runningJobIDs": runningJobs
        }),
        PREFIX=args.jobPrefix,
        WALLTIME=args.walltime,
        MEM=args.mem
    )
    callingJobID = qsub(CALLING_SCRIPT)
    runningJobs.append(callingJobID)
    
    # Write and qsub normalisation pipeline script
    make_normalise_script(Container({
            "numJobs": len(contigIDs),
            "workingDir": os.getcwd(),
            "genomeDir": os.path.dirname(args.fastaFile),
            "genome": os.path.basename(args.fastaFile),
            "contig": args.contigName,
            "outputFileName": NORMALISE_SCRIPT,
            "runningJobIDs": runningJobs
        }),
        PREFIX=args.jobPrefix
    )
    normaliseJobID = qsub(NORMALISE_SCRIPT)
    runningJobs.append(normaliseJobID)
    
    # Write and qsub concatenation script
    make_concatenation_script(Container({
            "numJobs": len(contigIDs),
            "workingDir": os.getcwd(),
            "outputFileName": CONCAT_SCRIPT,
            "runningJobIDs": runningJobs
        }),
        PREFIX=args.jobPrefix
    )
    concatenationJobID = qsub(CONCAT_SCRIPT)
    
    # Report job IDs
    print(f"Job IDs: calling={callingJobID}, normalisation={normaliseJobID}, " + 
          f"concatenation={concatenationJobID}")
    if args.afterok != None:
        print(f"Calling job will begin after {args.afterok} completes.")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

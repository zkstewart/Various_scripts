#! python3
# variant_calling_pipeline.py
# Script to call variants using the bcftools mpileup -> call pipeline.
# Has some extra optional features like phasing with WhatsHap and
# generating per-variant or per-sequence reports.

import os, argparse, sys, re, subprocess

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))) # 3 dirs up is where we find GFF3IO
from Function_packages import ZS_Utility, ZS_SeqIO

# Specify file names being produced by this script
CONTIG_LIST = "contiglist.txt"
BAM_LIST = "bamlist.txt"
CALLING_SCRIPT = "variant_call.sh"
NORMALISE_SCRIPT = "variant_normalise.sh"
CONCAT_SCRIPT = "variant_concatenation.sh"
OUTPUT_PREFIX = "merged"

####

# Define functions
def validate_args(args):
    def _not_discoverable_error(program):
        print(f"ERROR: {program} not discoverable in your system PATH.")
        quit()
    
    # Validate input file locations
    if not os.path.isfile(args.fastaFile):
        print('I am unable to locate the genome FASTA file (' + args.fastaFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    args.fastaFile = os.path.abspath(args.fastaFile)
    
    if not os.path.isdir(args.bamDirectory):
        print('I am unable to locate the BAM directory (' + args.bamDirectory + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    
    # Validate BAM suffix
    foundABAM = False
    for file in os.listdir(args.bamDirectory):
        if file.endswith(args.bamSuffix):
            foundABAM = True
            break
    if not foundABAM:
        print(f'No BAM files with suffix "{args.bamSuffix}" found in directory "{args.bamDirectory}"')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    
    # Validate program discoverability
    if ZS_Utility.wsl_which("samtools") is None:
        _not_discoverable_error("samtools")
    
    if ZS_Utility.wsl_which("bcftools") is None:
        _not_discoverable_error("bcftools")
    
    if ZS_Utility.wsl_which("tabix") is None:
        _not_discoverable_error("tabix")
    
    # Validate behavioural inputs
    if not re.match(r"^\d+:\d{2}:\d{2}$", args.walltime):
        print("--walltime should be in the format of 'HH:MM:SS'")
        quit()
    if not re.match(r"^\d+G$", args.mem):
        print("--mem should be in the format of '*G' where * is a number")
        quit()
    if args.jobPrefix != "":
        if not args.jobPrefix.endswith("_"):
            args.jobPrefix += "_" # Add underscore if not present
    
    # Validate output file names
    if os.path.exists(CALLING_SCRIPT):
        raise FileExistsError(f"ERROR: '{CALLING_SCRIPT}' already exists.")
    if os.path.exists(NORMALISE_SCRIPT):
        raise FileExistsError(f"ERROR: '{NORMALISE_SCRIPT}' already exists.")

def get_fasta_ids(fastaFile):
    '''
    Parameters:
        fastaFile -- a string indicating the location of a FASTA file
    Returns:
        contigIDs -- a list of strings representing the sequence IDs in the FASTA file
    '''
    contigIDs = []
    with open(fastaFile, "r") as fastaIn:
        for line in fastaIn:
            if line.startswith(">"):
                contigIDs.append(line[1:].rstrip("\r\n "))
    return contigIDs

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

def make_calling_script(argsDict, PREFIX="", WALLTIME="72:00:00", MEM="40G"):
    '''
    Function to format a script for qsub submission to the HPC cluster to perform
    variant calling using the bcftools mpileup->call pipeline in parallel.
    
    Parameters:
        argsDict -- a dictionary containing the following keys:
            numJobs -- the number of jobs to submit to the HPC cluster;
                       should be equal to the number of contigs in the genome
            workingDir -- the directory to run the script in; should be the CWD
            genomeDir -- the directory containing the genome FASTA file
            genome -- the name of the genome FASTA file
            outputFileName -- the name of the script file to write
        PREFIX -- OPTIONAL; a string to prepend to the job name; default == ""
        WALLTIME -- OPTIONAL; a string indicating how much walltime to give
                    the job; default == "72:00:00"
        MEM -- OPTIONAL; a string indicating how much memory to give the job;
               default == "40G"
    '''
    
    scriptText = \
"""#!/bin/bash -l
#PBS -N {PREFIX}call
#PBS -l walltime={WALLTIME}
#PBS -l mem={MEM}
#PBS -l ncpus=1
#PBS -J 1-{numJobs}

cd {workingDir}

####

# Specify the location of the genome FASTA
GENOMEDIR={genomeDir}
GENOME={genome}

# Specify the name of the contiglist and bamlist file
CONTIG_LIST={CONTIG_LIST}
BAMLIST={BAM_LIST}

####

# STEP 1: Get the contig to work on
CONTIG=$(eval $(cat ${{CONTIG_LIST}} | head -n ${{PBS_ARRAY_INDEX}} | tail -n 1))

# STEP 2: Run bcftools mpileup->call pipeline
bcftools mpileup -Ou \\
    -f ${{GENOMEDIR}}/${{GENOME}} \\
    -r ${{CONTIG}} \\
    --bam-list ${{BAM_LIST}} \\
    -q 10 -Q 20 -a AD | bcftools call -m -v -Oz -o ${{CONTIG}}.vcf.gz

# STEP 3: Index the VCF file
tabix ${{CONTIG}}.vcf.gz
tabix -C ${{CONTIG}}.vcf.gz

""".format(
    PREFIX=PREFIX,
    WALLTIME=WALLTIME,
    MEM=MEM,
    numJobs=argsDict.numJobs,
    workingDir=argsDict.workingDir,
    genomeDir=argsDict.genomeDir,
    genome=argsDict.genome,
    CONTIG_LIST=CONTIG_LIST, # global variable
    BAM_LIST=BAM_LIST # global variable
)

    # Write script to file
    with open(argsDict.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_normalise_script(argsDict, PREFIX="", WALLTIME="08:00:00", MEM="15G"):
    '''
    Function to format a script for qsub submission to the HPC cluster to perform
    variant normalisation.
    
    Parameters:
        argsDict -- a dictionary containing the following keys:
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
#PBS -N {PREFIX}norm
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

# Specify the name of the contiglist file
CONTIG_LIST={CONTIG_LIST}

####

# STEP 1: Get the contig to work on
CONTIG=$(eval $(cat ${{CONTIG_LIST}} | head -n ${{PBS_ARRAY_INDEX}} | tail -n 1))

# STEP 2: Split multiallelic records to biallelic
bcftools norm -m- -Oz -o ${{CONTIG}}.split.vcf.gz -N ${{CONTIG}}.vcf.gz

# STEP 3: Rejoin biallic sites into multiallelic sites
bcftools norm -m+ -Oz -o ${{CONTIG}}.rejoin.vcf.gz -N ${{CONTIG}}.split.vcf.gz

# STEP 4: Left-align and normalise everything
bcftools norm -f ${{GENOMEDIR}}/${{GENOME}} -Ov -o ${{CONTIG}}.normalised.vcf ${{CONTIG}}.rejoin.vcf.gz

# STEP 5: vt decompose SNPs
vt decompose_blocksub ${{CONTIG}}.normalised.vcf > ${{CONTIG}}.decomposed.vcf

# STEP 6: Make file ready for next steps
bgzip -c ${{CONTIG}}.decomposed.vcf > ${{CONTIG}}.decomposed.vcf.gz
bcftools index ${{CONTIG}}.decomposed.vcf.gz

""".format(
    PREFIX=PREFIX,
    WALLTIME=WALLTIME,
    MEM=MEM,
    numJobs=argsDict.numJobs,
    workingDir=argsDict.workingDir,
    genomeDir=argsDict.genomeDir,
    genome=argsDict.genome,
    CONTIG_LIST=CONTIG_LIST, # global variable
    afterokLine = "#PBS -W depend=afterok:{0}".format(":".join(argsDict.runningJobIDs)) if argsDict.runningJobIDs != [] else ""
)

    # Write script to file
    with open(argsDict.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_concatenation_script(argsDict, PREFIX="", WALLTIME="08:00:00", MEM="15G"):
    '''
    Function to format a script for qsub submission to the HPC cluster to perform
    concatenation of multiple VCF files produced in parallel for the same genome.
    
    Parameters:
        argsDict -- a dictionary containing the following keys:
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
#PBS -N {PREFIX}concat
#PBS -l walltime={WALLTIME}
#PBS -l mem={MEM}
#PBS -l ncpus=1
{afterokLine}

cd {workingDir}

####

# Specify output file name prefix
OUTPUT_PREFIX={OUTPUT_PREFIX}

####

# STEP 1: Get our file list
declare -a VCFFILES
i=0
for f in *.decomposed.vcf.gz; do
    VCFFILES[${{i}}]=$(echo "${{f}}");
    i=$((i+1));
done

# STEP 2: Get our input files argument
SEPARATOR=" "
VCFTOOLS_ARG="$( printf "${{SEPARATOR}}%s" "${{VCFFILES[@]}}" )"

# STEP 3: Merge individual VCFs
bcftools merge -Oz -o ${{OUTPUT_PREFIX}}.vcf.gz ${{VCFTOOLS_ARG}}

# STEP 4: Index VCF
tabix ${{OUTPUT_PREFIX}}.vcf.gz;
tabix -C ${{OUTPUT_PREFIX}}.vcf.gz;

""".format(
    PREFIX=PREFIX,
    WALLTIME=WALLTIME,
    MEM=MEM,
    numJobs=argsDict.numJobs,
    workingDir=argsDict.workingDir,
    CONTIG_LIST=CONTIG_LIST, # global variable
    OUTPUT_PREFIX=OUTPUT_PREFIX, # global variable
    afterokLine = "#PBS -W depend=afterok:{0}".format(":".join(argsDict.runningJobIDs)) if argsDict.runningJobIDs != [] else ""
)

    # Write script to file
    with open(argsDict.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

## Main
def main():
    # User input
    usage = """%(prog)s accepts a reference genome and a directory containing one or more BAM
    files, one for each sample mapped to that reference genome. It will produce scripts
    suitable for submission to a HPC cluster to perform variant calling using the bcftools
    mpileup->call pipeline in parallel. This is followed by a normalisation process and
    concatenation of the resulting VCF files. The final VCF file will be named """ + \
    f"{OUTPUT_PREFIX}.vcf.gz."
    
    p = argparse.ArgumentParser(description=usage)
    # Reqs
    p.add_argument("-f", dest="fastaFile",
                   required=True,
                   help="Input genome FASTA file")
    p.add_argument("-b", dest="bamDirectory",
                   required=True,
                   help="Input directory containing BAM file(s)")
    # Opts (behavioural)
    p.add_argument("--bamSuffix", dest="bamSuffix",
                   required=False,
                   help="""Indicate the suffix of the BAM files to look for;
                   default == '.sorted.bam'""",
                   default=".sorted.bam")
    p.add_argument("--walltime", dest="walltime",
                   required=False,
                   help="""Optionally, specify how much walltime you want the variant
                   prediction job to have on the HPC cluster; default == '72:00:00'""",
                   default="72:00:00")
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
        ZS_SeqIO.StandardProgramRunners.samtools_faidx(args.fastaFile, args.samtools)
    
    # Parse FASTA file for sequence IDs
    contigIDs = get_fasta_ids(args.fastaFile)
    
    # Write out FASTA contigs file
    if os.path.exists(CONTIG_LIST):
        if not os.path.isfile(CONTIG_LIST):
            raise FileExistsError(f"ERROR: '{CONTIG_LIST}' exists and is not a file.")
        else:
            print(f"contiglist file '{CONTIG_LIST}' already exists; skipping.")
    else:
        with open(CONTIG_LIST, "w") as fileOut:
            for contigID in contigIDs:
                fileOut.write(f"{contigID}\n")
    
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
    make_calling_script(
        {
            "numJobs": len(contigIDs),
            "workingDir": os.getcwd(),
            "genomeDir": os.path.dirname(args.fastaFile),
            "genome": os.path.basename(args.fastaFile),
            "outputFileName": CALLING_SCRIPT
        },
        PREFIX=args.jobPrefix,
        WALLTIME=args.walltime,
        MEM=args.mem
    )
    callingJobID = qsub(CALLING_SCRIPT)
    
    # Write and qsub normalisation pipeline script
    make_normalise_script(
        {
            "numJobs": len(contigIDs),
            "workingDir": os.getcwd(),
            "genomeDir": os.path.dirname(args.fastaFile),
            "genome": os.path.basename(args.fastaFile),
            "outputFileName": NORMALISE_SCRIPT,
            "runningJobIDs": [callingJobID]
        },
        PREFIX=args.jobPrefix
    )
    normaliseJobID = qsub(NORMALISE_SCRIPT)
    
    # Write and qsub concatenation script
    make_concatenation_script(
        {
            "numJobs": len(contigIDs),
            "workingDir": os.getcwd(),
            "outputFileName": CONCAT_SCRIPT,
            "runningJobIDs": [callingJobID, normaliseJobID]
        },
        PREFIX=args.jobPrefix
    )
    concatenationJobID = qsub(CONCAT_SCRIPT)
    
    # Report job IDs
    print(f"Job IDs: calling={callingJobID}, normalisation={normaliseJobID}, " + 
          f"concatenation={concatenationJobID}")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

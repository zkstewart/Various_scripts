#! python3
# dorado_correct_pipeline.py
# Script to automate the distributed processing of dorado correct
# of nanopore reads. Separates the CPU-bound from GPU-bound tasks.

import os, argparse, sys, re, shutil, subprocess

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
    # Validate input file locations
    if not os.path.isfile(args.fastqFile):
        raise FileNotFoundError(f"Unable to locate the FASTQ reads file '{args.fastqFile}'")
    args.fastqFile = os.path.abspath(args.fastqFile)
    
    # Validate program discoverability
    if args.dorado == None:
        args.dorado = shutil.which("dorado")
        if args.dorado == None:
            raise FileNotFoundError("dorado not discoverable in your system PATH.")
    else:
        if not os.path.isfile(args.dorado):
            raise FileNotFoundError(f"Unable to locate dorado at '{args.dorado}'")
    args.dorado = os.path.abspath(args.dorado)
    
    # Validate dorado lib location
    "dorado executable should be in <base>/bin and lib should be at <base>/lib"
    args.lib = os.path.join(os.path.dirname(os.path.dirname(args.dorado)), "lib")
    if not os.path.isdir(args.lib):
        raise FileNotFoundError(f"Unable to locate dorado lib directory where expected at '{args.lib}'")
    
    # Validate numeric inputs
    if args.numBlocks < 1:
        raise ValueError("--numBlocks should be a positive integer")
    if args.cpuCPU < 1:
        raise ValueError("--cpuCPU should be a positive integer")
    if args.cpuGPU < 1:
        raise ValueError("--cpuGPU should be a positive integer")
    
    # Validate other PBS inputs
    walltimeRegex = re.compile(r"^\d+:\d{2}:\d{2}$")
    memRegex = re.compile(r"^\d+G$")
    
    if not walltimeRegex.match(args.walltimeCPU):
        print("--walltimeCPU should be in the format of 'HH:MM:SS'")
        quit()
    if not walltimeRegex.match(args.walltimeGPU):
        print("--walltimeGPU should be in the format of 'HH:MM:SS'")
        quit()
    if not memRegex.match(args.memCPU):
        print("--memCPU should be in the format of '*G' where * is a number")
        quit()
    if not memRegex.match(args.memGPU):
        print("--memGPU should be in the format of '*G' where * is a number")
        quit()
    
    if args.jobPrefix != "":
        if not args.jobPrefix.endswith("_"):
            args.jobPrefix += "_" # Add underscore if not present
    
    # Validate output directory
    if not os.path.isdir(args.outputDirectory):
        os.makedirs(args.outputDirectory, exist_ok=True)
        print(f"Created output directory '{args.outputDirectory}' as part of argument validation")
    else:
        print(f"Output directory '{args.outputDirectory}' already exists; will overwrite scripts...")
    args.outputDirectory = os.path.abspath(args.outputDirectory)

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

def make_cpu_script(argsContainer, PREFIX="", WALLTIME="72:00:00", CPU=10, MEM="40G"):
    '''
    Function to format a script for qsub submission to the HPC cluster to perform
    dorado correct's CPU-bound step of the pipeline.
    
    Parameters:
        argsContainer -- an object containing the following keys:
            workingDir -- the directory to run the script in; should be the CWD
            outputDir -- the directory to write the step 1 PAF files to
            fastqDir -- the directory containing the FASTQ file to process
            fastqFile -- the filename of the FASTQ file to process
            numBlocks -- the number of blocks to process the FASTQ file in
            dorado -- the full path to the dorado executable
            lib -- the full path to the dorado lib directory
            outputFileName -- the name of the script file to write
        PREFIX -- OPTIONAL; a string to prepend to the job name; default == ""
        WALLTIME -- OPTIONAL; a string indicating how much walltime to give
                    the job; default == "72:00:00"
        MEM -- OPTIONAL; a string indicating how much memory to give the job;
               default == "40G"
    '''
    
    scriptText = \
"""#!/bin/bash -l
#PBS -N {PREFIX}s1
#PBS -l walltime={WALLTIME}
#PBS -l mem={MEM}
#PBS -l ncpus={CPU}
{PBSJ}

cd {workingDir}

####

# Specify the location of the dorado bin and lib folders
DORADOEXE={DORADOEXE}
DORADOLIB={DORADOLIB}

# Specify FASTQ file for correction
FQDIR={FQDIR}
FQFILE={FQFILE}

# Specify directory to write step 1 PAF files to
OUTDIR={OUTDIR}

# Specify computational parameters
CPUS={CPU}

####

# STEP 0: Make job index 0-based
declare -i INDEX
INDEX=${{PBS_ARRAY_INDEX}}-1

# STEP 0.5: Make sure LIB path is set
export LD_LIBRARY_PATH=${{DORADOLIB}}:$LD_LIBRARY_PATH

# STEP 1: Run dorado's CPU-bound step for this block
mkdir -p ${{OUTDIR}}
${{DORADOEXE}} correct --device cpu \\
    --threads ${{CPUS}} \\
    --run-block-id ${{INDEX}} \\
    --to-paf \\
    ${{FQDIR}}/${{FQFILE}} > ${{OUTDIR}}/${{FQFILE}}.block_${{INDEX}}.paf

""".format(
    PREFIX=PREFIX,
    WALLTIME=WALLTIME,
    CPU=CPU,
    MEM=MEM,
    PBSJ=f"#PBS -J 1-{argsContainer.numBlocks}" if argsContainer.numBlocks > 1 else "PBS_ARRAY_INDEX=1",
    workingDir=argsContainer.workingDir,
    FQDIR=argsContainer.fastqDir,
    FQFILE=argsContainer.fastqFile,
    OUTDIR=argsContainer.outputDir,
    DORADOEXE=argsContainer.dorado,
    DORADOLIB=argsContainer.lib
)

    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_gpu_script(argsContainer, PREFIX="", WALLTIME="72:00:00", CPU=10, MEM="40G"):
    '''
    Function to format a script for qsub submission to the HPC cluster to perform
    dorado correct's GPU-bound step of the pipeline.
    
    Parameters:
        argsContainer -- an object containing the following keys:
            workingDir -- the directory to run the script in; should be the CWD
            outputDir -- the directory to write the step 2 FASTA files to
            pafDir -- the directory containing the PAF files to process
            fastqFile -- the filename of the FASTQ file to process
            numBlocks -- the number of blocks to process the FASTQ file in
            dorado -- the full path to the dorado executable
            lib -- the full path to the dorado lib directory
            outputFileName -- the name of the script file to write
            step1JobID -- the job ID of the CPU-bound dorado correct job
        PREFIX -- OPTIONAL; a string to prepend to the job name; default == ""
        WALLTIME -- OPTIONAL; a string indicating how much walltime to give
                    the job; default == "72:00:00"
        MEM -- OPTIONAL; a string indicating how much memory to give the job;
               default == "40G"
    '''
    
    scriptText = \
"""#!/bin/bash -l
#PBS -N {PREFIX}s2
#PBS -l walltime={WALLTIME}
#PBS -l mem={MEM}
#PBS -l ncpus={CPU}
#PBS -l ngpus=1
#PBS -l gputype=A100
#PBS -W depend=afterok:{PREVJOB}
{PBSJ}

cd {workingDir}

####

# Set CUDA_VISIBLE_DEVICES appropriately
## See https://github.com/microsoft/DeepSpeed/issues/5278#issuecomment-2232782045

UUID=$(echo $CUDA_VISIBLE_DEVICES | cut -d',' -f1)
ID=$(nvidia-smi --id=$UUID --query-gpu=index --format=csv,noheader)
export CUDA_VISIBLE_DEVICES=$ID

####

# Specify the location of the dorado bin and lib folders
DORADOEXE={DORADOEXE}
DORADOLIB={DORADOLIB}

# Specify FASTQ file for correction
FQDIR={FQDIR}
FQFILE={FQFILE}

# Specify the location of step 1 PAF files
PAFDIR={PAFDIR}

# Specify directory to write step 2 PAF files to
OUTDIR={OUTDIR}

# Specify computational parameters
CPUS={CPU}

####

# STEP 0: Make job index 0-based
declare -i INDEX
INDEX=${{PBS_ARRAY_INDEX}}-1

# STEP 0.5: Make sure LIB path is set
export LD_LIBRARY_PATH=${{DORADOLIB}}:$LD_LIBRARY_PATH

# STEP 1: Run dorado's GPU-bound step for this block
mkdir -p ${{OUTDIR}}
${{DORADOEXE}} correct ${{PAFDIR}}/${{FQFILE}} \\
    --device cuda:all --threads ${{CPUS}} \\
    --from-paf ${{PAFDIR}}/${{FQFILE}}.block_${{INDEX}}.paf > ${{OUTDIR}}/${{FQFILE}}.block_${{INDEX}}.fasta

""".format(
    PREFIX=PREFIX,
    WALLTIME=WALLTIME,
    CPU=CPU,
    MEM=MEM,
    PREVJOB=argsContainer.step1JobID,
    PBSJ=f"#PBS -J 1-{argsContainer.numBlocks}" if argsContainer.numBlocks > 1 else "PBS_ARRAY_INDEX=1",
    workingDir=argsContainer.workingDir,
    PAFDIR=argsContainer.pafDir,
    FQDIR=argsContainer.fastqDir,
    FQFILE=argsContainer.fastqFile,
    OUTDIR=argsContainer.outputDir,
    DORADOEXE=argsContainer.dorado,
    DORADOLIB=argsContainer.lib
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
    usage = """%(prog)s automates the HPC cluster processing of dorado correct
    of nanopore reads. The script separates the CPU-bound from GPU-bound tasks to
    maximise resource utilisation.
    
    Note: To derive the value to use for --num-blocks, run: dorado correct <-f file>
    --compute-num-blocks
    """
    
    p = argparse.ArgumentParser(description=usage)
    # Reqs
    p.add_argument("-f", dest="fastqFile",
                   required=True,
                   help="Input FASTQ reads file")
    p.add_argument("-n", dest="numBlocks",
                   required=True,
                   type=int,
                   help="Specify the number of blocks the FASTQ file can be processed in")
    p.add_argument("-o", dest="outputDirectory",
                   required=True,
                   help="Output directory to write results to")
    # Opts (program locations)
    p.add_argument("--dorado", dest="dorado",
                   required=False,
                   help="""Optionally, specify the location of the dorado executable;
                   if unspecified, the program will search the system PATH""",
                   default=None)
    # Opts (PBS)
    p.add_argument("--walltimeCPU", dest="walltimeCPU",
                   required=False,
                   help="""Optionally, specify how much walltime you want the CPU-bound
                   job to have on the HPC cluster; default == '48:00:00'""",
                   default="48:00:00")
    p.add_argument("--walltimeGPU", dest="walltimeGPU",
                   required=False,
                   help="""Optionally, specify how much walltime you want the GPU-bound
                   job to have on the HPC cluster; default == '16:00:00'""",
                   default="16:00:00")
    p.add_argument("--cpuCPU", dest="cpuCPU",
                   required=False,
                   type=int,
                   help="""Optionally, specify how many cores you want each CPU-bound
                   subjob to have on the HPC cluster; default == 10""",
                   default=10)
    p.add_argument("--cpuGPU", dest="cpuGPU",
                   required=False,
                   type=int,
                   help="""Optionally, specify how many cores you want the GPU-bound
                   job to have on the HPC cluster; default == 4""",
                   default=4)
    p.add_argument("--memCPU", dest="memCPU",
                   required=False,
                   help="""Optionally, specify how much memory you want the CPU-bound
                   job to have on the HPC cluster; default == '450G'""",
                   default="450G")
    p.add_argument("--memGPU", dest="memGPU",
                   required=False,
                   help="""Optionally, specify how much memory you want the GPU-bound
                   job to have on the HPC cluster; default == '60G'""",
                   default="60G")
    p.add_argument("--jobPrefix", dest="jobPrefix",
                   required=False,
                   help="""Optionally, specify a prefix to append to each job;
                   default == '' (no prefix); suggest keeping this short""",
                   default="")
    
    args = p.parse_args()
    validate_args(args)
    
    # Write and qsub CPU-bound script
    step1Dir = os.path.join(args.outputDirectory, "step1")
    
    cpuScriptName = os.path.join(args.outputDirectory, "run_dorado_correct_step1.sh")
    make_cpu_script(Container({
            "workingDir": args.outputDirectory,
            "outputDir": step1Dir,
            "fastqDir": os.path.dirname(args.fastqFile),
            "fastqFile": os.path.basename(args.fastqFile),
            "numBlocks": args.numBlocks,
            "dorado": args.dorado,
            "lib": args.lib,
            "outputFileName": cpuScriptName
        }),
        PREFIX=args.jobPrefix,
        WALLTIME=args.walltimeCPU,
        CPU=args.cpuCPU,
        MEM=args.memCPU
    )
    cpuJobID = qsub(cpuScriptName)
    
    # Write and qsub GPU-bound script
    step2Dir = os.path.join(args.outputDirectory, "step2")
    
    gpuScriptName = os.path.join(args.outputDirectory, "run_dorado_correct_step2.sh")
    make_gpu_script(Container({
            "workingDir": args.outputDirectory,
            "outputDir": step2Dir,
            "pafDir": step1Dir,
            "fastqDir": os.path.dirname(args.fastqFile),
            "fastqFile": os.path.basename(args.fastqFile),
            "numBlocks": args.numBlocks,
            "dorado": args.dorado,
            "lib": args.lib,
            "outputFileName": gpuScriptName,
            "step1JobID": cpuJobID
        }),
        PREFIX=args.jobPrefix,
        WALLTIME=args.walltimeGPU,
        CPU=args.cpuGPU,
        MEM=args.memGPU
    )
    gpuJobID = qsub(gpuScriptName)
    
    # Report job IDs
    print(f"Job IDs: step1={cpuJobID}, step2={gpuJobID}")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

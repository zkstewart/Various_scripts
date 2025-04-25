#! python3
# Script to handle Nanopore single-end read mapping via minimap2
# It generally follows the process of https://www.mdpi.com/1422-0067/21/23/9177/htm

import os, sys, argparse, gzip, shutil

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from illumina.illumina_map import parse_metadata_csv, get_files_from_prefix, process_readgroups

# Define functions
def validate_args(args):
    # Validate input file location
    if not os.path.isfile(args.fastaFile):
        raise FileNotFoundError(f"I am unable to locate the genome FASTA file ({args.fastaFile})")
    if not os.path.isfile(args.metadataFile):
        raise FileNotFoundError(f"I am unable to locate the metadata file ({args.metadataFile})")
    if not os.path.isdir(args.fastaqDirectory):
        raise FileNotFoundError(f"I am unable to locate the FASTA/Q files directory ({args.fastaqDirectory})")
    
    # Validate program discoverability
    if args.minimap2 == None:
        args.minimap2 = shutil.which("minimap2")
    if not os.path.isfile(args.minimap2):
        raise FileNotFoundError(f"I am unable to locate the minimap2 executable file ({args.minimap2})")
    
    if args.runChopper: # only check if chopper is needed
        if args.chopper == None:
            args.chopper = shutil.which("chopper")
        if not os.path.isfile(args.chopper):
            raise FileNotFoundError(f"I am unable to locate the chopper executable file ({args.chopper})")
    
    # Validate numeric inputs
    if args.cpus < 1:
        raise ValueError(f"--cpus value '{args.cpus}' should be a positive integer (greater than zero)")
    if args.minLength > 0 and not args.runChopper:
        raise ValueError(f"Provide --minLength value only if --runChopper is set")
    if args.maxLength > 0 and not args.runChopper:
        raise ValueError(f"Provide --maxLength value only if --runChopper is set")

def is_gz(file):
    if file.endswith(".gz"):
        return True
    else:
        with gzip.open(file, "r") as fh:
            try:
                fh.read(1)
                return True
            except OSError:
                return False

def create_cmd_file(fastaqFiles, speciesIds, readgroups, genomeFile, minimap2Exe, preset, 
                    chopperExe, runChopper, minLength, maxLength, outputFileName,
                    readsSuffix=".fq.gz", cpus=1):    
    # Validations
    if os.path.isfile(outputFileName):
        raise FileExistsError(f"'{outputFileName}' already exists; delete, move, or rename then try again")
    assert len(fastaqFiles) == len(readgroups)
    
    # Write to file
    with open(outputFileName, "w") as fileOut:
        for i in range(len(fastaqFiles)):
            sid = speciesIds[i]
            fq = fqFiles[i]
            rg = readgroups[i]
            
            # Format the cmd with any pre-processing QC steps
            if runChopper:
                if is_gz(fq):
                    cmd = f"gunzip -c {fq} | {chopperExe}"
                else:
                    cmd = f"cat {fq} | {chopperExe}"
                
                if minLength > 0:
                    cmd += f" -l {minLength}"
                if maxLength > 0:
                    cmd += f" --maxlength {maxLength}"
                cmd += f"2> {sid}.qc | "
            else:
                if is_gz(fq):
                    cmd = f"gunzip -c {fq} | "
                else:
                    cmd = f"cat {fq} | "
            
            # Add and write the minimap2 command
            cmd += f"{minimap2Exe} -t {cpus} -R '{rg}' -ax {preset} {genomeFile} - > {sid}.sam\n"
            fileOut.write(cmd)

def create_shell_script(cmdFile, numJobs, outputFileName="run_nanopore_map.sh", cpus=1):
    # Specify hard-coded script features
    jobname = "nanopore_map"
    walltime = "12:00:00"
    mem = "20G"
    
    # Setup the script's contents
    formatStr = """#!/bin/bash -l
#PBS -N {jobname}
#PBS -l walltime={walltime}
#PBS -l mem={mem}
#PBS -l ncpus={cpus}
#PBS -J 1-{numJobs}

cd $PBS_O_WORKDIR

eval $(cat {cmdFile} | head -n ${{PBS_ARRAY_INDEX}} | tail -n 1)
"""
    # Write to file
    with open(outputFileName, "w") as fileOut:
        fileOut.write(formatStr.format(
            jobname=jobname,
            walltime=walltime,
            mem=mem,
            cpus=cpus,
            numJobs=numJobs,
            cmdFile=cmdFile
        ))

def main():
    # User input
    usage = """%(prog)s enables easier running of minimap2 with readgroup
    tags for downstream popgen analysis of Nanopore reads. It will output a shell script
    amenable to batched job submission via PBS.
    
    The --prefix column must contain a unique value that identifies each sample's read files
    (forward and reverse). This is used solely to locate the files in the -d directory.
    
    The --id column must contain unique values for each sample. This will dictate the output
    SAM files and their embedded read group identifiers, which will dictate what the samples
    are called in your VCF files.
    
    During variant calling, if two samples have the same --sample value, they may be pooled together.
    You probably want this to be unique for each sample and be the same as your --id value unless
    you have a good reason to pool samples.
    
    The species ID column should contain unique values for each sample. It should
    also be the prefix to each FASTA/Q file! The genotype columns should, when concatenated,
    render a unique value for each sample as well.
    """
    p = argparse.ArgumentParser(description=usage)
    ## Required
    p.add_argument("-d", dest="fastaqDirectory",
                   required=True,
                   help="Input directory where FASTA/Q files reside")
    p.add_argument("-f", dest="fastaFile",
                   required=True,
                   help="Input genome fasta file")
    p.add_argument("-csv", dest="metadataFile",
                   required=True,
                   help="Input Metadata CSV or TSV file")
    p.add_argument("-p", dest="preset",
                   required=True,
                   choices=["map-ont", "asm5", "asm10", "asm20",
                            "lr:hq", "splice", "splice:hq", "ava-ont"],
                   help="Specify the preset to use for minimap2")
    ## Optional (locations)
    p.add_argument("--minimap2", dest="minimap2",
                   required=False,
                   help="""If minimap2 is not available in your PATH variable, input
                   the full path to the minimap2 executable here""",
                   default=None)
    p.add_argument("--chopper", dest="chopper",
                   required=False,
                   help="""If chopper is not available in your PATH variable, input
                   the full path to the chopper executable here""",
                   default=None)
    p.add_argument("--readsSuffix", dest="readsSuffix",
                   required=False,
                   help="""Specify the file suffix which ALL your FASTQ files have
                   e.g., '.fq.gz' by default""",
                   default=".fq.gz")
    ## Optional (behavioural)
    p.add_argument("--runChopper", dest="runChopper",
                   required=False,
                   action="store_true",
                   help="""Set this flag if you want to run chopper on the FASTA/Q files""",
                   default=False)
    p.add_argument("--minLength", dest="minLength",
                   type=int,
                   required=False,
                   help="Optionally, specify a minimum read length if --runChopper is set",
                   default=-1)
    p.add_argument("--maxLength", dest="maxLength",
                   type=int,
                   required=False,
                   help="Optionally, specify a maximum read length if --runChopper is set",
                   default=-1)
    ## Optional (metadata)
    p.add_argument("--prefix", dest="prefixCol",
                   required=False,
                   help="Column name where file prefix is located (default == 'prefix')",
                   default="prefix")
    p.add_argument("--id", dest="idCol",
                   required=False,
                   help="Column name where species ID is located (default == 'id')",
                   default="id")
    p.add_argument("--sample", dest="sampleCol",
                   required=False,
                   help="Column name where sample ID is located (default == 'sm')",
                   default="sm")
    p.add_argument("--library", dest="library",
                   required=False,
                   help="String to use for library e.g., 'lib1' by default",
                   default="lib1")
    p.add_argument("--unit", dest="unit",
                   required=False,
                   help="String to use for unit e.g., 'unit1' by default",
                   default="unit1")
    args = p.parse_args()
    args.platform = "nanopore" # hard-coded as script is for nanopore data
    validate_args(args)
    
    # Parse CSV file columns
    prefixes, ids, samples = parse_metadata_csv(args.metadataFile, args.prefixCol, args.idCol, args.sampleCol)
    
    # Validate that files exist and get their locations
    fastqFiles = get_files_from_prefix(args.fastaqDirectory, prefixes)
    
    # Create readgroups
    readgroups = process_readgroups(ids, samples, args.platform, args.library, args.unit)
    
    # Create cmd file
    cmdFileName = "cmd_nanopore_map.txt"
    create_cmd_file(speciesIds, args.fastaqDirectory, readgroups, args.fastaFile,
                    args.minimap2, args.preset, args.chopper, args.runChopper,
                    args.minLength, args.maxLength, cmdFileName,
                    args.readsSuffix, args.cpus)
    
    # Create shell script
    create_shell_script(cmdFileName, len(readgroups), cpus=args.cpus)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

#! python3

import os, argparse

# Define functions for later use
## Argument validation
def validate_args(args):
    # Validate input location
    if not os.path.isdir(args.fastqsDir):
        print('I am unable to locate the input FASTQs directory (' + args.fastqsDir + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Handle file overwrites
    if os.path.isfile(args.outputFileName):
        print(args.outputFileName + ' already exists. Specify a different output file name or delete, move, or rename this file and run the program again.')
        quit()

def main():
    #### USER INPUT SECTION
    usage = """%(prog)s receives a directory containing FASTQs that have been trimmed
    by trimmomatic wherein samples have been split across several lanes. Assuming the
    file naming standard of Illumina (_L00#), this script will produce a shell script
    that can be qsubbed to concatenate the FASTQs into a single forward/reverse file
    per sample.
    """
    
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-d", dest="fastqsDir", required=True,
                help="Input directory containing trimmed FASTQs")
    # Opts
    p.add_argument("-o", dest="outputFileName", required=False,
                help="Optionally specify shell script name (default==run_read_prep.sh)", default="run_read_prep.sh")
    p.add_argument("--suffix", dest="fastqSuffix", required=False,
                help="Optionally specify FASTQ suffix(default==.fq.gz)", default=".fq.gz")
    
    args = p.parse_args()
    validate_args(args)
    
    # Locate files from the directory
    fileNames = []
    for file in os.listdir(args.fastqsDir):
        if file.endswith(args.fastqSuffix):
            fileNames.append(file)
    
    # Get the prefixes and suffixes for reads
    uniquePrefixes=set()
    uniqueSuffixes=set()
    for fileName in fileNames:
        prefix = fileName.split("_L00")[0]
        uniquePrefixes.add(prefix)
        
        suffix = fileName.split("_L00")[1].split(".")[0] 
        uniqueSuffixes.add(suffix)
    
    # Ensure that suffixes match for everything
    for prefix in uniquePrefixes:
        for suffix in uniqueSuffixes:
            laneSuffix = f"_L00{suffix}"
            assert any([fileName.startswith(prefix + laneSuffix) for fileName in fileNames])
    
    # Get ordered suffixes
    orderedSuffixes = list(uniqueSuffixes)
    orderedSuffixes.sort(key = lambda x: int(x))
    
    # Format cat commands to join lane files
    catCmds = []
    for prefix in uniquePrefixes:
        cmd1 = "cat " + " ".join([f"${{TRIMMEDDIR}}/{prefix}_L00{s}.trimmed_1P.fq.gz" for s in orderedSuffixes]) + f" > {prefix}_1.fq.gz"
        cmd2 = "cat " + " ".join([f"${{TRIMMEDDIR}}/{prefix}_L00{s}.trimmed_2P.fq.gz" for s in orderedSuffixes]) + f" > {prefix}_2.fq.gz"
        catCmds.append(cmd1)
        catCmds.append(cmd2)
    
    # Write the script file
    script = f'''#!/bin/bash -l
#PBS -N prep_reads
#PBS -l walltime=120:00:00
#PBS -l mem=5G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

TRIMMEDDIR={os.path.abspath(args.fastqsDir)}
'''
    script += "\n".join(catCmds)

    with open(args.outputFileName, "w") as fileOut:
        fileOut.write(script)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

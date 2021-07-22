#! python3
# setup_hced_shell_script.py
# Simple script to generate a shell script file
# for submission to the QUT HPC job manager which
# will rotate all sequences to have their start
# aligned approximately with a single specified
# reference genome.

import os, argparse

# Define functions
def validate_args(args):
    # Validate that all arguments have been provided
    for key, value in vars(args).items():
        if value == None:
            print(key + ' argument was not specified. Fix this and try again.')
            quit()
    # Validate input file locations
    if not os.path.isfile(args.referenceGenome):
        print('I am unable to locate the reference genome file (' + args.referenceGenome + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    for genomeFile in args.targetGenomes:
        if not os.path.isfile(genomeFile):
            print('I am unable to locate at least one of your target genome file\'s (' + genomeFile + ')')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    if not os.path.isfile(args.hcedExe):
        print('I am unable to locate the hCED executable file (' + args.hcedExe + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print('File already exists at output location (' + args.outputFileName + ')')
        print('Make sure you specify a unique file name and try again.')
        quit()

def format_hced_script(referenceGenome, targetGenomes, hcedExeLocation, outputFileName):
    # Setup default script head
    scriptLines = []
    scriptLines += [
        "#!/bin/bash -l",
        "#PBS -N testHced",
        "#PBS -l walltime=02:00:00",
        "#PBS -l mem=30G",
        "#PBS -l ncpus=1\n",
        "cd $PBS_O_WORKDIR\n",
        "mkdir -p tmp",
        "mkdir -p intermediate",
        "mkdir -p hCED_result\n"
    ]
    # Specify files as an array for shell script
    #targetGenomes = ["/path/to/genome1", "path/to/genome2"]
    scriptLines.append("TARGETS=({0})\n".format(" ".join(targetGenomes)))

    # Set up loop for hCED operations
    #referenceGenome = "/path/to/reference"
    #hcedExeLocation = "/home/n8942188/various_programs/hCED/hCED"
    scriptLines += [
        r"for t in ${TARGETS[@]}; do",
        "    PREFIX=$(basename $t);",
        "    cat {0} $t > tmp/tmp_hced.fasta;".format(referenceGenome),
        "    {0} -i tmp/tmp_hced.fasta -o intermediate/$PREFIX.hced.fasta;".format(hcedExeLocation),
        "done\n"
    ]

    # Concatenate hCED outputs into a single MSA
    scriptLines.append("cat {0} > hCED_result/hCED_result.fasta".format(referenceGenome))
    scriptLines.append("cd hCED_result")
    
    scriptLines += [
        r"for t in ${TARGETS[@]}; do",
        "    PREFIX=$(basename $t);",
        "    tail -n 2 ${PBS_O_WORKDIR}/intermediate/$PREFIX.hced.fasta >> hCED_result.fasta;",
        "done\n"
    ]

    # Write script to file
    with open(outputFileName, "w") as fileOut:
        fileOut.write("\n".join(scriptLines))

def main():
    # User input
    usage = """%(prog)s receives various arguments, including a number of
    target genomes and a single reference genome, and constructs a shell script
    good for submission to the QUT HPC
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-r", dest="referenceGenome",
        help="Input reference genome file around which all other sequences will be \"rotated\"")
    p.add_argument("-t", dest="targetGenomes", nargs="+", default=[],
        help="Input target genome files which will have their start/end adjusted according to the reference")
    p.add_argument("-hced", dest="hcedExe", 
        help="Specify the full path to the hCED executable file")
    p.add_argument("-o", dest="outputFileName",
        help="Output file name for the shell script")
    args = p.parse_args()
    validate_args(args)

    # Generate script file
    format_hced_script(args.referenceGenome, args.targetGenomes, args.hcedExe, args.outputFileName)

if __name__ == "__main__":
    main()

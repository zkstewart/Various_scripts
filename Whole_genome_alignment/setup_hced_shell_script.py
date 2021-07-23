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
    if not os.path.isfile(args.fastaHandlingCode):
        print('I am unable to locate the fasta handling master code file (' + args.fastaHandlingCode + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print('File already exists at output location (' + args.outputFileName + ')')
        print('Make sure you specify a unique file name and try again.')
        quit()

def format_hced_script(referenceGenome, targetGenomes, hcedExeLocation, fastaHandlingCode, suffix, outputFileName):
    # Setup default script head
    scriptLines = []
    scriptLines += [
        "#!/bin/bash -l",
        "#PBS -N testHced",
        "#PBS -l walltime=02:00:00",
        "#PBS -l mem=10G",
        "#PBS -l ncpus=1\n",
        "cd $PBS_O_WORKDIR\n",
        "mkdir -p tmp",
        "mkdir -p intermediate",
        "mkdir -p hCED_result\n"
    ]
    # Specify files as an array for shell script
    #targetGenomes = ["/path/to/genome1", "path/to/genome2"]
    scriptLines.append("TARGETS=({0})\n".format(" ".join(targetGenomes)))

    # Setup BLAST db for finding out if we need to reverse complement our sequence or not
    scriptLines.append("makeblastdb -in {0} -dbtype nucl".format(referenceGenome))

    # Set up loop for hCED operations
    #referenceGenome = "/path/to/reference"
    #hcedExeLocation = "/home/n8942188/various_programs/hCED/hCED"
    scriptLines += [
        r"for t in ${TARGETS[@]}; do",
        "    BASE=$(basename $t);",
        "    PREFIX=${BASE%%.fasta};",
        "    rm tmp/hced_tmp_target.*;",
        ## > Rename the contigs to be unique
        "    python {0} -f rename -s ${{PREFIX}}_{1}_seq{{}} -i $t -o tmp/hced_tmp_target.fasta;".format(fastaHandlingCode, suffix),
        "    rm tmp/hced_tmp_target.list;", # This prevents the if [[]] block below from erroring out when the .list file exists
        ## > Perform a BLAST to see if we need to reverse complement or not
        "    blastn -query tmp/hced_tmp_target.fasta -db {0} -outfmt \"6 qseqid sseqid qstart qend sstart send evalue bitscore qframe sframe\" -out tmp/hced_tmp_target.outfmt6;".format(referenceGenome),
        ## > Check the BLAST result
        "    SFRAME=$(head -n 1 tmp/hced_tmp_target.outfmt6 | awk '{print $10}');", # print $10 grabs the sframe field of the custom outfmt above
        ## > Reverse complement depending on sframe
        "    if [[ \"$SFRAME\" = \"-1\" ]]; then mv tmp/hced_tmp_target.fasta tmp/hced_tmp_target.fasta.tmp; python {0} -f reversecomplement2multi -n 60 -i tmp/hced_tmp_target.fasta.tmp -o tmp/hced_tmp_target.fasta; fi;".format(fastaHandlingCode),
        ## > Resume normal operation for hCED
        "    cat tmp/hced_tmp_target.fasta {0} > tmp/tmp_hced.fasta;".format(referenceGenome),
        "    {0} -i tmp/tmp_hced.fasta -o intermediate/$PREFIX.hced.fasta;".format(hcedExeLocation),
        "done\n"
    ]

    # Concatenate hCED outputs into a single MSA
    scriptLines.append("cat {0} > hCED_result/hCED_result.fasta".format(referenceGenome))
    scriptLines.append("cd hCED_result")
    scriptLines += [
        r"for t in ${TARGETS[@]}; do",
        "    BASE=$(basename $t);",
        "    PREFIX=${BASE%%.fasta}",
        "    tail -n 2 ${PBS_O_WORKDIR}/intermediate/$PREFIX.hced.fasta >> hCED_result.fasta;",
        "done\n"
    ]

    # Fix the MSA to have consistent multispacing
    scriptLines.append("python {0} -f single2multi -n 60 -i hCED_result.fasta -o hCED_result.fix.fasta".format(fastaHandlingCode))

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
    p.add_argument("-f", dest="fastaHandlingCode", 
        help="Specify the full path to the fasta handling master code python file")
    p.add_argument("-s", dest="suffix",
        help="Specify suffix to add to all sequence names (suffix proceeds a _ character; omit this _ here)")
    p.add_argument("-o", dest="outputFileName",
        help="Output file name for the shell script")
    args = p.parse_args()
    validate_args(args)

    # Generate script file
    format_hced_script(args.referenceGenome, args.targetGenomes, args.hcedExe, args.fastaHandlingCode, args.suffix, args.outputFileName)

if __name__ == "__main__":
    main()

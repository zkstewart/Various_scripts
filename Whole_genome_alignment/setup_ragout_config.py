#! python3
# setup_ragout_config.py
# Simple script to generate the configuration recipe
# file needed by ragout to perform whole-genome
# scaffolding with multiple incomplete references

import os, argparse

# Define functions
def validate_args(args):
    # Validate that all arguments have been provided
    for key, value in vars(args).items():
        if value == None:
            print(key + ' argument was not specified. Fix this and try again.')
            quit()
    # Validate input file locations
    #fullPathReferences = []
    for refFile in args.references:
        if not os.path.isfile(refFile):
            print('I am unable to locate the reference genome file (' + refFile + ')')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
        #else:
        #    fullPathReferences.append(os.path.abspath(refFile))
    if not os.path.isfile(args.target):
        print('I am unable to locate the target genome file (' + args.target + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    #else:
    #    fullPathTarget = os.path.abspath(args.target)
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print('File already exists at output location (' + args.outputFileName + ')')
        print('Make sure you specify a unique file name and try again.')
        quit()
    #return fullPathReferences, fullPathTarget

def format_ragout_config(referenceGenomeFiles, targetGenomeFile, outputConfigName, isDraft):
    # Format reference values
    refNames = []
    refPaths = []
    for refFile in referenceGenomeFiles:
        refName = os.path.basename(refFile.rsplit(".", maxsplit=1)[0]).lower()
        path = "{0}.fasta = {1}".format(refName, refFile)
        refNames.append(refName)
        refPaths.append(path)

    # Format target values
    targetName = os.path.basename(targetGenomeFile.rsplit(".", maxsplit=1)[0]).lower()
    targetPath = "{0}.fasta = {1}".format(targetName, targetGenomeFile)

    # Write to file
    with open(outputConfigName, "w") as fileOut:
        fileOut.write("#reference and target genome names (required)\n")
        fileOut.write(".references = {0}\n".format(",".join(refNames)))
        fileOut.write(".target = {0}\n\n".format(targetName)) # Leave a blank line after

        fileOut.write("#paths to genome fasta files (required for Sibelia)\n")
        for path in refPaths:
            fileOut.write("{0}\n".format(path))
        fileOut.write("{0}\n".format(targetPath))

        if isDraft:
            fileOut.write("\n*.fasta = true\n") # Specify that all genomes are in draft form, with separation between it and prior details

def main():
    # User input
    usage = """%(prog)s receives various arguments, including a number of
    reference genomes and a single target genome, and constructs a config file
    for use by Ragout to scaffold the target genome. It does not run Ragout itself.
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-r", dest="references", nargs="+",
        help="Input reference genome files")
    p.add_argument("-t", dest="target",
        help="Input target genome file")
    p.add_argument("-o", dest="outputFileName",
        help="Output file name for the config file")
    p.add_argument("-d", dest="drafts", action="store_true", default=False,
        help="Optionally specify that the references are drafts (unspecified == not draft)")
    args = p.parse_args()
    validate_args(args)

    # Generate config file
    format_ragout_config(args.references, args.target, args.outputFileName, args.drafts)

if __name__ == "__main__":
    main()

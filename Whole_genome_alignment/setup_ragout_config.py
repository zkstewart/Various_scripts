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
    for refFile in args.draftRefs:
        if not os.path.isfile(refFile):
            print('I am unable to locate the reference genome file (' + refFile + ')')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    for refFile in args.completeRefs:
        if not os.path.isfile(refFile):
            print('I am unable to locate the reference genome file (' + refFile + ')')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    if not os.path.isfile(args.target):
        print('I am unable to locate the target genome file (' + args.target + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print('File already exists at output location (' + args.outputFileName + ')')
        print('Make sure you specify a unique file name and try again.')
        quit()

def format_ragout_config(draftReferenceFiles, completeReferenceFiles, targetGenomeFile, outputConfigName):
    # Format reference values
    draftRefNames = []
    draftRefPaths = []
    for refFile in draftReferenceFiles:
        refName = os.path.basename(refFile.rsplit(".", maxsplit=1)[0]).lower()
        path = "{0}.fasta = {1}".format(refName, refFile)
        draftRefNames.append(refName)
        draftRefPaths.append(path)
    
    completeRefNames = []
    completeRefPaths = []
    for refFile in completeReferenceFiles:
        refName = os.path.basename(refFile.rsplit(".", maxsplit=1)[0]).lower()
        path = "{0}.fasta = {1}".format(refName, refFile)
        completeRefNames.append(refName)
        completeRefPaths.append(path)

    # Format target values
    targetName = os.path.basename(targetGenomeFile.rsplit(".", maxsplit=1)[0]).lower()
    targetPath = "{0}.fasta = {1}".format(targetName, targetGenomeFile)

    # Write to file
    with open(outputConfigName, "w") as fileOut:
        # Specify the names for genomes
        fileOut.write("#reference and target genome names (required)\n")
        fileOut.write(".references = {0}\n".format(",".join(draftRefNames + completeRefNames)))
        fileOut.write(".target = {0}\n\n".format(targetName)) # Leave a blank line after

        # Format the draft names, specifying that they are draft quality
        fileOut.write("#paths to genome fasta files (required for Sibelia)\n")
        for i in range(len(draftRefPaths)):
            path = draftRefPaths[i]
            name = draftRefNames[i]
            fileOut.write("{0}\n".format(path))
            fileOut.write("{0}.draft = true\n".format(name))
        
        # Format the complete names, leaving quality default to indicate they are complete
        for i in range(len(completeRefPaths)):
            path = completeRefPaths[i]
            fileOut.write("{0}\n".format(path))
        
        # Add the target name
        fileOut.write("{0}\n".format(targetPath))

def main():
    # User input
    usage = """%(prog)s receives various arguments, including a number of
    reference genomes and a single target genome, and constructs a config file
    for use by Ragout to scaffold the target genome. It does not run Ragout itself.
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-rd", dest="draftRefs", nargs="+", default=[],
        help="Input reference genome files which are DRAFT quality")
    p.add_argument("-rc", dest="completeRefs", nargs="+", default=[],
        help="Input reference genome files which are COMPLETE/CHROMOSOME quality")
    p.add_argument("-t", dest="target",
        help="Input target genome file")
    p.add_argument("-o", dest="outputFileName",
        help="Output file name for the config file")
    args = p.parse_args()
    validate_args(args)

    # Generate config file
    format_ragout_config(args.draftRefs, args.completeRefs, args.target, args.outputFileName)

if __name__ == "__main__":
    main()

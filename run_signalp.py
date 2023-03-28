#! python3
# run_signalp.py
# Script to take a FASTA file input and run SignalP
# on it for signal peptide predictions

import os, argparse, sys

sys.path.append(os.path.dirname(__file__))
from Function_packages import ZS_SignalPIO

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.fastaFile):
        print(f'I am unable to locate the input FASTA file ({args.fastaFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate program location
    if not os.path.isfile(args.signalpExe):
        print(f'I am unable to locate the SignalP exe file ({args.signalpExe})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print(f'File already exists at output location ({args.outputFileName})')
        print('Make sure you specify a unique file name and try again.')
        quit()

## Main
def main():
    # User input
    usage = """%(prog)s reads in a FASTA file and generates an output TSV indicating
    the signal peptide predictions made. It's assumed your FASTA contains protein
    sequences; nucleotide input may lead to unpredictable outcomes.
    
    This script relies on the ZS_SignalPIO package, which at this time attempts to
    handle SignalP versions 4, 5, and 6. However, version 4 is non-functional, and
    version 6 is only functional on Linux (not WSL!).
    """
    p = argparse.ArgumentParser(description=usage)
    # Required
    p.add_argument("-i", dest="fastaFile",
                   required=True,
                   help="Input FASTA file containing protein sequences.")
    p.add_argument("-s", dest="signalpExe",
                   required=True,
                   help="Input full path to the signalP executable")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file name for signalP predictions")
    # Optional
    p.add_argument("--org", dest="organism",
                   required=False,
                   choices=["euk", "eukarya", "gram+", "gram-", "arch", "other"],
                   help="""Optionally, specify which organism type to search for
                   (default == "euk"; make sure to use the right term for your
                   version of signalP)""")
    
    args = p.parse_args()
    validate_args(args)
    
    # Create and configure signalP handler object
    sigp = ZS_SignalPIO.SignalP(args.fastaFile, args.signalpExe)
    sigp.organism = args.organism
    
    # Run prediction
    sigpResultsDict = sigp.signalp()
    
    # Create output file
    with open(args.outputFileName, "w") as fileOut:
        # Write header
        fileOut.write("#sequence_id\tstart\tend\n")
        # Content lines
        for seqid, coords in sigpResultsDict.items():
            start, end = coords
            fileOut.write(f"{seqid}\t{start}\t{end}\n")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

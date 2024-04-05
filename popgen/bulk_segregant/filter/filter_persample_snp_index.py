#! python3
# filter_qtlseq_snp_index.py
# Script to enable easy filtering of the snp_index.p##.tsv
# file produced by QTL-seq

import os, argparse
from collections import OrderedDict

# Define functions
def validate_args(args):
    # Validate input file location
    if not os.path.isfile(args.snpIndexFile):
        print('I am unable to locate the SNP index p## file (' + args.snpIndexFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate numeric inputs
    if args.differenceRatio != None:
        if 0 > args.differenceRatio:
            print("differenceRatio argument only makes sense to be 0 or greater")
            quit()
        if 1 < args.differenceRatio:
            print("differenceRatio argument only makes sense to be less than 1")
            quit()
    # Handle file overwrites
    if os.path.isfile(args.outputFileName):
        print(args.outputFileName + ' already exists. Delete/move/rename this file and run the program again.')
        quit()

def main():
    # User input
    usage = """%(prog)s filters a QTLseq persample analysis.
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="snpIndexFile",
                   required=True,
                   help="Input SNP index p## file produced by QTL-seq")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Specify the output file name")
    # Opts
    p.add_argument("--differenceRatio", dest="differenceRatio",
                   required=True,
                   type=float,
                   help="Optionally, filter to retain the indicated variant type")
    args = p.parse_args()
    validate_args(args)
    
    # Write new output file
    firstLine = True
    with open(args.snpIndexFile, "r") as fileIn, open(args.outputFileName, "w") as fileOut:
        for line in fileIn:
            if firstLine == True:
                fileOut.write(line)
                firstLine = False
            else:
                sl = line.rstrip("\r\n ").split("\t")
                if float(sl[-1]) >= args.differenceRatio:
                    fileOut.write(line)
        
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

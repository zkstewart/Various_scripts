#! python3
# convert_orthofinder_to_cafe.py
# Script to take the OrthoFinder output file 
# (Orthogroups.GeneCount.tsv) and make it amenable
# to CAFE operations.

import os, argparse
from Bio import SeqIO

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.tsvFileName):
        print('I am unable to locate the Orthogroups.GeneCount.tsv file (' + args.tsvFileName + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print('File already exists at output location (' + args.outputFileName + ')')
        print('Make sure you specify a unique file name and try again.')
        quit()

def convert_orthofinder_genecount_to_cafe(genecountFileName):
    '''
    Returns:
        cafeLinesList -- a list of lines to be written to an output file
    '''
    firstLine = True
    cafeLinesList = []
    with open(genecountFileName, "r") as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n").split("\t")
            if firstLine == True:
                cafeLinesList.append("Desc\tFamily ID\t{0}".format("\t".join(sl[1:-1]))) # drop the Total column at end
                firstLine = False
            else:
                familyID, counts = sl[0], sl[1:-1] # drop the Total column at end
                desc = "(null)"
                cafeLinesList.append("{0}\t{1}\t{2}".format(desc, familyID, "\t".join(counts)))
    return cafeLinesList

def main():
    # User input
    usage = """%(prog)s reads in the Orthogroups.GeneCount.tsv file producted by OrthoFinder
    and reformats it for use by CAFE5.
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-t", dest="tsvFileName", required=True,
        help="Input Orthogroups.GeneCount.tsv file location")
    p.add_argument("-o", dest="outputFileName", required=True,
        help="Output file name for the modified tsv file")
    args = p.parse_args()
    validate_args(args)

    # Parse OrthoFinder file
    cafeLinesList = convert_orthofinder_genecount_to_cafe(args.tsvFileName)

    # Produce output file
    with open(args.outputFileName, "w") as fileOut:
        fileOut.write("\n".join(cafeLinesList))
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

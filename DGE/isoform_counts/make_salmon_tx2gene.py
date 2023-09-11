#! python3
# make_salmon_tx2gene.py
# Reads in a GFF3 and produces a tx2gene-amenable file

import os, argparse, sys

# Load functions from other scripts
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from Function_packages import ZS_GFF3IO

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.gff3File):
        print('I am unable to locate the GFF3 file (' + args.gff3File + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print('File already exists at output location (' + args.outputFileName + ')')
        print('Make sure you specify a unique file name and try again.')
        quit()

## Main
def main():
    # User input
    usage = """%(prog)s accepts a GFF3 file, and produces a tx2gene
    TSV file where the first column is the transcript ID, and the second
    is the gene ID. This file can be used, after formatting list in R akin
    to list(counts = ..., abundance = ..., length = ..., countsFromAbundance="no")
    by calling summarizeToGene(list(...), tx2gene).
    """
    p = argparse.ArgumentParser(description=usage)
    ## Required
    p.add_argument("-i", dest="gff3File",
                   required=True,
                   help="Input GFF3 file")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file name")
    
    args = p.parse_args()
    validate_args(args)
    
    # Load GFF3
    gff3Obj = ZS_GFF3IO.GFF3(args.gff3File, False) # non-strict parsing
    
    # Write output
    with open(args.outputFileName, "w") as fileOut:
        fileOut.write("transcript_id\tgene_id\n")
        for geneFeature in gff3Obj.types["gene"]:
            if hasattr(geneFeature, "mRNA"):
                geneID = geneFeature.ID
                for mrnaFeature in geneFeature.mRNA:
                    mrnaID = mrnaFeature.ID
                    fileOut.write(f"{mrnaID}\t{geneID}\n")
    
    # Alert user to program success
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

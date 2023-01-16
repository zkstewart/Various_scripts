#! python3
# clustered_annot_table.py
# A simple script to receive an annotation table produced by the
# ZKS series of scripts and reduce it to only the representative
# sequence of each cluster.

import os, argparse

def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.representativesTSV):
        print(f"I was unable to locate the Corset representatives file at '{args.representativesTSV}'")
        print("Make sure you specified the right file name and location, then try again.")
        quit()
    if not os.path.isfile(args.annotationTableFile):
        print(f"I was unable to locate the annotation table file at '{args.annotationTableFile}'")
        print("Make sure you specified the right file name and location, then try again.")
        quit()
    # Ensure that the output location is sensible
    if os.path.isfile(args.outputFile):
        print(f"The specified output file '{args.outputFile}' already exists. This script will not overwrite existing files.")
        print("Make sure to move this file or specify a different file name, then try again.")
        quit()

def main():
    #### USER INPUT SECTION
    usage = """%(prog)s will read in TSV file containing two columns linking
    1) the Corset cluster ID to 2) the representative transcript ID. The output
    annotation table will retain only rows which correspond to representatives,
    with the ID being replaced with the Corset cluster ID.
    """
    # Required
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-r", dest="representativesTSV",
                   required=True,
                   help="""Specify the location of the representatives .TSV file
                   produced by corset_cluster_to_fasta.py""")
    p.add_argument("-t", dest="annotationTableFile",
                   required=True,
                   help="Specify the annotation table file location")
    p.add_argument("-o", dest="outputFile",
                   required=True,
                   help="Specify the name of the output clustered annotation table file")
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse representatives
    representatives = {}
    with open(args.representativesTSV, "r") as fileIn:
        for line in fileIn:
            l = line.rstrip("\r\n ")
            if l == "":
                continue
            else:
                clustID, geneID = l.split("\t")
                representatives[clustID] = geneID
                representatives[geneID] = clustID
    
    # Parse table to get output
    with open(args.annotationTableFile, "r") as fileIn, open(args.outputFile, "w") as fileOut:
        for line in fileIn:
            if line.startswith("#"):
                fileOut.write(line)
                continue
            
            sl = line.rstrip("\r\n").split("\t")
            id = sl[0]
            
            # Figure out if we want to write this line to file
            if id in representatives:
                fileOut.write(line)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

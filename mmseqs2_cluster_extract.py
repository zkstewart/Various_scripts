#! python3
# mmseqs_cluster_extract.py
# Script to take in a cluster file from MMseqs2 alongside its FASTA
# and pull out just the representative sequences

import argparse, os
from Bio import SeqIO

# Various functions for program operations
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.clusterFile) and not os.path.islink(args.clusterFile):
        print(f'I am unable to locate the MMseqs2 cluster file ({args.clusterFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.fastaFile) and not os.path.islink(args.fastaFile):
        print(f'I am unable to locate the FASTA file ({args.fastaFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate output file location
    if os.path.exists(args.outputFileName):
        print('The specified output file name already exists (' + args.outputFileName + ')')
        print('Specify a different name or move/delete/rename the existing file and try again.')
        quit()

def parse_mms2_representatives(clusterFile):
    repIDs = []
    with open(clusterFile, "r") as fileIn:
        for line in fileIn:
            rep, member = line.rstrip("\r\n ").split("\t")
            repIDs.append(rep)
    return set(repIDs)

# Main call
def main():
    #### USER INPUT SECTION
    usage = """Script to receive an MMseqs2 cluster file in two-column TSV format
    where left column == representative ID, and right column == member sequence ID.
    Alongside the originally clustered FASTA, this script will produce an output
    FASTA subset to just representative sequences."""
    p = argparse.ArgumentParser(description=usage)
    # Required
    p.add_argument("-c", dest="clusterFile",
                   required=True,
                   help="Input cluster file")
    p.add_argument("-f", dest="fastaFile",
                   required=True,
                   help="Input FASTA file")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output FASTA file name for representative sequences")
    args = p.parse_args()
    validate_args(args)
    
    # Parse cluster file
    repIDs = parse_mms2_representatives(args.clusterFile)
    
    # Load in FASTA file
    records = SeqIO.to_dict(SeqIO.parse(args.fastaFile, "fasta"))
    
    # Produce output
    with open(args.outputFileName, "w") as fileOut:
        for repID in repIDs:
            fileOut.write(records[repID].format("fasta"))
    
    # Done!
    print('Program completed successfully!')

if __name__ == '__main__':
    main()

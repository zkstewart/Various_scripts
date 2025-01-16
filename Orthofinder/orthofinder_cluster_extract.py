#! python3
# orthofinder_cluster_extract.py
# Script to take in a Orthogroups.tsv file from OrthoFinder alongside its FASTA
# and pull out just the representative sequences

import argparse, os
from Bio import SeqIO

# Various functions for program operations
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.tsvFile) and not os.path.islink(args.tsvFile):
        print(f'I am unable to locate the OrthoFinder .tsv file ({args.tsvFile})')
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

def parse_orthofinder_representatives(tsvFile, recordsList, outputFileName, seqPrefix=""):
    repIDs = []
    with open(tsvFile, "r") as fileIn, open(outputFileName, "w") as fileOut:
        firstLine = True
        for line in fileIn:
            if firstLine:
                firstLine = False
                continue
            
            sl = line.rstrip("\r\n ").split("\t")
            
            seqIDs = [ x for column in sl[1:] for x in column.split(", ") if x != "" ]
            
            # Get longest sequence
            longestSeq = ""
            for i, column in enumerate(sl[1:]):
                if column != "":
                    for seqID in column.split(", "):
                        try:
                            seq = records[seqID]
                        except KeyError:
                            seq = records[seqPrefix + seqID]
                        if len(seq) > len(longestSeq):
                            longestSeq = seq
            
            # Write to file
            fileOut.write(longestSeq.format("fasta"))
    return set(repIDs)

# Main call
def main():
    #### USER INPUT SECTION
    usage = """Script to receive an OrthoFinder Orthogroups.tsv file
    alongside the originally clustered FASTA to produce an output
    FASTA subset to just representative sequences."""
    p = argparse.ArgumentParser(description=usage)
    # Required
    p.add_argument("-i", dest="tsvFile",
                   required=True,
                   help="Input cluster file")
    p.add_argument("-f", dest="fastaFiles",
                   required=True,
                   nargs="+",
                   help="""Input FASTA files for sequences; must provide in same
                   order as in the OrthoFinder .tsv file header""")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output FASTA file name for representative sequences")
    # Optional
    p.add_argument("--prefix", dest="seqPrefix",
                   required=False,
                   help="""Optionally, if your sequence IDs have a prefix that involves a pipe (|)
                   which OrthoFinder has removed, you can specify it here to allow for ID indexing
                   to work correctly.""",
                   default="")
    args = p.parse_args()
    validate_args(args)
    
    # Load in FASTA file
    recordsList = [ SeqIO.to_dict(SeqIO.parse(f, "fasta") for f in args.fastaFiles) ]
    
    # Parse cluster file and produce output
    write_orthofinder_representatives(args.tsvFile, recordsList, args.outputFileName, args.seqPrefix)
    
    # Done!
    print('Program completed successfully!')

if __name__ == '__main__':
    main()

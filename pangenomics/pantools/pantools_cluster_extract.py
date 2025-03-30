#! python3
# mmseqs_cluster_extract.py
# Script to take in a cluster file from MMseqs2 alongside its FASTA
# and pull out just the representative sequences

import argparse, os, sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from Function_packages import ZS_SeqIO

# Various functions for program operations
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.clusterFile) and not os.path.islink(args.clusterFile):
        raise FileNotFoundError(f"I am unable to locate the pantools cluster file '{args.clusterFile}'")
    for fastaFile in args.fastaFiles:
        if not os.path.isfile(fastaFile) and not os.path.islink(fastaFile):
            raise FileNotFoundError(f"I am unable to locate the FASTA file '{fastaFile}'")
    # Validate output file location
    if os.path.exists(args.outputFileName):
        raise FileExistsError(f"The specified output file name already exists ({args.outputFileName})")

def parse_pantools_clusters(fileName):
    '''
    Parameters:
        fileName -- a string indicating the location of the Corset TSV file for parsing.
    Returns:
        clusterDict -- a dictionary with structure like:
                       {
                           0: ['seqID1', 'seqID2'],
                           1: ['seqID3'],
                           ...
                       }
    '''
    clusterDict = {}
    clusterNum = 0
    with open(fileName, "r") as fileIn:
        firstLine = True
        for line in fileIn:
            l = line.rstrip("\r\n ")
            if firstLine:
                firstLine = False
            elif l != "":
                clustID, seqIDs = line.rstrip("\r\n ").split(": ")
                clusterDict[clusterNum] = [ x.rsplit("#", maxsplit=1)[0] for x in seqIDs.split(" ") if x != "" ]
                clusterNum += 1
    return clusterDict

def pick_pantools_representatives(clusterDict, fastaCollection):
    repIDs = []
    for clusterID, seqIDs in clusterDict.items():
        # Get longest sequence
        longestSeq = ""
        for seqID in seqIDs:
            seq = fastaCollection[seqID]
            if len(seq) > len(longestSeq):
                longestSeq = seq
        # Store representative ID
        repIDs.append(longestSeq.name)
    return set(repIDs)

# Main call
def main():
    #### USER INPUT SECTION
    usage = """Script to receive an MMseqs2 cluster file in two-column TSV format
    where left column == representative ID, and right column == member sequence ID.
    Alongside the originally clustered FASTA(s), this script will produce an output
    FASTA subset to just representative sequences."""
    p = argparse.ArgumentParser(description=usage)
    # Required
    p.add_argument("-c", dest="clusterFile",
                   required=True,
                   help="Input cluster file")
    p.add_argument("-f", dest="fastaFiles",
                   required=True,
                   nargs="+",
                   help="Input FASTA file(s)")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output FASTA file name for representative sequences")
    args = p.parse_args()
    validate_args(args)
    
    # Parse cluster file
    clusterDict = parse_pantools_clusters(args.clusterFile)
    
    # Load in FASTA file(s)
    fastaCollection = ZS_SeqIO.FastaCollection(args.fastaFiles)
    
    # Pick representatives
    repIDs = pick_pantools_representatives(clusterDict, fastaCollection)
    
    # Produce output
    with open(args.outputFileName, "w") as fileOut:
        for repID in repIDs:
            record = fastaCollection[repID]
            fileOut.write(f">{record.name}\n{str(record)}\n")
    
    # Done!
    print('Program completed successfully!')

if __name__ == '__main__':
    main()

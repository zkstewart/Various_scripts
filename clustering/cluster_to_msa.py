#! python3
# cluster_to_msa.py
# Script to take in a cluster file from MMseqs2 alongside its FASTA
# and pull out MSAs for each

import os, sys, argparse
from Bio import SeqIO

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Function_packages import ZS_AlignIO, ZS_SeqIO

class DirectoryNotFoundError(Exception):
    pass

# Various functions for program operations
def validate_args(args):
    # Validate input file locations
    args.clusterFile = os.path.abspath(args.clusterFile)
    if not os.path.isfile(args.clusterFile):
        raise FileNotFoundError(f"-i file '{args.clusterFile}' does not exist!")
    args.fastaFile = os.path.abspath(args.fastaFile)
    if not os.path.isfile(args.fastaFile):
        raise FileNotFoundError(f"-f file '{args.fastaFile}' does not exist!")
    
    # Validate numeric arguments
    if args.threads < 1:
        raise ValueError("--threads must be >= 1")
    
    # Validate output file location
    args.outputFileLocation = os.path.abspath(args.outputFileLocation)
    if not os.path.exists(args.outputFileLocation):
        parentDir = os.path.dirname(args.outputFileLocation)
        if os.path.isdir(parentDir):
            os.mkdir(args.outputFileLocation)
            print(f"# Created '{args.outputFileLocation}' as part of argument validation")
        else:
            raise DirectoryNotFoundError(f"Cannot create '{args.outputFileLocation}' since its parent "
                                         f"dir '{parentDir}' also does not exist")

def parse_mms2_clusters(clusterFile):
    clusterDict = {}
    with open(clusterFile, "r") as fileIn:
        for line in fileIn:
            rep, member = line.rstrip("\r\n ").split("\t")
            clusterDict.setdefault(rep, [])
            clusterDict[rep].append(member)
    return clusterDict

# Main call
def main():
    #### USER INPUT SECTION
    usage = """Script to receive an MMseqs2 cluster file in two-column TSV format
    where left column == representative ID, and right column == member sequence ID.
    Alongside the originally clustered FASTA, this script will produce MSA outputs
    for each cluster."""
    p = argparse.ArgumentParser(description=usage)
    # Required
    p.add_argument("-i", dest="clusterFile",
                   required=True,
                   help="Input cluster file")
    p.add_argument("-f", dest="fastaFile",
                   required=True,
                   help="Input FASTA file")
    p.add_argument("-m", dest="molecule",
                   required=True,
                   choices=["nucleotide", "protein"],
                   help="Type of molecule to align")
    p.add_argument("-o", dest="outputFileLocation",
                   required=True,
                   help="Output location for cluster MSAs")
    p.add_argument("--threads", dest="threads",
                   required=False,
                   type=int,
                   help="Number of threads for MAFFT alignment (default==1)",
                   default=1)
    p.add_argument("--mafft", dest="mafft",
                   required=False,
                   help="Full path to mafft executable (default=find in PATH)",
                   default="mafft")
    args = p.parse_args()
    validate_args(args)
    
    # Parse cluster file
    clusterDict = parse_mms2_clusters(args.clusterFile)
    
    # Load in FASTA file
    records = ZS_SeqIO.FASTA(args.fastaFile)
    
    # Produce MSA outputs
    aligner = ZS_AlignIO.MAFFT(args.mafft, algorithm="einsi", thread=args.threads, molecule=args.molecule)
    for cluster, members in clusterDict.items():
        # Format a FASTA object with each member
        FASTA_obj = ZS_SeqIO.FASTA(fastaFile=None, isAligned=True)
        for memberID in members:
            fastaSeqObj = records[memberID]
            FASTA_obj.add(fastaSeqObj)
        
        # Align the FASTA
        alignedFASTA_obj = aligner.align(FASTA_obj)
        
        # Write to file
        outputFileName = os.path.join(args.outputFileLocation, f"{cluster}.fasta")
        alignedFASTA_obj.write(outputFileName, asAligned=True)
    
    # Done!
    print('Program completed successfully!')

if __name__ == '__main__':
    main()

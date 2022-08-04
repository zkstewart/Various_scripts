#! python3
# msa_concat.py
# Script for working with multiple MSA files (pre-aligned)
# to concatenate them into a single file prior to further
# analysis e.g., phylogeny.

import sys, argparse, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))) # 3 dirs up is where we find dependencies
from Function_packages import ZS_SeqIO

def validate_args(args):
    # Validate input data location
    if not os.path.isdir(args.msaDir):
        print('I am unable to locate the directory where the MSA FASTA files are (' + args.msaDir + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Handle file output
    if os.path.isfile(args.outputFileName):
        print(args.outputFileName + ' already exists. Specify a new file name or move/rename the existing file.')
        quit()

def main():
    #### USER INPUT SECTION
    usage = """%(prog)s will take all the files inside a specified directory and concatenate them
    into a single MSA suitable for subsequent analyses e.g., phylogenetics.
    """
    p = argparse.ArgumentParser(description=usage)
    # Required
    p.add_argument("-i", dest="msaDir", required=True,
                help="Specify the location of the directory containing one or more MSA FASTA files")
    p.add_argument("-o", dest="outputFileName", required=True,
                help="Specify the file name for the output MSA FASTA")
    args = p.parse_args()
    validate_args(args)
    
    # Locate all FASTA files
    files = [os.path.join(args.msaDir, file) for file in os.listdir(args.msaDir)]
    assert len(files) > 1, \
        "There's not enough files to concat in '{0}'?".format(args.msaDir)
    
    # Load FASTA files
    fastaObjs = []
    for file in files:
        f = ZS_SeqIO.FASTA(file, isAligned=True)
        f.make_uppercase()
        fastaObjs.append(f)
    
    # Ensure that all FASTA files contain the same sequences
    prevIDs = None
    prevFile = None
    for FASTA_obj in fastaObjs:
        if prevIDs == None:
            prevIDs = set(FASTA_obj.ids)
            prevFile = FASTA_obj.fileOrder[0][0]
        elif prevIDs != set(FASTA_obj.ids):
            print("FASTA files '{0}' and '{1}' are not compatible".format(prevFile, FASTA_obj.fileOrder[0][0]))
            quit()
    
    # Ensure the internal ordering of all MSAs are equivalent
    prevIDs = None
    isOrdered = True
    for FASTA_obj in fastaObjs:
        if prevIDs == None:
            prevIDs = FASTA_obj.ids
        elif prevIDs != FASTA_obj.ids:
            isOrdered = False
            break
    
    if not isOrdered:
        for FASTA_obj in fastaObjs:
            FASTA_obj.seqs.sort()
    
    # Concatenate FASTA objects
    for x in range(1, len(fastaObjs)):
        fastaObjs[0].concat(fastaObjs[x]) # it all gets concatenated into the first FASTA object
    
    # Write output
    fastaObjs[0].write(args.outputFileName, asAligned=True)
    
    # Done!
    print('Program completed successfully!')

if __name__ == "__main__":
    main()

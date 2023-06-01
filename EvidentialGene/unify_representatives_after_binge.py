#! python3
# unify_representatives_after_binge.py

# Script to take in the two other EvidentialGene files that weren't processed
# with BINge and unify them to have the same sequences, same IDs, in the same order.

import os, argparse
from Bio import SeqIO
from pyfaidx import Fasta

# Define classes
class FastaCollection:
    '''
    Wrapper for pyfaidx Fasta objects which allows multiple to be combined
    and queried as one logical entity.
    
    Parameters:
        fastaFiles -- a list of strings pointing to the locations of FASTA files
                      which are to be loaded in using pyfaidx.Fasta
    '''
    def __init__(self, fastaFiles):
        self.fastaFiles = fastaFiles
        self.records = []
        
        self._parse_fastas()
    
    def _parse_fastas(self):
        for fastaFile in self.fastaFiles:
            self.records.append(Fasta(fastaFile))
    
    def __getitem__(self, key):
        for records in self.records:
            try:
                return records[key]
            except:
                pass
        raise KeyError(f"'{key}' not found in collection")
    
    def __repr__(self):
        return (f"<FastaCollection object;num_records='{len(self.records)}';" +
                f"fastaFiles={self.fastaFiles}"
        )

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.representativesFasta):
        print(f'I am unable to locate the representatives FASTA file ({args.representativesFasta})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    for fileArg in [args.files1, args.files2]:
        for file in fileArg:
            if not os.path.isfile(file):
                print(f"I am unable to locate file at '{file}'")
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
    # Validate output file location
    suffixes = [ f.split(".")[-1].rstrip("._") for f in [args.files1[0], args.files2[0]] ]
    args.outputFileNames = []
    for suffix in suffixes:
        outFileName = args.outputPrefix + "." + suffix
        if os.path.isfile(outFileName):
            print(f'File already exists at output location ({outFileName})')
            print('Make sure you specify a unique file name and try again.')
            quit()
        args.outputFileNames.append(outFileName)

## Main
def main():
    # User input
    usage = """%(prog)s reads in the output FASTA of BINge_representatives.py which
    may be in transcript, CDS, or protein format. Here, you can provide one or more
    files from the other two formats so this program can obtain the corresponding
    sequences from these other formats, producing output files with the representative
    naming system.
    
    Note that this assumes all sequences can be found in the EvidentialGene files; you
    might need to run fix_evg_missing_orfs.py first.
    """
    p = argparse.ArgumentParser(description=usage)
    # Required
    p.add_argument("-r", dest="representativesFasta",
                   required=True,
                   help="Input representatives FASTA file")
    p.add_argument("-f1", dest="files1",
                   required=True,
                   nargs="+",
                   help="Input one or more files of the first format")
    p.add_argument("-f2", dest="files2",
                   required=True,
                   nargs="+",
                   help="Input one or more files of the second format")
    p.add_argument("-o", dest="outputPrefix",
                   required=True,
                   help="Output prefix for filtered file(s)")
    
    args = p.parse_args()
    validate_args(args)
    
    # Associate representatives with their cluster IDs
    representativesDict = {}
    with open(args.representativesFasta, "r") as fileIn:
        fastaRecords = SeqIO.parse(fileIn, "fasta")
        for record in fastaRecords:
            clusterID = record.id
            representativeID = record.description.split("representative=")[1].split(" ")[0]
            representativesDict[clusterID] = representativeID
    
    # Load transcripts from files1 and files2 into memory for quick access
    files1Records = FastaCollection(args.files1)
    files2Records = FastaCollection(args.files2)
    
    # Write output files
    for i in range(2):
        # Obtain values for this iteration
        outputFileName = args.outputFileNames[i]
        
        if i == 0:
            records = files1Records
        else:
            records = files2Records
        
        # Get sequences for each cluster, writing in order
        with open(outputFileName, "w") as fileOut:
            for clusterID, representativeID in representativesDict.items():
                seq = str(records[representativeID])
                fileOut.write(f">{clusterID} representative={representativeID}\n{seq}\n")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

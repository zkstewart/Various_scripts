#! python3
# fix_evg_missing_orfs.py

# Script to take in the aa and cds output files from EvidentialGene
# and fix any issues with missing main ORFs for sequences.

import os, argparse
from Bio import SeqIO

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isdir(args.dropsetDir):
        print(f'I am unable to locate the dropset directory ({args.dropsetDir})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.fastaFile):
        print(f'I am unable to locate the .fasta file ({args.fastaFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.aaFile):
        print(f'I am unable to locate the .aa file ({args.aaFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.cdsFile):
        print(f'I am unable to locate the .cds file ({args.cdsFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate output file location
    args.outputPrefix = args.outputPrefix.rstrip("._ ")
    args.outputFileNames = []
    for suffix in [".aa", ".cds", ".fixed_ids"]:
        outFileName = os.path.join(args.outputPrefix, suffix)
        if os.path.isfile(outFileName):
            print(f'File already exists at output location ({outFileName})')
            print('Make sure you specify a unique file name and try again.')
            quit()
        args.outputFileNames.append(outFileName)

## Main
def main():
    # User input
    usage = """%(prog)s reads in the output files of EvidentialGene and locates any
    instances where a sequence indicated in your .fasta files are missing from the 
    .aa and .cds files.
    """
    p = argparse.ArgumentParser(description=usage)
    # Required
    p.add_argument("-drop", dest="dropsetDir",
                   required=True,
                   help="Input EvidentialGene dropset directory")
    p.add_argument("-fasta", dest="fastaFile",
                   required=True,
                   help="Input EvidentialGene .fasta file")
    p.add_argument("-aa", dest="aaFile",
                   required=True,
                   help="Input EvidentialGene .aa file")
    p.add_argument("-cds", dest="cdsFile",
                   required=True,
                   help="Input EvidentialGene .cds file")
    p.add_argument("-o", dest="outputPrefix",
                   required=True,
                   help="Output prefix for filtered file(s)")
    
    args = p.parse_args()
    validate_args(args)
    
    # Store all transcript IDs
    transIDs = set()
    with open(args.fastaFile, "r") as fileIn:
        fastaRecords = SeqIO.parse(fileIn, "fasta")
        for record in fastaRecords:
            transIDs.add(record.id)
    
    # Locate all coding IDs which have a default ORF (not utrorf)
    seqIDs = set()
    with open(args.aaFile, "r") as fileIn:
        aaRecords = SeqIO.parse(fileIn, "fasta")
        for record in aaRecords:
            if not "utrorf" in record.id:
                seqIDs.add(record.id)
    
    # Locate sequences which have a utrorf but not a default ORF
    missingIDs = set()
    with open(args.aaFile, "r") as fileIn:
        aaRecords = SeqIO.parse(fileIn, "fasta")
        for record in aaRecords:
            baseID = record.id.split("utrorf")[0]
            if ("utrorf" in record.id and baseID not in seqIDs) or baseID not in transIDs:
                missingIDs.add(baseID)
    
    print(f"Note: Identified {missingIDs} missing sequences.")
    
    # Create new .aa and .cds files containing the missing ORFs from the dropset
    suffixes = [".aa", ".cds"]
    for suffix in suffixes:
        foundSequences = 0
        # Locate the relevant dropset file
        dropsetFile = [
            os.path.join(args.dropsetDir, f)
            for f in os.listdir(args.dropsetDir)
            if f.endswith(suffix)
        ]
        assert len(dropsetFile) == 1, \
            print(f"Could not locate the dropset '{suffix}' file!")
        dropsetFile = dropsetFile[0]
        
        # Load it in and write a revised output file containing any missing sequences
        fastaFile = args.aaFile if suffix == ".aa" else args.cdsFile
        outputFile = args.outputFileNames[suffixes.index(suffix)]
        
        with open(fastaFile, "r") as fastaIn, \
        open(dropsetFile, "r") as dropsetIn, \
        open(outputFile, "w") as fileOut:
            # Write a copy of the fasta file first
            for line in fastaIn:
                fileOut.write(line.rstrip("\r\n ") + "\n")
            
            # Parse in the dropset sequences, ...
            records = SeqIO.parse(dropsetIn, "fasta")
            
            # ... and write any of its records to file if needed
            for record in records:
                if record.id in missingIDs:
                    foundSequences += 1
                    fileOut.write(f">{record.description}\n{str(record.seq)}\n")
    
    print(f"Note: Restored {foundSequences} missing sequences.")
    
    # Write out the .fixed_ids file to allow any other potential fixes to occur
    with open(args.outputFileNames[-1], "w") as fileOut:
        for seqID in missingIDs:
            fileOut.write(f"{seqID}\n")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

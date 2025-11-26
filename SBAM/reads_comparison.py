#! python3
# Script to receive two BAM files and compare the
# reads which 1) occur in both files and 2) occur
# in only one or the other file.

import os, argparse
from contextlib import contextmanager

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.file1):
        raise FileNotFoundError(f"-1 file '{args.file1}' is not a file or does not exist!")
    if not os.path.isfile(args.file2):
        raise FileNotFoundError(f"-2 file '{args.file2}' is not a file or does not exist!")
    
    # Validate output file location
    if os.path.exists(args.outputFileName):
        raise FileExistsError(f"-o file '{args.outputFileName}' already exists!")

@contextmanager
def read_sbam_file(fileName):
    if fileName.endswith("bam"):
        with open(fileName, "rb") as f:
            yield f
    else:
        with open(fileName, "r") as f:
            yield f

def parse_for_readids(sbamFileName):
    readIDs = []
    with read_sbam_file(sbamFileName) as fileIn:
        for line in fileIn:
            if line.startswith("@"):
                continue
            else:
                sl = line.split("\t")
                readIDs.append(sl[0])
    return set(readIDs)

def main():
    # Establish parser
    usage = """%(prog)s will compare two SAM/BAM files to count the number of
    reads which 1) occur in both files and 2) occur in only one
    or the other file.
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-1", dest="file1",
                   required=True,
                   help="Input (S/B)AM file #1")
    p.add_argument("-2", dest="file2",
                   required=True,
                   help="Input (S/B)AM file #2")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Location to write output statistics")
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse each file
    file1IDs = parse_for_readids(args.file1)
    file2IDs = parse_for_readids(args.file2)
    
    # Compare read ID sets
    both = file1IDs.intersection(file2IDs)
    
    # Format statistics
    file1BothPct = (len(both) / len(file1IDs)) * 100
    file2BothPct = (len(both) / len(file2IDs)) * 100
    stats = [
        f"Files share {len(both)} reads in common",
        f"File 1 has {len(file1IDs)} reads; common reads are {file1BothPct}%",
        f"File 2 has {len(file2IDs)} reads; common reads are {file2BothPct}%"
    ]
    
    # Output statistics
    with open(args.outputFileName, "w") as fileOut:
        for stat in stats:
            print(stat)
            fileOut.write(stat + "\n")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

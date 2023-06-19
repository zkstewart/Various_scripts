#! python3
# reads_n_count
# Simply calculate the amount of gap ('N') characters in a fastq file to
# estimate the quality of a sequencing run

import os, argparse, locale, gzip
from Bio import SeqIO
from statistics import median, mean
from contextlib import contextmanager

locale.setlocale(locale.LC_ALL, '')

# Define functions for later use
## Validate arguments
def validate_args(args):
        # Validate input file locations
        if not os.path.isdir(args.inputDir):
            print('I am unable to locate the input dir containing FASTQ files (' + args.inputDir + ')')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
        # Handle file overwrites
        if os.path.isfile(args.output):
            print(args.output + ' already exists. Delete/move/rename this file and run the program again.')
            quit()

@contextmanager
def open_fq_file(filename):
    if filename.endswith(".gz"):
        with gzip.open(filename, "rt") as f:
            yield f
    else:
        with open(filename) as f:
            yield f

def main():
    usage = """%(prog)s reads in the FASTQ files contained within a provided directory
    and calculates the proportion of sequencing reads that are N's. It can handle
    plain FASTQ and gzip'd files, and expects files to have the .fq or .fastq file
    suffix.
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="inputDir",
                   help="input directory name")
    p.add_argument("-o", dest="output",
                   help="output file name")
    
    args = p.parse_args()
    validate_args(args)
    
    # Get fastq files list
    RECOGNISED_SUFFIXES = [".fq", ".fastq"]
    fastqFiles = []
    for file in os.listdir(args.inputDir):
        nogzFileName = file.lower().rstrip(".gz")
        if any([ nogzFileName.endswith(suffix) for suffix in RECOGNISED_SUFFIXES ]):
            fastqFiles.append(os.path.join(args.inputDir, file))
    
    # Exit if no files found
    if len(fastqFiles) == 0:
        print(f"No files were found with recognised FASTQ suffixes at '{args.inputDir}'")
        print("Program will exit now")
        quit()
    
    # Provide info if files were found
    else:
        print(f"{len(fastqFiles)} files were found; will begin processing now")
    
    # Count the number of N characters
    nchar = 0
    total = 0
    
    with open(args.output, "w") as fileOut:
        # Write header to file
        fileOut.write("filename\tnum_ns\ttotal_len\tn_pct\n")
        
        # Iterate through FQ files
        for fqFile in fastqFiles:
            # Get info from this file
            with open_fq_file(fqFile) as inFile:
                while True:
                    # Get this four-line block for a sequence read
                    seqID = inFile.readline()
                    if len(seqID) == 0: # exit condition if we reach EOF
                        break
                    
                    seq = inFile.readline().strip("\r\n ").lower()
                    qualID =  inFile.readline()
                    qual = inFile.readline()
                    assert seqID.startswith("@"), \
                        f"{fqFile} doesn't appear to be a standard four-line FASTQ file!! Bad seq ID == {seqID}"
                    
                    # Store relevant info
                    nchar += seq.count("n")
                    total += len(seq)
            nPct = nchar / total
            
            # Print and write to file
            fqDetails = f"{os.path.basename(fqFile)}\t{nchar}\t{total}\t{nPct}"
            print(fqDetails)
            fileOut.write(f"{fqDetails}\n")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

#! python3
# gc_count
# Simply calculate the amount of G and C characters in a fastq or fasta file

import os, argparse, locale, gzip
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
def open_gz_file(filename):
    if filename.endswith(".gz"):
        with gzip.open(filename, "rt") as f:
            yield f
    else:
        with open(filename) as f:
            yield f

def gc_fa_file(faFile):
    gcchar = 0
    total = 0
    with open_gz_file(faFile) as fileIn:
        for line in fileIn:
            if not line.startswith(">"):
                l = line.rstrip("\r\n ").lower()
                gcchar += (l.count("g") + l.count("c"))
                total += len(l)
    
    return gcchar, total

def gc_fq_file(fqFile):
    gcchar = 0
    total = 0
    with open_gz_file(fqFile) as fileIn:
        while True:
            # Get this four-line block for a sequence read
            seqID = fileIn.readline()
            if len(seqID) == 0: # exit condition if we reach EOF
                break
            
            seq = fileIn.readline().rstrip("\r\n ").lower()
            qualID = fileIn.readline()
            qual = fileIn.readline()
            assert seqID.startswith("@"), \
                f"{fqFile} doesn't appear to be a standard four-line FASTQ file!! Bad seq ID == {seqID}"
            
            # Store relevant info
            gcchar += (seq.count("g") + seq.count("c"))
            total += len(seq)
    
    return gcchar, total

def main():
    usage = """%(prog)s reads in the FASTA/Q files contained within a provided directory
    and calculates the proportion of sequencing reads that are N's. It can handle
    plain FASTA/Q and gzip'd files, and expects files to have the .fq/.fastq file
    or .fa/.fasta/.fna suffix.
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="inputDir",
                   help="input directory name")
    p.add_argument("-o", dest="output",
                   help="output file name")
    
    args = p.parse_args()
    validate_args(args)
    
    # Get fastq files list
    RECOGNISED_SUFFIXES = [".fq", ".fastq", ".fa", ".fasta", ".fna"]
    fastFiles = []
    for file in os.listdir(args.inputDir):
        nogzFileName = file.lower().rstrip(".gz")
        if any([ nogzFileName.endswith(suffix) for suffix in RECOGNISED_SUFFIXES ]):
            fastFiles.append(os.path.join(args.inputDir, file))
    
    # Exit if no files found
    if len(fastFiles) == 0:
        print(f"No files were found with recognised suffixes at '{args.inputDir}'")
        print("Program will exit now")
        quit()
    
    # Provide info if files were found
    else:
        print(f"{len(fastFiles)} files were found; will begin processing now")
    
    # Count the number of GC characters
    gcchar = 0
    total = 0
    
    with open(args.output, "w") as fileOut:
        # Write header to file
        fileOut.write("filename\tnum_gc\ttotal_len\tgc_pct\n")
        
        # Iterate through FQ files
        for fFile in fastFiles:
            # Figure out if A or Q file
            with open_gz_file(fFile) as fileIn:
                firstLine = fileIn.readline()
                if firstLine.startswith(">"):
                    isA = True
                elif firstLine.startswith("@"):
                    isA = False
                else:
                    raise Exception(
                        f"{fFile} not recognised as FASTA or FASTQ format!"
                    )
            
            # Handle based on file type
            if isA:
                gcchar, total = gc_fa_file(fFile)
            else:
                gcchar, total = gc_fq_file(fFile)
            gcPct = gcchar / total
            
            # Print and write to file
            fqDetails = f"{os.path.basename(fFile)}\t{gcchar}\t{total}\t{gcPct}"
            print(fqDetails)
            fileOut.write(f"{fqDetails}\n")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

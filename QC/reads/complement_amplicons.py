#! python 3
# complement_amplicons.py
# A script to receive a multiple sequence alignment of reference sequences
# alongside a FASTQ file of amplicon reads. This script will use a simple
# Levenshtein ratio to determine whether the read should be reverse complemented
# to match the reference sequence or not.

import os, argparse, gzip

from contextlib import contextmanager
from Bio import SeqIO
from Bio.Seq import Seq
from Levenshtein import ratio

def validate_args(args):
    # Validate input data
    if not os.path.isfile(args.fastqFile):
        print(f'I am unable to locate the reads FASTQ file ({args.fastqFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    
    args.referenceAmplicon = args.referenceAmplicon.replace("-", "")
    if len(args.referenceAmplicon) < 50:
        print(f"You've provided a reference amplicon that is very short ({len(args.referenceAmplicon)} bp sans any gaps)")
        print("I am assuming you've done something wrong, so I'm going to exit now.")
        quit()
    # Validate numeric arguments
    if args.minimumLevenshtein < 0:
        print("minimumLevenshtein value must be no less than 0...")
        quit()
    if args.minimumLevenshtein > 1:
        print("minimumLevenshtein value must be no greater than 1...")
        quit()
    # Handle file output
    if os.path.isdir(args.outputFileName):
        print('The specified output file name already exists as a directory!')
        print('You should specify a file name that does not exist.')
        print("Program will exit now.")
        quit()
    if os.path.isfile(args.outputFileName):
        print('The specified output file already exists! This program will not allowing overwriting.')
        print('You should specify a name that does not exist, or move the existing file elsewhere.')
        print("Program will exit now.")
        quit()

@contextmanager
def open_gz_file(filename):
    if filename.endswith(".gz"):
        with gzip.open(filename, "rt") as f:
            yield f
    else:
        with open(filename) as f:
            yield f

def main():
    usage = """%(prog)s will receive a FASTQ file of amplicon reads and a reference amplicon
    sequence. It will then compare each read to the reference amplicon, and if the reverse
    complement of the read is a better match, it will reverse complement the read and write
    it to a new FASTQ file.
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="fastqFile",
                   required=True,
                   help="Specify the location of the reads FASTQ file")
    p.add_argument("-a", dest="referenceAmplicon",
                   required=True,
                   help="""Specify an amplicon sequence to use as a guide for identifying
                   when reverse complementation is needed.""")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file name for the reverse complemented amplicons")
    # Optional
    p.add_argument("--minimum_lev", dest="minimumLevenshtein",
                   required=False,
                   type=float,
                   help="""Optionally, specify a minimum value for the Levenshtein
                   ratio between the read and the reference amplicon. Should be
                   a float between 0->1 (inclusive); default == 0 for no filtration,
                   but a value of 0.85-0.90 might be reasonable""",
                   default=0)
    
    args = p.parse_args()
    validate_args(args)
    
    numReads = 0
    numComplemented = 0
    numSkipped = 0
    with open_gz_file(args.fastqFile) as fastqIn, open(args.outputFileName, "w") as fileOut:
        records = SeqIO.parse(fastqIn, "fastq")
        
        for record in records:
            numReads += 1
            
            # Get sequence and its reverse complement
            seq = str(record.seq)
            rcSeq = str(Seq(seq).reverse_complement())
            outputFormat = record.format("fastq")
            
            # Calculate how closely the two complements match the reference
            seqMetric = ratio(seq, args.referenceAmplicon)
            rcSeqMetric = ratio(rcSeq, args.referenceAmplicon)
            
            # Enforce minimum Levenshtein ratio if relevant
            if seqMetric < args.minimumLevenshtein and rcSeqMetric < args.minimumLevenshtein:
                numSkipped += 1
                continue
            
            # If the reverse complement is better, update the FASTQ format for output
            if rcSeqMetric > seqMetric:
                qual = outputFormat.split("+")[-1].strip("\r\n")
                rcQual = qual[::-1]
                
                outputFormat = outputFormat.replace(seq, rcSeq).replace(qual, rcQual)
                numComplemented += 1
            
            # Write the complement that best matched the exemplar amplicon
            fileOut.write(outputFormat)
    
    # Print statistics then end program
    print("# complement_amplicons output statistics")
    print(f"# Processing of file: '{args.fastqFile}'")
    print(f"# > Of {numReads} amplicons, {numComplemented} were reverse complemented")
    print(f"# > This equates to {(numComplemented / numReads) * 100}% of amplicons")
    print(f"# > Additionally, {numSkipped} were skipped due to Levenshtein cut-off ({args.minimumLevenshtein})")
    print(f"# > This equates to {(numSkipped / numReads) * 100}% of amplicons")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

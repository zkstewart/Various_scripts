#! python3
# rotate_fasta_seq.py
# Simple script to perform minor rotations
# of specified sequences to start at new positions 
# in the original bit

import os, argparse
from Bio import SeqIO
from Bio.Seq import Seq

# Define functions
def validate_args(args):
    # Validate that all arguments have been provided
    for key, value in vars(args).items():
        if value == None:
            print(key + ' argument was not specified. Fix this and try again.')
            quit()
    # Validate input file location
    if not os.path.isfile(args.fastaFile):
        print('I am unable to locate the FASTA file (' + args.fastaFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate integer inputs
    if args.seqIndex < 0:
        print("seqIndex needs to be a positive integer")
        quit()
    if args.seqStart < 0:
        print("seqStart needs to be a positive integer")
        quit()
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print('File already exists at output location (' + args.outputFileName + ')')
        print('Make sure you specify a unique file name and try again.')
        quit()

def rotate_seq(fastaFileName, seqIndex, seqStart, outputFileName):
    with open(fastaFileName, "r") as fileIn, open(outputFileName, "w") as fileOut:
        records = SeqIO.parse(fileIn, "fasta")
        i = 0
        for record in records:
            i += 1
            if i != seqIndex:
                fileOut.write(record.format("fasta"))
            else:
                description = record.description
                seq = str(record.seq)
                newSeq = Seq(seq[seqStart-1: ] + seq[: seqStart-1])
                record.seq = newSeq
                fileOut.write(record.format("fasta"))
                #fileOut.write(">{0}\n{1}\n".format(description, newSeq))

def main():
    # User input
    usage = """%(prog)s receives a FASTA file and integer index for the sequence
    to rotate in the file (counting from 1). It will make the sequence start at the specified base
    and put the previous starting "bit" to the end of the sequence.
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-f", dest="fastaFile",
        help="Input FASTA file")
    p.add_argument("-si", dest="seqIndex", type=int,
        help="Specify the 1-based index of the sequence to modify")
    p.add_argument("-ss", dest="seqStart", type=int,
        help="Specify the 1-based position for the sequence to start at")
    p.add_argument("-o", dest="outputFileName",
        help="Output file name for the config file")
    args = p.parse_args()
    validate_args(args)

    # Convert MAF file
    rotate_seq(args.fastaFile, args.seqIndex, args.seqStart, args.outputFileName)

    print("Program completed successfully!")

if __name__ == "__main__":
    main()

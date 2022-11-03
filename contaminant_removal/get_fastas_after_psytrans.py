#! python3
# get_fastas_after_psytrans.py
# PsyTrans will have run with a CDS nucleotide file;
# afterwards, this script will enable obtaining the
# target sequences from .aa and .fasta files for example
# derived from EvidentialGene

import argparse, os
from Bio import SeqIO

# Various functions for program operations
def validate_args(args):
    # Validate file locations
    if not os.path.isfile(args.psyTransFile):
        print('I am unable to locate the PsyTrans FASTA file (' + args.psyTransFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    for fastaFile in args.inputFastaFiles:
        if not os.path.isfile(fastaFile):
            print('I am unable to locate the input FASTA file (' + fastaFile + ')')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    # Validate output file locations
    for fastaFile in args.inputFastaFiles:
        outFile = os.path.basename(fastaFile)
        if os.path.isfile(outFile):
            print('The intended output file name already exists in the current location (' + outFile + ')')
            print('Remove this file somehow or move to a new directory and try again.')
            quit()

# Main call
def main():
        #### USER INPUT SECTION
        usage = """After running PsyTrans, an output FASTA file will be produced. This script will
        get any corresponding amino acid or mRNA transcript FASTA files containing the corresponding
        sequences output from PsyTrans. Files will be produced with identical file names at the current
        working directory."""
        # Required
        p = argparse.ArgumentParser(description=usage)
        p.add_argument("-p", dest="psyTransFile",
                       required = True,
                       help="Specify the file output from PsyTrans")
        p.add_argument("-i", dest="inputFastaFiles", nargs = "+",
                       required = True,
                       help="Specify the file(s) you want to get corresponding sequences from.")
        args = p.parse_args()
        
        validate_args(args)
        
        # Parse PsyTrans FASTA sequence IDs
        with open(args.psyTransFile, "r") as fileIn:
            psyIDs = {line[1:].rstrip("\r\n "): None for line in fileIn if line.startswith(">")}
        
        # Parse all input FASTAs and get the appropriate sequences as output
        for fastaFile in args.inputFastaFiles:
            records = SeqIO.parse(fastaFile, "fasta")
            foundSeqs = []
            with open(os.path.basename(fastaFile), "w") as fileOut:
                for record in records:
                    if record.id in psyIDs:
                        fileOut.write(f">{record.description}\n{str(record.seq)}\n")
                        foundSeqs.append(record.id)
            if len(foundSeqs) == len(psyIDs):
                print(f"Found all sequences (n={len(psyIDs)}) in '{fastaFile}'")
            else:
                print(f"FAILED: Did not find all sequences (missing n={len(psyIDs) - len(foundSeqs)}) in '{fastaFile}' !!")
                missingSeqs = set(psyIDs.keys()).difference(set(foundSeqs))
                print("An excerpt of missed IDs is below:\n{0}".format(
                    "; ".join(list(missingSeqs)[0:10])
                ))
        
        # Done!
        print('Program completed successfully!')

if __name__ == '__main__':
        main()

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
        working directory. Note that this script anticipates file names that will come out of my own
        EvidentialGene pipeline, and so utrorf sequences are handled specifically in that light."""
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
        psyIDs = {}
        psyIDsNoUTR = {}
        with open(args.psyTransFile, "r") as fileIn:
            for line in fileIn:
                if line.startswith(">"):
                    seqid = line[1:].rstrip("\r\n ")
                    psyIDs[seqid] = None
                    psyIDsNoUTR[seqid.split("utrorf")[0]] = None
        
        # Parse all input FASTAs and get the appropriate sequences as output
        for fastaFile in args.inputFastaFiles:
            # Quickly check to see if we need to account for utrorf sequences
            hasUtrOrfs = any([r.id.endswith("utrorf") for r in SeqIO.parse(fastaFile, "fasta")])
            idsToCheckAgainst = psyIDs if hasUtrOrfs else psyIDsNoUTR
            
            # Parse all records for ID matches
            foundSeqs = []
            records = SeqIO.parse(fastaFile, "fasta")
            with open(os.path.basename(fastaFile), "w") as fileOut:
                for record in records:
                    if hasUtrOrfs:
                        if record.id in psyIDs:
                            fileOut.write(f">{record.description}\n{str(record.seq)}\n")
                            foundSeqs.append(record.id)
                    else:
                        if record.id in psyIDsNoUTR:
                            fileOut.write(f">{record.description}\n{str(record.seq)}\n")
                            foundSeqs.append(record.id)
            
            # Check to see if all sequences were found
            if len(foundSeqs) == len(idsToCheckAgainst):
                print(f"Found all sequences (n={len(foundSeqs)}) in '{fastaFile}'")
            else:
                print(f"FAILED: Did not find all sequences (missing n={len(idsToCheckAgainst) - len(foundSeqs)}) in '{fastaFile}' !!")
                missingSeqs = set(foundSeqs).difference(set(idsToCheckAgainst.keys()))
                print("An excerpt of missed IDs is below:\n{0}".format(
                "; ".join(list(missingSeqs)[0:10])
                ))
        
        # Done!
        print('Program completed successfully!')

if __name__ == '__main__':
        main()

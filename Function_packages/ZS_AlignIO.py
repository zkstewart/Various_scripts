#! python3
# ZS_AlignIO.py
# Contains various Classes to perform manipulations involving
# FASTA objects that are to be aligned or have been aligned.

# TBD
## Aligner class that can
## 1) align a FASTA and update its FastASeq.gap_seq values [-]
## 2) align nucleotides as peptide then convert back to nucleotide [-]
## 3) detect outliers in a MSA by sequence conservation [-]
## 4) detect outliers in a MSA by phylogeny mismatch [-]
## FASTA class (incl. key:value map function) [+]
## Stats functions that apply to FASTA class [-]

import os
from Bio.SeqIO.FastaIO import SimpleFastaParser
from collections import Counter

class Aligner:
    '''
    Relevant attributes include:
        TBD
    '''
    def __init__(self):
        # Validate input types
        # assert type(id).__name__ == "str"
        raise NotImplementedError()

class MSA:
    '''
    Relevant attributes include:
        TBD
    '''
    def __init__(self):
        raise NotImplementedError()

if __name__ == "__main__":
    pass

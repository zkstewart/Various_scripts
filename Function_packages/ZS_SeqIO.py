#! python3
# ZS_SeqIO.py
# Contains various Classes to perform manipulations involing
# FASTA sequences and MSAs.

# TBD
## Argument validation
## FASTA sequence class (name, altname, raw seq, gapped seq) [+]
## MSA class (incl. key:value map function) [WIP]
## Stats functions that apply to MSA class [-]

import os
from Bio.SeqIO.FastaIO import SimpleFastaParser

class FastASeq:
    def __init__(self, id, seq=None, alt=None, gapSeq=None):
        # Validate input types
        assert type(id).__name__ == "str"
        if seq != None:
            assert type(seq).__name__ == "str"
        if alt != None:
            assert type(alt).__name__ == "str"
        if gapSeq != None:
            assert type(gapSeq).__name__ == "str"
        assert seq != None or gapSeq != None
        
        # Set alt and/or ID
        self.id = id.replace(">", "").strip(" \t\r\n")
        self.alt = alt
        if alt != None:
            self.alt = alt.replace(">", "").strip(" \t\r\n")
        
        # Set seq and/or gap_seq
        if seq != None:
            self.seq = seq.replace("-", "").replace("\r", "").replace("\n", "").replace(" ", "")
        else: # gap_seq is guaranteed to exist is seq is None
            self.seq = gapSeq.replace("-", "").replace("\r", "").replace("\n", "").replace(" ", "")
        self.gap_seq = gapSeq
        if gapSeq != None:
            self.gap_seq = gapSeq.replace("\r", "").replace("\n", "").replace(" ", "")
    
    def get_str(self, withAlt=False, withGap=False):
        # Validate that object allows optional parameters
        if withAlt and type(self.alt).__name__ != "str":
            raise Exception("Alt ID not set on FastASeq with ID {0}".format(self.id))
        if withGap and type(self.gap_seq).__name__ != "str":
            raise Exception("Gap seq not set on FastASeq with ID {0}".format(self.id))
        
        # Produce expected output
        id = self.id if not withAlt else self.alt
        seq = self.seq if not withGap else self.gap_seq
        return ">{0}\n{1}".format(id, seq)
    
    def get_list(self, withAlt=False, withGap=False):
        # Validate that object allows optional parameters
        if withAlt and type(self.alt).__name__ != "str":
            raise Exception("Alt ID not set on FastASeq with ID {0}".format(self.id))
        if withGap and type(self.gap_seq).__name__ != "str":
            raise Exception("Gap seq not set on FastASeq with ID {0}".format(self.id))
        
        # Produce expected output
        id = self.id if not withAlt else self.alt
        seq = self.seq if not withGap else self.gap_seq
        return [id, seq]
    
    def extend(self, seq):
        '''
        This method assumes the user knows what they're doing.
        If this object already has a non-null gap_seq value, we
        will extend that with the method param as-is and extend the
        .seq value with the appropriately hyphen-neutralised value.
        Otherwise, we just modify the .seq attribute.
        '''
        # Validate input type
        assert type(seq).__name__ == "str"
        
        # Extend gap_seq and/or seq
        self.seq += seq.replace("-", "").replace("\r", "").replace("\n", "").replace(" ", "") # always exists
        if self.gap_seq != None:
            self.gap_seq += seq.replace("\r", "").replace("\n", "").replace(" ", "")
    
    def __str__(self):
        seq = self.seq if self.seq != None else self.gap_seq
        return ">{0}\n{1}".format(self.id, seq)
    
    def __repr__(self):
        return "FastASeq(id='{0}',seq='{1}',alt={4}{2}{4},gap_seq={5}{3}{5})".format(
            self.id, self.seq, self.alt, self.gap_seq,
            "'" if self.alt != None else "", "'" if self.gap_seq != None else "",
        )

class MSA:
    '''
    The MSA class is intended to serve as a container for FastASeq objects.
    It enables the FastASeq objects to be treated as a whole and have
    manipulations applied to them as a group e.g., setting their alt
    ID values, running an alignment and setting their gap_seq values, etc.
    
    To produce a MSA, you must provide it with a FASTA file. Although
    annoying (potentially), this allows us the ability to add more sequences
    to the MSA either through 1) concatenation (i.e., horizontally
    extending all sequences in order), or 2) addition (i.e., vertically
    expanding the alignment with more sequences).
    '''
    def __init__(self, fastaFile, isAligned=False):
        assert type(fastaFile).__name__ == "str"
        assert type(isAligned).__name__ == "bool"
        
        self.seqs = []
        self.fileOrder = []
        self.add(fastaFile, isAligned)
    
    def add(self, fastaFile, isAligned=False):
        # Validate value type and file existence
        assert type(fastaFile).__name__ == "str"
        assert type(isAligned).__name__ == "bool"
        if not os.path.isfile(fastaFile):
            raise Exception("{0} does not exist; can't load".format(fastaFile))
    
        # Load in file & add sequences to alignment
        with open(fastaFile, "r") as fileIn:
            for id, seq in SimpleFastaParser(fileIn):
                if isAligned:
                    seqObj = FastASeq(id=id, gapSeq=seq)
                else:
                    seqObj = FastASeq(id=id, seq=seq)
                self.seqs.append(seqObj)
        
        self.fileOrder.append([fastaFile, "add"])
    
    def concat(self, fastaFile):
        # Validate value type and file existence
        assert type(fastaFile).__name__ == "str"
        if not os.path.isfile(fastaFile):
            raise Exception("{0} does not exist; can't load".format(fastaFile))

        # Validate that concatenating file is compatible
        '''
        This is a crude estimator but it's quicker than
        using SimpleFastaParser and if the FASTA file format
        is broken we'll end up erroring out anyway
        '''
        seqCount = 0
        with open(fastaFile, "r") as fileIn:
            for line in fileIn:
                if line.startswith(">"):
                    seqCount += 1
        if seqCount != len(self.seqs):
            raise Exception("Concatenation not possible as sequence count differs")
        
        # Load in file & perform concatenation
        with open(fastaFile, "r") as fileIn:
            ongoingCount = 0
            for id, seq in SimpleFastaParser(fileIn):
                self.seqs[ongoingCount].extend(seq)
                ongoingCount += 1
        
        self.fileOrder.append([fastaFile, "concat"])
    
    def set_alt_ids_via_list(self, altList):
        raise NotImplementedError()
    
    def set_alt_ids_via_file(self, altFile):
        raise NotImplementedError()
    
    def __getitem__(self, key):
        return self.seqs[key]
    
    def __str__(self):
        addCount = len([f[0] for f in self.fileOrder if f[1] == "add"])
        concatCount = len([f[0] for f in self.fileOrder if f[1] == "concat"])
        
        return "MSA; contains {0} seqs; added {1} file{3}; concat {2} file{4}".format(
            len(self.seqs), addCount, concatCount, "s" if addCount > 1 else "",
            "s" if concatCount > 1 else ""
        )

if __name__ == "__main__":
    pass

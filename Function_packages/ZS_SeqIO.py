#! python3
# ZS_SeqIO.py
# Contains various Classes to perform manipulations involing
# FASTA sequences and MSAs.

# TBD
## Argument validation
## FASTA sequence class (name, altname, raw seq, gapped seq) [+]
## FASTA class (incl. key:value map function) [+]
## Stats functions that apply to FASTA class [-]

import os
from Bio.SeqIO.FastaIO import SimpleFastaParser

class FastASeq:
    '''
    Relevant attributes include:
        id -- A string indicating the name of this FASTA sequence
        alt -- A string representing an alternative way of naming this FASTA sequence
        seq -- A string of the nucleotide or protein sequence
        gap_seq -- A string of the nucleotide of protein sequence inclusive of gaps
                   introduced via sequence alignment
    '''
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
        '''
        We're not going to support long FASTA descriptions. It's rarely
        used by other programs anyway, and it makes a whole hassle.
        '''
        self.id = id.split(" ")[0].replace(">", "").strip(" \t\r\n")
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
        '''
        Params:
            withAlt -- A Boolean indicating whether the return should be the ID or alt
            withGap -- A Boolean indicating whether the return should be the seq or gap_seq
        Returns a string representation of this object with FASTA formatting
        '''
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
        '''
        Params:
            withAlt -- A Boolean indicating whether the return should be the ID or alt
            withGap -- A Boolean indicating whether the return should be the seq or gap_seq
        Returns a list containing:
            [0] -- the ID in accordance with withAlt param
            [1] -- the sequence in accordance with withGap param
        '''
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

class FASTA:
    '''
    The FASTA class is intended to serve as a container for FastASeq objects.
    It enables the FastASeq objects to be treated as a whole and have
    manipulations applied to them as a group e.g., setting their alt
    ID values, running an alignment and setting their gap_seq values, etc.
    
    To produce a FASTA, you must provide it with a FASTA file. Although
    annoying (potentially), this allows us the ability to add more sequences
    to the FASTA either through 1) concatenation (i.e., horizontally
    extending all sequences in order), or 2) addition (i.e., vertically
    expanding the alignment with more sequences).
    
    Relevant attributes include:
        seqs -- A list containing FastASeq objects
        fileOrder -- A list containing strings indicating the FASTA files 
                     that form part of this FASTA object
        consensus -- A FastaASeq object containing the consensus representation
                     of this FASTA object
        isAligned -- A Boolean indicating whether the FASTA file has been aligned
        
    '''
    def __init__(self, fastaFile, isAligned=False):
        assert type(fastaFile).__name__ == "str"
        assert type(isAligned).__name__ == "bool"
        
        self.seqs = []
        self.fileOrder = []
        self.consensus = None
        self.isAligned = isAligned
        self.add(fastaFile, isAligned)
    
    def add(self, fastaFile, isAligned=False):
        '''
        This method will read in the fastaFile and add the sequences vertically into
        the FASTA object (via appending to the internal .seqs list).
        
        Params:
            fastaFile -- a string providing the file name with or without the path to
                         find this file; will be interpreted by os.path
            isAligned -- a Boolean indicating whether the file contents have been aligned
        '''
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
        '''
        This method will read in the fastaFile and add the sequences horizontally into
        the FASTA object. This occurs via concatenating each sequence in fastaFile to the
        internal .seqs list assuming equivalent ordering. Note that we do not need to
        provide a isAligned Boolean for this method since we can figure out whether the
        file is expected to be aligned via an implicit assumption that presence of .gap_seq
        attributes in the underlying FastASeq objects means that it is aligned (see method
        header of the extend() method in FastASeq class).
        
        Params:
            fastaFile -- a string providing the file name with or without the path to
                         find this file; will be interpreted by os.path
        '''
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
            for _, seq in SimpleFastaParser(fileIn):
                self.seqs[ongoingCount].extend(seq)
                ongoingCount += 1
        
        self.fileOrder.append([fastaFile, "concat"])
    
    def set_alt_ids_via_list(self, altList):
        '''
        The number of values in the list must match that of the FASTA or an
        error will be raised. As mentioned in self_alt_ids_via_file, the order
        is also assumed to match that of the FASTA file(s) that are stored
        in self.seqs.
        '''
        # Validate value type & compatibility
        assert type(altList).__name__ == "list"
        if len(altList) != len(self.seqs):
            raise Exception(f"""
                Alt ID setting not possible as sequence count differs\n
                Num alt IDs={len(altList)}, Num FASTA seqs={len(self.seqs)}
            """)
        
        # Update FastASeq objects with IDs list
        for i in range(len(altList)):
            self.seqs[i].alt = altList[i]
        
    def set_alt_ids_via_dict(self, altDict):
        '''
        When using this function, all values in the FASTA must be discovered as
        a key within altDict or an error will be raised.
        '''
        # Validate value type & compatibility
        assert type(altDict).__name__ == "dict"
        numFound = sum([1 for FastASeq_obj in self.seqs if FastASeq_obj.id in altDict])            
        if numFound != len(self.seqs):
            raise Exception(f"""
                Alt ID setting not possible as not all sequences were discovered\n
                Num found IDs={numFound}, Num FASTA seqs={len(self.seqs)}
            """)
        
        # Update FastASeq objects with IDs list
        for i in range(len(self.seqs)):
            alt = altDict[self.seqs.id]
            self.seqs[i].alt = alt
    
    def set_alt_ids_via_file(self, altFile):
        '''
        This method accepts two kinds of alt file. The first is a simple list with one
        value per line. In this case, we assume the list is ordered equivalently to the
        original FASTA file that was loaded in. The second case is a key:value
        list separated by tabs. In this case the order can differ, and the list can
        contain more entries than present in the FASTA.
        '''
        # Validate value type and file existence
        assert type(altFile).__name__ == "str"
        if not os.path.isfile(altFile):
            raise Exception("{0} does not exist; can't update alt IDs".format(altFile))
        
        # Check what format the alt file is in
        isDict = False # isList is an implied value opposite to isDict
        with open(altFile, "r") as fileIn:
            for line in fileIn:
                if "\t" in line:
                    isDict = True
                break
               
        # Parse altFile
        altList, altDict = [], {}
        with open(altFile, "r") as fileIn:
            for line in fileIn:
                line = line.rstrip("\r\n ")
                if isDict:
                    l = line.split("\t")
                    altDict[l[0]] = l[1] # bidirectional indexing
                    altDict[l[1]] = l[0] # allows ID file to be orig:alt or alt:orig
                else:
                    altList.append(line)
        
        # Update FastASeq objects with our appropriate alt ID structure
        if isDict:
            self.set_alt_ids_via_dict(altDict)
        else:
            self.set_alt_ids_via_list(altList)
    
    def generate_consensus(self):
        # Validate sensibility and possibility of running this method
        if not self.isAligned:
            raise Exception("Consensus generation only relevant and/or possible if FASTA is aligned")
        ## Check if all FastASeq objects have the same length of their .gap_seq attribute
        raise NotImplementedError()
    
    def gc_content(self):
        raise NotImplementedError()
    
    def variable_sites(self):
        if self.consensus == None:
            raise Exception("Can't calculate variable sites since consensus generation has not occurred")
        raise NotImplementedError()
    
    def __getitem__(self, key):
        return self.seqs[key]
    
    def __str__(self):
        addCount = len([f[0] for f in self.fileOrder if f[1] == "add"])
        concatCount = len([f[0] for f in self.fileOrder if f[1] == "concat"])
        
        return "FASTA; contains {0} seqs; added {1} file{3}; concat {2} file{4}".format(
            len(self.seqs), addCount, concatCount, "s" if addCount > 1 else "",
            "s" if concatCount > 1 else ""
        )

if __name__ == "__main__":
    pass

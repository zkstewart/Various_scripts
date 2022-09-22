#! python3
# ZS_SeqIO.py
# Contains various Classes to perform manipulations involving
# FASTA sequences and MSAs.

import os, inspect
from Bio.SeqIO.FastaIO import SimpleFastaParser
from collections import Counter
from copy import deepcopy

class FastASeq:
    '''
    Relevant attributes include:
        id -- A string indicating the name of this FASTA sequence
        alt -- A string representing an alternative way of naming this FASTA sequence
        seq -- A string of the nucleotide or protein sequence
        gap_seq -- A string of the nucleotide of protein sequence inclusive of gaps
                   introduced via sequence alignment
    '''
    TRANSLATION_TABLE = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 'ACN': 'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CTN': 'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 'CCN': 'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 'CGN': 'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GTN': 'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 'GCN': 'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 'GGN': 'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TCN': 'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
        # Ambiguity handling
        'ATM':'I', 'ATW':'I', 'ATY':'I', 'ATH':'I',
        'ACM':'T', 'ACR':'T', 'ACW':'T', 'ACS':'T', 'ACY':'T', 'ACK':'T', 'ACV':'T', 'ACH':'T', 'ACD':'T', 'ACB':'T',
        'AAY':'N', 'AAR':'K',
        'AGY':'S', 'AGR':'R',
        'CTM':'L', 'CTR':'L', 'CTW':'L', 'CTS':'L', 'CTY':'L', 'CTK':'L', 'CTV':'L', 'CTH':'L', 'CTD':'L', 'CTB':'L',
        'CCM':'P', 'CCR':'P', 'CCW':'P', 'CCS':'P', 'CCY':'P', 'CCK':'P', 'CCV':'P', 'CCH':'P', 'CCD':'P', 'CCB':'P', 
        'CAY':'H', 'CAR':'Q',
        'CGM':'R', 'CGR':'R', 'CGW':'R', 'CGS':'R', 'CGY':'R', 'CGK':'R', 'CGV':'R', 'CGH':'R', 'CGD':'R', 'CGB':'R', 
        'GTM':'V', 'GTR':'V', 'GTW':'V', 'GTS':'V', 'GTY':'V', 'GTK':'V', 'GTV':'V', 'GTH':'V', 'GTD':'V', 'GTB':'V', 
        'GCM':'A', 'GCR':'A', 'GCW':'A', 'GCS':'A', 'GCY':'A', 'GCK':'A', 'GCV':'A', 'GCH':'A', 'GCD':'A', 'GCB':'A', 
        'GAY':'D', 'GAR':'E',
        'GGM':'G', 'GGR':'G', 'GGW':'G', 'GGS':'G', 'GGY':'G', 'GGK':'G', 'GGV':'G', 'GGH':'G', 'GGD':'G', 'GGB':'G', 
        'TCM':'S', 'TCR':'S', 'TCW':'S', 'TCS':'S', 'TCY':'S', 'TCK':'S', 'TCV':'S', 'TCH':'S', 'TCD':'S', 'TCB':'S', 
        'TTY':'F', 'TTR':'L',
        'TAY':'Y', 'TAR':'*',
        'TGY':'C', 
    }
    
    def __init__(self, id, seq=None, alt=None, gapSeq=None):
        # Validate input types
        assert isinstance(id, str)
        if seq != None:
            assert isinstance(seq, str)
        if alt != None:
            assert isinstance(alt, str)
        if gapSeq != None:
            assert isinstance(gapSeq, str)
        assert seq != None or gapSeq != None
        
        # Set alt and/or ID
        self.id = id.split(" ")[0].replace(">", "").strip(" \t\r\n")
        self.alt = alt
        if alt != None:
            self.alt = alt.replace(">", "").strip(" \t\r\n")
        
        # Set description line
        '''
        We're not supporting the use of descriptions directly in this Class, but we'll
        store the values so that users can directly interact with the attributes and
        expand upon behaviours not otherwise supported directly.
        '''
        self.description = id.replace(">", "").strip(" \t\r\n")
        
        # Set seq and/or gap_seq
        if seq != None:
            self.seq = seq.replace("-", "").replace("\r", "").replace("\n", "").replace(" ", "")
        else: # gap_seq is guaranteed to exist is seq is None
            self.seq = gapSeq.replace("-", "").replace("\r", "").replace("\n", "").replace(" ", "")
        self.gap_seq = gapSeq
        if gapSeq != None:
            self.gap_seq = gapSeq.replace("\r", "").replace("\n", "").replace(" ", "")
        
        # Helpful flag for checking data type
        self.isFastASeq = True
    
    def get_str(self, withAlt=False, withDescription=False, withGap=False):
        '''
        Params:
            withAlt -- A Boolean indicating whether the return should be alt ID
            withDescription -- A Boolean indicating whether the return should be the full description
            withGap -- A Boolean indicating whether the return should be the seq or gap_seq
        Returns:
            fastaStr -- a string representation of this object with FASTA formatting
        '''
        # Validate that object allows optional parameters
        if withAlt and not isinstance(self.alt, str):
            raise Exception("Alt ID not set on FastASeq with ID {0}".format(self.id))
        if withGap and not isinstance(self.gap_seq, str):
            raise Exception("Gap seq not set on FastASeq with ID {0}".format(self.id))
        assert isinstance(withDescription, bool)
        assert not (withDescription and withAlt), "Can't specify withDescription and withAlt at same time!"
        
        # Produce expected output
        if not withAlt and not withDescription:
            id = self.id
        elif withAlt:
            id = self.alt
        else:
            id = self.description
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
        assert isinstance(seq, str)
        
        # Extend gap_seq and/or seq
        self.seq += seq.replace("-", "").replace("\r", "").replace("\n", "").replace(" ", "") # always exists
        if self.gap_seq != None:
            self.gap_seq += seq.replace("\r", "").replace("\n", "").replace(" ", "")
    
    def get_translation(self, findBestFrame=False, strand=None, frame=None):
        '''
        This method will translate a nucleotide sequence into its peptide sequence
        with a standard translation table. Since this class is not yet aware of whether
        it is a nucleotide or protein already, the onus is on you as a user to not be
        dumb. An error will be raised if you try to translate a protein into a protein.
        
        This function allows you to constrain the finding of a best frame by specifying
        findBestFrame as True, and setting one of the strand or frame values. If you set
        strand=1, we'll find the best translation frame within the +ve strand. If you set
        frame=0, we'll find the best strand (+ve or -ve) with an ORF starting at frame 0.
        
        Params:
            findBestFrame -- a Boolean indicating whether you want this program to find
                             the longest translation for the sequence. This is relevant
                             when you don't know what frame the sequence is in. On short
                             sequences without any stop codons this might not work well.
                             If set, strand and frame arguments will be ignored.
            strand -- an integer indicating the strand the translation should occur on;
                      a positive 1 indicates forward (+ve) strand, and negative 1 is reverse
                      (-ve) strand. Only set this if you know you need a specific strand.
            frame -- an integer indicating the frame the translation has occurred in;
                     value ranges from 0 (first position in the sequence) to 2 (at the third
                     position in the sequence). Only set this if you know you need a specific
                     frame.
        Returns:
            protein -- a string of amino acid residues translated from this instance's .seq.
            strand -- an integer indicating the strand the translation has occurred on;
                      a positive 1 indicates forward (+ve) strand, and negative 1 is reverse
                      (-ve) strand.
            frame -- an integer indicating the frame the translation has occurred in;
                     value ranges from 0 (first position in the sequence) to 2 (at the third
                     position in the sequence).
        '''
        # Validate input type & values
        assert isinstance(findBestFrame, bool)
        assert isinstance(strand, int) or strand == None
        if isinstance(strand, int):
            assert strand in [1, -1], "Strand was provided but was not an appropriate value"
        assert isinstance(frame, int) or frame == None
        if isinstance(frame, int):
            assert frame in range(0, 3), "Frame was provided but was not an appropriate value"
        
        # Prevent incompatible parameter combination
        if findBestFrame == False:
            assert strand != None and frame != None, "If findBestFrame is False, you must set strand and frame values"
        if findBestFrame == True:
            "It doesn't make sense to findBestFrame if we've told it the strand and frame to search in"
            assert strand == None or frame == None, "If findBestFrame is True, strand and/or frame values must NOT be set"
        
        # Handle normal cases (findBestFrame is False)
        if not findBestFrame:
            if strand == 1 and frame == 0:
                return self._dna_to_protein(self.seq), strand, frame
            else:
                nuc = self.seq if strand == 1 else self.get_reverse_complement()
                length = 3 * ((len(nuc)-frame) // 3)
                frameNuc = nuc[frame:frame+length]
                frameProt = self._dna_to_protein(frameNuc)
                return frameProt, strand, frame
        # Handle other cases
        else:
            longest = [0, "", strand, frame] # [length, sequence, strand, frame]
            strandsToCheck = [1, -1] if strand == None else [strand]
            framesToCheck = range(3) if frame == None else range(frame, frame+1)
            for strand in strandsToCheck:
                nuc = self.seq if strand == 1 else self.get_reverse_complement()
                for frame in framesToCheck:
                    length = 3 * ((len(nuc)-frame) // 3)
                    frameNuc = nuc[frame:frame+length]
                    frameProt = self._dna_to_protein(frameNuc)
                    # Add bias to ORF if its translation ends in a stop codon
                    """If strand==1 has an ORF ending in a stop codon, and strand==-1
                    has an unbroken ORF not ending in a stop codon, without this
                    biaser, we'd end up picking the strand==-1 ORF. Ideally, we'd like
                    to prioritise strand==1, and also not penalise a sequence for
                    ending in a stop codon. Any other stop codon can go jump, however."""
                    if frameProt.count("*") == 1 and frameProt[-1] == "*":
                        orfs = [frameProt] # this prevents us splitting on "*" and counts its length
                    else:
                        orfs = frameProt.split("*")
                    for orf in orfs:
                        l = len(orf)
                        if l > longest[0]:
                            longest = [l, orf, strand, frame]
            _, protein, strand, frame = longest
            return protein, strand, frame
    
    def _dna_to_protein(self, dnaString):
        '''
        Hidden method for FastASeq to convert a string DNA sequence into
        its protein sequence. Peforms a simple codon substitution based on
        the static TRANSLATION_TABLE associated with this Class.
        
        Params:
            dnaString -- a string of nucleotides. Try to make sure this only
                         includes recognised nucleotides (no ambiguous characters)
                         or you'll end up with a bunch of 'X's in your translation.
        Returns:
            protein -- a string of amino acid residues translated from dnaString.
        '''
        protein = ""
        for i in range(0, len(dnaString), 3):
            codon = dnaString[i:i+3]
            if len(codon) < 3:
                continue
            protein += self.TRANSLATION_TABLE[codon.upper()] if codon.upper() in self.TRANSLATION_TABLE else "X"
        return protein
    
    def get_reverse_complement(self, staticSeq=None):
        '''
        Converts a nucleotide sequence into its reverse complement. Since the Class
        is unaware of whether it is a nucleotide or protein sequence, the onus is on you
        to not be dumb. You'll get a jumbled up protein if you run this method on a protein.
        This does not work in-place, it only returns a string as output.
        
        Can be run statically if the optional parameter staticSeq is provided
        
        Params:
            staticSeq -- optional only! If you specify this, this method acts like
                         a static function rather than being instance-based. It should
                         be a nucleotide as a string sequence.
        
        Returns:
            nucleotide -- a string of this .seq (or the staticSeq) after being
                          reverse complemented (e.g., by saying they're smarter
                          than they look)
        '''
        if staticSeq != None:
            assert isinstance(staticSeq, str)
            sequence = staticSeq
        else:
            sequence = self.seq

        reverseComplement = sequence[::-1].lower()
        reverseComplement = reverseComplement.replace('a', 'T')
        reverseComplement = reverseComplement.replace('t', 'A')
        reverseComplement = reverseComplement.replace('c', 'G')
        reverseComplement = reverseComplement.replace('g', 'C')
        # Ambiguity handling
        reverseComplement = reverseComplement.replace('m', 'K')
        reverseComplement = reverseComplement.replace('r', 'Y')
        reverseComplement = reverseComplement.replace('w', 'W')
        reverseComplement = reverseComplement.replace('s', 'S')
        reverseComplement = reverseComplement.replace('y', 'R')
        reverseComplement = reverseComplement.replace('k', 'M')
        reverseComplement = reverseComplement.replace('v', 'B')
        reverseComplement = reverseComplement.replace('h', 'D')
        reverseComplement = reverseComplement.replace('d', 'H')
        reverseComplement = reverseComplement.replace('b', 'V')
        return reverseComplement.upper()
    
    def trim_left(self, length, asAligned=False):
        '''
        This method trims the left side of this sequence to the length
        specified. Trimming the base sequence (asAligned=False) will behave as
        expected. Sequence in .gap_seq will be trimmed such that any thing in .seq
        that was trimmed will be considered a gap position in the .gap_seq attribute
        (if it exists).
        
        Trimming the aligned sequence (asAligned=True) will consider gap positions
        as part of the length. Any sequence it trims off will be updated in the .seq
        value as well.
        
        Params:
            length -- an integer indicating how many positions to trim.
            asAligned -- a Boolean indicating whether the raw sequence or aligned sequence
                         including gap characters should be trimmed.
        '''
        # Validate value type
        assert isinstance(length, int)
        assert length >= 0, "Length doesn't make sense as a negative integer!"
        if length == 0: # Trimming with length 0 has no impact, and this prevents [:-0] issue
            return

        # Validate aligned validity and possibility
        assert isinstance(asAligned, bool)
        if asAligned:
            if self.gap_seq == None:
                raise Exception(inspect.cleandoc("""
                                Sequence with ID {0} lacks a gap seq value; 
                                can't trim asAligned""".format(self.id)))
        
        # Perform .seq trimming
        if asAligned == False:
            trimmedSeq = self.seq[0:length]
            self.seq = self.seq[length:]
            
            # Update .gap_seq if relevant
            if self.gap_seq != None:
                trimLength = len(trimmedSeq)
                ongoingCount = 0
                
                for i in range(len(self.gap_seq)):
                    if ongoingCount >= trimLength:
                        break
                    
                    letter = self.gap_seq[i]
                    if letter == "-":
                        continue
                    else:
                        self.gap_seq = self.gap_seq[0:i] + "-" + self.gap_seq[i+1:]
                        ongoingCount += 1
                assert trimLength == ongoingCount, "Discovered that .seq and .gap_seq attributes are out of sync while trimming!"
        
        # Perform .gap_seq trimming
        else:
            trimmedSeq = self.gap_seq[0:length]
            self.gap_seq = self.gap_seq[length:]
            
            # Update .seq if relevant
            if trimmedSeq != "-"*length: # i.e., if we didn't just trim off gap sequence
                trimLength = len(trimmedSeq.replace("-",""))
                self.seq = self.seq[trimLength:]
    
    def trim_right(self, length, asAligned=False):
        '''
        This method trims the right side of all FASTA sequences to the length
        specified. Trimming the base sequence (asAligned=False) will behave as
        expected for the .seq value. Sequence in .gap_seq will be trimmed such
        that any thing in .seq that was trimmed will be considered a gap position
        in the .gap_seq attribute (if it exists).
        
        Trimming the aligned sequence (asAligned=True) will consider gap positions
        as part of the length. Any sequence it trims off will be updated in the .seq
        value as well.
        
        Params:
            length -- an integer indicating how many positions to trim.
            asAligned -- a Boolean indicating whether the raw sequences or aligned sequences
                         including gap characters should be trimmed.
        '''
        # Validate length type and sensibility
        assert isinstance(length, int)
        assert length >= 0, "Length doesn't make sense as a negative integer!"
        if length == 0: # Trimming with length 0 has no impact, and this prevents [:-0] issue
            return

        # Validate aligned validity and possibility
        assert isinstance(asAligned, bool)
        if asAligned:
            if self.gap_seq == None:
                raise Exception(inspect.cleandoc("""
                                Sequence with ID {0} lacks a gap seq value; 
                                can't trim asAligned""".format(self.id)))
    
        # Perform .seq trimming
        if asAligned == False:
            trimmedSeq = self.seq[-length:]
            self.seq = self.seq[:-length]
            
            # Update .gap_seq if relevant
            if self.gap_seq != None:
                trimLength = len(trimmedSeq)
                ongoingCount = 0
                
                for i in range(len(self.gap_seq)-1, -1, -1): # Loop through from right -> left
                    if ongoingCount >= trimLength:
                        break
                    
                    letter = self.gap_seq[i]
                    if letter == "-":
                        continue
                    else:
                        self.gap_seq = self.gap_seq[0:i] + "-" + self.gap_seq[i+1:]
                        ongoingCount += 1
                assert trimLength == ongoingCount, "Discovered that .seq and .gap_seq attributes are out of sync while trimming!"
        
        # Perform .gap_seq trimming
        else:
            trimmedSeq = self.gap_seq[-length:]
            self.gap_seq = self.gap_seq[:-length]
            
            # Update .seq if relevant
            if trimmedSeq != "-"*length: # i.e., if we didn't just trim off gap sequence
                trimLength = len(trimmedSeq.replace("-",""))
                self.seq = self.seq[:-trimLength]
    
    def make_uppercase(self):
        '''
        This method simply makes sure all .seq and .gap_seq letters are in upper case
        (if applicable).
        '''
        self.seq = self.seq.upper()
        if self.gap_seq != None:
            self.gap_seq = self.gap_seq.upper()
    
    def __str__(self):
        seq = self.seq if self.seq != None else self.gap_seq
        return ">{0}\n{1}".format(self.id, seq)
    
    def __repr__(self):
        return "FastASeq(id='{0}',seq='{1}',alt={4}{2}{4},gap_seq={5}{3}{5})".format(
            self.id, self.seq if len(self.seq) < 200 else "{0} ... {1}".format(self.seq[0:100], self.seq[-100:]),
            self.alt,
            "<None>" if self.gap_seq == None else self.gap_seq if len(self.gap_seq) < 200 else "{0} ... {1}".format(self.gap_seq[0:100], self.gap_seq[-100:]),
            "'" if self.alt != None else "", "'" if self.gap_seq != None else ""
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
    
    Alternatively, if you want a blank FASTA which you progressively add
    sequences to, specify fastaFile as None.
    
    Relevant attributes include:
        seqs -- A list containing FastASeq objects
        fileOrder -- A list containing strings indicating the FASTA files 
                     that form part of this FASTA object
        consensus -- A FastaASeq object containing the consensus representation
                     of this FASTA object
        isAligned -- A Boolean indicating whether the FASTA file has been aligned
    '''
    def __init__(self, fastaFile, isAligned=False):
        assert isinstance(fastaFile, str) or fastaFile == None
        assert isinstance(isAligned, bool)
        
        self.seqs = []
        self.fileOrder = []
        self.consensus = None
        self.isAligned = isAligned
        if fastaFile != None:
            self.add(fastaFile, isAligned)
        
        # Helpful flag for checking data type
        self.isFASTA = True
    
    def add(self, fasta, isAligned=False):
        '''
        This method will read in the fastaFile and add the sequences vertically into
        the FASTA object (via appending to the internal .seqs list).
        
        This method acts like an overloaded method i.e., you can provide either a string
        for a file location, or a ZS_SeqIO.FASTA object.
        
        Params:
            fasta -- a string providing the file name with or without the path to
                     find this file; will be interpreted by os.path... OR, a ZS_SeqIO.FASTA
                     object.
            isAligned -- a Boolean indicating whether the fasta contents have been aligned
        '''
        # Validate parameter type
        assert isinstance(isAligned, bool)
        
        # Diverge into one of the three overload methods
        if isinstance(fasta, str):
            if not os.path.isfile(fasta):
                raise Exception("{0} does not exist; can't add".format(fasta))
            else:
                self._add_from_file(fasta, isAligned)
        elif hasattr(fasta, "isFASTA") and fasta.isFASTA is True:
            self._add_from_FASTA_object(fasta, isAligned)
        elif hasattr(fasta, "isFastASeq") and fasta.isFastASeq is True:
            self._add_from_FastASeq_obj(fasta, isAligned)
        else:
            raise Exception("Unknown type '{0}' provided to FASTA.add()".format(type(fasta).__name__))
    
    def _add_from_FASTA_object(self, FASTA_obj, isAligned):
        # Stop potential infinite looping bug
        "If FASTA_obj == self, then it'll loop infinitely. I don't entirely understand why, I only have a hunch."
        FASTA_obj = deepcopy(FASTA_obj)
        # Validate that adding sequences will work
        "We do this check first so we don't modify the object unless we're sure it will work"
        for FastASeq_obj in FASTA_obj:
            if isAligned and FastASeq_obj.gap_seq == None:
                raise Exception("FastASeq object with id '{0}' should be aligned but has no .gap_seq value; adding sequences failed".format(FastASeq_obj.id))
        # Add sequences to alignment
        for FastASeq_obj in FASTA_obj:
            self.seqs.append(FastASeq_obj)
        
        self.fileOrder.append(["FASTA_object", "add"])
    
    def _add_from_FastASeq_obj(self, FastASeq_obj, isAligned):
        # Validate that adding sequences will work
        "We do this check first so we don't modify the object unless we're sure it will work"
        if isAligned and FastASeq_obj.gap_seq == None:
            raise Exception("FastASeq object with id '{0}' should be aligned but has no .gap_seq value; adding sequences failed".format(FastASeq_obj.id))
        # Add sequences to alignment
        self.seqs.append(FastASeq_obj)
        
        self.fileOrder.append(["FastASeq_object", "add"])
    
    def _add_from_file(self, fastaFile, isAligned):
        # Load in file & add sequences to alignment
        with open(fastaFile, "r") as fileIn:
            for id, seq in SimpleFastaParser(fileIn):
                if isAligned:
                    FastASeq_obj = FastASeq(id=id, gapSeq=seq)
                else:
                    FastASeq_obj = FastASeq(id=id, seq=seq)
                self.seqs.append(FastASeq_obj)
        
        self.fileOrder.append([fastaFile, "add"])
    
    def concat(self, fasta):
        '''
        This method will read in the fasta parameter and add the sequences horizontally into
        this (self) FASTA object. This occurs via concatenating each sequence in fasta to the
        internal .seqs list assuming equivalent ordering. Note that we do not need to
        provide a isAligned Boolean for this method since we can figure out whether the
        file is expected to be aligned via an implicit assumption that presence of .gap_seq
        attributes in the underlying FastASeq objects means that it is aligned (see method
        header of the extend() method in FastASeq class).
        
        This method acts like an overloaded method i.e., you can provide either a string
        for a file location, or a ZS_SeqIO.FASTA object.
        
        Params:
            fasta -- a string providing the file name with or without the path to
                     find this file; will be interpreted by os.path... OR, a ZS_SeqIO.FASTA
                     object.
        '''
        # Diverge into one of the three overload methods
        if isinstance(fasta, str):
            if not os.path.isfile(fasta):
                raise Exception("{0} does not exist; can't concat".format(fasta))
            else:
                self._concat_from_file(fasta)
        elif hasattr(fasta, "isFASTA") and fasta.isFASTA is True:
            self._concat_from_FASTA_obj(fasta)
        elif hasattr(fasta, "isFastASeq") and fasta.isFastASeq is True:
            self._concat_from_FastASeq_obj(fasta)
        else:
            raise Exception("Unknown type '{0}' provided to FASTA.concat()".format(type(fasta).__name__))
    
    def _concat_from_FASTA_obj(self, FASTA_obj):
        # Validate that concatenating object is compatible
        if len(FASTA_obj) != len(self.seqs):
            raise Exception("Concatenation not possible as sequence count differs")
        
        # Perform concatenation
        for i in range(len(FASTA_obj)):
            FastASeq_obj = FASTA_obj[i]
            if FastASeq_obj.gap_seq != None:
                self.seqs[i].extend(FastASeq_obj.gap_seq)
            else:
                self.seqs[i].extend(FastASeq_obj.seq)
        
        self.fileOrder.append(["FASTA_object", "concat"])
    
    def _concat_from_FastASeq_obj(self, FastASeq_obj):
        # Validate that concatenating object is compatible
        if len(self.seqs) != 1:
            raise Exception("Concatenation not possible as sequence count differs")
        
        # Perform concatenation
        if FastASeq_obj.gap_seq != None:
            self.seqs[0].extend(FastASeq_obj.gap_seq)
        else:
            self.seqs[0].extend(FastASeq_obj.seq)
        
        self.fileOrder.append(["FastASeq_object", "concat"])
    
    def _concat_from_file(self, fastaFile):
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
    
    def insert(self, index, FastASeq_obj):
        '''
        This method will insert a FastASeq object into this FASTA instance. Behaviourally
        this is the same as a list .insert() since our underlying datastructure of the
        .seq attribute is a list.
        
        Params:
            index -- any positive or negative integer.
            FastASeq_obj -- a ZS_SeqIO FastASeq object.
        '''
        # Validate value type and file existence
        #assert isinstance(FastASeq_obj, FastASeq) # Fails, lame
        assert type(FastASeq_obj).__name__ == "FastASeq" or type(FastASeq_obj).__name__ == "ZS_SeqIO.FastASeq"
        assert isinstance(index, int)
        
        self.seqs.insert(index, FastASeq_obj)
    
    def trim_left(self, length, asAligned=False):
        '''
        This method trims the left side of all FASTA sequences to the length
        specified. Trimming the base sequence (asAligned=False) will behave as
        expected for the .seq value. Sequence in .gap_seq will be trimmed such
        that any thing in .seq that was trimmed will be considered a gap position
        in the .gap_seq attribute (if it exists).
        
        Trimming the aligned sequence (asAligned=True) will consider gap positions
        as part of the length. Any sequence it trims off will be updated in the .seq
        value as well.
        
        Params:
            length -- an integer indicating how many positions to trim.
            asAligned -- a Boolean indicating whether the raw sequences or aligned sequences
                         including gap characters should be trimmed.
        '''
        # Validate value type
        assert isinstance(length, int)
        assert length >= 0, "Length doesn't make sense as a negative integer!"

        # Validate aligned validity and possibility
        assert isinstance(asAligned, bool)
        if asAligned:
            if not self.isAligned:
                raise Exception("FASTA object isn't flagged as being aligned; cant trim asAligned")
            for FastASeq_obj in self.seqs:
                if FastASeq_obj.gap_seq == None:
                    raise Exception(inspect.cleandoc("""
                                    Sequence with ID {0} lacks a gap seq value; 
                                    can't trim asAligned""".format(FastASeq_obj.id)))
        
        # Pass along trim method command to all contained FastASeq objects
        for FastASeq_obj in self.seqs:
            FastASeq_obj.trim_left(length, asAligned)
    
    def trim_right(self, length, asAligned=False):
        '''
        This method trims the right side of all FASTA sequences to the length
        specified. Trimming the base sequence (asAligned=False) will behave as
        expected for the .seq value. Sequence in .gap_seq will be trimmed such
        that any thing in .seq that was trimmed will be considered a gap position
        in the .gap_seq attribute (if it exists).
        
        Trimming the aligned sequence (asAligned=True) will consider gap positions
        as part of the length. Any sequence it trims off will be updated in the .seq
        value as well.
        
        Params:
            length -- an integer indicating how many positions to trim.
            asAligned -- a Boolean indicating whether the raw sequences or aligned sequences
                         including gap characters should be trimmed.
        '''
        # Validate length type and sensibility
        assert isinstance(length, int)
        assert length >= 0, "Length doesn't make sense as a negative integer!"

        # Validate aligned validity and possibility
        assert isinstance(asAligned, bool)
        if asAligned:
            if not self.isAligned:
                raise Exception("FASTA object isn't flagged as being aligned; cant trim asAligned")
            for FastASeq_obj in self.seqs:
                if FastASeq_obj.gap_seq == None:
                    raise Exception(inspect.cleandoc("""
                                    Sequence with ID {0} lacks a gap seq value; 
                                    can't trim asAligned""".format(FastASeq_obj.id)))
        
        # Pass along trim method command to all contained FastASeq objects
        for FastASeq_obj in self.seqs:
            FastASeq_obj.trim_right(length, asAligned)
    
    def slice_cols(self, start, end):
        '''
        This method returns a slice of the FASTA object's columns within the specified
        range (end non-inclusive same as range() function). It returns this as a NEW object
        rather than modifying the existing object.
        
        Note that this method is only relevant for aligned FASTA objects. It will return
        an error if .isAligned is not True.
        
        Params:
            start -- an integer indicating which column to start from (inclusive, 0-based)
            end -- an integer indicating which column to stop at (non-inclusive, 0-based)
        Returns:
            FASTA_obj -- a new FASTA object.
        '''
        # Validate that this object is aligned
        if not self.isAligned:
            raise Exception("FASTA object isn't flagged as being aligned; cant slice columns")
        for FastASeq_obj in self.seqs:
            if FastASeq_obj.gap_seq == None:
                raise Exception(inspect.cleandoc("""
                                Sequence with ID {0} lacks a gap seq value; 
                                can't slice columns""".format(FastASeq_obj.id)))
        
        # Validate value types and sensibility
        assert isinstance(start, int)
        assert isinstance(end, int)
        assert start >= 0, "start doesn't make sense as a negative integer!"
        assert end <= len(self.seqs[0].gap_seq), "end doesn't make sense if it's longer than the aligned length ({0})!".format(len(self.seqs[0].gap_seq))
        
        # Perform slice operation
        '''
        We can consider this slice operation as a kind of "trimming" where we want to trim
        from the right and left any sequence surrounding our chosen slice.
        '''
        result_FASTA_obj = deepcopy(self) # make a copy of this object for making changes to
        startTrim = start
        endTrim = len(self.seqs[0].gap_seq) - end
        for FastASeq_obj in result_FASTA_obj:
            FastASeq_obj.trim_left(startTrim, asAligned=True)
            FastASeq_obj.trim_right(endTrim, asAligned=True)
        return result_FASTA_obj
    
    def cut_cols(self, start, end):
        '''
        This method returns the FASTA object minus the columns within the specified
        range (end non-inclusive same as range() function). It returns this as a NEW object
        rather than modifying the existing object.
        
        Note that this method is only relevant for aligned FASTA objects. It will return
        an error if .isAligned is not True.
        
        Params:
            start -- an integer indicating which column to start from (inclusive, 0-based)
            end -- an integer indicating which column to stop at (non-inclusive, 0-based)
        Returns:
            FASTA_obj -- a new FASTA object.
        '''
        # Validate that this object is aligned
        if not self.isAligned:
            raise Exception("FASTA object isn't flagged as being aligned; cant slice columns")
        for FastASeq_obj in self.seqs:
            if FastASeq_obj.gap_seq == None:
                raise Exception(inspect.cleandoc("""
                                Sequence with ID {0} lacks a gap seq value; 
                                can't slice columns""".format(FastASeq_obj.id)))
        
        # Validate value types and sensibility
        assert isinstance(start, int)
        assert isinstance(end, int)
        assert start >= 0, "start doesn't make sense as a negative integer!"
        assert end <= len(self.seqs[0].gap_seq), "end doesn't make sense if it's longer than the aligned length ({0})!".format(len(self.seqs[0].gap_seq))
        
        # Perform cut operation
        '''
        We can consider this cut operation as two separate operations. First, we
        trim from the right of our left-most fragment, and we also trim from the
        left of our right-most fragment. Second, we concatenate the fragments
        together.
        '''
        left_FASTA_obj = deepcopy(self)
        right_FASTA_obj = deepcopy(self)
        
        leftTrim = len(self.seqs[0].gap_seq) - start
        for left_FastASeq_obj in left_FASTA_obj:
            left_FastASeq_obj.trim_right(leftTrim, asAligned=True)
        
        rightTrim = end
        for right_FastASeq_obj in right_FASTA_obj:
            right_FastASeq_obj.trim_left(rightTrim, asAligned=True)
        
        # Merge the left and right segments back together & return
        left_FASTA_obj.concat(right_FASTA_obj)
        
        return left_FASTA_obj
    
    def set_alt_ids_via_list(self, altList):
        '''
        The number of values in the list must match that of the FASTA or an
        error will be raised. As mentioned in self_alt_ids_via_file, the order
        is also assumed to match that of the FASTA file(s) that are stored
        in self.seqs.
        '''
        # Validate value type & compatibility
        assert isinstance(altList, list)
        if len(altList) != len(self.seqs):
            raise Exception(inspect.cleandoc(f"""
                Alt ID setting not possible as sequence count differs.
                Num alt IDs={len(altList)}, Num FASTA seqs={len(self.seqs)}
            """))
        
        # Update FastASeq objects with IDs list
        for i in range(len(altList)):
            self.seqs[i].alt = altList[i]
    
    def set_ids_via_list(self, idsList):
        '''
        The number of values in the list must match that of the FASTA or an
        error will be raised. The order is assumed to match that of the FASTA
        file(s) that are stored in self.seqs.
        '''
        # Validate value type & compatibility
        assert isinstance(idsList, list)
        if len(idsList) != len(self.seqs):
            raise Exception(inspect.cleandoc(f"""
                ID setting not possible as sequence count differs.
                Num IDs={len(idsList)}, Num FASTA seqs={len(self.seqs)}
            """))
        
        # Update FastASeq objects with IDs list
        for i in range(len(idsList)):
            self.seqs[i].id = idsList[i]
        
    def set_alt_ids_via_dict(self, altDict):
        '''
        When using this function, all values in the FASTA must be discovered as
        a key within altDict or an error will be raised.
        '''
        # Validate value type & compatibility
        assert isinstance(altDict, dict)
        numFound = sum([1 for FastASeq_obj in self.seqs if FastASeq_obj.id in altDict])            
        if numFound != len(self.seqs):
            raise Exception(inspect.cleandoc(f"""
                Alt ID setting not possible as not all sequences were discovered
                Num found IDs={numFound}, Num FASTA seqs={len(self.seqs)}
            """))
        
        # Update FastASeq objects with IDs list
        for i in range(len(self.seqs)):
            alt = altDict[self.seqs.id]
            self.seqs[i].alt = alt
    
    def set_ids_via_dict(self, idsDict):
        '''
        When using this function, all values in the FASTA must be discovered as
        a key within idsDict or an error will be raised.
        '''
        # Validate value type & compatibility
        assert isinstance(idsDict, dict)
        numFound = sum([1 for FastASeq_obj in self.seqs if FastASeq_obj.id in idsDict])            
        if numFound != len(self.seqs):
            raise Exception(inspect.cleandoc(f"""
                ID setting not possible as not all sequences were discovered
                Num found IDs={numFound}, Num FASTA seqs={len(self.seqs)}
            """))
        
        # Update FastASeq objects with IDs list
        for i in range(len(self.seqs)):
            id = idsDict[self.seqs.id]
            self.seqs[i].id = id
    
    def set_alt_ids_via_file(self, altFile):
        '''
        This method accepts two kinds of alt file. The first is a simple list with one
        value per line. In this case, we assume the list is ordered equivalently to the
        original FASTA file that was loaded in. The second case is a key:value
        list separated by tabs. In this case the order can differ, and the list can
        contain more entries than present in the FASTA.
        '''
        # Validate value type and file existence
        assert isinstance(altFile, str)
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
    
    def set_ids_via_file(self, idsFile):
        '''
        This method accepts two kinds of IDs file. The first is a simple list with one
        value per line. In this case, we assume the list is ordered equivalently to the
        original FASTA file that was loaded in. The second case is a key:value
        list separated by tabs. In this case the order can differ, and the list can
        contain more entries than present in the FASTA.
        '''
        # Validate value type and file existence
        assert isinstance(idsFile, str)
        if not os.path.isfile(idsFile):
            raise Exception("{0} does not exist; can't update IDs".format(idsFile))
        
        # Check what format the alt file is in
        isDict = False # isList is an implied value opposite to isDict
        with open(idsFile, "r") as fileIn:
            for line in fileIn:
                if "\t" in line:
                    isDict = True
                break
               
        # Parse altFile
        idsList, idsDict = [], {}
        with open(idsFile, "r") as fileIn:
            for line in fileIn:
                line = line.rstrip("\r\n ")
                if isDict:
                    l = line.split("\t")
                    idsDict[l[0]] = l[1] # bidirectional indexing
                    idsDict[l[1]] = l[0] # allows ID file to be orig:alt or alt:orig
                else:
                    idsList.append(line)
        
        # Update FastASeq objects with our appropriate alt ID structure
        if isDict:
            self.set_ids_via_dict(idsDict)
        else:
            self.set_ids_via_list(idsList)
    
    def make_uppercase(self):
        '''
        This method simply makes sure all FastASeq letters are in upper case.
        '''
        for FastASeq_obj in self.seqs:
            FastASeq_obj.make_uppercase()
    
    def generate_consensus(self):
        '''
        This method produces a very basic consensus sequence through a voting mechanism
        whereby the most common nucleotide/amino acid will be the representative. This
        might not be biologically valid, but it provides a metric by which human investigations
        can occur especially via calculation of variable site statistics.
        
        Method sets self.consensus, but also returns:
            consensus -- A string of the consensus sequence including gaps
        '''
        # Validate sensibility and possibility of running this method
        if not self.isAligned:
            raise Exception("Consensus generation only relevant and/or possible if FASTA is aligned")
        
        prevLength = None
        for FastASeq_obj in self.seqs:
            if FastASeq_obj.gap_seq == None:
                raise Exception("Can't generate consensus if FastaASeq .gap_seq attribute not set")
            
            thisLength = len(FastASeq_obj.gap_seq)
            if prevLength != None and prevLength != thisLength:
                raise Exception(".gap_seq attributes are set but alignment length differs")
            
            prevLength = thisLength
        
        # Generate consensus sequence
        self.consensus = ""
        for i in range(prevLength):
            positionList = []
            for x in range(len(self.seqs)):
                positionList.append(self.seqs[x].gap_seq[i])
            positionCount = Counter(positionList)
            self.consensus += positionCount.most_common(1)[0][0] # most_common(1) gives a list with (residue, count) tuple
        
        return self.consensus
        
    def gc_content(self):
        '''
        Calculates the proportion of sequence that is G/C out of all nucleotides in all sequences.
        i.e., in a FASTA with 3 sequences of length 10bp each, if each sequence has 3 G's or C's
        in them, our GC calculation will be:
            (3*3) / (3*10) == 0.3
        
        Note that this function only has biological relevance if the sequences are nucleotides.
        This function makes no checks to assert that this is the case. It's on the user to be
        sensible.
        
        Returns:
            proportion -- A float value from 0->1 representing the proportion of sequences that
                          are G or C.
        '''
        gcCount = 0
        totalCount = 0
        for FastASeq_obj in self.seqs:
            for nucleotide in FastASeq_obj.seq.lower():
                if nucleotide in ["g", "c"]:
                    gcCount += 1
                totalCount += 1
        
        return gcCount / totalCount
    
    def conserved_proportion(self):
        '''
        Calculates the proportion of sequence positions that matches the consensus sequence.
        In an alignment like this:
        
            CONSENSUS: ABCDEF
            SEQ1:      ABCDEF
            SEQ2:      ABCABC
        
        Iterating through each position, we'd determine how many share the same sequence as
        the consensus. We'd get a result like this:
        
            POSITION:  0  |  1  |  2  |  3  |  4  |  5
            SHARED:    1  |  1  |  1  | .5  | .5  | .5
        
        We then average the shared values i.e.,: 
        
            summed shared / length of consensus == (1 + 1 + 1 + .5 + .5 + .5) / 6 == 0.75
        
        Hence, returning a value of 0.75 indicating that 75 percent of the sequence positions
        are shared in common or, informally, are "conserved". Note that we ignore gaps in
        the consensus sequence to avoid bias from alignments where a single sequence is introducing
        large gaps, artificially inflating the proportion.
        '''
        # Validate that we can run this method successfully
        if self.consensus == None:
            raise Exception("Can't calculate variable sites since consensus generation has not occurred")
        
        # Perform calculation
        shared = []
        for i in range(len(self.consensus)):
            thisPosition = self.consensus[i].lower()
            if thisPosition == "-":
                continue
            
            positionProportion = []
            for x in range(len(self.seqs)):
                positionProportion.append(1 if self.seqs[x].gap_seq[i].lower() == thisPosition else 0)
            shared.append(sum(positionProportion) / len(positionProportion))
        
        if sum(shared) == 0:
            return 0
        else:
            return sum(shared) / len(shared)
    
    def write(self, outputFileName, withAlt=False, withDescription=False, asAligned=False, withConsensus=False):
        '''
        Writes the FASTA object out to a file in appropriate FASTA formatting. Various
        method params allow the customisation of how the file is produced.
        
        Params:
            outputFileName -- a string indicating the file location to write to. This file
                              must not exist or an error will be raised.
            withAlt -- a Boolean indicating whether sequences should be labelled with their
                       ID value or with the alternative IDs provided. An error will be raised
                       if no alternative IDs exist.
            withDescription -- a Boolean indicating whether sequences should be labelled with
                               their ID value or with the full description name. An error will
                               be raised if you specify this and withAlt simultaneously.
            asAligned -- a Boolean indicating whether the raw sequences or aligned sequences
                         including gap characters should be output. An error will be raised
                         if .isAligned == False or if any sequences lack a .gap_seq value.
            withConsensus -- a Boolean indicating whether the output file should be inclusive
                             of an additional sequence at the top with ">consensus" ID and the
                             generated consensus sequence as its value.
        '''
        # Validate alternative ID validity and possibility
        assert isinstance(withAlt, bool)
        if withAlt:
            for FastASeq_obj in self.seqs:
                if FastASeq_obj.alt == None:
                    raise Exception(inspect.cleandoc("""
                                    Sequence with ID {0} lacks an alt ID; 
                                    can't write withAlt""".format(FastASeq_obj.id)))
        
        # Validate description validity
        assert isinstance(withDescription, bool)
        assert not (withDescription and withAlt), "Can't specify withDescription and withAlt at same time!"
        
        # Validate aligned validity and possibility
        assert isinstance(asAligned, bool)
        if asAligned:
            if not self.isAligned:
                raise Exception("FASTA object isn't flagged as being aligned; cant write asAligned")
            for FastASeq_obj in self.seqs:
                if FastASeq_obj.gap_seq == None:
                    raise Exception(inspect.cleandoc("""
                                    Sequence with ID {0} lacks a gap seq value; 
                                    can't write asAligned""".format(FastASeq_obj.id)))
        
        # Validate consensus validity and possibility
        assert isinstance(withConsensus, bool)
        if withConsensus:
            if self.consensus == None:
                raise Exception("FASTA object doesn't have a consensus sequence; cant write withConsensus")
        
        # Validate output value types and file non-existence
        assert isinstance(outputFileName, str)
        if os.path.isfile(outputFileName):
            raise Exception("{0} already exists; can't write output file".format(outputFileName))

        # Actually write the output file
        with open(outputFileName, "w") as fileOut:
            if withConsensus:
                fileOut.write(">consensus\n{0}\n".format(self.consensus))
            for FastASeq_obj in self.seqs:
                s = FastASeq_obj.get_str(withAlt=withAlt, withGap=asAligned, withDescription=withDescription)
                fileOut.write("{0}\n".format(s))
    
    @property
    def ids(self):
        return [seq.id for seq in self.seqs]
    
    @ids.setter
    def ids(self, value):
        if isinstance(value, str):
            self.set_ids_via_file(value)
        elif isinstance(value, dict):
            self.set_ids_via_dict(value)
        elif isinstance(value, list):
            self.set_ids_via_list(value)
    
    @property
    def alts(self):
        return [seq.alt for seq in self.seqs]
    
    @alts.setter
    def alts(self, value):
        if isinstance(value, str):
            self.set_alt_ids_via_file(value)
        elif isinstance(value, dict):
            self.set_alt_ids_via_dict(value)
        elif isinstance(value, list):
            self.set_alt_ids_via_list(value)
    
    def __iter__(self):
        return iter(self.seqs)
    
    def __len__(self):
        return len(self.seqs)
    
    def __getitem__(self, key):
        if isinstance(key, int):
            return self.seqs[key]
        else:
            for value in self.seqs:
                if value.id == key or value.alt == key or value.description == key:
                    return value
    
    def __contains__(self, item):
        return True if self[item] is not None else False
    
    def __str__(self):
        addCount = len([f[0] for f in self.fileOrder if f[1] == "add"])
        concatCount = len([f[0] for f in self.fileOrder if f[1] == "concat"])
        
        return "FASTA; contains {0} seqs; added {1} file{3}; concat {2} file{4}".format(
            len(self.seqs), addCount, concatCount, "s" if addCount > 1 else "",
            "s" if concatCount > 1 else ""
        )
    
    def __repr__(self):
        return "FASTA(seqs={0} <FastASeq> objs, consensus={1}, isAligned={2}, fileOrder={3})".format(
            len(self.seqs), self.consensus, self.isAligned,
            str(self.fileOrder) if len(str(self.fileOrder)) < 200 else "{0} ... {1}".format(str(self.fileOrder)[0:100], str(self.fileOrder)[-100:])
        )

if __name__ == "__main__":
    pass

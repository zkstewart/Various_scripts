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

import os, platform, sys, subprocess, hashlib, time, random
from copy import deepcopy
from Bio.Align.Applications import MafftCommandline

sys.path.append(os.path.dirname(__file__))
from ZS_SeqIO import FASTA

class MAFFT:
    '''
    The MAFFT Class is intended to be used alongside the ZS_SeqIO.FASTA Class to
    handle the alignment of FASTA sequences with OOP magics. In short, it provides
    easy access to MAFFT's E-INSi and L-INSi algorithms with default parameters
    otherwise.
    
    It's primary purpose is the .run() method which receives a ZS_SeqIO.FASTA instance
    and updates the .gap_seq attributes of its underlying ZS_SeqIO.FastASeq objects.
    '''
    def __init__(self, mafftDir):
        # Validate input type and location
        assert isinstance(mafftDir, str)
        assert os.path.isdir(mafftDir)
        if platform.system() == "Windows":
            if not os.path.isfile(os.path.join(mafftDir, "mafft.bat")):
                raise Exception("{0} does not exist".format(os.path.join(mafftDir, "mafft.bat")))
        else:
            if not os.path.isfile(os.path.join(mafftDir, "mafft")) and not os.path.isfile(os.path.join(mafftDir, "mafft.exe")):
                raise Exception("mafft or mafft.exe does not exist at {0}".format(mafftDir))
        
        # Establish commandline function
        if platform.system() == "Windows":
            self.cline = MafftCommandline(os.path.join(mafftDir, "mafft.bat"))
        else:
            if os.path.isfile(os.path.join(mafftDir, "mafft")):
                self.cline = MafftCommandline(os.path.join(mafftDir, "mafft"))
            else:
                self.cline = MafftCommandline(os.path.join(mafftDir, "mafft.exe"))
        
        # Set default attributes
        self.mafftDir = mafftDir
        self.cline.genafpair = True   # E-INSi
        self.cline.localpair = False  # L-INSi
        self.cline.thread = 1
    
    def use_einsi(self):
        '''
        This is a method to toggle the use in favour of E-INSi which is equivalent
        to setting the --genafpair value for MAFFT's command line.
        '''
        self.cline.genafpair = True
        self.cline.localpair = False
    
    def use_linsi(self):
        '''
        This is a method to toggle the use in favour of L-INSi which is equivalent
        to setting the --localpair value for MAFFT's command line.
        '''
        self.cline.localpair = True
        self.cline.genafpair = False
    
    def set_threads(self, num):
        '''
        This method allows the use of multithreading by MAFFT. num should be a valid
        integer greater than 0. Be sensible since there's no upper limit checking.
        '''
        assert isinstance(num, int)
        if num < 0:
            raise Exception("Number of threads must be more than 0")
        
        self.cline.thread = num
        
    def run(self, FASTA_obj):
        '''
        Handles the execution of MAFFT alignment from start to finish. It must:
        
            1) Create a temporary file for MAFFT to read from the .FASTA object
            2) Perform the alignment with MAFFT
            3) Parse the output file and modify the .FASTA object appropriately
            4) Clean up temporary and output files
        
        Before calling this, you might consider the other methods of this class e.g.,
        set_threads() and use_linsi() if you don't want the default E-INSi behaviour.
        
        Params:
            FASTA_obj -- an object of ZS_SeqIO.FASTA class.
        '''
        # Validate input value type
        assert type(FASTA_obj).__name__ == "FASTA" or type(FASTA_obj).__name__ == "ZS_SeqIO.FASTA"
        #assert isinstance(FASTA_obj, FASTA)
        
        # Create temporary file
        tmpHash = hashlib.sha256(bytes(str(FASTA_obj.fileOrder[0][0]) + str(time.time()) + str(random.randint(0, 100000)), 'utf-8') ).hexdigest()
        tmpFileName = self._tmp_file_name_gen("mafft_tmp" + tmpHash[0:20], "fasta")
        FASTA_obj.write(tmpFileName)
        
        # Run
        self.cline.input = tmpFileName
        stdout, stderr = self.cline()
        if stdout == '':
            raise Exception("MAFFT error text below" + str(stderr))
        
        # Clean MAFFT output
        stdout = stdout.split("\n")
        while stdout[-1] == "\n" or stdout[-1] == "" or stdout[-1].startswith("Terminate batch job"):   # Remove junk
            del stdout[-1]
        
        # Parse output back into the FASTA object
        thisSeq = None
        ongoingCount = 0
        for line in stdout:
            if line.startswith(">"):
                if thisSeq == None:
                    thisSeq = []
                else:
                    FASTA_obj[ongoingCount].gap_seq = "".join(thisSeq)
                    thisSeq = []
                    ongoingCount += 1
            else:
                thisSeq.append(line)
        FASTA_obj[ongoingCount].gap_seq = "".join(thisSeq) # handle last sequence in the iteration
        
        # Set flag that this FASTA object has been aligned
        FASTA_obj.isAligned = True
        
        # Clean up temporary file
        os.unlink(tmpFileName)
    
    def run_nucleotide_as_protein(self, FASTA_obj, strand=None, frame=None):
        '''
        Refer to method header of run() for the fundamentals of this method.
        Refer to FastASeq method header of get_translation() for descriptions
        of the other parameters.
        
        What this method does differently to run() is that it will:
        
            1) Translate nucleotides into their corresponding protein sequence
            2) Align the protein sequences
            3) Map the gaps in the protein alignment back to the nucleotide sequence
        
        In effect, it allows nucleotide sequences to be aligned as codons rather than
        as individual base pairs, which might be more biologically relevant when dealing
        with closely related sequences.
        
        Note that findBestFrame is an implicit parameter. If strand is None or the indexed list
        value is None, and/or the frame value is similar, then findBestFrame will be set to True
        since we'll need to derive the appropriate strand and/or frame depending on which is None.
        If you set both strand and frame, or the indexed list value is set, then we obviously
        don't need to find how to translate it now do we?
        
        Parameters:
            findBestFrame -- a boolean indicating whether we should get the FastASeq.get_translation()
                             function figure out the best translation.
            strand -- one of three types; 1) an integer of 1 to get a +ve strand translation, or -1 for 
                      a -ve strand translation, 2) None to not constrain strandedness of the translation
                      if findBestFrame is True, or 3) a list of integers and/or None.
            frame -- one of three types; 1) an integer in the range(0,3) corresponding to 0-based
                     frame number, 2) None None to not constrain the frame of the translation
                     if findBestFrame is True, or 3) a list of integers and/or None.
        '''
        # Validate input FASTA parameter
        assert type(FASTA_obj).__name__ == "FASTA" or type(FASTA_obj).__name__ == "ZS_SeqIO.FASTA"
        
        # Validate strand parameter
        assert isinstance(strand, int) or strand == None or isinstance(strand, list)
        if isinstance(strand, int):
            assert strand in [1, -1], "Strand was provided but was not an appropriate value"
        if isinstance(strand, list):
            for _strand in strand:
                if isinstance(_strand, int):
                    assert _strand in [1, -1], "Strand was provided as a list, but contains inappropriate value(s)"
                else:
                    assert _strand == None, "Strand was provided as a list, but contains unknown types"
        
        # Validate frame parameter
        assert isinstance(frame, int) or frame == None or isinstance(frame, list)
        if isinstance(frame, int):
            assert frame in range(0, 3), "Frame was provided but was not an appropriate value"
        if isinstance(frame, list):
            for _frame in frame:
                if isinstance(_frame, int):
                    assert _frame in range(0,3), "Frame was provided as a list, but contains inappropriate value(s)"
                else:
                    assert _frame == None, "Frame was provided as a list, but contains unknown types"
        
        # Create a dummy FASTA object with our translations
        dummy = deepcopy(FASTA_obj)
        strands, frames = [], []
        for i in range(len(dummy)):
            # Get the appropriate strand/frame values
            if isinstance(strand, list):
                _strand = strand[i]
            else:
                _strand = strand
            if isinstance(frame, list):
                _frame = frame[i]
            else:
                _frame = frame
            
            # Run the translation
            if _strand != None and _frame != None:
                protein, _strand, _frame = dummy.seqs[i].get_translation(False, _strand, _frame)
            else:
                protein, _strand, _frame = dummy.seqs[i].get_translation(True, _strand, _frame)
            # Update and store results
            strands.append(_strand)
            frames.append(_frame)
            dummy.seqs[i].seq = protein
        
        # Align it
        self.run(dummy)
        
        # Map back aligned proteins to their original nucleotide counterparts
        longestTrimAmount = 0
        trimAmounts = []
        for i in range(len(dummy)):
            # if dummy[i].id == "Phascogale_calura_AM_M44886":
            #     stophere
            strand = strands[i]
            frame = frames[i]
            nuc = FASTA_obj[i].seq if strand == 1 else FASTA_obj[i].get_reverse_complement()
            alignedProt = dummy[i].gap_seq
            
            # Perform mapping procedure
            trimAmount = 0 if frame == 0 else 1 if frame == 2 else 2
            trimAmounts.append(trimAmount)
            longestTrimAmount = max(trimAmount, longestTrimAmount)
            alignedNuc = nuc[0: trimAmount] # put any extra bits before frame start in now
            codonIndex = trimAmount
            for proteinIndex in range(0, len(alignedProt)):
                if alignedProt[proteinIndex] != "-":
                    alignedNuc += nuc[codonIndex: codonIndex+3] # codon length == 3
                    codonIndex += 3 # iterate our nucleotide position marker
                else:
                    alignedNuc += "---"
            alignedNuc += nuc[codonIndex:] # add any potential stop codons back
            
            # Update original FASTA values
            FASTA_obj[i].gap_seq = alignedNuc
        
        # Add in any gap sequence necessary to accommodate trimmed bits
        for i in range(len(dummy)):
            trimAmount = trimAmounts[i]
            FASTA_obj[i].gap_seq = "-"*(longestTrimAmount-trimAmount) + FASTA_obj[i].gap_seq
        
        # Briefly fix up any weirdness
        "I think my codon system kinda fucks up the last codon by missing 1-2 gaps at times"
        maxLen = max([len(FastASeq_obj.gap_seq) for FastASeq_obj in FASTA_obj])
        for FastASeq_obj in FASTA_obj:
            if len(FastASeq_obj.gap_seq) != maxLen:
                FastASeq_obj.gap_seq += "-"*(maxLen - len(FastASeq_obj.gap_seq))
        
        # Set flag that this FASTA object has been aligned
        FASTA_obj.isAligned = True
    
    def add(self, aligned_FASTA_obj, add_FASTA_obj):
        '''
        Handles the execution of MAFFT alignment, running in the mode where new sequences
        are added into an existing alignment. This method otherwise behaves similarly to
        the .run() method, but it:
        
            1) Creates TWO temporarys file for MAFFT to read from the TWO FASTA objects provided
            2) Performs the alignment with MAFFT
            3) Parses the output file and without modifying input FASTA objects
            4) Cleans up temporary and output files
        
        Before calling this, you might consider the other methods of this class e.g.,
        set_threads() and use_linsi() if you don't want the default E-INSi behaviour.
        
        Params:
            aligned_FASTA_obj -- an object of ZS_SeqIO.FASTA class which has already been
                                 aligned and is to have new sequences added into it.
            add_FASTA_obj -- an object of ZS_SeqIO.FASTA class which has not been aligned
                               and will be added into the existing alignment
        Returns:
            result_FASTA_obj -- an object of ZS_SeqIO.FASTA class which results from MAFFT
                                adding the add_FASTA_obj into aligned_FASTA_obj
        '''
        # Validate input value type
        # assert isinstance(aligned_FASTA_obj, FASTA) ## This fails -- annoying!
        assert type(aligned_FASTA_obj).__name__ == "FASTA" or type(aligned_FASTA_obj).__name__ == "ZS_SeqIO.FASTA"
        assert aligned_FASTA_obj.isAligned, "aligned_FASTA_obj must be aligned first!"
        # assert isinstance(add_FASTA_obj, FASTA)
        assert type(add_FASTA_obj).__name__ == "FASTA" or type(add_FASTA_obj).__name__ == "ZS_SeqIO.FASTA"
        
        # Create temporary files
        tmpHash = hashlib.sha256(bytes(str(aligned_FASTA_obj.fileOrder[0][0]) + str(time.time()) + str(random.randint(0, 100000)), 'utf-8') ).hexdigest()
        
        tmpAlignedFileName = self._tmp_file_name_gen("mafft_aligned_tmp" + tmpHash[0:20], "fasta")
        aligned_FASTA_obj.write(tmpAlignedFileName, asAligned=True, withDescription=True)
        
        tmpAddedFileName = self._tmp_file_name_gen("mafft_added_tmp" + tmpHash[0:20], "fasta")
        add_FASTA_obj.write(tmpAddedFileName, withDescription=True)
        
        tmpOutputFileName = self._tmp_file_name_gen("mafft_output_tmp" + tmpHash[0:20], "fasta")
        
        # Format command for running
        cmd = "{0} --thread {1} --add \"{2}\" \"{3}\" > \"{4}\"".format(
            os.path.join(self.mafftDir, "mafft"), self.cline.thread,
            tmpAddedFileName, tmpAlignedFileName, tmpOutputFileName
        )
        
        # Run
        run_mafft = subprocess.Popen(cmd, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
        stdout, stderr = run_mafft.communicate()
        if stdout == '':
            raise Exception("MAFFT error text below" + str(stderr))
        
        # Parse output file into new FASTA object
        result_FASTA_obj = FASTA(tmpOutputFileName, isAligned=True)
        
        # Clean up temporary files
        os.unlink(tmpAlignedFileName)
        os.unlink(tmpAddedFileName)
        os.unlink(tmpOutputFileName)
        
        # Return new result
        return result_FASTA_obj
    
    def _tmp_file_name_gen(self, prefix, suffix):
        '''
        Hidden function for use by Class methods.
        Params:
            prefix -- a string for a file prefix e.g., "tmp"
            suffix -- a string for a file suffix e.g., "fasta". Note that we don't
                      use a "." in this, since it's inserted between prefix and suffix
                      automatically.
        Returns:
            tmpName -- a string for a file name which does not exist in the current dir.
        '''
        ongoingCount = 1
        while True:
            if not os.path.isfile("{0}.{1}".format(prefix, suffix)):
                return "{0}.{1}".format(prefix, suffix)
            elif os.path.isfile("{0}.{1}.{2}".format(prefix, ongoingCount, suffix)):
                ongoingCount += 1
            else:
                return "{0}.{1}.{2}".format(prefix, ongoingCount, suffix)

if __name__ == "__main__":
    pass

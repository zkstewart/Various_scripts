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

import os, platform
from Bio.Align.Applications import MafftCommandline

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
        assert type(mafftDir).__name__ == "str"
        assert os.path.isdir(mafftDir)
        if platform.system() == "Windows":
            if not os.path.isfile(os.path.join(mafftDir, "mafft.bat")):
                raise Exception("{0} does not exist".format(os.path.join(mafftDir, "mafft.bat")))
        else:
            if not os.path.isfile(os.path.join(mafftDir, "mafft")) or not os.path.isfile(os.path.join(mafftDir, "mafft.exe")):
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
        assert type(num).__name__ == "int"
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
        assert type(FASTA_obj).__name__ == "FASTA"
        
        # Create temporary file
        tmpFileName = self._tmp_file_name_gen("mafft_tmp", "fasta")
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
        
        # Clean up temporary file
        os.unlink(tmpFileName)
    
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

class MSA:
    '''
    Relevant attributes include:
        TBD
    '''
    def __init__(self):
        raise NotImplementedError()

if __name__ == "__main__":
    pass

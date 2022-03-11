#! python3
# ZS_HmmIO.py
# Contains various Classes to perform manipulations involving
# MSAs and HMMs

# TBD
## HMMER class that can
## 1) search a HMM against a FASTA and retrieve results [-]
## 2) ???

import os, subprocess
from ZS_SeqIO import FASTA

class HMM:
    '''
    This Class provides methods for creating a HMM file from FASTA files or ZS_SeqIO.FASTA
    objects. Alternatively, you can point it to the location of an existing HMM file.
    The HMMER Class then makes use of HMM instances to perform various HMMER operations.
    
    If you build a HMM using this class from a FASTA file, it MUST already be aligned. This
    Class is dumb and assumes any FASTA file is aligned. If you provide it a ZS_SeqIO.FASTA
    instance, it will check to see if the .isAligned flag was set and raise errors if not.
    '''
    def __init__(self, hmmerDir):
        # Validate input type and location
        assert isinstance(hmmerDir, str)
        assert os.path.isdir(hmmerDir)
        
        # Validate that all necessary HMMER executables can be located
        for exe in ["hmmpress", "hmmbuild", "hmmsearch"]:
            if not os.path.isfile(os.path.join(hmmerDir, exe)) and not os.path.isfile(os.path.join(hmmerDir, exe + ".exe")):
                raise Exception("{0} does not exist at {1}".format(exe, hmmerDir))
        
        # Set attributes
        self.hmmerDir = hmmerDir
        self.FASTA = None # Stores a string or ZS_SeqIO.FASTA object
        self.FASTA_is_file = False # Flag so we know if .FASTA is a string file location
        self.FASTA_is_obj = False # Flag so we know if .FASTA is a ZS_SeqIO.FASTA object
        self.useAlts = False # Flag relevant when FASTA_is_obj for temp file creation
        self.hmmFile = None
    
    def load_FASTA_from_file(self, fastaFile):
        '''
        This method receives a FASTA file name and, if it's locateable, will store
        the location of our FASTA file for future use.
        
        Params:
            fastaFile -- a string indicating the location of a FASTA file relative
                         to the current working dir (or just give the full path).
        '''
        # Validate input type and location
        assert isinstance(fastaFile, str)
        if not os.path.isfile(fastaFile):
            raise Exception("{0} is not a file or does not exist".format(fastaFile))
        
        self.FASTA = fastaFile
        self.FASTA_is_file = True
    
    def load_FASTA_from_object(self, FASTA_obj, useAlts=False):
        '''
        This method receives a ZS_SeqIO.FASTA object for storage prior to HMM conversion.
        
        Params:
            FASTA_obj -- a ZS_SeqIO.FASTA instance of an aligned file.
            useAlts -- a boolean to indicate whether you want the FASTA object to rely upon
                       its alt IDs or use the default IDs.
        '''
        # Validate input type and location
        assert isinstance(FASTA_obj, FASTA)
        if not FASTA_obj.isAligned:
            raise Exception("FASTA object .isAligned flag is not set")
        
        self.FASTA = FASTA_obj
        self.FASTA_is_object = True
        self.useAlts = useAlts
    
    def convert_to_HMM(self, hmmName):
        '''
        Once the .FASTA attribute is set via loading a FASTA, this method will handle the
        conversion of the FASTA into a HMM file using hmmbuild and hmmpress.
        
        Params:
            hmmName -- a string providing the file name for our created HMM. This can include
                       the path to where you want the file to be written
        '''
        # Validate input type and location
        assert isinstance(hmmName, str)
        if os.path.isfile(hmmName):
            raise Exception("""{0} already exists so it won't be overwritten. Maybe try using
                            .load_HMM_file() if you want to make use of an already existing 
                            HMM""".format(hmmName))
            
        if os.path.dirname(hmmName) != "":
            if not os.path.isdir(os.path.dirname(hmmName)):
                raise Exception("""{0} does not exist, and therefore we won't write a file into
                                a non-existing directory. Maybe try creating this directory first
                                or specifying a location that already exists
                                """.format(os.path.dirname(hmmName)))
        
        # Create a temporary file if .FASTA is an object
        fileName = self.FASTA # If not self.FASTA_is_obj, then this remains our default value
        if self.FASTA_is_obj:
            fileName = self._tmp_file_name_gen("hmmbuild_tmp", "fasta") # Overwrite the default fileName here
            self.FASTA.write(fileName, withAlt = self.useAlts, asAligned = True) # Always write asAligned since it should be an MSA
        
        # Run hmmbuild & hmmpress
        self.hmmbuild(fileName, hmmName)
        self.hmmpress(hmmName)
        self.hmmFile = hmmName # Sets our instance attribute so we know the HMM file exists
        
        # Clean up temporary file if relevant
        if self.FASTA_is_obj:
            os.unlink(fileName)
    
    def load_HMM_file(self, hmmName):
        '''
        This method allows us to bypass the HMM creation stage if it's already been done before.
        In these cases, you'll want to run this method first before using a HMM instance by the
        HMMER Class.
        
        Params:
            hmmName -- a string providing the file name for our created HMM. This can include
                       the path to where you want the file to be written
        '''
        # Validate input type and location
        assert isinstance(hmmName, str)
        if not os.path.isfile(hmmName):
            raise Exception("{0} is not a file or does not exist".format(hmmName))
        
        self.hmmFile = hmmName # Sets our instance attribute so we know the HMM file exists
    
    def hmmbuild(self, fastaFile, outputFileName):
        '''
        This method isn't intended to be called directly by users, but is written static-like in case
        you want access to this method without going through the fuss of working with this Class "properly".
        By writing it like this, you at least need to specify the location of the HMMER executables first.
        
        Params:
            fastaFile -- a string indicating the location of a FASTA file to be hmmbuild-ed.
            outputFileName -- a string indicating the name and, optionally, the path of the HMM file to be created.
        '''
        # Validate input type and location
        assert isinstance(fastaFile, str)
        if not os.path.isfile(fastaFile):
            raise Exception("{0} is not a file or does not exist".format(fastaFile))
        
        assert isinstance(outputFileName, str)
        if os.path.isfile(outputFileName):
            raise Exception("""{0} already exists so it won't be overwritten. Maybe try using
                            .load_HMM_file() if you want to make use of an already existing 
                            HMM""".format(outputFileName))
        
        cmd = "{0} \"{1}\" \"{2}\"".format(os.path.join(self.hmmerDir, "hmmbuild"),  outputFileName, fastaFile)
        run_hmmbuild = subprocess.Popen(cmd, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
        _, hmmerr = run_hmmbuild.communicate()
        if hmmerr.decode("utf-8") != '':
            raise Exception("""hmmbuild error text below\n{0}\nProgram crashed when processing {1}
                            """.format(str(hmmerr.decode("utf-8")), outputFileName))
    
    def hmmpress(self, hmmFile):
        '''
        This method isn't intended to be called directly by users, but is written static-like in case
        you want access to this method without going through the fuss of working with this Class "properly".
        By writing it like this, you at least need to specify the location of the HMMER executables first.
        
        Params:
            hmmFile -- a string indicating the location of a HMM file to be hmmpress-ed.
        '''
        # Validate input type and location
        assert isinstance(hmmFile, str)
        if not os.path.isfile(hmmFile):
            raise Exception("{0} is not a file or does not exist".format(hmmFile))
        
        cmd = "{0} -f \"{1}\"".format(os.path.join(self.hmmerDir, "hmmpress"),  hmmFile)
        run_hmmpress = subprocess.Popen(cmd, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
        _, hmmerr = run_hmmpress.communicate()
        if hmmerr.decode("utf-8") != '':
            raise Exception("""hmmpress error text below\n{0}\nProgram crashed when processing {1}
                            """.format(str(hmmerr.decode("utf-8")), hmmFile))
    
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

#! python3
# ZS_HmmIO.py
# Contains various Classes to perform manipulations involving
# HMMs and FASTA files using HMMER

import os, subprocess, inspect, sys
sys.path.append(os.path.dirname(__file__))

from ZS_SeqIO import FASTA
from domtblout_handling import hmmer_parse, nhmmer_parse # Make these available to things loading HmmIO

class HMM:
    '''
    This Class provides methods for creating a HMM file from FASTA files or ZS_SeqIO.FASTA
    objects. Alternatively, you can point it to the location of an existing HMM file.
    The HMMER Class then makes use of HMM instances to perform various HMMER operations.
    
    If you build a HMM using this class from a FASTA file, it MUST already be aligned. This
    Class is dumb and assumes any FASTA file is aligned. If you provide it a ZS_SeqIO.FASTA
    instance, it will check to see if the .isAligned flag was set and raise errors if not.
    
    Pay attention to .isNucleotide! It defaults to False since most HMMs are proteins, but
    occasionally you'll have a nucleotide HMM and this field is important when used alongside
    the HMMER Class!
    '''
    def __init__(self, hmmerDir, isNucleotide=False):
        # Validate input type and location
        assert isinstance(hmmerDir, str)
        assert os.path.isdir(hmmerDir)
        
        # Validate that all necessary HMMER executables can be located
        for exe in ["hmmpress", "hmmbuild", "hmmsearch"]:
            if not os.path.isfile(os.path.join(hmmerDir, exe)) and not os.path.isfile(os.path.join(hmmerDir, exe + ".exe")):
                raise Exception("{0} does not exist at {1}".format(exe, hmmerDir))
        
        # Set attributes
        self.hmmerDir = hmmerDir
        self.isNucleotide = isNucleotide # Flag to specify whether the HMM is based on nucleotide or protein sequence
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
        self.FASTA_is_obj = True
        self.useAlts = useAlts
    
    def create_HMM(self, hmmName, hmmBuildExtraArgs=""):
        '''
        Once the .FASTA attribute is set via loading a FASTA, this method will handle the
        creation of a HMM file using hmmbuild and hmmpress.
        
        Params:
            hmmName -- a string providing the file name for our created HMM. This can include
                       the path to where you want the file to be written
            hmmBuildExtraArgs -- a string that optionally allows you to add extra arguments 
                                 to the command e.g., "--dna" may be needed for hmmbuild to 
                                 work successfully.
        '''
        # Validate input type and location
        assert isinstance(hmmBuildExtraArgs, str)
        assert isinstance(hmmName, str)
        if os.path.isfile(hmmName):
            raise Exception(inspect.cleandoc("""
                            {0} already exists so it won't be overwritten. Maybe try using
                            .load_HMM_file() if you want to make use of an already existing
                            HMM""".format(hmmName)))
            
        if os.path.dirname(hmmName) != "":
            if not os.path.isdir(os.path.dirname(hmmName)):
                raise Exception(inspect.cleandoc("""
                                {0} does not exist, and therefore we won't write a file into
                                a non-existing directory. Maybe try creating this directory first
                                or specifying a location that already exists
                                """.format(os.path.dirname(hmmName))))
        
        # Create a temporary file if .FASTA is an object
        fileName = self.FASTA # If not self.FASTA_is_obj, then this remains our default value
        if self.FASTA_is_obj:
            fileName = self._tmp_file_name_gen("hmmbuild_tmp", "fasta") # Overwrite the default fileName here
            self.FASTA.write(fileName, withAlt = self.useAlts, asAligned = True) # Always write asAligned since it should be an MSA
        
        # Run hmmbuild & hmmpress
        self.hmmbuild(fileName, hmmName, extraArgs=hmmBuildExtraArgs)
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
        
        # Check if we need to run hmmpress
        if not os.path.isfile(hmmName + '.h3f') and not os.path.isfile(hmmName + '.h3i') and not os.path.isfile(hmmName + '.h3m') and not os.path.isfile(hmmName + '.h3p'):
            self.hmmpress(hmmName)
        
        # Store prepared file
        self.hmmFile = hmmName # Sets our instance attribute so we know the HMM file exists
    
    def hmmbuild(self, fastaFile, outputFileName, extraArgs=""):
        '''
        This method isn't intended to be called directly by users, but is written static-like in case
        you want access to this method without going through the fuss of working with this Class "properly".
        By writing it like this, you at least need to specify the location of the HMMER executables first.
        
        Params:
            fastaFile -- a string indicating the location of a FASTA file to be hmmbuild-ed.
            outputFileName -- a string indicating the name and, optionally, the path of the HMM file to be created.
            extraArgs -- a string that optionally allows you to add extra arguments to the command
                         e.g., "--dna" may be needed for hmmbuild to work successfully.
        '''
        # Validate input type and location
        assert isinstance(extraArgs, str)
        assert isinstance(fastaFile, str)
        if not os.path.isfile(fastaFile):
            raise Exception("{0} is not a file or does not exist".format(fastaFile))
        
        assert isinstance(outputFileName, str)
        if os.path.isfile(outputFileName):
            raise Exception(inspect.cleandoc("""
                            {0} already exists so it won't be overwritten. Maybe try using
                            .load_HMM_file() if you want to make use of an already existing 
                            HMM""".format(outputFileName)))
        
        cmd = "{0} {3} \"{1}\" \"{2}\"".format(os.path.join(self.hmmerDir, "hmmbuild"),  outputFileName, fastaFile, extraArgs)
        run_hmmbuild = subprocess.Popen(cmd, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
        _, hmmerr = run_hmmbuild.communicate()
        if hmmerr.decode("utf-8") != '':
            raise Exception(inspect.cleandoc("""
                            hmmbuild error text below\n{0}\nProgram crashed when processing {1}
                            """.format(str(hmmerr.decode("utf-8")), outputFileName)))
    
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
            raise Exception(inspect.cleandoc("""
                            hmmpress error text below\n{0}\nProgram crashed when processing {1}
                            """.format(str(hmmerr.decode("utf-8")), hmmFile)))
    
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

class HMMER:
    '''
    This Class provides methods that make use of a ZS_HmmIO.HMM instance to run
    various HMMER modules.
    
    An important but understated part of this Class is its reliance on the
    HMM_obj's .isNucleotide attribute. This controls whether we run hmmsearch
    or nhmmer, which has minor differences in how we parse the output. It's 
    important you get this right, though.
    
    Params:
        HMM_obj -- a ZS_HmmIO.HMM instance which has been fully setup to point
                   to a prepared HMM file.
    '''
    def __init__(self, HMM_obj):
        # Validate input type and set-up
        assert isinstance(HMM_obj, HMM)
        if HMM_obj.hmmFile == None:
            raise Exception("HMM object not set up properly; lacks .hmmFile value")
        
        # Set attributes
        self.HMM = HMM_obj
        self.threads = 1 # Default to single-threaded
        self.Evalue = 0.01 # Default E-value of hmmsearch
        
        self.outputFileName = None # Default to automatically generate a suitable name
        self.domDict = None # This will be set to a dictionary containing parsed results
        
        self.FASTA = None # Stores a string or ZS_SeqIO.FASTA object
        self.FASTA_is_file = False # Flag so we know if .FASTA is a string file location
        self.FASTA_is_obj = False # Flag so we know if .FASTA is a ZS_SeqIO.FASTA object
        self.useAlts = False # Flag relevant when FASTA_is_obj for temp file creation
    
    def set_threads(self, num):
        '''
        This method allows the use of multithreading by HMMER. num should be a valid
        integer greater than 0. Be sensible since there's no upper limit checking.
        '''
        assert isinstance(num, int)
        if num < 0:
            raise Exception("Number of threads must be more than 0")
        
        self.threads = num
    
    def set_Evalue(self, num):
        '''
        This method allows the significance threshold to be altered. num should be a valid
        float or integer greater than 0. Be sensible since there's no upper limit checking.
        '''
        assert isinstance(num, float) or isinstance(num, int)
        if num < 0:
            raise Exception("Evalue must be more than 0")
        
        self.Evalue = num
    
    def set_output_name(self, name):
        '''
        This method allows the output file name to be explicitly specified rather than
        generated based on the target and HMM files.
        '''
        assert isinstance(name, str)
        if name == "":
            raise Exception("Name must not be empty")
        if os.path.isfile(name):
            raise Exception("{0} already exists; we will not allow HMMER to overwrite it".format(name))
        if os.path.dirname(name) != "":
            if not os.path.isdir(os.path.dirname(name)):
                raise Exception(inspect.cleandoc("""
                                {0} does not exist, and therefore we won't write a file into
                                a non-existing directory. Maybe try creating this directory first
                                or specifying a location that already exists
                                """.format(os.path.dirname(name))))
        
        self.outputFileName = name
    
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
        This method receives a ZS_SeqIO.FASTA object for storage prior to HMMER execution.
        
        Params:
            FASTA_obj -- a ZS_SeqIO.FASTA instance of an aligned file.
            useAlts -- a boolean to indicate whether you want the FASTA object to rely upon
                       its alt IDs or use the default IDs.
        '''
        # Validate input type and location
        assert isinstance(FASTA_obj, FASTA)
        
        self.FASTA = FASTA_obj
        self.FASTA_is_obj = True
        self.useAlts = useAlts
    
    def run_search(self):
        '''
        Once the .FASTA attribute is set via loading a FASTA, this method will handle
        the execution of HMMER's hmmsearch using this instance's HMM as a database and
        the FASTA as the target.
        
        Params:
            hmmName -- a string providing the file name for our created HMM. This can include
                       the path to where you want the file to be written
        '''
        # Validate possibility of running this function
        if self.FASTA == None:
            raise Exception("Can't run hmmsearch; HMMER instance not configured with a FASTA target.")
        if self.outputFileName == None:
            raise Exception("Can't run hmmsearch; HMMER instance not configured with an output file name.")
        
        # Create a temporary file if .FASTA is an object
        fileName = self.FASTA # If not self.FASTA_is_obj, then this remains our default value
        if self.FASTA_is_obj:
            fileName = self._tmp_file_name_gen("hmmsearch_tmp", "fasta") # Overwrite the default fileName here
            self.FASTA.write(fileName, withAlt = self.useAlts, asAligned = False) # As a target file we don't need gaps
        
        # Run hmmsearch
        self.hmmsearch(threads=self.threads, evalue=self.Evalue, outputFileName=self.outputFileName, hmmFile=self.HMM.hmmFile, fastaFile=fileName, isNucleotide=self.HMM.isNucleotide)
        
        # Parse search results
        self.domDict = nhmmer_parse(self.outputFileName, self.Evalue, extendedDetails=True) if self.HMM.isNucleotide else hmmer_parse(self.outputFileName, self.Evalue)
        
        # Clean up temporary file if relevant
        if self.FASTA_is_obj:
            os.unlink(fileName)
    
    def hmmsearch(self, threads=1, evalue=0.01, outputFileName="", hmmFile="", fastaFile="", isNucleotide=False):
        '''
        This method isn't intended to be called directly by users, but is written static-like in case
        you want access to this method without going through the fuss of working with this Class "properly".
        It lets you run hmmsearch without any post-processing or parsing.
        
        Params:
            hmmFile -- a string indicating the location of a HMM file to be hmmpress-ed.
            threads -- an integer controlling the number of HMMER threads to specify.
            evalue -- a float indicating the significance threshold for results to be returned.
            outputFileName -- a string value for a file to write to; should not already exist!
            hmmFile -- a string value for a HMM file.
            fastaFile -- a string value for a target FASTA file to search.
            isNucleotide -- a boolean indicating whether the we're using a nucleotide model or not.
        '''
        # Validate input types and locations
        assert isinstance(threads, int)
        assert isinstance(evalue, int) or isinstance(evalue, float)
        
        assert isinstance(outputFileName, str), "Try .set_output_name() first"
        if outputFileName == "":
            raise Exception("Name must not be empty")
        if os.path.isfile(outputFileName):
            raise Exception("{0} already exists; we will not allow HMMER to overwrite it".format(outputFileName))
        if os.path.dirname(outputFileName) != "":
            if not os.path.isdir(os.path.dirname(outputFileName)):
                raise Exception(inspect.cleandoc("""
                                {0} does not exist, and therefore we won't write a file into
                                a non-existing directory. Maybe try creating this directory first
                                or specifying a location that already exists
                                """.format(os.path.dirname(outputFileName))))
        
        if not os.path.isfile(hmmFile):
            raise Exception("{0} is not a file or does not exist".format(hmmFile))
        if not os.path.isfile(fastaFile):
            raise Exception("{0} is not a file or does not exist".format(fastaFile))
        
        # Perform the hmmsearch operation
        program = "nhmmer" if isNucleotide else "hmmsearch"
        outputArg = "tblout" if isNucleotide else "domtblout"
        cmd = "{0} --cpu {1} -E {2} --{3} {4} \"{5}\" \"{6}\"".format(
            os.path.join(self.HMM.hmmerDir, program), threads, evalue, outputArg, outputFileName, hmmFile, fastaFile
        )
        run_hmmsearch = subprocess.Popen(cmd, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
        _, hmmerr = run_hmmsearch.communicate()
        if hmmerr.decode("utf-8") != '':
            raise Exception(inspect.cleandoc("""
                            hmmsearch error text below\n{0}\nProgram crashed when processing {1}
                            """.format(str(hmmerr.decode("utf-8")), hmmFile)))
    
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

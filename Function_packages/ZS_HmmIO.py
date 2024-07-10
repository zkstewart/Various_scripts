#! python3
# ZS_HmmIO.py
# Contains various Classes to perform manipulations involving
# HMMs and FASTA files using HMMER

import os, subprocess, inspect, sys, hashlib, time, random, platform

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import ZS_Utility, ZS_SeqIO
from domtblout_handling import hmmer_parse, nhmmer_parse # Make these available to things loading HmmIO

class HMM:
    '''
    This Class provides methods for creating a HMM file from a FASTA file or ZS_SeqIO.FASTA
    object. It is intended for Linux and WSL use. The native Windows executables
    are not supported; you should compile through WSL and use the executables from there.
    
    Parameters:
        hmmerDir -- a string indicating the location of the HMMER executables including 'hmmpress'
                    and 'hmmbuild'.
    '''
    def __init__(self, hmmerDir, isNucleotide=False):
        self.hmmerDir = hmmerDir
        self.isNucleotide = isNucleotide
    
    @property
    def hmmerDir(self):
        return self._hmmerDir
    
    @hmmerDir.setter
    def hmmerDir(self, value):
        convertedValue = ZS_Utility.convert_to_wsl_if_not_unix(value)
        assert ZS_Utility.wsl_exists(convertedValue, isFolder=True), \
            f"hmmer folder not found at '{convertedValue}' after WSL compatibility conversion"
        self._hmmerDir = convertedValue
        self.hmmpress = value + "/hmmpress"
        self.hmmbuild = value + "/hmmbuild"
    
    @property
    def hmmpress(self):
        return self._hmmpress
    
    @hmmpress.setter
    def hmmpress(self, value):
        convertedValue = ZS_Utility.convert_to_wsl_if_not_unix(value)
        assert ZS_Utility.wsl_exists(convertedValue), \
            f"hmmpress executable not found at '{convertedValue}' after WSL compatibility conversion"
        self._hmmpress = convertedValue
    
    @property
    def hmmbuild(self):
        return self._hmmbuild
    
    @hmmbuild.setter
    def hmmbuild(self, value):
        convertedValue = ZS_Utility.convert_to_wsl_if_not_unix(value)
        assert ZS_Utility.wsl_exists(convertedValue), \
            f"hmmbuild executable not found at '{convertedValue}' after WSL compatibility conversion"
        self._hmmbuild = convertedValue
    
    def create(self, inputFasta, outputFileName, isNucleotide=False):
        '''
        This method will handle the creation of a HMM file using hmmbuild and hmmpress.
        
        Params:
            inputFasta -- a string indicating the location of a FASTA file OR a ZS_SeqIO.FASTA object
                          OR a string pointing to an existing HMM file.
            outputFileName -- a string providing the file name for our created HMM. This can include
                              the path to where you want the file to be written
            isNucleotide -- OPTIONAL; a boolean to indicate whether the provided FASTA contains nucleotide
                            sequences or not. Default == False i.e., FASTA contains protein sequences.
        '''
        # Validate input value type
        if isinstance(inputFasta, str):
            assert os.path.isfile(inputFasta), f"ERROR: HMM.create() could not find the FASTA file '{inputFasta}'"
            fastaFileName, fastaIsTemporary = inputFasta, False
        elif hasattr(inputFasta, "isFASTA") and inputFasta.isFASTA is True:
            fastaFileName, fastaIsTemporary = ZS_SeqIO.Conversion.get_filename_for_input_sequences(inputFasta)
        else:
            raise Exception(f"ERROR: HMM.create() requires a FASTA file or FASTA object as input; did not understand '{inputFasta}'")
        
        # Validate that output file does not already exist
        assert not os.path.exists(outputFileName), f"ERROR: HMM.create() will not overwrite existing file '{outputFileName}'"
        
        # Run hmmbuild & hmmpress
        HMM.hmmbuild(self.hmmbuild, fastaFileName, outputFileName, isNucleotide)
        HMM.hmmpress(self.hmmpress, outputFileName)
        
        # Clean up temporary file if relevant
        if fastaIsTemporary:
            os.unlink(fastaFileName)
    
    @staticmethod
    def hmmbuild(hmmbuildExe, fastaFileName, outputFileName, isNucleotide=False):
        '''
        Parameters:
            hmmbuildExe -- a string indicating the location of the hmmbuild executable.
            fastaFileName -- a string indicating the location of a FASTA file to be hmmbuild-ed.
            outputFileName -- a string indicating the name and, optionally, the path of the HMM file to be created.
            isNucleotide -- OPTIONAL; a boolean to indicate whether the provided FASTA contains nucleotide
                            sequences or not. Default == False i.e., FASTA contains protein sequences.
        '''
        assert os.path.isfile(fastaFileName), \
            f"ERROR: HMM.hmmbuild() could not find FASTA file '{fastaFileName}'"
        
        # Construct the cmd for subprocess
        cmd = ZS_Utility.base_subprocess_cmd(hmmbuildExe)
        
        # Handle nucleotide sequences
        if isNucleotide:
            cmd += ["--dna"]
        
        # Set input and output file arguments
        cmd += [
            ZS_Utility.convert_to_wsl_if_not_unix(outputFileName), # output comes first, it's weird
            ZS_Utility.convert_to_wsl_if_not_unix(fastaFileName)
        ]
        
        # Format commands for Linux
        if platform.system() != "Windows":
            cmd = " ".join(cmd)
        
        # Run the command
        run_build = subprocess.Popen(cmd, shell = True,
                                     stdout = subprocess.PIPE,
                                     stderr = subprocess.PIPE)
        buildout, builderr = run_build.communicate()
        
        # Check to see if there was an error
        if builderr.decode("utf-8") != "":
            raise Exception(("ERROR: HMM.hmmbuild() encountered an error; have a look " +
                            f'at the stdout ({buildout.decode("utf-8")}) and stderr ' + 
                            f'({builderr.decode("utf-8")}) to make sense of this.'))
    
    @staticmethod
    def hmmpress(hmmpressExe, hmmFileName):
        '''
        Parameters:
            hmmpressExe -- a string indicating the location of the hmmpress executable.
            hmmFileName -- a string indicating the location of a HMM file to be hmmpress-ed.
        '''
        assert os.path.isfile(hmmFileName), \
            f"ERROR: HMM.hmmpress() could not find HMM file '{hmmFileName}'"
        
        # Construct the cmd for subprocess
        cmd = ZS_Utility.base_subprocess_cmd(hmmpressExe)
        
        # Set input file argument
        cmd += ["-f", ZS_Utility.convert_to_wsl_if_not_unix(hmmFileName)] # -f flag will force overwrite
        
        # Format commands for Linux
        if platform.system() != "Windows":
            cmd = " ".join(cmd)
        
        # Run the command
        run_press = subprocess.Popen(cmd, shell = True,
                                     stdout = subprocess.PIPE,
                                     stderr = subprocess.PIPE)
        pressout, presserr = run_press.communicate()
        
        # Check to see if there was an error
        if presserr.decode("utf-8") != "":
            raise Exception(("ERROR: HMM.hmmpress() encountered an error; have a look " +
                            f'at the stdout ({pressout.decode("utf-8")}) and stderr ' + 
                            f'({presserr.decode("utf-8")}) to make sense of this.'))

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
            tmpHash = hashlib.sha256(bytes(str(self.FASTA.fileOrder[0][0]) + str(time.time()) + str(random.randint(0, 100000)), 'utf-8') ).hexdigest()
            fileName = ZS_Utility.tmp_file_name_gen("hmmsearch_tmp" + tmpHash[0:20], "fasta") # Overwrite the default fileName here
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

if __name__ == "__main__":
    pass

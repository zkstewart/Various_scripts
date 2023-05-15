#! python3
# ZS_SignalPIO.py
# Contains Class(es) to manipulate FastASeq and FASTA
# objects for the prediction of signal peptides.

import os, sys, subprocess, platform, shutil

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from ZS_SeqIO import Conversion
import ZS_Utility

class SignalP:
    '''
    The SignalP Class provides easy access to the signalP program to perform queries
    using ZS_FastASeq or ZS_FASTA objects. SignalP versions supported include
    (TBD...).
    
    If this Class is being used on Windows, it's assumed you have a functioning
    WSL setup; not meeting this expectation will result in errors. If you're running
    this on Linux, then no worries!
    
    Its primary purpose is the .get_sigp_results() method.
    
    Attributes:
        query (REQUIRED) -- a string, FASTA, or FastASeq object from the ZS_SeqIO suite.
        signalpExe (REQUIRED) -- a string indicating the location of the SignalP executable.
        clean (OPTIONAL) -- a Boolean indicating whether to clean up after BLAST search is
                            complete (i.e., delete the result file) or to keep it. Defaults
                            to True, which means we will DELETE the results file after it
                            has been parsed.
        _version (PRIVATE) -- an integer indicating which version of SignalP is being operated;
                             is set intrinsically by the signalpExe property.
    '''
    def __init__(self, query, signalpExe):
        # Validate input types
        assert type(query).__name__ == "str" \
            or type(query).__name__ == "FASTA" \
            or type(query).__name__ == "ZS_SeqIO.FASTA" \
            or type(query).__name__ == "FastASeq" \
            or type(query).__name__ == "ZS_SeqIO.FastASeq"
        
        assert type(signalpExe).__name__ == "str"
        
        # Validate that inputs exist if specified as a string (file location)
        if type(query).__name__ == "str" and not os.path.isfile(query):
            raise FileNotFoundError("Query parameter is a string, but does not point to an existing file location")
        self.query = query
        
        # Set attributes
        self.clean = True
        self.signalpExe = signalpExe # internally validates
        self.organism = "euk" # default for Eukaryote, versions 4/5/6 accept "euk" which is nice
        self.results = None # default to empty, will be set when running signalP
    
    @property
    def clean(self):
        return self._clean
    
    @clean.setter
    def clean(self, value):
        '''
        Setter for the clean attribute, which controls
        whether signalP output files are kept or not.
        
        Parameters:
            value -- a Boolean of True or False
        '''
        assert isinstance(value, bool)
        self._clean = value
    
    @property
    def signalpExe(self):
        return self._signalpExe
    
    @signalpExe.setter
    def signalpExe(self, value):
        '''
        Setter for the signalpExe attribute, which validates
        that the signalP executable file exists and is valid.
        
        Parameters:
            value -- a Boolean of True or False
        '''
        sigpIsValid, notValidMessage = SignalP.signalP_is_valid(value)
        if not sigpIsValid:
            raise Exception(f"SignalP executable failed to validate because {notValidMessage}")
        
        self._signalpExe = value
        
        # Intrinsic action: set _version attribute
        self.version = SignalP.get_signalP_version(value)
    
    @property
    def version(self):
        return self._version
    
    @version.setter
    def version(self, value):
        '''
        Setter for the version attribute, which controls
        signalP execution behaviour to work with versions
        4, 5, or 6 of signalP
        
        Parameters:
            value -- an integer in the list:
                     [4, 5, 6]
        '''
        VALID_VERSIONS = [4, 5, 6]
        assert value in VALID_VERSIONS
        self._version = value
    
    @property
    def organism(self):
        return self._organism
    
    @organism.setter
    def organism(self, value):
        '''
        Setter for the organism attribute, which controls
        which organism type signalP is expecting to handle.
        
        Parameters:
            value -- a string with supported options depending on SignalP
                     version i.e.,:
                        4 = ["euk", "gram+", "gram-"]
                        5 = ["arch", "euk", "gram+", "gram-"]
                        6 = ["eukarya", "other", "euk"]
        '''
        VALID_4 = ["euk", "gram+", "gram-"]
        VALID_5 = ["arch", "euk", "gram+", "gram-"]
        VALID_6 = ["eukarya", "other", "euk"]
        
        if self.version == 4:
            assert value.lower() in VALID_4, \
                f"'{value}' organism not supported; must be in list {VALID_4}"
        elif self.version == 5:
            assert value.lower() in VALID_5, \
                f"'{value}' organism not supported; must be in list {VALID_5}"
        elif self.version == 6:
            assert value.lower() in VALID_6, \
                f"'{value}' organism not supported; must be in list {VALID_6}"
        
        self._organism = value.lower()
    
    @property
    def results(self):
        return self._results
    
    @results.setter
    def results(self, value):
        '''
        Setter for the results attribute, which stores the most recent
        signalP run from this object. Will be None prior to run, after
        which it's expected to be set with a dictionary value.
        
        Parameters:
            value -- None, OR a dictionary with structure like:
                     {
                        'seqid1': [start, end],
                        'seqid2': [start, end],
                        ...
                    }
        '''
        assert value == None or isinstance(value, dict), \
            "SignalP results attribute can only be None or a dictionary"
        
        self._results = value
    
    @staticmethod
    def signalP_is_valid(signalpExe):
        '''
        Static helper to ascertain whether the signalP executable provided
        is compatible with the ZS_SignalPIO class.
        
        Parameters:
            signalpExe -- a string indicating the location of the executable itself
        Returns:
            sigpIsValid -- a boolean; True if signalP is validated, False otherwise
            errorMessage -- a string indicating what went wrong if validation failed
        '''
        EXPECTED_STRINGS = ["usage", "version"]
        
        # Check that executable file exists
        if not os.path.isfile(signalpExe):
            return False, f"'{signalpExe}' doesn't exist!"
        
        # Check that executable file returns default help message
        cmd = ZS_Utility.base_subprocess_cmd(signalpExe) + ["-h"]
        run_sigp = subprocess.Popen(" ".join(cmd), shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        sigpout, sigperr = run_sigp.communicate()
        out, err = sigpout.decode("utf-8").lower(), sigperr.decode("utf-8").lower()
        if not any([eStr in out or eStr in err for eStr in EXPECTED_STRINGS]):
            "Different SignalP versions output on stdout or stderr which is stupid"
            return False, sigperr.decode("utf-8")
        else:
            return True, ""
    
    @staticmethod
    def get_signalP_version(signalpExe):
        '''
        Static helper to check for the signalP executable version. The executable
        itself will be validated beforehand by the signalP_is_valid() static method.
        
        Parameters:
            signalpExe -- a string indicating the location of the executable itself
        Returns:
            version -- an integer of 4, 5, or 6 depending on the version; any other
                       version will result in an error since we do not handle them.
        '''
        # Pre-validate
        sigpIsValid, notValidMessage = SignalP.signalP_is_valid(signalpExe)
        if not sigpIsValid:
            raise Exception("signalP_is_valid indicates invalid signalP executable prior to signalP version checking!")
        
        # Check which instruction we need to obtain version
        cmd = ZS_Utility.base_subprocess_cmd(signalpExe) + ["-h"]
        run_sigp = subprocess.Popen(" ".join(cmd), shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        sigpout, sigperr = run_sigp.communicate()
        
        if "-version" in sigperr.decode("utf-8").lower():
            "Need to have it as a plain string to work??"
            cmd = ZS_Utility.base_subprocess_cmd(signalpExe) + ["--version"] # cmd for version 5 or 6
        else:
            cmd = ZS_Utility.base_subprocess_cmd(signalpExe) + ["-V"] # cmd for version 4
        
        # Get version stdout
        run_version = subprocess.Popen(" ".join(cmd), shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        versionout, versionerr = run_version.communicate()
        if versionerr.decode("utf-8") != "":
            raise Exception("get_signalP_version ran into unexpected failure with stderr of {0}".format(
                versionerr.decode("utf-8")
            ))
        versionMessage = versionout.decode("utf-8").lower()
        
        # Extract version from stdout
        if "signalp 4" in versionMessage:
            raise NotImplementedError("SignalP 4 refuses to work properly; hence, not implemented. Sorry.")
            return 4
        elif "signalp version 5" in versionMessage:
            return 5
        elif "signalP 6" in versionMessage:
            if platform.system() == "Windows":
                raise NotImplementedError("Unfortunately, I haven't implemented SignalP 6 on Windows yet")
            return 6
        else:
            raise ValueError("get_signalP_version could not understand or does not support {0}".format(
                versionMessage
            ))
    
    def signalp(self):
        '''
        Performs the signalP operation.
        
        This method does not use .query because it may be in the FASTA or
        FastASeq object type, rather than a string which we require here.
        
        KNOWN ERRORS:
            - Segmentation faults in WSL for signalP 4; I tried to make it compatible
              but it just won't work.
            - Doesn't handle signalP 6 on Windows; I just haven't installed and tested
              it is all. It takes up a lot of file storage to install this version...
        '''
        # Get query as a FASTA file
        fastaQuery, isTemporary = Conversion.get_filename_for_input_sequences(self.query)
        fastaQuery = os.path.abspath(fastaQuery) # need abspath for later
        
        # Get tmp dir (2 dirs back)
        if platform.system() == "Windows":
            tmpDir = ZS_Utility.convert_windows_to_wsl_path(
                os.path.join(os.path.dirname(os.path.dirname(self.signalpExe)), "tmp")
            )
        else:
            tmpDir = os.path.join(os.path.dirname(os.path.dirname(self.signalpExe)), "tmp")
        
        # Get prefix for output files
        prefixHash = Conversion.get_hash_for_input_sequences(self.query)
        
        # Format cmd depending on SignalP version
        if self.version == 4:
            # Get our tmp output file
            sigp4TmpFile = ZS_Utility.tmp_file_name_gen(f"{prefixHash}_tmp", "gff3")
            
            # Format the cmd now
            cmd = ZS_Utility.base_subprocess_cmd(self.signalpExe)
            cmd += [
                "-t", self.organism,
                "-f", "short",
                "-n", f'{sigp4TmpFile}',
                f'{ZS_Utility.convert_windows_to_wsl_path(fastaQuery)}' if platform.system() == "Windows"
                    else fastaQuery
            ]
        elif self.version == 5:
            # cd into the SignalP bin dir before doing anything else
            if platform.system() == "Windows":
                cmd = [
                    "wsl", "--cd",
                    os.path.dirname(ZS_Utility.convert_windows_to_wsl_path(self.signalpExe)),
                    "-e"
                ]
            else:
                cmd = [
                    "cd", os.path.dirname(self.signalpExe), ";"
                ]
            
            # Add signalP instructions
            cmd += [
                f"./{os.path.basename(self.signalpExe)}",
                "-org", self.organism,
                "-format", "short",
                "-gff3",
                "-prefix", prefixHash,
                "-tmp", tmpDir,
                "-fasta", f'{ZS_Utility.convert_windows_to_wsl_path(fastaQuery)}' 
                    if platform.system() == "Windows" else fastaQuery
            ]
        elif self.version == 6:
            # Make a tmp dir just for signalP 6
            sigp6TmpDir = os.path.join(os.getcwd(), f"{prefixHash}_tmp")
            
            # Format the cmd now
            cmd = ZS_Utility.base_subprocess_cmd(self.signalpExe)
            cmd += [
                "-org", self.organism,
                "-format", "txt",
                "--output_dir", f'{ZS_Utility.convert_windows_to_wsl_path(sigp6TmpDir)}' 
                    if platform.system() == "Windows" else sigp6TmpDir,
                "-fasta", f'{ZS_Utility.convert_windows_to_wsl_path(fastaQuery)}' 
                    if platform.system() == "Windows" else fastaQuery
            ]
        
        # Run signalP with the formatted cmd
        run_sigp = subprocess.Popen(" ".join(cmd), shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        sigpout, sigperr = run_sigp.communicate()
        if sigperr.decode("utf-8") != '':
            raise Exception('SignalP error text below\n' + sigperr.decode("utf-8"))
        
        # Parse output file
        if self.version == 4:
            raise NotImplementedError("Sigp4 not handled yet, sorry!")
            sigpGff3 = sigp4TmpFile
        elif self.version == 5:
            sigpGff3 = os.path.join(os.path.dirname(self.signalpExe), prefixHash + ".gff3")
        elif self.version == 6:
            sigpGff3 = os.path.join(sigp6TmpDir, "output.gff3")
        
        self.results = SignalP.parse_sigp_gff3(sigpGff3) # store for safe keeping
        
        # Clean up output and temporary files
        if isTemporary:
            os.unlink(fastaQuery)
        
        if self.version == 4 and self.clean == True:
            raise NotImplementedError("Sigp4 not handled yet, sorry!")
            os.unlink(sigp4TmpFile)
        elif self.version == 5 and self.clean == True:
            summaryFile = os.path.join(os.path.dirname(self.signalpExe), prefixHash + "_summary.signalp5")
            os.unlink(sigpGff3)
            os.unlink(summaryFile) # signalp5 also creates an extra summary file; we need to clean it, too
        elif self.version == 6 and self.clean == True:
            shutil.rmtree(sigp6TmpDir)
        
        return self.results # return for calling function
    
    @staticmethod
    def parse_sigp_gff3(sigpGff3):
        '''
        Parameters:
            sigpGff3 -- a string indicating the location of a GFF3 produced by
                        signalP
        Returns:
            sigpDict -- a dictionary with structure like:
                        {
                            'seqid1': [start, end],
                            'seqid2': [start, end],
                            ...
                        }
        '''
        sigpDict = {}
        with open(sigpGff3, "r") as fileIn:
            for line in fileIn:
                if line.startswith("#"):
                    continue
                else:
                    seqid, sigpVersion, predictionType, \
                        start, end, score, _, _, _ = line.rstrip("\r\n ").split("\t")
                    sigpDict[seqid] = [int(start), int(end)]
        return sigpDict
    
    def __str__(self):
        return "SignalP; query file ='{0}', running SignalP version {1} with '{2}' prediction, {3})".format(
            self.query, self.version, self.organism,
            "contains .results value" if self.results != None else "no .results value"
        )
    
    def __repr__(self):
        return "SignalP(query='{0}',sigpExe='{1}',version={2},organism={3},clean={4},results={5})".format(
            self.query, self.signalpExe, self.version, self.organism, self.clean,
            "<dict>" if self.results != None else "<None>"
        )

if __name__ == "__main__":
    pass

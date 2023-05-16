#! python3
# ZS_SignalPIO.py
# Contains Class(es) to manipulate FastASeq and FASTA
# objects for the prediction of signal peptides.

import os, sys, subprocess, platform, shutil
from pathlib import Path

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
        cygwinDir (OPTIONAL) -- a string indicating the location of the Cygwin bin directory;
                                only relevant if planning to run SignalP 4 on Windows.
        clean (OPTIONAL) -- a Boolean indicating whether to clean up after BLAST search is
                            complete (i.e., delete the result file) or to keep it. Defaults
                            to True, which means we will DELETE the results file after it
                            has been parsed.
        _version (PRIVATE) -- an integer indicating which version of SignalP is being operated;
                             is set intrinsically by the signalpExe property.
    '''
    def __init__(self, query, signalpExe, cygwinDir=None):
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
        self.cygwinDir = cygwinDir # defaults to None in init since only needed for Windows signalP 4
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
        sigpIsValid, notValidMessage = SignalP.signalP_is_valid(value, self.cygwinDir) # setting this comes AFTER cygwinDir
        if not sigpIsValid:
            raise Exception((f"SignalP executable failed to validate because {notValidMessage}\n\n" +
                             "If this above error message doesn't really help, a possible explanation " +
                             "is that you're trying to run SignalP 4 on Windows without having provided " +
                             "the cygwinDir value to the SignalP class on instantiation. Or maybe " +
                             "your SignalP 4 file has the wrong configuration. Hopefully that helps."))
        
        self._signalpExe = value
        
        # Intrinsic action: set _version attribute
        self.version = SignalP.get_signalP_version(value, self.cygwinDir)
    
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
    def cygwinDir(self):
        return self._cygwinDir
    
    @cygwinDir.setter
    def cygwinDir(self, value):
        '''
        Setter for the cygwinDir attribute, which stores
        the location of the Cygwin executable files
        (i.e., bin dir where bash.exe exists).
        
        Should be a valid string pointing to the Cygwin
        directory, or None.
        
        Parameters:
            value -- a string pointing to the Cygwin
                     bin directory
        '''
        if value != None:
            assert os.path.isdir(value)
            assert os.path.isfile(os.path.join(value, "bash")) \
                or os.path.isfile(os.path.join(value, "bash.exe"))
            
            if platform.system() != "Windows":
                print(("SignalP will allow cygwinDir to be set, but you should know that " + 
                    "it is irrelevant to do so in a non-Windows environment as I have " +
                    f"detected you are running '{platform.system()}' instead."))
        
        self._cygwinDir = value
    
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
    def cygwin_program_execution_check(cygwinDir, signalpExe, signalpCmd):
        '''
        Static helper to check whether the signalP executable can be run
        using Cygwin. This is necessary for signalP 4 to run on Windows.
        
        Parameters:
            cygwinDir -- the directory where Cygwin executables (including bash.exe)
                         are located
            signalpExe -- the full path to the signalp executable file
            signalpCmd -- a string containing the command to feed in to signalP;
                          since validation is finnicky, only two commands are supported
                          herein i.e., "-h" for help, and "--
        '''
        # Validate that input parameters have any chance of success
        assert os.path.isdir(cygwinDir), \
            f"cygwin_program_execution_check won't even try to run because '{cygwinDir}' does not exist as a directory!"
        assert os.path.isfile(signalpExe), \
            f"cygwin_program_execution_check won't even try to run because '{signalpExe}' does not exist as a file!"
        assert os.path.isfile(os.path.join(cygwinDir, "bash")) or os.path.isfile(os.path.join(cygwinDir, "bash.exe")), \
            f"cygwin_program_execution_check won't even try to run because the bash.exe file does not exist at '{cygwinDir}'!"
        
        # Create script for testing signalP execution
        prefixHash = Conversion.get_hash_for_input_sequences(" ".join([cygwinDir, signalpExe, signalpCmd]))
        scriptFileName = ZS_Utility.tmp_file_name_gen(f"cygwin_test_{prefixHash}", "sh")
        scriptContents = Path(signalpExe).as_posix() + " {0}".format(signalpCmd.strip(" "))
        
        with open(scriptFileName, "w") as fileOut:
            fileOut.write(scriptContents)
        
        # Format cmd for execution
        cmd = os.path.join(cygwinDir, "bash") + " -l -c " + os.path.abspath(scriptFileName).replace("\\", "/")
        run_cmd = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
        cmdout, cmderr = run_cmd.communicate()
        os.remove(scriptFileName) # clean up temporary file
        
        # Check if it ran successfully
        if cmderr.decode("utf-8") != "" and not "perl: warning: falling back to the standard locale" in cmderr.decode("utf-8").lower():
            '''Need the above extra check for signalP since, on Windows at least, you can receive perl warnings which don't impact
            program operations. I think if that 'falling back' line is in stderr, nothing more serious will be present in stderr -
            this isn't completely tested, however.'''
            
            raise ValueError(("cygwin_program_execution_check has indicated a failure to execute " +
                              "signalP using Cygwin.\n stderr is below for debugging purposes.\n " +
                              cmderr.decode("utf-8")))
    
    @staticmethod
    def _popen_handler(cmd, cygwinDir=None, keepOut=True, keepErr=True):
        '''
        Hidden static function for helping to run a cmd list depending on system
        and whether Cygwin is being used or not.
        
        Parameters:
            cmd -- a list containing strings for use with Popen.
            cygwinDir -- None if Cygwin is not needed, or the string location of
                         the Cygwin bin directory
        Returns:
            stdout -- utf-8 decoded string of the stdout value
            stderr -- utf-8 decoded string of the stderr value
        '''
        assert isinstance(cmd, list), \
            "_format_subprocess_cmd failed because arguments provided are not in a list"
        
        if cygwinDir == None or platform.system() != "Windows":
            """This extra platform.system() check prevents spurious setting of cygwinDir from
            initiating a Cygwin execution on a non-Windows OS"""
            run = subprocess.Popen(" ".join(cmd),
                                   shell = True,
                                   stdout = subprocess.PIPE if keepOut else subprocess.DEVNULL,
                                   stderr = subprocess.PIPE if keepErr else subprocess.DEVNULL)
            stdout, stderr = run.communicate()
        else:
            # Strip any wsl values if needed
            if " ".join(cmd[0: 3]) == "wsl ~ -e":
                cmd = cmd[3:]
            
            # Create script for execution
            prefixHash = Conversion.get_hash_for_input_sequences(" ".join(cmd))
            scriptFileName = ZS_Utility.tmp_file_name_gen(f"cygwin_sigp_{prefixHash}", "sh")
            scriptContents = Path(ZS_Utility.convert_wsl_to_windows_path(cmd[0])).as_posix() + " {0}".format(" ".join(cmd[1: ])).replace('\\', '/')
            
            with open(scriptFileName, "w") as fileOut:
                fileOut.write(scriptContents)
            
            # Format cmd for execution
            cmd = os.path.join(cygwinDir, "bash") + " -l -c " + os.path.abspath(scriptFileName).replace("\\", "/")
            run = subprocess.Popen(cmd,
                                   shell = True,
                                   stdout = subprocess.PIPE if keepOut else subprocess.DEVNULL,
                                   stderr = subprocess.PIPE if keepErr else subprocess.DEVNULL)
            stdout, stderr = run.communicate()
            os.remove(scriptFileName) # clean up temporary file
        
        stdout = stdout.decode("utf-8").lower() if stdout != None else stdout
        stderr = stderr.decode("utf-8").lower() if stderr != None else stderr
        return stdout, stderr
    
    @staticmethod
    def signalP_is_valid(signalpExe, cygwinDir=None):
        '''
        Static helper to ascertain whether the signalP executable provided
        is compatible with the ZS_SignalPIO class.
        
        Parameters:
            signalpExe -- a string indicating the location of the executable itself
            cygwinDir -- either None if Cygwin is not needed for running, or the
                         location of the Cygwin bin directory where the bash.exe file
                         is located
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
        sigpout, sigperr = SignalP._popen_handler(cmd, cygwinDir)
        if not any([eStr in sigpout or eStr in sigperr for eStr in EXPECTED_STRINGS]):
            "Different SignalP versions output on stdout or stderr which is stupid"
            return False, sigperr
        else:
            return True, ""
    
    @staticmethod
    def get_signalP_version(signalpExe, cygwinDir=None):
        '''
        Static helper to check for the signalP executable version. The executable
        itself will be validated beforehand by the signalP_is_valid() static method.
        
        Parameters:
            signalpExe -- a string indicating the location of the executable itself
            cygwinDir -- either None if Cygwin is not needed for running, or the
                         location of the Cygwin bin directory where the bash.exe file
                         is located
        Returns:
            version -- an integer of 4, 5, or 6 depending on the version; any other
                       version will result in an error since we do not handle them.
        '''
        # Pre-validate
        sigpIsValid, notValidMessage = SignalP.signalP_is_valid(signalpExe, cygwinDir)
        if not sigpIsValid:
            raise Exception(("signalP_is_valid indicates invalid signalP executable prior to signalP " +
                             f"version checking! (error = '{notValidMessage}')"))
        
        # Check which instruction we need to obtain version
        cmd = ZS_Utility.base_subprocess_cmd(signalpExe) + ["-h"]
        sigpout, sigperr = SignalP._popen_handler(cmd, cygwinDir, keepOut=False)
        
        if "-version" in sigperr:
            "Need to have it as a plain string to work??"
            cmd = ZS_Utility.base_subprocess_cmd(signalpExe) + ["--version"] # cmd for version 5 or 6
        else:
            cmd = ZS_Utility.base_subprocess_cmd(signalpExe) + ["-V"] # cmd for version 4
        
        # Get version stdout
        versionout, versionerr = SignalP._popen_handler(cmd, cygwinDir)
        if versionerr != "":
            raise Exception("get_signalP_version ran into unexpected failure with stderr of {0}".format(
                versionerr
            ))
        
        # Extract version from stdout
        if "signalp 4" in versionout:
            return 4
        elif "signalp version 5" in versionout:
            return 5
        elif "signalP 6" in versionout:
            return 6
        else:
            raise ValueError("get_signalP_version could not understand or does not support {0}".format(
                versionout
            ))
    
    def signalp(self, withScore=False):
        '''
        Performs the signalP operation, automatically producing the output, parsing it,
        and cleaning up afterwards (if .clean == True).
        
        Parameters:
            withScore -- a bool indicating whether the score of the signal peptide
                         prediction should be returned (True) or not (False)
        Returns:
            resultsDict -- a dictionary containing signalP results as arising from
                           SignalP.parse_sigp_gff3()
        
        KNOWN ERRORS:
            - Segmentation faults in WSL for signalP 4; I tried to make it compatible
              but it just won't work.
              - Hence, this version is only supported when running through Cygwin.
            - Doesn't handle signalP 6 on Windows; I just haven't installed and tested
              it is all. It takes up a lot of file storage to install this version...
        '''
        # Exit if unhandled parameter / system combinations exist
        if platform.system() == "Windows":
            if self.version == 4 and self.cygwinDir == None:
                raise AttributeError(("SignalP 4 can only run on Windows with the use of Cygwin. " +
                                      "You must set the .cygwinDir attribute of this object first " +
                                      "prior to calling the .signalp() method"))
            elif self.version == 4:
                SignalP.cygwin_program_execution_check(self.cygwinDir, self.signalpExe, "-h")
            
            if self.version == 6:
                raise NotImplementedError("Unfortunately, I haven't implemented SignalP 6 on Windows yet")
        
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
            sigp4TmpFile = os.path.abspath(ZS_Utility.tmp_file_name_gen(f"{prefixHash}_tmp", "gff3"))
            
            # Format the cmd now
            cmd = ZS_Utility.base_subprocess_cmd(self.signalpExe)
            cmd += [
                "-t", self.organism,
                "-f", "short",
                "-n", f'{sigp4TmpFile}',
                f'{ZS_Utility.convert_windows_to_wsl_path(fastaQuery)}'
                    if platform.system() == "Windows" and self.version != 4
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
        sigpout, sigperr = SignalP._popen_handler(cmd, self.cygwinDir, keepOut=False)
        if sigperr != '' and not "no sequences predicted with a signal peptide" in sigperr:
            raise Exception('SignalP error text below\n' + sigperr)
        
        # Parse output file
        if self.version == 4:
            sigpGff3 = sigp4TmpFile
        elif self.version == 5:
            sigpGff3 = os.path.join(os.path.dirname(self.signalpExe), prefixHash + ".gff3")
        elif self.version == 6:
            sigpGff3 = os.path.join(sigp6TmpDir, "output.gff3")
        
        if os.path.isfile(sigpGff3):
            self.results = SignalP.parse_sigp_gff3(sigpGff3, withScore) # store for safe keeping
        else:
            self.results = {} # this can occur with signalP 4 when no sequences are predicted
        
        # Clean up output and temporary files
        if isTemporary:
            os.unlink(fastaQuery)
        
        if self.version == 4 and self.clean == True and os.path.isfile(sigpGff3):
            os.unlink(sigp4TmpFile)
        elif self.version == 5 and self.clean == True and os.path.isfile(sigpGff3):
            summaryFile = os.path.join(os.path.dirname(self.signalpExe), prefixHash + "_summary.signalp5")
            os.unlink(sigpGff3)
            os.unlink(summaryFile) # signalp5 also creates an extra summary file; we need to clean it, too
        elif self.version == 6 and self.clean == True and os.path.isdir(sigp6TmpDir):
            shutil.rmtree(sigp6TmpDir)
        
        return self.results # return for calling function
    
    @staticmethod
    def parse_sigp_gff3(sigpGff3, withScore=False):
        '''
        Parameters:
            sigpGff3 -- a string indicating the location of a GFF3 produced by
                        signalP
            withScore -- a bool indicating whether the score of the signal peptide
                         prediction should be returned (True) or not (False)
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
                    if withScore:
                        sigpDict[seqid] = [int(start), int(end), float(score)]
                    else:
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

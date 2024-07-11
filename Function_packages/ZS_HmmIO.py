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
    
    def create(self, inputFasta, outputFileName, hmmName=None, isNucleotide=False):
        '''
        This method will handle the creation of a HMM file using hmmbuild and hmmpress.
        
        Params:
            inputFasta -- a string indicating the location of a FASTA file OR a ZS_SeqIO.FASTA object
                          OR a string pointing to an existing HMM file.
            outputFileName -- a string providing the file name for our created HMM. This can include
                              the path to where you want the file to be written
            hmmName -- OPTIONAL; a string indicating the name of the HMM to be created, or None to have the name
                       be obtained from the FASTA file name. Default == None.
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
        HMM.hmmbuild(self.hmmbuild, fastaFileName, outputFileName, hmmName, isNucleotide)
        HMM.hmmpress(self.hmmpress, outputFileName)
        
        # Clean up temporary file if relevant
        if fastaIsTemporary:
            os.unlink(fastaFileName)
    
    @staticmethod
    def hmmbuild(hmmbuildExe, fastaFileName, outputFileName, hmmName=None, isNucleotide=False):
        '''
        Parameters:
            hmmbuildExe -- a string indicating the location of the hmmbuild executable.
            fastaFileName -- a string indicating the location of a FASTA file to be hmmbuild-ed.
            outputFileName -- a string indicating the name and, optionally, the path of the HMM file to be created.
            hmmName -- OPTIONAL; a string indicating the name of the HMM to be created, or None to have the name
                       be obtained from the FASTA file name. Default == None.
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
        
        # Handle hmmName
        if not hmmName == None:
            cmd += ["-n", hmmName]
        
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
    This Class provides methods that make use of a HMM file to run
    various HMMER modules.
    
    Parameters:
        hmmFile -- a string indicating the location of a HMM file that has been
                   built and pressed using HMMER.
        threads -- OPTONAL; an integer controlling the number of HMMER threads to specify.
        evalue -- OPTIONAL; a float indicating the significance threshold for results to be returned.
    '''
    def __init__(self, hmmerDir, hmmFile, threads=1, evalue=10):
        self.hmmerDir = hmmerDir
        self.hmmFile = hmmFile
        self.threads = threads
        self.evalue = 10
        self.domDict = None
    
    @property
    def hmmerDir(self):
        return self._hmmerDir
    
    @hmmerDir.setter
    def hmmerDir(self, value):
        convertedValue = ZS_Utility.convert_to_wsl_if_not_unix(value)
        assert ZS_Utility.wsl_exists(convertedValue, isFolder=True), \
            f"hmmer folder not found at '{convertedValue}' after WSL compatibility conversion"
        self._hmmerDir = convertedValue
        self.hmmsearch = value + "/hmmsearch"
        self.nhmmer = value + "/nhmmer"
    
    @property
    def hmmsearch(self):
        return self._hmmsearch
    
    @hmmsearch.setter
    def hmmsearch(self, value):
        convertedValue = ZS_Utility.convert_to_wsl_if_not_unix(value)
        assert ZS_Utility.wsl_exists(convertedValue), \
            f"hmmsearch executable not found at '{convertedValue}' after WSL compatibility conversion"
        self._hmmsearch = convertedValue
    
    @property
    def nhmmer(self):
        return self._nhmmer
    
    @nhmmer.setter
    def nhmmer(self, value):
        convertedValue = ZS_Utility.convert_to_wsl_if_not_unix(value)
        assert ZS_Utility.wsl_exists(convertedValue), \
            f"nhmmer executable not found at '{convertedValue}' after WSL compatibility conversion"
        self._nhmmer = convertedValue
    
    @property
    def hmmFile(self):
        return self._hmmFile
    
    @hmmFile.setter
    def hmmFile(self, value):
        convertedValue = ZS_Utility.convert_to_wsl_if_not_unix(value)
        assert ZS_Utility.wsl_exists(convertedValue), \
            f"HMM file not found at '{convertedValue}' after WSL compatibility conversion"
        self._hmmFile = convertedValue
    
    @property
    def threads(self):
        return self._threads
    
    @threads.setter
    def threads(self, value):
        assert isinstance(value, int), "threads must be an integer"
        assert value > 0, "threads must be greater than 0"
        self._threads = value
    
    @property
    def evalue(self):
        return self._evalue
    
    @evalue.setter
    def evalue(self, value):
        assert isinstance(value, float) or isinstance(value, int), "evalue must be an integer or float"
        assert value > 0, "evalue must be greater than 0"
        self._evalue = value
    
    def run(self, inputFasta, outputFileName, isNucleotide=False):
        '''
        Handles the running of hmmsearch using the HMM file (set on object initialisation)
        and a FASTA file provided herein.
        
        Parameters:
            inputFasta -- a string indicating the location of a FASTA file OR a ZS_SeqIO.FASTA object
                          OR a string pointing to an existing HMM file.
            outputFileName -- a string providing the file name for our created HMM. This can include
                              the path to where you want the file to be written
            isNucleotide -- OPTIONAL; a boolean to indicate whether the provided FASTA contains nucleotide
                            sequences or not. Default == False i.e., FASTA contains protein sequences.
        Returns:
            domDict -- a dictionary containing the parsed results of the search. This value is also
                       stored on the object for later access at self.domDict.
        '''
        # Validate input value type
        if isinstance(inputFasta, str):
            assert os.path.isfile(inputFasta), f"ERROR: HMMER.search() could not find the FASTA file '{inputFasta}'"
            fastaFileName, fastaIsTemporary = inputFasta, False
        elif hasattr(inputFasta, "isFASTA") and inputFasta.isFASTA is True:
            fastaFileName, fastaIsTemporary = ZS_SeqIO.Conversion.get_filename_for_input_sequences(inputFasta)
        else:
            raise Exception(f"ERROR: HMMER.search() requires a FASTA file or FASTA object as input; did not understand '{inputFasta}'")
        
        # Validate that output file does not already exist
        assert not os.path.exists(outputFileName), f"ERROR: HMMER.search() will not overwrite existing file '{outputFileName}'"
        
        # Run hmmsearch or nhmmer
        HMMER.search(self.nhmmer if isNucleotide else self.hmmsearch, self.hmmFile,
                     fastaFileName, outputFileName, self.threads, self.evalue)
        
        # Clean up temporary file if relevant
        if fastaIsTemporary:
            os.unlink(fastaFileName)
        
        # Parse search results
        self.domDict = nhmmer_parse(outputFileName, self.evalue, extendedDetails=True) \
                       if isNucleotide else hmmer_parse(outputFileName, self.evalue)
        
        # Return result
        return self.domDict
    
    @staticmethod
    def search(queryExe, hmmFileName, fastaFileName, outputFileName, threads=1, evalue=10):
        '''
        Parameters:
            queryExe -- a string indicating the location of the hmmsearch or nhmmer executable;
                        use nhmmer for nucleotide searches and hmmsearch for protein searches.
            hmmFileName -- a string indicating the location of a HMM file to be hmmpress-ed.
            fastaFileName -- a string indicating the location of a FASTA file to be searched.
            outputFileName -- a string indicating the location to write output to.
            threads -- OPTIONAL; an integer controlling the number of HMMER threads to specify;
                       default is 1.
            evalue -- OPTIONAL; a float indicating the significance threshold for results to be returned;
                      default is 10.
        Returns:
            domDict -- a dictionary containing the parsed results of the search.
        '''
        assert ZS_Utility.wsl_exists(hmmFileName), \
            f"ERROR: HMMER.search() could not find HMM file '{hmmFileName}'"
        assert not ZS_Utility.wsl_exists(outputFileName), \
            f"ERROR: HMMER.search() will not allow overwriting '{outputFileName}'"
        assert isinstance(threads, int), "threads must be an integer"
        assert threads > 0, "threads must be greater than 0"
        assert isinstance(evalue, float) or isinstance(evalue, int), "evalue must be an integer or float"
        assert evalue > 0, "evalue must be greater than 0"
        
        # Figure out which program we are running
        assert queryExe.endswith("hmmsearch") or queryExe.endswith("nhmmer"), \
            "ERROR: HMMER.search() requires either hmmsearch or nhmmer as the queryExe"
        isNucleotide = queryExe.endswith("nhmmer")
        
        # Construct the cmd for subprocess
        cmd = ZS_Utility.base_subprocess_cmd(queryExe)
        
        # Set other arguments
        cmd += [
            "--cpu", str(threads), "-E", str(evalue),
            "--tblout" if isNucleotide else "--domtblout",
            ZS_Utility.convert_to_wsl_if_not_unix(outputFileName),
            ZS_Utility.convert_to_wsl_if_not_unix(hmmFileName),
            ZS_Utility.convert_to_wsl_if_not_unix(fastaFileName)
        ]
        
        # Format commands for Linux
        if platform.system() != "Windows":
            cmd = " ".join(cmd)
        
        # Run the command
        run_search = subprocess.Popen(cmd, shell = True,
                                      stdout = subprocess.PIPE,
                                      stderr = subprocess.PIPE)
        searchout, searcherr = run_search.communicate()
        
        # Check to see if there was an error
        if searcherr.decode("utf-8") != "":
            raise Exception(("ERROR: HMMER.search() encountered an error; have a look " +
                            f'at the stdout ({searchout.decode("utf-8")}) and stderr ' + 
                            f'({searcherr.decode("utf-8")}) to make sense of this.'))

if __name__ == "__main__":
    pass

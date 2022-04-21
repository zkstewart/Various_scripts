#! python3
# ZS_ClustIO.py
# Specifies the CDHIT Class for performing CD-HIT reduction of
# string fasta files, ZS_SeqIO.FASTA and ZS_SeqIO.FastASeq objects.

import os, sys, subprocess, hashlib, time, random, shutil

sys.path.append(os.path.dirname(__file__))
from ZS_SeqIO import FASTA

class CDHIT:
    '''
    The CDHIT Class provides easy access to CD-HIT functionality to perform
    redundancy reduction of FASTA files. These can be in the form of providing
    the location of a FASTA file as string, or with ZS_SeqIO modules i.e.,
    FASTA and FastASeq.
    
    Attributes:
        fasta (REQUIRED) -- a string, FASTA, or FastASeq object from the ZS_SeqIO suite.
        molecule (REQUIRED) -- a string, indicating whether the FASTA file contains protein
                               or nucleotide sequences. Value must be in the list:
                               [n, nucl, nucleotide, p, prot, protein]
        cdhitDir (OPTIONAL) -- a string indicating the location of the CD-HIT executables,
                               if they're not available from the PATH variable.
    '''
    def __init__(self, fasta, molecule, cdhitDir=None):
        # Validate input types and values
        assert type(fasta).__name__ == "str" \
            or type(fasta).__name__ == "FASTA" \
            or type(fasta).__name__ == "ZS_SeqIO.FASTA" \
            or type(fasta).__name__ == "FastASeq" \
            or type(fasta).__name__ == "ZS_SeqIO.FastASeq"
        
        assert isinstance(molecule, str)
        assert molecule.lower() in ["n", "nucl", "nucleotide", "p", "prot", "protein"]
        
        assert isinstance(cdhitDir, str) or cdhitDir == None
        
        # Validate that fasta exists if specified as a string (file location)
        if type(fasta).__name__ == "str" and not os.path.isfile(fasta):
            raise Exception("fasta parameter is a string, but does not point to an existing file location")
        self.fasta = fasta
        
        # Validate that cdhitDir exists and contains the relevant executables if specified as a string
        if isinstance(cdhitDir, str):
            assert os.path.isdir(cdhitDir), "cdhitDir does not exist"
            assert os.path.isfile(os.path.join(cdhitDir, "cd-hit")) or os.path.isfile(os.path.join(cdhitDir, "cd-hit.exe")), "cd-hit executable not found at \"{0}\"".format(cdhitDir)
            assert os.path.isfile(os.path.join(cdhitDir, "cd-hit-est")) or os.path.isfile(os.path.join(cdhitDir, "cd-hit-est.exe")), "cd-hit-est executable not found at \"{0}\"".format(cdhitDir)
        else:
            assert shutil.which("cd-hit") != None, "cd-hit executable not found in path"
            assert shutil.which("cd-hit-est") != None, "cd-hit-est executable not found in path"
        self.cdhitDir = cdhitDir
        
        # Coerce molecule into a consistent value
        self.molecule = "nucleotide" if molecule.lower() in ["n", "nucl", "nucleotide"] else "protein"
        
        # Set default attributes
        self.identity = 0.9
        self.word_length = 5 if self.molecule == "protein" else 8
        self.local = False
        self.shorter_cov_pct = 0.0
        self.longer_cov_pct = 0.0
        self.mem = 1000 # 1GB by default
        self.threads = 1
        self.clean = True
    
    def set_identity(self, num):
        '''
        Relates to the -c parameter for CD-HIT. Controls the sequence identity
        threshold for clustering sequences together. This method also controls
        the word length setting, since the word length should always correspond
        to the identity threshold anyway.
        
        Parameters:
            num -- should be a valid float in the range of 0 to 1 (inclusive)
        '''
        assert isinstance(num, float)
        if not 0<=num<=1:
            raise Exception("Identity (-c) must be between 0 and 1 (inclusive)")
        
        self.identity = num
        
        # Set word length intelligently, depending on molecule type
        if self.molecule == "nucleotide":
            if 0.9<=num<=1.0:
                self.word_length = 8
            elif 0.88<=num<=0.9:
                self.word_length = 7
            elif 0.85<=num<=0.88:
                self.word_length = 6
            elif 0.80<=num<=0.85:
                self.word_length = 5
            elif 0.75<=num<=0.80:
                self.word_length = 4
            else:
                self.word_length = 3 # is this okay? I don't know
        else:
            if 0.7<=num<=1.0:
                self.word_length = 5
            elif 0.6<=num<=0.7:
                self.word_length = 4
            elif 0.5<=num<=0.6:
                self.word_length = 3
            else:
                self.word_length = 2 # assuming this is also okay
    
    def set_local(self):
        '''
        Relates to the -G parameter for CD-HIT. Changes how identity is scored
        to only consider the local region in which an alignment was found.
        
        This is NOT the default CD-HIT behaviour.
        '''
        self.local = True
    
    def set_global(self):
        '''
        Relates to the -G parameter for CD-HIT. Changes how identity is scored
        to consider the entire sequence length as part of the calculation.
        
        This is the default CD-HIT behaviour
        '''
        self.local = False
    
    def set_shorter_cov_pct(self, num):
        '''
        Relates to the -aS parameter for CD-HIT. Controls minimum threshold
        of sequence that must align well from the shorter sequence. For example,
        a value of 0.5 means at least half of the shorter sequence in a pairwise
        alignment must align against the longer sequence. If it's less, it
        will be discarded.
        
        0.0 means no threshold will be enforced.
        
        Parameters:
            num -- should be a valid float in the range of 0 to 1 (inclusive)
        '''
        assert isinstance(num, float)
        if not 0<=num<=1:
            raise Exception("shorter_cov_pct (-aS) must be between 0 and 1 (inclusive)")
        
        self.shorter_cov_pct = num
    
    def set_longer_cov_pct(self, num):
        '''
        Relates to the -aL parameter for CD-HIT. Controls minimum threshold
        of sequence that must align well from the longer sequence. For example,
        a value of 0.5 means at least half of the longest sequence in a pairwise
        alignment must align against the shorter sequence. If it's less, it
        will be discarded.
        
        0.0 means no threshold will be enforced.
        
        Parameters:
            num -- should be a valid float in the range of 0 to 1 (inclusive)
        '''
        assert isinstance(num, float)
        if not 0<=num<=1:
            raise Exception("longer_cov_pct (-aL) must be between 0 and 1 (inclusive)")
        
        self.longer_cov_pct = num
    
    def set_mem(self, num):
        '''
        Relates to the -M parameter for CD-HIT. Controls how much memory is allowed
        to be used by the algorithm. Number is in megabytes.
        
        Parameters:
            num -- should be a valid integer that is minimally bounded at 100.
        '''
        assert isinstance(num, int)
        if not 100<=num:
            raise Exception("mem (-M) must be at least 100 megabytes")
        
        self.mem = num
    
    def set_threads(self, num):
        '''
        Relates to the -T parameter for CD-HIT. Controls how many threads CD-HIT
        will use.
        
        0.0 means it will use as many threads as there are CPU cores.
        
        Parameters:
            num -- should be a valid integer that is minimally bounded at 0.
        '''
        assert isinstance(num, int)
        if not 0<=num:
            raise Exception("threads (-T) must be greater than or equal to 0")
        
        self.threads = num
    
    def set_clean(self, clean):
        '''
        This method allows the clean attribute to be set, which controls
        whether CD-HIT output files are kept or not.
        
        Parameters:
            clean -- a Boolean of True or False
        '''
        assert isinstance(clean, bool)
        
        self.clean = clean

    def cdhit(self, fasta, outputDir, outputFasta):
        '''
        Performs the CD-HIT operation.
        
        This method does not use .fasta because it may be in the FASTA or
        FastASeq object types, rather than a string which we require here.
        
        Parameters:
            fasta -- a string indicating the location of a FASTA file to cluster.
            outputDir -- a string indicating the location to write the output FASTA and .clstr files to.
            outFile -- a string indicating the name for the output FASTA file. File must not already exist!
        '''
        assert isinstance(fasta, str)
        assert isinstance(outputDir, str)
        assert isinstance(outputFasta, str)
        
        assert os.path.isfile(fasta), "fasta file does not exist"
        assert os.path.isdir(outputDir), "output directory does not exist"
        assert os.path.basename(outputFasta) == outputFasta, "output fasta file needs to be just the file name; its location is specified in the outputDir method parameter"
        assert not os.path.isfile(os.path.join(outputDir, outputFasta)), "\"{0}\" already exists; cdhit method won't overwrite it".format(os.path.join(outputDir, outputFasta))
        
        if self.molecule == "nucleotide":
            program = os.path.join(self.cdhitDir if self.cdhitDir != None else "", 'cd-hit-est')
        else:
            program = os.path.join(self.cdhitDir if self.cdhitDir != None else "", 'cd-hit')
        
        cmd = f"{program} -i {fasta} -o {os.path.join(outputDir, outputFasta)} -c {self.identity} -n {self.word_length} -G {0 if self.local else 1} -aS {self.shorter_cov_pct} -aL {self.longer_cov_pct} -M {self.mem} -T {self.threads}"
        run_cdhit = subprocess.Popen(cmd, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
        cdout, cderr = run_cdhit.communicate()
        if cderr.decode("utf-8") != '':
            raise Exception('CD-HIT Error text below' + str(cderr.decode("utf-8")))
    
    def get_cdhit_results(self, workingDir="."):
        '''
        This function pipelines the process of obtaining CD-HIT results. Intermediate files are
        deleted automatically, and hence this function will only result in the return of the
        clustered FASTA object.
        
        Parameters:
            workingDir -- a string indicating the location to write CD-HIT results to.
        Returns:
            FASTA_obj -- a ZS_SeqIO.FASTA object of the clustered CD-HIT results.
            cdhitResultFile -- a string indicating the file name of the results file. If self.clean is True,
                               this will instead return None.
        '''
        
        assert os.path.isdir(workingDir), "workingDir must already exist, or just leave it as default to write to current working directory"
        
        # Get file names for query and target after data type coercion
        f, fIsTemporary = self._get_filename_for_fasta()
        
        # Get hash for temporary file creation
        tmpHash = hashlib.sha256(bytes(f + str(time.time()) + str(random.randint(0, 100000)), 'utf-8') ).hexdigest()
        
        # Run CD-HIT
        tmpResultName = self._tmp_file_name_gen("cdhit_result_tmp" + tmpHash[0:20], "fasta")
        self.cdhit(f, workingDir, tmpResultName) # "." for working directory being the current one
        
        # Parse CD-HIT results
        result_FASTA_obj = FASTA(tmpResultName)
        
        # Clean up f temporary file
        if fIsTemporary:
            os.unlink(f)
        
        # Clean up results and return (if relevant)
        if self.clean:
            os.unlink(tmpResultName)
            os.unlink(tmpResultName + ".clstr")
            return result_FASTA_obj, None
        # Or just return results
        else:
            return result_FASTA_obj, tmpResultName
    
    def _get_filename_for_fasta(self):
        '''
        Hidden method for use when boiling down one of the three data types (FASTA, FastASeq, and string)
        into a string representing the fasta file name.
        
        If it's already a string, this does nothing. Otherwise, it will make sure a file exists
        with the FASTA data contents for use by CD-HIT.
        
        Returns:
            fasta -- a string indicating the file name for the .fasta value. If it is
                     already a string, this will be the same. Otherwise, it will be a file
                     name containing the contents of the FASTA or FastASeq object.
            isTemporary -- a Boolean indicating whether the returned value has been
                           created by this method as a temporary file (True) or if it was
                           already existing (False, i.e., .fasta was already a string)
        '''
        
        f = self.fasta
        
        # Get a hash for temporary file creation
        hashForTmp = f if isinstance(f, str) \
                     else f.fileOrder[0][0] if type(f).__name__ == "ZS_SeqIO.FASTA" or type(f).__name__ == "FASTA" \
                     else f.id
        tmpHash = hashlib.sha256(bytes(hashForTmp + str(time.time()) + str(random.randint(0, 100000)), 'utf-8') ).hexdigest()
        
        # If f is a FastASeq, get it as a FASTA object
        if type(f).__name__ == "FastASeq" or type(f).__name__ == "ZS_SeqIO.FastASeq":
            tmpFastaName = self._tmp_file_name_gen("f_fasta_tmp" + tmpHash[0:20], "fasta")
            with open(tmpFastaName, "w") as fileOut:
                fileOut.write(">{0}\n{1}\n".format(f.id, f.seq))
            f = FASTA(tmpFastaName)
            os.unlink(tmpFastaName)
        
        # If f is a FASTA, make it into a file
        isTemporary = False
        if type(f).__name__ == "FASTA" or type(f).__name__ == "ZS_SeqIO.FASTA":
            tmpFName = self._tmp_file_name_gen("cdhit_tmp" + tmpHash[0:20], "fasta")
            f.write(tmpFName)
            f = tmpFName # after this point, f will be a string indicating a FASTA file name
            isTemporary = True # if we set this, f was not originally a string
        
        return f, isTemporary
    
    def _tmp_file_name_gen(self, prefix, suffix):
        '''
        Hidden function for use by this Class when creating temporary files.
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

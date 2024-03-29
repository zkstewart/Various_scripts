#! python3
# ZS_ClustIO.py
# Specifies the CDHIT Class for performing CD-HIT reduction of
# string fasta files, ZS_SeqIO.FASTA and ZS_SeqIO.FastASeq objects.

import os, sys, subprocess, shutil, platform
from pathlib import Path
from hashlib import sha256

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from ZS_SeqIO import FASTA, Conversion
from ZS_Utility import tmp_file_name_gen, base_subprocess_cmd, convert_windows_to_wsl_path

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
        self.identity = 0.9 # implicitly sets self.word_length
        self.local = False
        self.shorter_cov_pct = 0.0
        self.longer_cov_pct = 0.0
        self.mem = 1000 # 1GB by default
        self.threads = 1
        self.clean = True
        self.description_length = 0
        
        # Results storage values
        self.resultFasta = None
        self.resultClusters = None
    
    @property
    def identity(self):
        return self._identity
    
    @identity.setter
    def identity(self, num):
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
        
        self._identity = num
        
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
    
    @property
    def word_length(self):
        return self._word_length
    
    @word_length.setter
    def word_length(self, num):
        '''
        Relates to the -n parameter for CD-HIT. Controls algorithmic behaviour
        of CD-HIT, and should be set to a value according to the identity parameter.
        
        Usually, you should allow the identity setter to automatically specify the
        word length, but this setter allows you to override that behaviour.
        
        Parameters:
            num -- should be a valid integer that is minimally bounded at 1 and
                   maximally bounded at 11.
        '''
        assert isinstance(num, int)
        if not (1 <= num <= 11):
            raise Exception("word length (-n) must be in the range of 1 to 11 (inclusive)")
        
        self._word_length = num
    
    @property
    def mem(self):
        return self._mem
    
    @mem.setter
    def mem(self, num):
        '''
        Relates to the -M parameter for CD-HIT. Controls how much memory is allowed
        to be used by the algorithm. Number is in megabytes.
        
        Parameters:
            num -- should be a valid integer that is minimally bounded at 100.
        '''
        assert isinstance(num, int)
        if not num >= 100:
            raise Exception("mem (-M) must be at least 100 megabytes")
        
        self._mem = num
    
    @property
    def threads(self):
        return self._threads
    
    @threads.setter
    def threads(self, num):
        '''
        Relates to the -T parameter for CD-HIT. Controls how many threads CD-HIT
        will use.
        
        0 means it will use as many threads as there are CPU cores.
        
        Parameters:
            num -- should be a valid integer that is minimally bounded at 0.
        '''
        assert isinstance(num, int)
        if not num >= 0:
            raise Exception("threads (-T) must be greater than or equal to 0")
        
        self._threads = num
    
    @property
    def clean(self):
        return self._clean
    
    @clean.setter
    def clean(self, value):
        '''
        This method allows the clean attribute to be set, which controls
        whether CD-HIT output files are kept or not.
        
        Parameters:
            clean -- a Boolean of True or False
        '''
        assert isinstance(value, bool)
        
        self._clean = value
    
    @property
    def description_length(self):
        return self._description_length
    
    @description_length.setter
    def description_length(self, num):
        '''
        Relates to the -d parameter for CD-HIT. Controls the presentation of
        the output cluster file.
        
        Parameters:
            num -- should be a valid integer that is minimally bounded at 0.
        '''
        assert isinstance(num, int)
        if not num >= 0:
            raise Exception("description length (-d) must be >= 0")
        
        self._description_length = num
    
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
        if not (0<=num<=1):
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
        if not (0<=num<=1):
            raise Exception("longer_cov_pct (-aL) must be between 0 and 1 (inclusive)")
        
        self.longer_cov_pct = num

    def cdhit(self, fasta, outputDir, outputFasta):
        '''
        Performs the CD-HIT operation.
        
        This method does not use .fasta because it may be in the FASTA or
        FastASeq object types, rather than a string which we require here.
        
        Parameters:
            fasta -- a string indicating the location of a FASTA file to cluster.
            outputDir -- a string indicating the location to write the output FASTA and .clstr files to.
            outFile -- a string indicating the name for the output FASTA file. File must not already exist!
        Returns:
            cmd -- a string indicating the command executed for CD-HIT clustering.
        '''
        # Validate parameters
        assert isinstance(fasta, str)
        assert isinstance(outputDir, str)
        assert isinstance(outputFasta, str)
        
        assert os.path.isfile(fasta), "fasta file does not exist"
        assert os.path.isdir(outputDir), "output directory does not exist"
        assert os.path.basename(outputFasta) == outputFasta, \
            "output fasta file needs to be just the file name; its location is specified in the outputDir method parameter"
        assert not os.path.isfile(os.path.join(outputDir, outputFasta)), \
            f"'{os.path.join(outputDir, outputFasta)}' already exists; cdhit method won't overwrite it"
        
        # Figure out which CD-HIT executable we're using
        if self.molecule == "nucleotide":
            program = os.path.join(self.cdhitDir if self.cdhitDir != None else "", 'cd-hit-est')
        else:
            program = os.path.join(self.cdhitDir if self.cdhitDir != None else "", 'cd-hit')
        
        # Begin formatting cmd, converting to WSL paths where needed
        if platform.system() == "Windows":
            fasta = convert_windows_to_wsl_path(fasta)
            outputFile = convert_windows_to_wsl_path(os.path.join(outputDir, outputFasta))
        else:
            outputFile = os.path.join(outputDir, outputFasta)
        cmd = base_subprocess_cmd(program)
        
        # Format cmd and run it
        cmd += list(map(str, ["-i", fasta, "-o", outputFile,
               "-c", self.identity, "-n", self.word_length, "-G", "0" if self.local else "1", 
               "-aS", self.shorter_cov_pct, "-aL", self.longer_cov_pct, "-M", self.mem,
               "-T", self.threads, "-d", self.description_length]))
        
        if platform.system() != "Windows":
            run_cdhit = subprocess.Popen(" ".join(cmd), shell = True,
                                         stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        else:
            run_cdhit = subprocess.Popen(cmd, shell = True,
                                         stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        
        cdout, cderr = run_cdhit.communicate()
        if cderr.decode("utf-8") != '':
            raise Exception('CD-HIT Error text below' + str(cderr.decode("utf-8")))
        
        return cmd
    
    @staticmethod
    def parse_clstr_file(clstrFile):
        '''
        For this to be effective, you should make sure CD-HIT was run with -d 0
        so as to give the sequence ID in a format expected here.
        
        Parameters:
            clstrFile -- a string pointing to the location of a CD-HIT output cluster
                         file.
        Returns:
            clstrDict -- a dictionary with structure like:
                         {
                             0: [seqid1, seqid2, ...],
                             1: [ ... ],
                             ...
                         }
        '''
        clstrDict = {}
        with open(clstrFile, "r") as fileIn:
            for line in fileIn:
                sl = line.rstrip("\r\n ").split()
                
                # Handle cluster ID lines
                if line.startswith(">"):
                    thisCluster = int(sl[1])
                    clstrDict[thisCluster] = []
                
                # Handle content lines
                else:
                    seqID = sl[2].strip(">.")
                    clstrDict[thisCluster].append(seqID)
        return clstrDict
    
    def get_cdhit_results(self, workingDir=".", returnFASTA=True, returnClusters=False):
        '''
        This function pipelines the process of obtaining CD-HIT results. Intermediate files are
        deleted automatically, and hence this function will only result in the return of the
        clustered FASTA object.
        
        Parameters:
            workingDir -- a string indicating the location to write CD-HIT results to.
            returnFASTA -- a boolean indicating whether we should parse the output FASTA
                           file at the end of this, storing its value in .resultFasta as
                           a ZS_SeqIO.FASTA object.
            returnClusters -- a boolean indicating whether we should parse the .clstr
                              file at the end of this, storing its value in .resultClusters
                              as a dictionary.
        Returns:
            FASTA_obj -- a ZS_SeqIO.FASTA object of the clustered CD-HIT results.
            cdhitResultFile -- a string indicating the file name of the results file. If self.clean is True,
                               this will instead return None.
        '''
        assert os.path.isdir(workingDir), \
            "workingDir must already exist, or just leave it as default to write to current working directory"
        
        # Validate that running this function will result in some sort of output
        if self.clean == True:
            if returnFASTA == False and returnClusters == False:
                raise ValueError(("get_cdhit_results will not return you anything " +
                                  "at the end of running, since .clean is set to True " +
                                  "and you set both returnFASTA and returnClusters to False. " +
                                  "Since there's no point running, I'm not going to."))
        
        # Get file name after data type coercion
        f, fIsTemporary = Conversion.get_filename_for_input_sequences(self.fasta)
        
        # Get hash for temporary file creation
        tmpHash = Conversion.get_hash_for_input_sequences(f)
        
        # Run CD-HIT
        tmpResultName = tmp_file_name_gen("cdhit_result_tmp" + tmpHash, "fasta")
        self.cdhit(f, workingDir, tmpResultName) # "." for working directory being the current one
        
        # Parse CD-HIT results if desired
        if returnFASTA == True:
            result_FASTA_obj = FASTA(tmpResultName)
        else:
            result_FASTA_obj = None
        
        # Clean up f temporary file
        if fIsTemporary:
            os.unlink(f)
        
        # Parse clstr file if desired
        if returnClusters == True:
            clstrDict = CDHIT.parse_clstr_file(tmpResultName + ".clstr")
        else:
            clstrDict = None
        
        # Store results
        self.resultFasta = result_FASTA_obj
        self.resultClusters = clstrDict
        
        # Clean up results (if relevant) and return the output file name
        if self.clean:
            os.unlink(tmpResultName)
            os.unlink(tmpResultName + ".clstr")
            return None
        # Or just return results
        else:
            return tmpResultName
    
    def __repr__(self):
        return (
            f"<CDHIT object;identity={self.identity};local={self.local};" +
            f"shorter_cov_pct={self.shorter_cov_pct};longer_cov_pct={self.longer_cov_pct};" +
            f"mem={self.mem};threads={self.threads};clean={self.clean};" +
            f"description_length={self.description_length};" +
            f"resultFASTA contains data={'NO' if self.resultFasta == None else 'YES'};" +
            f"resultClusters contains data={'NO' if self.resultClusters == None else 'YES'}"
        )

class MM_Clust:
    '''
    The MM_Clust Class behaves as an abstract class for inheritance of basic
    validation logic of parameters shared in common.
    
    Attributes:
        mmDB (REQUIRED) -- a MM_DB object from ZS_BlastIO.
        evalue (OPTIONAL) -- a positive float with a minimum of 0.0 controlling the
                             E-value threshold for clustering; default == 1e-3.
        identity (OPTIONAL) -- a positive float in the range 0.0 -> 1.0 controlling the
                               sequence identity threshold for clustering; default == 0.9.
        cov_pct (OPTIONAL) -- a positive float in the range 0.0 -> 1.0 controlling the
                              amount of aligned residues in both shorter and longer
                              sequences; default == 0.8.
        clust_mode (OPTIONAL) -- a string in the list ["set-cover", "connected-component",
                                 "greedy"], corresponding to modes 0, 1, and 2,3;
                                 default == "set-cover", but optimal is probably
                                 "connected-component".
        threads (OPTIONAL) -- a positive integer for how many threads to use when running
                              MMseqs2 clustering (default==1).
        tmpDir (OPTIONAL) -- a string location for where MMseqs2 should keep temp files;
                             if unspecified, it will use the same location as
                             mmDB.tmpDir.
    '''
    def __init__(self, mmDB, evalue=1e-3, identity=0.9,
                 cov_pct=0.8, clust_mode="set-cover", threads=1, tmpDir=None):
        self.mmDB = mmDB
        self.evalue = evalue
        self.identity = identity
        self.cov_pct = cov_pct
        self.clust_mode = clust_mode
        self.threads = threads
        self.tmpDir = tmpDir
    
    @property
    def mmDB(self):
        return self._mmDB
    
    @mmDB.setter
    def mmDB(self, value):
        assert type(value).__name__ == "MM_DB" \
            or type(value).__name__ == "ZS_BlastIO.MM_DB" \
            or hasattr(value, "isMM_DB") and value.isMM_DB == True
        
        self._mmDB = value
    
    @property
    def evalue(self):
        return self._evalue
    
    @evalue.setter
    def evalue(self, value):
        assert isinstance(value, float) or isinstance(value, int)
        assert 0 <= value, "evalue must be >= 0"
        
        self._evalue = value
    
    @property
    def identity(self):
        return self._identity
    
    @identity.setter
    def identity(self, value):
        assert isinstance(value, float) or isinstance(value, int)
        assert 0 <= value <= 1.0, "identity must be in the range of 0.0 -> 1.0 (inclusive)"
        
        self._identity = value
    
    @property
    def cov_pct(self):
        return self._cov_pct
    
    @cov_pct.setter
    def cov_pct(self, value):
        assert isinstance(value, float) or isinstance(value, int)
        assert 0 <= value <= 1.0, "cov_pct must be in the range of 0.0 -> 1.0 (inclusive)"
        
        self._cov_pct = value
    
    @property
    def clust_mode(self):
        return self._clust_mode
    
    @clust_mode.setter
    def clust_mode(self, value):
        assert isinstance(value, str)
        assert value in ["set-cover", "connected-component", "greedy"], \
            "clust_mode must be one of the supported Linclust algorithms"
        
        self._clust_mode = ["set-cover", "connected-component", "greedy"].index(value)
    
    @property
    def threads(self):
        return self._threads
    
    @threads.setter
    def threads(self, value):
        assert isinstance(value, int)
        assert 0 < value, "threads must be a positive integer"
        
        self._threads = value
    
    @property
    def tmpDir(self):
        return self._tmpDir
    
    @tmpDir.setter
    def tmpDir(self, value):
        if value == None:
            value = self.mmDB.tmpDir
        assert type(value).__name__ == "str" \
            or isinstance(value, Path)
        value = os.path.abspath(value)
        
        if not os.path.isdir(os.path.dirname(value)):
            raise Exception((f"tmpDir's parent location ('{os.path.dirname(value)}') " +
                            "does not exist"))
        
        if not os.path.isdir(value):
            os.mkdir(value)
        
        self._tmpDir = value
    
    def cluster(self):
        raise NotImplementedError("Abstract method must be overridden")
    
    def tabulate(self):
        raise NotImplementedError("Abstract method must be overridden")
    
    def clean_all(self):
        raise NotImplementedError("Abstract method must be overridden")
    
    def parse_tsv(self, tabulatedFile):
        '''
        Parses the output of .tabulate() (aka mmseqs createtsv) and produces
        a cluster dictionary.
        
        Parameters:
            tabulatedFile -- a string indicating the file location of the
                             output from createtsv.
        Returns:
            clusterDict -- a dictionary with structure like:
                           {
                               0: [seqid1, seqid2, ...],
                               1: [ ... ],
                               ...
                           }
        '''
        clusterDict = {}
        clusterNum = -1
        lastCluster = None
        with open(tabulatedFile, "r") as fileIn:
            for line in fileIn:
                l = line.rstrip("\r\n ")
                if l != "":
                    clustID, seqID = l.split("\t")
                    if clustID != lastCluster:
                        clusterNum += 1
                    lastCluster = clustID
                    
                    clusterDict.setdefault(clusterNum, [])
                    clusterDict[clusterNum].append(seqID)
        return clusterDict
    
    def _header_parse_tsv(self, tabulatedFile):
        '''
        Parses the output of .tabulate() (aka mmseqs createtsv) and produces
        a cluster dictionary.
        
        This function is now deprecated since the method of creating a tsv
        has since changed. But, I'd like to keep the code around just in case.
        
        Parameters:
            tabulatedFile -- a string indicating the file location of the
                             output from createtsv.
        Returns:
            clusterDict -- a dictionary with structure like:
                           {
                               0: [seqid1, seqid2, ...],
                               1: [ ... ],
                               ...
                           }
        '''
        # Format the location of relevant files
        dbname = os.path.abspath(f"{self.mmDB.fasta}_seqDB")
        dbHeader = os.path.abspath(f"{dbname}_h")
        assert os.path.isfile(dbHeader), \
            f"MM_Clust object can't find file at expected location '{dbHeader}'!"
        
        # Parse the db header
        headDict = {}
        with open(dbHeader, "r") as fileIn:
            ongoingCount = 0
            for line in fileIn:
                seqID = line.strip("\x00").strip("\r\n ").split(" ")[0]
                headDict[ongoingCount] = seqID
                ongoingCount += 1
        
        # Parse the table file to a cluster dictionary
        clusterDict = {}
        with open(tabulatedFile, "r") as fileIn:
            lastID = None
            ongoingCount = -1 # first loop will iterate this to 0
            for line in fileIn:
                # Get details from two column TSV line
                clusterRep, seqIndex = line.rstrip("\r\n ").split("\t")
                seqIndex = int(seqIndex)
                
                # If we haven't seen this representative before, it's a new cluster
                if not clusterRep == lastID:
                    ongoingCount += 1
                
                # Store this in our dict
                clusterDict.setdefault(ongoingCount, [])
                clusterDict[ongoingCount].append(headDict[seqIndex])
                
                # Remember this line's representative for next loop
                lastID = clusterRep
        
        return clusterDict

class MM_Linclust(MM_Clust):
    '''
    The MM_Linclust Class provides the logic for running Linclust with an MM_DB object
    as input. The MMseqs executable location and tmpDir will be pulled from the
    mmDB input; a new tmpDir for running linclust can be specified if desired.
    
    Attributes:
        super class attributes -- see docstring for MM_Clust.
    '''
    def __init__(self, mmDB, evalue=1e-3, identity=0.9,
                 cov_pct=0.8, clust_mode="set-cover", threads=1, tmpDir=None):
        
        super().__init__(mmDB, evalue, identity, cov_pct,
                        clust_mode, threads, tmpDir)
    
    def cluster(self):
        '''
        Runs the Linclust process on the given mmDB.
        '''
        
        # Run DB generation & indexing if relevant
        self.mmDB.generate()
        self.mmDB.index()
        
        # Get hash for parameters combination
        strForHash = str(self.evalue) + str(self.identity) + str(self.cov_pct) + \
            str(self.clust_mode)
        strHash = sha256(bytes(strForHash, 'utf-8')).hexdigest()
        
        # Specify file locations
        fasta = self.mmDB.fasta
        dbname = os.path.abspath(f"{self.mmDB.fasta}_{strHash}_linclustDB")
        tmpDir = self.tmpDir
        
        # Skip if db already exists
        if os.path.isfile(dbname) or os.path.isfile(dbname + ".dbtype"):
            logString = f"# Skipping '{dbname}' linclust clustering..."
            return logString
        
        # Convert to WSL paths where needed
        if platform.system() == "Windows":
            fasta = convert_windows_to_wsl_path(fasta)
            dbname = convert_windows_to_wsl_path(dbname)
            tmpDir = convert_windows_to_wsl_path(tmpDir)
        
        # Format command
        cmd = base_subprocess_cmd(self.mmDB.mmseqsExe)
        cmd += ["linclust", f"{fasta}_seqDB", dbname, tmpDir,
               "--min-seq-id", str(self.identity), "-c", str(self.cov_pct),
               "-e", str(self.evalue), "--cluster-mode", str(self.clust_mode),
               "--threads", str(self.threads)]
        
        # Clustering
        logString = "# Running linclust with: " + " ".join(cmd)
        if platform.system() != "Windows":
            run_linclust = subprocess.Popen(" ".join(cmd), shell = True,
                                         stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        else:
            run_linclust = subprocess.Popen(cmd, shell = True,
                                         stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        linclustout, linclusterr = run_linclust.communicate()
        if linclusterr.decode("utf-8") != '':
            raise Exception('Linclust error text below\n' +
                            linclusterr.decode("utf-8"))
        
        return logString
    
    def tabulate(self, outputFileName):
        '''
        Tabulates a linclust database file.
        
        Parameters:
            outputFileName -- a string indicating the file name to write TSV formatted
                              clustering results to.
        '''
        # Get hash for parameters combination
        strForHash = str(self.evalue) + str(self.identity) + str(self.cov_pct) + \
            str(self.clust_mode)
        strHash = sha256(bytes(strForHash, 'utf-8')).hexdigest()
        
        # Specify file locations
        seqdbname = os.path.abspath(f"{self.mmDB.fasta}_seqDB")
        clustdbname = os.path.abspath(f"{self.mmDB.fasta}_{strHash}_linclustDB")
        
        # Skip if table already exists
        if os.path.isfile(outputFileName):
            logString = f"# Skipping '{outputFileName}' linclust table generation..."
            return logString
        
        # Convert to WSL paths where needed
        if platform.system() == "Windows":
            seqdbname = convert_windows_to_wsl_path(seqdbname)
            clustdbname = convert_windows_to_wsl_path(clustdbname)
            outputFileName = convert_windows_to_wsl_path(outputFileName)
        
        # Format command
        cmd = base_subprocess_cmd(self.mmDB.mmseqsExe)
        cmd += ["createtsv", seqdbname, seqdbname, clustdbname, outputFileName]
        
        # Tabulation
        logString = "# Running table generation with: " + " ".join(cmd)
        if platform.system() != "Windows":
            run_tabulate = subprocess.Popen(" ".join(cmd), shell = True,
                                         stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        else:
            run_tabulate = subprocess.Popen(cmd, shell = True,
                                         stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        tableout, tableerr = run_tabulate.communicate()
        if tableerr.decode("utf-8") != '':
            raise Exception('Linclust tabulation text below\n' +
                            tableerr.decode("utf-8"))
        
        return logString
    
    def clean_all(self):
        '''
        Function to invoke after performing Linclust, the results of which
        are no longerwanted. It should clean up all files with _linclustDB* suffix.
        '''
        # Get hash for parameters combination
        strForHash = str(self.evalue) + str(self.identity) + str(self.cov_pct) + \
            str(self.clust_mode)
        strHash = sha256(bytes(strForHash, 'utf-8')).hexdigest()
        
        # Specify file prefixes
        dbPrefix = os.path.basename(f"{self.mmDB.fasta}_{strHash}_linclustDB")
        dbDir = os.path.dirname(os.path.abspath(self.mmDB.fasta))
        
        # Locate and delete files
        for file in os.listdir(dbDir):
            if file.startswith(dbPrefix):
                os.unlink(os.path.join(dbDir, file))

class MM_Cascade(MM_Clust):
    '''
    The MM_Cascade Class provides the logic for running cascaded MMseqs2 clustering
    with an MM_DB object as input. The MMseqs executable location and tmpDir will be pulled
    from the mmDB input; a new tmpDir for running cascaded clusering can be specified if desired.
    
    Attributes:
        super class attributes -- see docstring for MM_Clust.
        sensitivity -- a float in the list [1,2,3,4,5,5.7,6,7,7.5]; default == 4.0,
                       but recommended to use 5.7 or 7.5.
        cluster_steps -- an int minimally bounded at 1 for how many cascaded clustering
                         steps to run; default == 3, recommended to keep that.
    '''
    def __init__(self, mmDB, evalue=1e-3, identity=0.9,
                 cov_pct=0.8, clust_mode="set-cover",
                 threads=1, tmpDir=None,
                 sensitivity=4.0, cluster_steps=3):
        
        super().__init__(mmDB, evalue, identity, cov_pct,
                        clust_mode, threads, tmpDir)
        
        self.sensitivity = sensitivity
        self.cluster_steps = cluster_steps
    
    @property
    def sensitivity(self):
        return self._sensitivity
    
    @sensitivity.setter
    def sensitivity(self, value):
        assert value in [1,2,3,4,5,5.7,6,7,7.5], \
            "sensitivity must be a value in the list [1,2,3,4,5,5.7,6,7,7.5]"
        
        self._sensitivity = value
    
    @property
    def cluster_steps(self):
        return self._cluster_steps
    
    @cluster_steps.setter
    def cluster_steps(self, value):
        assert isinstance(value, int)
        assert 0 < value, "cluster_steps must be a positive integer"
        
        self._cluster_steps = value
    
    def cluster(self):
        '''
        Runs the cascaded clustering process on the given mmDB.
        '''
        
        # Run DB generation & indexing if relevant
        self.mmDB.generate()
        self.mmDB.index()
        
        # Get hash for parameters combination
        strForHash = str(self.evalue) + str(self.identity) + str(self.cov_pct) + \
            str(self.clust_mode) + str(self.sensitivity)
        strHash = sha256(bytes(strForHash, 'utf-8')).hexdigest()
        
        # Specify file locations
        fasta = self.mmDB.fasta
        dbname = os.path.abspath(f"{self.mmDB.fasta}_{strHash}_clustDB")
        tmpDir = self.tmpDir
        
        # Skip if db already exists
        if os.path.isfile(dbname) or os.path.isfile(dbname + ".dbtype"):
            logString = f"# Skipping '{dbname}' cascaded clustering..."
            return logString
        
        # Convert to WSL paths where needed
        if platform.system() == "Windows":
            fasta = convert_windows_to_wsl_path(fasta)
            dbname = convert_windows_to_wsl_path(dbname)
            tmpDir = convert_windows_to_wsl_path(tmpDir)
        
        # Format command
        cmd = base_subprocess_cmd(self.mmDB.mmseqsExe)
        cmd += ["cluster", f"{fasta}_seqDB", dbname, tmpDir,
               "--min-seq-id", str(self.identity), "-c", str(self.cov_pct),
               "-e", str(self.evalue), "--cluster-mode", str(self.clust_mode),
               "-s", str(self.sensitivity), "--cluster-steps", str(self.cluster_steps),
               "--threads", str(self.threads)]
        
        # Clustering
        logString = "# Running cascaded clustering with: " + " ".join(cmd)
        if platform.system() != "Windows":
            run_cluster = subprocess.Popen(" ".join(cmd), shell = True,
                                         stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        else:
            run_cluster = subprocess.Popen(cmd, shell = True,
                                         stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        clustout, clusterr = run_cluster.communicate()
        if clusterr.decode("utf-8") != '':
            raise Exception('Cascaded clustering error text below\n' +
                            clusterr.decode("utf-8"))
        
        return logString
    
    def tabulate(self, outputFileName):
        '''
        Tabulates a cascaded database file.
        
        Parameters:
            outputFileName -- a string indicating the file name to write TSV formatted
                              clustering results to.
        '''
        # Get hash for parameters combination
        strForHash = str(self.evalue) + str(self.identity) + str(self.cov_pct) + \
            str(self.clust_mode) + str(self.sensitivity)
        strHash = sha256(bytes(strForHash, 'utf-8')).hexdigest()
        
        # Specify file locations
        seqdbname = os.path.abspath(f"{self.mmDB.fasta}_seqDB")
        clustdbname = os.path.abspath(f"{self.mmDB.fasta}_{strHash}_clustDB")
        
        # Skip if table already exists
        if os.path.isfile(outputFileName):
            logString = f"# Skipping '{outputFileName}' cascaded table generation..."
            return logString
        
        # Convert to WSL paths where needed
        if platform.system() == "Windows":
            seqdbname = convert_windows_to_wsl_path(seqdbname)
            clustdbname = convert_windows_to_wsl_path(clustdbname)
            outputFileName = convert_windows_to_wsl_path(outputFileName)
        
        # Format command
        cmd = base_subprocess_cmd(self.mmDB.mmseqsExe)
        cmd += ["createtsv", seqdbname, seqdbname, clustdbname, outputFileName]
        
        # Tabulation
        logString = "# Running table generation with: " + " ".join(cmd)
        if platform.system() != "Windows":
            run_tabulate = subprocess.Popen(" ".join(cmd), shell = True,
                                         stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        else:
            run_tabulate = subprocess.Popen(cmd, shell = True,
                                         stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        tableout, tableerr = run_tabulate.communicate()
        if tableerr.decode("utf-8") != '':
            raise Exception('Cascaded clustering tabulation text below\n' +
                            tableerr.decode("utf-8"))
        
        return logString
    
    def clean_all(self):
        '''
        Function to invoke after performing Linclust, the results of which
        are no longerwanted. It should clean up all files with _clustDB* suffix.
        '''
        # Get hash for parameters combination
        strForHash = str(self.evalue) + str(self.identity) + str(self.cov_pct) + \
            str(self.clust_mode) + str(self.sensitivity)
        strHash = sha256(bytes(strForHash, 'utf-8')).hexdigest()
        
        # Specify file prefixes
        dbPrefix = os.path.basename(f"{self.mmDB.fasta}_{strHash}_clustDB")
        dbDir = os.path.dirname(os.path.abspath(self.mmDB.fasta))
        
        # Locate and delete files
        for file in os.listdir(dbDir):
            if file.startswith(dbPrefix):
                os.unlink(os.path.join(dbDir, file))

if __name__ == "__main__":
    pass

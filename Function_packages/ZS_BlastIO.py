#! python3
# ZS_BlastIO.py
# Specifies the BLAST Class for performing BLAST searches
# using ZS_SeqIO.FASTA and ZS_SeqIO.FastASeq objects, as well
# as for parsing outfmt6 results.

import os, sys, subprocess, time, random, platform
from pathlib import Path
from hashlib import sha256

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from ZS_Utility import base_subprocess_cmd, convert_windows_to_wsl_path, \
    tmp_file_name_gen, get_codec
from ZS_SeqIO import Conversion

class BLAST:
    '''
    The BLAST Class provides easy access to several BLAST functions to perform
    searches of FASTAs against each other. Object values can be set to dictate
    what type of BLAST to perform, including any E-value cut-offs amongst other
    things. It's primary purpose is the .get_blast_results() method.
    
    Attributes:
        query (REQUIRED) -- a string, FASTA, or FastASeq object from the ZS_SeqIO suite.
        target (REQUIRED) -- a string, FASTA, or FastASeq object from the ZS_SeqIO suite.
        blastAlgorithm (REQUIRED) -- a string indicating what BLAST algorithm to run. Valid
                                     values include blastp, blastn, tblastn, and tblastx.
        evalue (OPTIONAL) -- a float or integer value specifying the E-value cut-off to enforce
                             when obtaining BLAST results. Defaults to 1e-5.
        threads (OPTIONAL) -- an integer value specifying the number of threads to run when
                              performing BLAST search. Defaults to 1.
        clean (OPTIONAL) -- a Boolean indicating whether to clean up after BLAST search is
                            complete (i.e., delete the result file) or to keep it. Defaults
                            to True, which means we will DELETE the results file after it
                            has been parsed.
    '''
    def __init__(self, query, target, blastAlgorithm):
        # Validate input types
        assert type(query).__name__ == "str" \
            or type(query).__name__ == "FASTA" \
            or type(query).__name__ == "ZS_SeqIO.FASTA" \
            or type(query).__name__ == "FastASeq" \
            or type(query).__name__ == "ZS_SeqIO.FastASeq"
        assert type(target).__name__ == "str" \
            or type(target).__name__ == "FASTA" \
            or type(target).__name__ == "ZS_SeqIO.FASTA" \
            or type(target).__name__ == "FastASeq" \
            or type(target).__name__ == "ZS_SeqIO.FastASeq"
        assert isinstance(blastAlgorithm, str)
        
        # Validate that inputs exist if specified as a string (file location)
        if type(query).__name__ == "str" and not os.path.isfile(query):
            raise Exception("Query parameter is a string, but does not point to an existing file location")
        if type(target).__name__ == "str" and not os.path.isfile(target):
            raise Exception("Target parameter is a string, but does not point to an existing file location")
        self.query = query
        self.target = target
        
        # Set default attributes
        self.set_blastAlgorithm(blastAlgorithm)
        self.evalue = 1e-5
        self.threads = 1
        self.clean = True
    
    def set_blastAlgorithm(self, algorithm):
        '''
        This method allows the blastAlgorithm parameter to be set, which is
        used to dictate what type of BLAST to run. It will also be interpretted
        when creating a BLAST database from the target (if needed).
        
        Parameters:
            algorithm -- a string with value contained in the list
                         ['blastp', 'blastn', 'tblastn', 'tblastx']
        '''
        assert algorithm.lower() in ['blastp', 'blastn', 'tblastn', 'tblastx'], "Algorithm value \"{0}\" is not recognised".format(algorithm.lower())
        
        self.blastAlgorithm = algorithm.lower()
    
    def set_evalue(self, num):
        '''
        This method allows the E-value threshold to be set.
        
        Parameters:
            num -- should be a valid float or integer greater than 0.
        '''
        assert isinstance(num, int) or isinstance(num, float)
        if num < 0:
            raise Exception("E-value cannot be less than 0")
        
        self.evalue = num
    
    def set_threads(self, num):
        '''
        This method allows the use of multithreading by BLAST
        
        Parameters:
            num -- should be a valid integer greater than 0. Be sensible
                   since there's no upper limit checking.
        '''
        assert isinstance(num, int)
        if num < 0:
            raise Exception("Number of threads must be more than 0")
        
        self.threads = num
    
    def set_clean(self, clean):
        '''
        This method allows the clean attribute to be set, which controls
        whether BLAST output files are kept or not.
        
        Parameters:
            clean -- a Boolean of True or False
        '''
        assert isinstance(clean, bool)
        
        self.clean = clean

    def blastdb_exists(self, fastaFile):
        '''
        Relies on a simple assumption that a .nsq or .psq file's presence
        indicates that a BLAST database was successfully created
        from the FASTA file.
        
        This method behaves statically to allow the .target value
        to come in multiple different data types.
        
        Parameters:
            fastaFile -- a string indicating the full path to the
                         FASTA file whose database existence we are
                         querying.
        Returns:
            dbExists -- a Boolean where True means the database exists,
                        and False means it does not exist.
        '''
        return os.path.isfile("{0}.nsq".format(fastaFile)) or os.path.isfile("{0}.psq".format(fastaFile))
    
    def makeblastdb(self, fastaFile):
        '''
        Makes a database out of the provided value for use in BLAST search.
        
        This method behaves statically to allow the .target value
        to come in multiple different data types.
        
        Parameters:
            fastaFile -- a string indicating the full path to the
                         FASTA file we want to create a database for.
        '''
        assert isinstance(fastaFile, str)
        
        # Coerce algorithm into molecule type
        if self.blastAlgorithm.lower() == "blastp":
            dbType = "prot"
        elif self.blastAlgorithm.lower() == "blastn":
            dbType = "nucl"
        elif self.blastAlgorithm.lower() == "tblastn":
            dbType = "nucl"
        elif self.blastAlgorithm.lower() == "tblastx": # I think?
            dbType = "prot"
        else:
            raise Exception("blastAlgorithm attribute \"{0}\" is not recognised by makeblastdb method".format(self.blastAlgorithm.lower()))

        # Format and run command
        cmd = 'makeblastdb -in "{0}" -dbtype {1} -out "{0}"'.format(fastaFile, dbType)
        run_makedb = subprocess.Popen(cmd, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        makedbout, makedberr = run_makedb.communicate()
        if makedberr.decode("utf-8") != '':
                raise Exception('Makeblastdb error text below\n' + makedberr.decode("utf-8")) 
    
    def blast(self, queryFasta, dbFastaFile, outFile):
        '''
        Performs the BLAST operation.
        
        This method does not use .query and .target because they may be in the FASTA or
        FastASeq object types, rather than a string which we require here.
        
        Parameters:
            queryFasta -- a string indicating the location of a FASTA file to use as query.
            dbFastaFile -- a string indicating the location of a FASTA file to use as the database.
            outFile -- a string indicating the location to write the output file to. File must not already exist!
        '''
        assert isinstance(queryFasta, str)
        assert isinstance(dbFastaFile, str)
        assert self.blastdb_exists(dbFastaFile), "BLAST database must already exist for \"{0}\"!".format(dbFastaFile)
        assert isinstance(outFile, str)
        assert not os.path.isfile(outFile), "\"{0}\" already exists; blast method won't overwrite it".format(outFile)
        
        cmd = '{0} -query "{1}" -db "{2}" -num_threads {3} -evalue {4} -out "{5}" -outfmt 6'.format(
            self.blastAlgorithm.lower(), queryFasta, dbFastaFile, self.threads, self.evalue, outFile
        )
        
        run_blast = subprocess.Popen(cmd, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        blastout, blasterr = run_blast.communicate()
        if blasterr.decode("utf-8") != '':
            raise Exception('BLAST error text below\n' + blasterr.decode("utf-8"))
    
    def get_blast_results(self):
        '''
        This function pipelines the process of obtaining BLAST results. Intermediate files are
        deleted automatically, and hence this function will only result in the return of the
        blastDict object. However, if self.clean is False, we will leave the output BLAST file
        (but still clean up any temporary FASTA files).
        
        Returns:
            blastDict -- a dict with structure:
                query_id: [[target_id, identity_pct, query_start, query_end, target_start, target_end, evalue], ...]
            blastResultFile -- a string indicating the file name of the results file. If self.clean is True,
                               this will instead return None.
        '''
        # Get file names for query and target after data type coercion
        q, qIsTemporary = Conversion.get_filename_for_input_sequences(self.query)
        t, tIsTemporary = Conversion.get_filename_for_input_sequences(self.target)
        
        # Make sure target file is ready for BLAST
        if not self.blastdb_exists(t):
            self.makeblastdb(t)
        
        # Get hash for temporary file creation
        tmpHash = sha256(bytes(q + t + str(time.time()) + str(random.randint(0, 100000)), 'utf-8') ).hexdigest()
        
        # Run BLAST
        tmpResultName = tmp_file_name_gen("{0}.vs.{1}_{2}".format(os.path.basename(q).rsplit(".", maxsplit=1)[0], os.path.basename(t).rsplit(".", maxsplit=1)[0], tmpHash[0:20]), "outfmt6")
        self.blast(q, t, tmpResultName)
        
        # Parse BLAST results
        blastResults = BLAST_Results(tmpResultName)
        blastResults.evalue = self.evalue
        blastResults.parse_blast_hit_coords()
        blastDict = blastResults.results
        
        # Clean up q and t temporary files
        if qIsTemporary:
            os.unlink(q)
            qDir = os.path.dirname(q)
            for file in os.listdir(qDir if qDir != "" else "."): # kill off any temporary database files too
                file = os.path.join(qDir, file)
                if file.startswith(q):
                    os.unlink(file)
        if tIsTemporary:
            os.unlink(t)
            tDir = os.path.dirname(t)
            for file in os.listdir(tDir if tDir != "" else "."): # kill off any temporary database files too
                file = os.path.join(tDir, file)
                if file.startswith(t):
                    os.unlink(file)
        
        # Clean up results and return (if relevant)
        if self.clean:
            os.unlink(tmpResultName)
            return blastDict, None
        # Or just return results
        else:
            return blastDict, tmpResultName

class BLAST_Results:
    '''
    The BLAST_Results Class provides parsing capability for outfmt6 files.
    
    Attributes:
        file -- a string indicating the location of the outfmt6 file that was
                parsed to generate this object instance
        evalue -- an int or float value indicating the cutoff to enforce for
                  retrieving hits for a query
        num_hits -- an int value which can be -1 (to retrieve all hits) or a value
                    > 0 to get that many hits per query
    '''
    def __init__(self, outfmt6File):
        # Validate input types
        assert isinstance(outfmt6File, str)
        
        # Validate that input file location exists
        if not os.path.isfile(outfmt6File):
            raise Exception(f"'{outfmt6File} does not point to an existing file location")
        self.file = outfmt6File
        
        # Set default property values
        self.evalue = 1e-5
        self.num_hits = -1
        
        # Also set helper attribute
        self.isBLAST_Results = True
    
    @property
    def evalue(self):
        return self._evalue
    
    @evalue.setter
    def evalue(self, value):
        assert isinstance(value, float) or isinstance(value, int), \
            "evalue can only be set to a float or integer value"
        
        if value < 0:
            raise ValueError("evalue must be a value >= 0")
        
        self._evalue = value
    
    @property
    def num_hits(self):
        return self._num_hits
    
    @num_hits.setter
    def num_hits(self, value):
        assert isinstance(value, int), \
            "num_hits must be an integer value"
        
        if value == 0:
            raise ValueError("num_hits cannot be 0 (which would retrieve no hits)")
        if value < -1:
            raise ValueError("num_hits cannot be less than -1 (which retrieves all hits)")
        
        self._num_hits = value
    
    def parse_blast_hit_coords(self):
        '''
        Parameters:
            self.file -- a string indicating the location of a BLAST results file in
                         outfmt6 format.
        Sets:
            self.results -- a dict with structure like:
                                query_id1: [
                                    [target_id, identity_pct, query_start, query_end, target_start, target_end, evalue],
                                    ...
                                ],
                                query_id2: [ ... ],
                                ...
        '''
        blastDict = {}
        
        with open(self.file, "r", encoding=get_codec(self.file)) as fileIn:
            for line in fileIn:
                # Extract details
                sl = line.rstrip("\r\n").split('\t')
                qid = sl[0]
                tid = sl[1]
                identityPct = float(sl[2])
                qstart = int(sl[6])
                qend = int(sl[7])
                tstart = int(sl[8])
                tend = int(sl[9])
                evalue = float(sl[10])
                bitscore = float(sl[11])
                
                # Skip if evalue isn't significant
                if evalue > self.evalue: # filter here since self.evalue might differ between BLAST run and parsing now
                    continue
                
                # Store result
                if qid not in blastDict:
                    blastDict[qid] = [[tid, identityPct, qstart, qend, tstart, tend, evalue, bitscore]]
                else:
                    blastDict[qid].append([tid, identityPct, qstart, qend, tstart, tend, evalue, bitscore])
        
        # Sort individual entries in blastDict
        for value in blastDict.values():
            value.sort(key = lambda x: (x[6], x[7])) # sort by evalue and bitscore
            for v in value:
                del v[7] # drop the bitscore since we don't need it after sorting
        
        # Enforce num_hits threshold
        if self.num_hits != -1:
            for key in blastDict.keys():
                blastDict[key] = blastDict[key][0:self.num_hits]
        
        self.results = blastDict
    
    def __iter__(self):
        return iter(self.results)
    
    def __len__(self):
        return len(self.results)
    
    def __getitem__(self, key):
        if key in self.results:
            return self.results[key]
    
    def __contains__(self, item):
        return True if self[item] is not None else False
    
    def __str__(self):
        return "BLAST_Results; parsed '{0}' file, contains results for {1} queries".format(
            self.file, len(self)
        )
    
    def __repr__(self):
        return "BLAST_Results(outfmt6='{0}', queries={1}, evalue={2}, num_hits={3})".format(
            self.file, len(self), self.evalue, self.num_hits
        )

class MM_DB:
    '''
    The MM_DB Class provides the logic for creating and handling MMSeqs2 databases.
    
    Attributes:
        fasta (REQUIRED) -- a string, FASTA, or FastASeq object from the ZS_SeqIO suite.
        mmseqsDir (REQUIRED) -- a string pointing to the location where the mmseqs executable
                                is found.
        tmpDir (REQUIRED) -- a string location for where MMSeqs2 should keep temp files.
        threads (OPTIONAL) -- a positive integer for how many threads to use when indexing
                              a database file.
    '''
    def __init__(self, fasta, mmseqsDir, tmpDir, threads=1):
        self.fasta = fasta
        self.mmseqsDir = mmseqsDir
        self.tmpDir = tmpDir
        self.threads = threads
        
        self.isMM_DB = True # flag to check object type
    
    @property
    def fasta(self):
        return self._fasta
    
    @fasta.setter
    def fasta(self, value):
        assert type(value).__name__ == "str" \
            or type(value).__name__ == "FASTA" \
            or type(value).__name__ == "ZS_SeqIO.FASTA" \
            or type(value).__name__ == "FastASeq" \
            or type(value).__name__ == "ZS_SeqIO.FastASeq"
        if type(value).__name__ == "str" and not os.path.isfile(value):
            raise Exception(("Fasta parameter is a string, but does not point " + 
                            "to an existing file location"))
        
        self._fasta = value
    
    @property
    def mmseqsDir(self):
        return self._mmseqsDir
    
    @mmseqsDir.setter
    def mmseqsDir(self, value):
        assert type(value).__name__ == "str" \
            or isinstance(value, Path)
        if not os.path.isdir(value):
            raise Exception(("mmseqsDir does not point to an existing directory"))
        
        self._mmseqsDir = value
        self.mmseqsExe = os.path.join(value, "mmseqs")
    
    @property
    def mmseqsExe(self):
        return self._mmseqsExe
    
    @mmseqsExe.setter
    def mmseqsExe(self, value):
        assert type(value).__name__ == "str" \
            or isinstance(value, Path)
        if not os.path.isfile(value):
            raise Exception(("mmseqsExe does not point to an existing file"))
        
        self._mmseqsExe = value
    
    @property
    def tmpDir(self):
        return self._tmpDir
    
    @tmpDir.setter
    def tmpDir(self, value):
        assert type(value).__name__ == "str" \
            or isinstance(value, Path)
        value = os.path.abspath(value)
        
        if not os.path.isdir(os.path.dirname(value)):
            raise Exception((f"tmpDir's parent location ('{os.path.dirname(value)}') " +
                             "does not exist"))
        
        if not os.path.isdir(value):
            os.mkdir(value)
        
        self._tmpDir = value
    
    @property
    def threads(self):
        return self._threads
    
    @threads.setter
    def threads(self, value):
        assert isinstance(value, int)
        assert 0 < value, "threads must be a positive integer"
        
        self._threads = value
    
    def clean_all(self):
        '''
        Function to invoke after using a MMSeqs2 database which is no longer
        wanted. It should clean up all files with _seqDB* suffix.
        '''
        dbPrefix = os.path.basename(f"{self.fasta}_seqDB")
        dbDir = os.path.dirname(os.path.abspath(self.fasta))
        
        for file in os.listdir(dbDir):
            if file.startswith(dbPrefix):
                os.unlink(os.path.join(dbDir, file))
    
    def generate(self):
        '''
        Creates a sequence DB for MMSeqs2 if it does not already exist.
        '''
        # Specify file locations
        dbname = os.path.abspath(f"{self.fasta}_seqDB")
        fasta = self.fasta
        
        # Skip if db already exists
        if os.path.isfile(dbname):
            logString = f"# Skipping '{dbname}' DB generation..."
            return logString
        
        # Convert to WSL paths where needed
        if platform.system() == "Windows":
            fasta = convert_windows_to_wsl_path(fasta)
            dbname = convert_windows_to_wsl_path(dbname)
        
        # Format command
        cmd = base_subprocess_cmd(self.mmseqsExe)
        cmd += ["createdb", fasta, dbname]
        
        # DB generation
        logString = "# DB generation with: " + " ".join(cmd)
        if platform.system() != "Windows":
            run_makedb = subprocess.Popen(" ".join(cmd), shell = True,
                                         stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        else:
            run_makedb = subprocess.Popen(cmd, shell = True,
                                         stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        makedbout, makedberr = run_makedb.communicate()
        if makedberr.decode("utf-8") != '':
            raise Exception('Make MMseqs2 query db error text below\n' +
                            makedberr.decode("utf-8"))
        
        return logString
    
    def index(self):
        '''
        Indexes as sequence DB for MMSeqs2 if it has not already been done.
        '''
        # Specify file locations
        dbname = os.path.abspath(f"{self.fasta}_seqDB")
        tmpDir = self.tmpDir
        
        # Skip if index already exists
        if MM_DB.mms2_index_exists(os.path.basename(dbname), os.path.dirname(dbname)):
            logString = f"# Skipping '{dbname}' DB indexing..."
            return logString
        
        # Convert to WSL paths where needed
        origDbname = dbname # retain the original value for output flag creation later
        if platform.system() == "Windows":
            dbname = convert_windows_to_wsl_path(dbname)
            tmpDir = convert_windows_to_wsl_path(tmpDir)
        
        # Format command
        cmd = base_subprocess_cmd(self.mmseqsExe)
        cmd += ["createindex", dbname, tmpDir, "--threads", str(self.threads)]
        
        # Run query index
        logString = "# DB indexing with: " + " ".join(cmd)
        if platform.system() != "Windows":
            run_index = subprocess.Popen(" ".join(cmd), shell = True,
                                         stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        else:
            run_index = subprocess.Popen(cmd, shell = True,
                                         stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        indexout, indexerr = run_index.communicate()
        if indexerr.decode("utf-8") != '':
                raise Exception('Indexing MMseqs2 query db error text below\n' + indexerr.decode("utf-8"))
        
        # Create output flag
        with open(origDbname + ".idxComplete", "w") as fileOut: # need orig name if WSL change occurred
            pass # create empty file
        
        return logString
    
    @staticmethod
    def mms2_index_exists(fileNamePrefix, directory):
        indexName = fileNamePrefix + '.idxComplete'
        if indexName in os.listdir(directory):
            return True
        else:
            return False

if __name__ == "__main__":
    pass

#! python3
# ZS_BlastIO.py
# Specifies the BLAST Class for performing BLAST searches
# using ZS_SeqIO.FASTA and ZS_SeqIO.FastASeq objects.

import os, sys, subprocess, hashlib, time, random

sys.path.append(os.path.dirname(__file__))
from ZS_SeqIO import FASTA

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

    def parse_blast_hit_coords(self, outfmt6File):
        '''
        Parameters:
            outfmt6File -- a string indicating the location of a BLAST results file in
                           outfmt6 format.
        Returns:
            blastDict -- a dict with structure:
                query_id: [[target_id, identity_pct, query_start, query_end, target_start, target_end, evalue], ...]
        '''
        blastDict = {}
        with open(outfmt6File, 'r') as fileIn:
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
        
        return blastDict
    
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
        q, qIsTemporary = self._get_filename_for_query_or_target(self.query)
        t, tIsTemporary = self._get_filename_for_query_or_target(self.target)
        
        # Make sure target file is ready for BLAST
        if not self.blastdb_exists(t):
            self.makeblastdb(t)
        
        # Get hash for temporary file creation
        tmpHash = hashlib.sha256(bytes(q + t + str(time.time()) + str(random.randint(0, 100000)), 'utf-8') ).hexdigest()
        
        # Run BLAST
        tmpResultName = self._tmp_file_name_gen("{0}.vs.{1}_{2}".format(os.path.basename(q).rsplit(".", maxsplit=1)[0], os.path.basename(t).rsplit(".", maxsplit=1)[0], tmpHash[0:20]), "outfmt6")
        self.blast(q, t, tmpResultName)
        
        # Parse BLAST results
        blastDict = self.parse_blast_hit_coords(tmpResultName)
        
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
    
    def _get_filename_for_query_or_target(self, qt):
        '''
        Hidden method for use when boiling down one of the three data types (FASTA, FastASeq, and string)
        into a string representing the query file name.
        
        If it's already a string, this does nothing. Otherwise, it will make sure a file exists
        with the FASTA data contents for use by BLAST.
        
        Parameters:
            qt -- the value of .query or .target, provided to this method. We won't
                  call this from self since this method (to reduce code repetition)
                  needs to work for both.
        Returns:
            qt -- a string indicating the file name for the input value of qt. If it was
                  already a string, this will be the same. Otherwise, it will be a file
                  name containing the contents of the FASTA or FastASeq object.
            isTemporary -- a Boolean indicating whether the returned qt value has been
                           created by this method as a temporary file (True) or if it was
                           already existing (False, i.e., qt was already a string)
        '''
        # Get a hash for temporary file creation
        hashForTmp = qt if isinstance(qt, str) \
                       else qt.fileOrder[0][0] if type(qt).__name__ == "ZS_SeqIO.FASTA" or type(qt).__name__ == "FASTA" \
                       else qt.id
        tmpHash = hashlib.sha256(bytes(hashForTmp + str(time.time()) + str(random.randint(0, 100000)), 'utf-8') ).hexdigest()
        
        # If qt is a FastASeq, make it a FASTA object
        if type(qt).__name__ == "FastASeq" or type(qt).__name__ == "ZS_SeqIO.FastASeq":
            tmpFastaName = self._tmp_file_name_gen("qt_fasta_tmp" + tmpHash[0:20], "fasta")
            with open(tmpFastaName, "w") as fileOut:
                fileOut.write(">{0}\n{1}\n".format(qt.id, qt.seq))
            qt = FASTA(tmpFastaName)
            os.unlink(tmpFastaName)
        
        # If qt is a FASTA, make it into a file
        isTemporary = False
        if type(qt).__name__ == "FASTA" or type(qt).__name__ == "ZS_SeqIO.FASTA":
            tmpQtName = self._tmp_file_name_gen("blast_tmp" + tmpHash[0:20], "fasta")
            qt.write(tmpQtName)
            qt = tmpQtName # after this point, qt will be a string indicating a FASTA file name
            isTemporary = True # if we set this, qt was not originally a string
        
        return qt, isTemporary
    
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

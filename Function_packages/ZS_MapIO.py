#! python3
# ZS_MapIO.py
# Specifies the GMAP Class for performing GMAP and read alignments.

import os, sys, subprocess, platform
from pathlib import Path

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from ZS_Utility import base_subprocess_cmd, convert_windows_to_wsl_path

class GMAP_DB:
    '''
    The GMAP_DB Class encapsulates the logic of indexing a FASTA file using GMAP.
    
    Attributes:
        fasta (REQUIRED) -- a string indicating the location of a FASTA file.
        gmapDir (REQUIRED) -- a string indicating the location of GMAP binaries.
    '''
    def __init__(self, fasta, gmapDir):
        self.fasta = fasta
        self.gmapDir = gmapDir
        
        self.isGMAP_DB = True # flag to check object type
    
    @property
    def fasta(self):
        return self._fasta
    
    @fasta.setter
    def fasta(self, value):
        assert type(value).__name__ == "str"
        if not os.path.isfile(value):
            raise Exception(("Fasta parameter is a string, but does not point " + 
                            "to an existing file location"))
        
        self._fasta = value
    
    @property
    def gmapDir(self):
        return self._gmapDir
    
    @gmapDir.setter
    def gmapDir(self, value):
        assert type(value).__name__ == "str" \
            or isinstance(value, Path)
        if not os.path.isdir(value):
            raise Exception(("gmapDir does not point to an existing directory"))
        
        self._gmapDir = value
        self.buildExe = os.path.join(value, "gmap_build")
    
    @property
    def buildExe(self):
        return self._buildExe
    
    @buildExe.setter
    def buildExe(self, value):
        assert type(value).__name__ == "str" \
            or isinstance(value, Path)
        if not os.path.isfile(value):
            raise Exception(("buildExe does not point to an existing file"))
        
        self._buildExe = value
    
    def index_exists(self):
        '''
        Relies on a simple assumption that a .gmap directory's presence
        indicates that a GMAP database was successfully created
        from the FASTA file.
        
        Returns:
            dbExists -- a Boolean where True means the database exists,
                        and False means it does not exist.
        '''
        return os.path.isdir(f"{self.fasta}.gmap")
    
    def index(self):
        '''
        Makes a database out of the .fasta value for use in BLAST search.
        '''
        # Skip if index exists
        if self.index_exists():
            raise FileExistsError("GMAP index already exists!")
        
        # Convert to WSL paths where needed
        if platform.system() == "Windows":
            fasta = convert_windows_to_wsl_path(self.fasta)
        else:
            fasta = self.fasta
        
        # Format command
        cmd = base_subprocess_cmd(self.buildExe)
        cmd += ["-D", os.path.dirname(fasta), "-d", f"{os.path.basename(fasta)}.gmap", fasta]
        
        # Run indexing
        if platform.system() != "Windows":
            run_index = subprocess.Popen(" ".join(cmd), shell = True,
                                         stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        else:
            run_index = subprocess.Popen(cmd, shell = True,
                                         stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        indexout, indexerr = run_index.communicate()
        
        # Raise exception if final part of stderr from successful run isn't found
        if "Writing localdb sarrays" not in indexerr.decode("utf-8"):
            raise Exception('GMAP indexing error text below\n' + indexerr.decode("utf-8"))

class GMAP:
    '''
    The GMAP Class provides easy access to the GMAP search function to perform
    mapping of FASTAs against each other. Object values can be set to control
    the parameters specified for a GMAP search.
    
    Attributes:
        query (REQUIRED) -- a string indicating the location of a FASTA file.
        target (REQUIRED) -- a string indicating the location of a FASTA file.
        gmapDir (REQUIRED) -- a string indicating the location of GMAP binaries.
        threads (OPTIONAL) -- an integer value specifying the number of threads to run when
                              performing BLAST search. Defaults to 1.
    '''
    def __init__(self, query, target, gmapDir, threads=1):
        self.query = query
        self.gmapDir = gmapDir
        self.target = target # sets self.db
        self.threads = threads
        
        self.isGMAP = True # flag to check object type
        
        # Set default values
        self.outputFormat = 2 # my default, not GMAPs
        self.npaths = 5
        self.chimeraMargin = 30
        self.batch = 2
        self.maxIntronLengthMiddle = 500000
        self.maxIntronLengthEnds = 10000
    
    @property
    def query(self):
        return self._query
    
    @query.setter
    def query(self, value):
        assert type(value).__name__ == "str"
        if not os.path.isfile(value):
            raise Exception(("Query parameter is a string, but does not point " + 
                            "to an existing file location"))
        
        self._query = value
    
    @property
    def target(self):
        return self._target
    
    @target.setter
    def target(self, value):
        assert type(value).__name__ == "str"
        if not os.path.isfile(value):
            raise Exception(("Target parameter is a string, but does not point " + 
                            "to an existing file location"))
        
        self._target = value
        self.db = GMAP_DB(value, self.gmapDir)
    
    @property
    def gmapDir(self):
        return self._gmapDir
    
    @gmapDir.setter
    def gmapDir(self, value):
        assert type(value).__name__ == "str" \
            or isinstance(value, Path)
        if not os.path.isdir(value):
            raise Exception(("gmapDir does not point to an existing directory"))
        
        self._gmapDir = value
        self.gmapExe = os.path.join(value, "gmap")
    
    @property
    def gmapExe(self):
        return self._gmapExe
    
    @gmapExe.setter
    def gmapExe(self, value):
        assert type(value).__name__ == "str" \
            or isinstance(value, Path)
        if not os.path.isfile(value):
            raise Exception(("gmapExe does not point to an existing file"))
        
        self._gmapExe = value
    
    @property
    def threads(self):
        return self._threads
    
    @threads.setter
    def threads(self, value):
        assert isinstance(value, int)
        if value < 0:
            raise Exception("Number of threads must be more than 0")
        
        self._threads = value
    
    @property
    def outputFormat(self):
        return self._outputFormat
    
    @outputFormat.setter
    def outputFormat(self, value):
        assert isinstance(value, int)
        if value < 1 or value > 9:
            raise Exception("GMAP --format only supports output formats 1-9")
        
        self._outputFormat = value
    
    @property
    def npaths(self):
        return self._npaths
    
    @npaths.setter
    def npaths(self, value):
        assert isinstance(value, int)
        if value < 0:
            raise Exception("GMAP --npaths must be >= 0")
        
        self._npaths = value
    
    @property
    def chimeraMargin(self):
        return self._chimeraMargin
    
    @chimeraMargin.setter
    def chimeraMargin(self, value):
        assert isinstance(value, int)
        if value < 0:
            raise Exception("GMAP --chimera-margin must be >= 0")
        
        self._chimeraMargin = value
    
    @property
    def batch(self):
        return self._batch
    
    @batch.setter
    def batch(self, value):
        assert isinstance(value, int)
        if value < 0 or value > 5:
            raise Exception("GMAP --batch must be in the range 0-5")
        
        self._batch = value
    
    @property
    def maxIntronLengthMiddle(self):
        return self._maxIntronLengthMiddle
    
    @maxIntronLengthMiddle.setter
    def maxIntronLengthMiddle(self, value):
        assert isinstance(value, int)
        if value < 1:
            raise Exception("GMAP --max-intronlength-middle must be a positive integer")
        
        self._maxIntronLengthMiddle = value
    
    @property
    def maxIntronLengthEnds(self):
        return self._maxIntronLengthEnds
    
    @maxIntronLengthEnds.setter
    def maxIntronLengthEnds(self, value):
        assert isinstance(value, int)
        if value < 1:
            raise Exception("GMAP --max-intronlength-ends must be a positive integer")
        
        self._maxIntronLengthEnds = value
    
    def index_exists(self):
        '''
        Checks if the .target value is in a GMAP index/db yet.
        '''
        return self.db.index_exists()
    
    def index(self):
        '''
        Makes a database out of the .target value for use in BLAST search.
        '''
        # Skip if index exists
        if self.db.index_exists():
            raise FileExistsError("GMAP index already exists!")
        
        # Create index otherwise
        self.db.index()
    
    def gmap(self, outFile, force=False):
        '''
        Performs the GMAP search operation.
        
        This method does not use .query and .target because they may be in the FASTA or
        FastASeq object types, rather than a string which we require here.
        
        Parameters:
            outFile -- a string indicating the location to write the output file to. File must not already exist!
            force -- a boolean indicating whether to overwrite the output file if it already exists.
        '''
        # Skip if index does not exist
        if not self.db.index_exists():
            raise FileNotFoundError("GMAP index does not exist! Cannot perform search.")
        
        # Error if output file already exists
        if os.path.isfile(outFile) and force == False:
            raise FileExistsError(f"GMAP output file '{outFile}' already exists! Will not overwrite.")
        
        # Convert to WSL paths where needed
        if platform.system() == "Windows":
            query = convert_windows_to_wsl_path(self.query)
            target = convert_windows_to_wsl_path(self.target)
            #outFile = convert_windows_to_wsl_path(outFile)
            "WSL redirects to a Windows formatted path... it's strange to me"
        else:
            query = self.query
            target = self.target
        
        # Format command
        cmd = base_subprocess_cmd(self.gmapExe)
        cmd += [
            "-D", os.path.dirname(target), "-d", f"{os.path.basename(target)}.gmap",
            "-f", str(self.outputFormat), "-n", str(self.npaths), "-x", str(self.chimeraMargin),
            "-B", str(self.batch), "-t", str(self.threads),
            f"--max-intronlength-middle={self.maxIntronLengthMiddle}",
            f"--max-intronlength-ends={self.maxIntronLengthEnds}",
            query, ">", outFile
        ]
        
        # Perform GMAP search
        if platform.system() != "Windows":
            run_gmap = subprocess.Popen(" ".join(cmd), shell = True,
                                        stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        else:
            run_gmap = subprocess.Popen(cmd, shell = True,
                                        stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        gmapout, gmaperr = run_gmap.communicate()
        if "queries/sec)" not in gmaperr.decode("utf-8"):
            raise Exception('GMAP searching error text below\n' +
                            gmaperr.decode("utf-8"))

class Salmon_DB:
    '''
    The Salmon_DB Class encapsulates the logic of indexing a FASTA file using salmon.
    
    Attributes:
        fasta (REQUIRED) -- a string indicating the location of a FASTA file.
        salmonDir (REQUIRED) -- a string indicating the location of salmon binaries.
    '''
    def __init__(self, fasta, salmonDir, threads=1):
        self.fasta = fasta
        self.salmonDir = salmonDir
        self.threads = threads
        
        self.isSalmon_DB = True # flag to check object type
    
    @property
    def fasta(self):
        return self._fasta
    
    @fasta.setter
    def fasta(self, value):
        assert type(value).__name__ == "str"
        if not os.path.isfile(value):
            raise Exception(("Fasta parameter is a string, but does not point " + 
                            "to an existing file location"))
        
        self._fasta = value
    
    @property
    def salmonDir(self):
        return self._salmonDir
    
    @salmonDir.setter
    def salmonDir(self, value):
        assert type(value).__name__ == "str" \
            or isinstance(value, Path)
        if not os.path.isdir(value):
            raise Exception(("salmonDir does not point to an existing directory"))
        
        self._salmonDir = value
        self.salmonExe = os.path.join(value, "salmon")
    
    @property
    def salmonExe(self):
        return self._salmonExe
    
    @salmonExe.setter
    def salmonExe(self, value):
        assert type(value).__name__ == "str" \
            or isinstance(value, Path)
        if not os.path.isfile(value):
            raise Exception(("salmonExe does not point to an existing file"))
        
        self._salmonExe = value
    
    @property
    def threads(self):
        return self._threads
    
    @threads.setter
    def threads(self, value):
        assert isinstance(value, int)
        if value < 0:
            raise Exception("Number of threads must be more than 0")
        
        self._threads = value
    
    def index_exists(self):
        '''
        Checks for the expected index directory.
        
        Returns:
            dbExists -- a Boolean where True means the database exists,
                        and False means it does not exist.
        '''
        return os.path.isdir(f"{self.fasta}.salmonDB")
    
    def index(self):
        '''
        Makes a database out of the .fasta value for use in salmon mapping.
        '''
        # Skip if index exists
        if self.index_exists():
            raise FileExistsError("salmon index already exists!")
        
        # Convert to WSL paths where needed
        if platform.system() == "Windows":
            fasta = convert_windows_to_wsl_path(self.fasta)
        else:
            fasta = self.fasta
        
        # Format command
        cmd = base_subprocess_cmd(self.salmonExe)
        cmd += ["index", "--threads", str(self.threads), "--transcripts", fasta,
                "--index", f"{fasta}.salmonDB"]
        
        # Run indexing
        if platform.system() != "Windows":
            run_index = subprocess.Popen(" ".join(cmd), shell = True,
                                         stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        else:
            run_index = subprocess.Popen(cmd, shell = True,
                                         stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        indexout, indexerr = run_index.communicate()
        
        # Raise exception if final part of stderr from successful run isn't found
        if "done building index" not in indexerr.decode("utf-8"):
            raise Exception('salmon indexing error text below\n' + indexerr.decode("utf-8"))

class Salmon:
    '''
    The Salmon Class provides access to the salmon read mapping functionality. It cooperates
    with a Salmon_DB class object to perform the mapping.
    
    Attributes:
        reads (REQUIRED) -- a list containing one string (single-end) or two strings (paired-end)
                            indicating the location of FASTQ file(s) to map.
        target (REQUIRED) -- a Salmon_DB object.
        salmonDir (REQUIRED) -- a string pointing to the location where the mmseqs executable
                                is found.
        threads (OPTIONAL) -- an integer value specifying the number of threads to run when
                              performing MMseqs2 search. Defaults to 1.
    '''
    def __init__(self, reads, target, salmonDir, threads=1):
        self.reads = reads
        self.target = target
        self.salmonDir = salmonDir
        self.threads = threads
        
        # Set helper attributes
        self.isSetup = False
        self.searchResult = None
        self.isSalmon = True
    
    @property
    def reads(self):
        return self._reads
    
    @reads.setter
    def reads(self, value):
        assert isinstance(value, list)
        assert len(value) == 1 or len(value) == 2, \
            "reads must be a list of one or two strings"
        for read in value:
            assert isinstance(read, str)
            if not os.path.isfile(read):
                raise FileNotFoundError(f"Read file does not exist: {read}")
        
        self._reads = value
    
    @property
    def target(self):
        return self._target
    
    @target.setter
    def target(self, value):
        assert type(value).__name__ == "Salmon_DB" \
            or type(value).__name__ == "ZS_MapIO.Salmon_DB" \
            or hasattr(value, "isSalmon_DB") and value.isSalmon_DB == True
        
        self._target = value
    
    @property
    def salmonDir(self):
        return self._salmonDir
    
    @salmonDir.setter
    def salmonDir(self, value):
        assert type(value).__name__ == "str" \
            or isinstance(value, Path)
        if not os.path.isdir(value):
            raise Exception(("salmonDir does not point to an existing directory"))
        
        self._salmonDir = value
        self.salmonExe = os.path.join(value, "salmon")
    
    @property
    def salmonExe(self):
        return self._salmonExe
    
    @salmonExe.setter
    def salmonExe(self, value):
        assert type(value).__name__ == "str" \
            or isinstance(value, Path)
        if not os.path.isfile(value):
            raise Exception(("salmonExe does not point to an existing file"))
        
        self._salmonExe = value
    
    @property
    def threads(self):
        return self._threads
    
    @threads.setter
    def threads(self, value):
        assert isinstance(value, int)
        assert 0 < value, "threads must be a positive integer"
        
        self._threads = value
    
    def setup(self):
        '''
        Indexes database (if relevant)
        '''
        if not self.target.index_exists():
            self.target.index()
        
        self.isSetup = True
    
    def quant(self, outputDirectory):
        '''
        Performs the salmon quant operation.
        
        Parameters:
            outputDirectory -- a string indicating the location to
                               write outputs to.
        '''
        if not self.isSetup:
            self.setup()
        
        # Specify file locations
        reads = [os.path.abspath(read) for read in self.reads]
        target = os.path.abspath(f"{self.target.fasta}.salmonDB")
        
        # Format command
        cmd = [
            self.salmonExe, "quant", "--threads", str(self.threads),
            "--index", target, "--libType", "A", "-o", outputDirectory
        ]
        
        # Add read files
        if len(reads) == 1:
            cmd += ["-r", reads[0]]
        elif len(reads) == 2:
            cmd += ["-1", reads[0], "-2", reads[1]]
        else:
            raise Exception("Incorrect number of read files provided")
        
        # Perform salmon search
        if platform.system() != "Windows":
            run_salmon = subprocess.Popen(" ".join(cmd), shell = True,
                                        stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        else:
            run_salmon = subprocess.Popen(cmd, shell = True,
                                        stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        
        salmonout, salmonerr = run_salmon.communicate()
        if "writing output" not in salmonerr.decode("utf-8"):
            raise Exception('salmon searching error text below\n' +
                            salmonerr.decode("utf-8"))

if __name__ == "__main__":
    pass

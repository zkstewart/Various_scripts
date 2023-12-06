#! python3
# ZS_MapIO.py
# Specifies the GMAP Class for performing GMAP alignments
# using ZS_SeqIO.FASTA and ZS_SeqIO.FastASeq objects.

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
    
    def gmap(self, outFile):
        '''
        Performs the GMAP search operation.
        
        This method does not use .query and .target because they may be in the FASTA or
        FastASeq object types, rather than a string which we require here.
        
        Parameters:
            queryFasta -- a string indicating the location of a FASTA file to use as query.
            dbFastaFile -- a string indicating the location of a FASTA file to use as the database.
            outFile -- a string indicating the location to write the output file to. File must not already exist!
        '''
        # Skip if index does not exist
        if not self.db.index_exists():
            raise FileNotFoundError("GMAP index does not exist! Cannot perform search.")
        
        # Error if output file already exists
        if os.path.isfile(outFile):
            raise FileExistsError(f"GMAP output file '{outFile}' already exists! Will not overwrite.")
        
        # Convert to WSL paths where needed
        if platform.system() == "Windows":
            query = convert_windows_to_wsl_path(self.query)
            target = convert_windows_to_wsl_path(self.target)
            #outFile = convert_windows_to_wsl_path(outFile)
            "WSL redirects to a Windows formatted path... it's strange to me"
        
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

if __name__ == "__main__":
    pass

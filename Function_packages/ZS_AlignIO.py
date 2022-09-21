#! python3
# ZS_AlignIO.py
# Contains various Classes to perform manipulations involving
# FASTA objects that are to be aligned or have been aligned.

# TBD
## Aligner class that can
## 1) align a FASTA and update its FastASeq.gap_seq values [-]
## 2) align nucleotides as peptide then convert back to nucleotide [-]
## 3) detect outliers in a MSA by sequence conservation [-]
## 4) detect outliers in a MSA by phylogeny mismatch [-]

import os, platform, sys, subprocess, hashlib, time, random, re, subprocess
from copy import deepcopy
from Bio.Align.Applications import MafftCommandline

sys.path.append(os.path.dirname(__file__))
from ZS_SeqIO import FASTA, FastASeq
from ZS_GFF3IO import Feature

import parasail
if platform.system() != 'Windows':
    from skbio.alignment import StripedSmithWaterman

def _tmp_file_name_gen(prefix, suffix):
    '''
    Hidden function for use by Class methods.
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

class MAFFT:
    '''
    The MAFFT Class is intended to be used alongside the ZS_SeqIO.FASTA Class to
    handle the alignment of FASTA sequences with OOP magics. In short, it provides
    easy access to MAFFT's E-INSi and L-INSi algorithms with default parameters
    otherwise.
    
    It's primary purpose is the .run() method which receives a ZS_SeqIO.FASTA instance
    and updates the .gap_seq attributes of its underlying ZS_SeqIO.FastASeq objects.
    '''
    def __init__(self, mafftDir):
        # Validate input type and location
        assert isinstance(mafftDir, str)
        assert os.path.isdir(mafftDir)
        if platform.system() == "Windows":
            if not os.path.isfile(os.path.join(mafftDir, "mafft.bat")):
                raise Exception("{0} does not exist".format(os.path.join(mafftDir, "mafft.bat")))
        else:
            if not os.path.isfile(os.path.join(mafftDir, "mafft")) and not os.path.isfile(os.path.join(mafftDir, "mafft.exe")):
                raise Exception("mafft or mafft.exe does not exist at {0}".format(mafftDir))
        
        # Establish commandline function
        if platform.system() == "Windows":
            self.cline = MafftCommandline(os.path.join(mafftDir, "mafft.bat"))
        else:
            if os.path.isfile(os.path.join(mafftDir, "mafft")):
                self.cline = MafftCommandline(os.path.join(mafftDir, "mafft"))
            else:
                self.cline = MafftCommandline(os.path.join(mafftDir, "mafft.exe"))
        
        # Set default attributes
        self.mafftDir = mafftDir
        self.cline.genafpair = True   # E-INSi
        self.cline.localpair = False  # L-INSi
        self.cline.thread = 1
    
    def use_einsi(self):
        '''
        This is a method to toggle the use in favour of E-INSi which is equivalent
        to setting the --genafpair value for MAFFT's command line.
        '''
        self.cline.genafpair = True
        self.cline.localpair = False
    
    def use_linsi(self):
        '''
        This is a method to toggle the use in favour of L-INSi which is equivalent
        to setting the --localpair value for MAFFT's command line.
        '''
        self.cline.localpair = True
        self.cline.genafpair = False
    
    def set_threads(self, num):
        '''
        This method allows the use of multithreading by MAFFT. num should be a valid
        integer greater than 0. Be sensible since there's no upper limit checking.
        '''
        assert isinstance(num, int)
        if num < 0:
            raise Exception("Number of threads must be more than 0")
        
        self.cline.thread = num
        
    def run(self, FASTA_obj):
        '''
        Handles the execution of MAFFT alignment from start to finish. It must:
        
            1) Create a temporary file for MAFFT to read from the .FASTA object
            2) Perform the alignment with MAFFT
            3) Parse the output file and modify the .FASTA object appropriately
            4) Clean up temporary and output files
        
        Before calling this, you might consider the other methods of this class e.g.,
        set_threads() and use_linsi() if you don't want the default E-INSi behaviour.
        
        Params:
            FASTA_obj -- an object of ZS_SeqIO.FASTA class.
        '''
        # Validate input value type
        assert type(FASTA_obj).__name__ == "FASTA" or type(FASTA_obj).__name__ == "ZS_SeqIO.FASTA"
                
        # Create temporary file
        tmpHash = hashlib.sha256(bytes(str(FASTA_obj.fileOrder[0][0]) + str(time.time()) + str(random.randint(0, 100000)), 'utf-8') ).hexdigest()
        tmpFileName = _tmp_file_name_gen("mafft_tmp" + tmpHash[0:20], "fasta")
        FASTA_obj.write(tmpFileName)
        
        # Run
        self.cline.input = tmpFileName
        stdout, stderr = self.cline()
        if stdout == '':
            raise Exception("MAFFT error text below" + str(stderr))
        
        # Clean MAFFT output
        stdout = stdout.split("\n")
        while stdout[-1] == "\n" or stdout[-1] == "" or stdout[-1].startswith("Terminate batch job"):   # Remove junk
            del stdout[-1]
        
        # Parse output back into the FASTA object
        thisSeq = None
        ongoingCount = 0
        for line in stdout:
            if line.startswith(">"):
                if thisSeq == None:
                    thisSeq = []
                else:
                    FASTA_obj[ongoingCount].gap_seq = "".join(thisSeq)
                    thisSeq = []
                    ongoingCount += 1
            else:
                thisSeq.append(line)
        FASTA_obj[ongoingCount].gap_seq = "".join(thisSeq) # handle last sequence in the iteration
        
        # Set flag that this FASTA object has been aligned
        FASTA_obj.isAligned = True
        
        # Clean up temporary file
        os.unlink(tmpFileName)
    
    def run_nucleotide_as_protein(self, FASTA_obj, strand=None, frame=None):
        '''
        Refer to method header of run() for the fundamentals of this method.
        Refer to FastASeq method header of get_translation() for descriptions
        of the other parameters.
        
        What this method does differently to run() is that it will:
        
            1) Translate nucleotides into their corresponding protein sequence
            2) Align the protein sequences
            3) Map the gaps in the protein alignment back to the nucleotide sequence
        
        In effect, it allows nucleotide sequences to be aligned as codons rather than
        as individual base pairs, which might be more biologically relevant when dealing
        with closely related sequences.
        
        Note that findBestFrame is an implicit parameter. If strand is None or the indexed list
        value is None, and/or the frame value is similar, then findBestFrame will be set to True
        since we'll need to derive the appropriate strand and/or frame depending on which is None.
        If you set both strand and frame, or the indexed list value is set, then we obviously
        don't need to find how to translate it now do we?
        
        Parameters:
            findBestFrame -- a boolean indicating whether we should get the FastASeq.get_translation()
                             function figure out the best translation.
            strand -- one of three types; 1) an integer of 1 to get a +ve strand translation, or -1 for 
                      a -ve strand translation, 2) None to not constrain strandedness of the translation
                      if findBestFrame is True, or 3) a list of integers and/or None.
            frame -- one of three types; 1) an integer in the range(0,3) corresponding to 0-based
                     frame number, 2) None None to not constrain the frame of the translation
                     if findBestFrame is True, or 3) a list of integers and/or None.
        '''
        # Validate input FASTA parameter
        assert type(FASTA_obj).__name__ == "FASTA" or type(FASTA_obj).__name__ == "ZS_SeqIO.FASTA"
        
        # Validate strand parameter
        assert isinstance(strand, int) or strand == None or isinstance(strand, list)
        if isinstance(strand, int):
            assert strand in [1, -1], "Strand was provided but was not an appropriate value"
        if isinstance(strand, list):
            for _strand in strand:
                if isinstance(_strand, int):
                    assert _strand in [1, -1], "Strand was provided as a list, but contains inappropriate value(s)"
                else:
                    assert _strand == None, "Strand was provided as a list, but contains unknown types"
        
        # Validate frame parameter
        assert isinstance(frame, int) or frame == None or isinstance(frame, list)
        if isinstance(frame, int):
            assert frame in range(0, 3), "Frame was provided but was not an appropriate value"
        if isinstance(frame, list):
            for _frame in frame:
                if isinstance(_frame, int):
                    assert _frame in range(0,3), "Frame was provided as a list, but contains inappropriate value(s)"
                else:
                    assert _frame == None, "Frame was provided as a list, but contains unknown types"
        
        # Create a dummy FASTA object with our translations
        dummy = deepcopy(FASTA_obj)
        strands, frames = [], []
        for i in range(len(dummy)):
            # Get the appropriate strand/frame values
            if isinstance(strand, list):
                _strand = strand[i]
            else:
                _strand = strand
            if isinstance(frame, list):
                _frame = frame[i]
            else:
                _frame = frame
            
            # Run the translation
            if _strand != None and _frame != None:
                protein, _strand, _frame = dummy.seqs[i].get_translation(False, _strand, _frame)
            else:
                protein, _strand, _frame = dummy.seqs[i].get_translation(True, _strand, _frame)
            
            # Modify stop codons to be "X"s
            '''
            We need to do this since MAFFT will NOT report stop codons in the output sequence,
            which poses problems for us when mapping codon positions back. However, we don't
            need to actively track these instances since our untranslation works by just reading
            the character and mapping it to the next 3 nucleotides - regardless of what the
            character is. Hence, it's easy to just sub out *'s for X's and it'll just work!
            '''
            protein = protein.replace("*", "X")
            
            # Update and store results
            strands.append(_strand)
            frames.append(_frame)
            dummy.seqs[i].seq = protein
        
        # Align it
        self.run(dummy)
        
        # Map back aligned proteins to their original nucleotide counterparts
        trimmedBits = []
        for i in range(len(dummy)):
            _strand = strands[i]
            _frame = frames[i]
            nuc = FASTA_obj[i].seq if _strand == 1 else FASTA_obj[i].get_reverse_complement()
            alignedProt = dummy[i].gap_seq
            
            # Perform mapping procedure
            alignedNuc = ""
            codonIndex = _frame
            for proteinIndex in range(0, len(alignedProt)):
                if alignedProt[proteinIndex] != "-":
                    alignedNuc += nuc[codonIndex: codonIndex+3] # codon length == 3
                    codonIndex += 3 # iterate our nucleotide position marker
                else:
                    alignedNuc += "---"
            
            # Add in any leftover bits e.g., stop codons wherever they fall
            '''
            This means, any leftover bit not included in the translation goes as
            close to the nearest sequence position as possible
            '''
            alignedNuc = alignedNuc.rstrip("-") # we'll pad gaps below
            alignedNuc += nuc[codonIndex:] # add any potential stop codons back
            
            # Update original FASTA values
            FASTA_obj[i].gap_seq = alignedNuc
            
            # Retain any bits we trimmed off from the start
            trimmedBits.append(nuc[0:_frame])
        
        # Pad out the end of the sequence for anything we rstripped earlier
        maxLen = max([len(FastASeq_obj.gap_seq) for FastASeq_obj in FASTA_obj])
        for FastASeq_obj in FASTA_obj:
            if len(FastASeq_obj.gap_seq) != maxLen:
                FastASeq_obj.gap_seq += "-"*(maxLen - len(FastASeq_obj.gap_seq))
        
        # Add in any bits we trimmed off from the start now
        for i in range(len(FASTA_obj)):
            FastASeq_obj = FASTA_obj[i]
            trimmedBit = trimmedBits[i]
            FastASeq_obj.gap_seq = trimmedBit + FastASeq_obj.gap_seq.lstrip("-")
        
        # Pad out the start of the sequence for anything we lstripped just above
        maxLen = max([len(FastASeq_obj.gap_seq) for FastASeq_obj in FASTA_obj])
        for FastASeq_obj in FASTA_obj:
            if len(FastASeq_obj.gap_seq) != maxLen:
                FastASeq_obj.gap_seq = "-"*(maxLen - len(FastASeq_obj.gap_seq)) + FastASeq_obj.gap_seq
        
        # Set flag that this FASTA object has been aligned
        FASTA_obj.isAligned = True
    
    def add(self, aligned_FASTA_obj, add_FASTA_obj):
        '''
        Handles the execution of MAFFT alignment, running in the mode where new sequences
        are added into an existing alignment. This method otherwise behaves similarly to
        the .run() method, but it:
        
            1) Creates TWO temporarys file for MAFFT to read from the TWO FASTA objects provided
            2) Performs the alignment with MAFFT
            3) Parses the output file and without modifying input FASTA objects
            4) Cleans up temporary and output files
        
        Before calling this, you might consider the other methods of this class e.g.,
        set_threads() and use_linsi() if you don't want the default E-INSi behaviour.
        
        Params:
            aligned_FASTA_obj -- an object of ZS_SeqIO.FASTA class which has already been
                                 aligned and is to have new sequences added into it.
            add_FASTA_obj -- an object of ZS_SeqIO.FASTA class which has not been aligned
                               and will be added into the existing alignment
        Returns:
            result_FASTA_obj -- an object of ZS_SeqIO.FASTA class which results from MAFFT
                                adding the add_FASTA_obj into aligned_FASTA_obj
        '''
        # Validate input value type
        # assert isinstance(aligned_FASTA_obj, FASTA) ## This fails -- annoying!
        assert type(aligned_FASTA_obj).__name__ == "FASTA" or type(aligned_FASTA_obj).__name__ == "ZS_SeqIO.FASTA"
        assert aligned_FASTA_obj.isAligned, "aligned_FASTA_obj must be aligned first!"
        # assert isinstance(add_FASTA_obj, FASTA)
        assert type(add_FASTA_obj).__name__ == "FASTA" or type(add_FASTA_obj).__name__ == "ZS_SeqIO.FASTA"
        
        # Create temporary files
        tmpHash = hashlib.sha256(bytes(str(aligned_FASTA_obj.fileOrder[0][0]) + str(time.time()) + str(random.randint(0, 100000)), 'utf-8') ).hexdigest()
        
        tmpAlignedFileName = _tmp_file_name_gen("mafft_aligned_tmp" + tmpHash[0:20], "fasta")
        aligned_FASTA_obj.write(tmpAlignedFileName, asAligned=True, withDescription=True)
        
        tmpAddedFileName = _tmp_file_name_gen("mafft_added_tmp" + tmpHash[0:20], "fasta")
        add_FASTA_obj.write(tmpAddedFileName, withDescription=True)
        
        tmpOutputFileName = _tmp_file_name_gen("mafft_output_tmp" + tmpHash[0:20], "fasta")
        
        # Format command for running
        cmd = "{0} --thread {1} --add \"{2}\" \"{3}\" > \"{4}\"".format(
            os.path.join(self.mafftDir, "mafft"), self.cline.thread,
            tmpAddedFileName, tmpAlignedFileName, tmpOutputFileName
        )
        
        # Run
        run_mafft = subprocess.Popen(cmd, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
        stdout, stderr = run_mafft.communicate()
        if stdout == '':
            raise Exception("MAFFT error text below" + str(stderr))
        
        # Parse output file into new FASTA object
        result_FASTA_obj = FASTA(tmpOutputFileName, isAligned=True)
        
        # Clean up temporary files
        os.unlink(tmpAlignedFileName)
        os.unlink(tmpAddedFileName)
        os.unlink(tmpOutputFileName)
        
        # Return new result
        return result_FASTA_obj

class SSW_Result:
    '''
    Simple object to act as a container for the results from SSW alignment.
    '''
    def __init__(self, queryAlign, targetAlign, score, queryStartIndex, targetStartIndex):
        self.queryAlign = queryAlign
        self.targetAlign = targetAlign
        self.score = score
        self.queryStartIndex = queryStartIndex
        self.targetStartIndex = targetStartIndex

class SSW:
    '''
    Class to encapsulate static methods used for performing alignments using SSW implementations.
    '''
    @staticmethod
    def ssw_parasail(queryString, targetString):
        '''
        Special implementation of striped Smith Waterman alignment for exon liftover
        project.
        '''
        # Perform SSW with parasail implementation
        profile = parasail.profile_create_sat(targetString, parasail.blosum62)
        alignment = parasail.sw_trace_striped_profile_sat(profile, queryString, 10, 1)
        queryAlign = alignment.traceback.ref # this ssw implementation sees things differently than I
        targetAlign = alignment.traceback.query
        
        # Figure out where we're starting for the alignments
        queryStartIndex = queryString.find(queryAlign.replace('-', ''))
        targetStartIndex = targetString.find(targetAlign.replace('-', ''))
        
        # Make and return an object containing all our results
        result = SSW_Result(queryAlign, targetAlign, alignment.score, queryStartIndex, targetStartIndex)
        return result

    @staticmethod
    def ssw_skbio(queryString, targetString):
        if platform.system() == 'Windows':
            print("skbio is not supported on Windows yet (as of last time this code was touched); won't proceed")
            return
        
        # Perform SSW with scikit.bio implementation
        query = StripedSmithWaterman(targetString)
        alignment = query(queryString)
        targetAlign = alignment.aligned_query_sequence
        queryAlign = alignment.aligned_target_sequence
        
        # Figure out where we're starting for the alignments
        queryStartIndex = queryString.find(queryAlign.replace('-', ''))
        targetStartIndex = targetString.find(targetAlign.replace('-', ''))
        
        # Make and return an object containing all our results
        result = SSW_Result(queryAlign, targetAlign, alignment.optimal_alignment_score, queryStartIndex, targetStartIndex)
        return result

class Exonerate:
    '''
    The Exonerate Class is intended to be used alongside the ZS_SeqIO.FASTA Class to
    handle the alignment of FASTA sequences with OOP magics. In short, it provides
    easy access to exonerate to perform various forms of alignment, returning a 
    parsed result.
    
    As of now, it's fairly limited and has only been tested with protein2genome
    and "--showtargetgff yes". In fact, the only exonerate output parser I have
    is limited to handling results where this is true. So, don't mess with showtargetgff.
    '''
    def __init__(self, exonerateExe, query, target):
        self.exonerateExe = exonerateExe
        self.query = query
        self.target = target
        
        # Set default attributes
        self.model = "ungapped" # default of exonerate v2.4.0
        self.exhaustive = False
        self.score = 100
        self.showtargetgff = True # exonerate is False by default, but this class only handles GFF
    
    @property
    def exonerateExe(self):
        return self._exonerateExe
    
    @exonerateExe.setter
    def exonerateExe(self, value):
        assert isinstance(value, str) and os.path.isfile(value), \
            f"'{value}' is not a string or does not point to a file"
        
        if platform.system() == "Windows":
            print("Exonerate Class on Windows assumes exonerate was built using WSL; carrying on...")
            self._exonerateExe = Exonerate.convert_windows_to_wsl_path(value)
        else:
            self._exonerateExe = value
    
    @property
    def model(self):
        return self._model
    
    @model.setter
    def model(self, value):
        options = [
            "ungapped", "ungapped:trans",
            "affine:global", "affine:bestfit", "affine:local", "affine:overlap",
            "est2genome", "ner",
            "protein2genome", "protein2genome:bestfit",
            "protein2dna", "protein2dna:bestfit",
            "coding2coding", "coding2genome", "cdna2genome", "genome2genome"
        ]
        if value not in options:
            raise ValueError(f"{value} is not recognised as an option within {options}")
        self._model = value
    
    @property
    def query(self):
        return self._query
    
    @query.setter
    def query(self, value):
        recognisedTypes = ["str", "FASTA", "ZS_SeqIO.FASTA", "FastASeq", "ZS_SeqIO.FastASeq"]
        assert type(value).__name__ in recognisedTypes, \
            "Query value type not handled by Exonerate class"
        self._query = value
    
    @property
    def target(self):
        return self._target
    
    @target.setter
    def target(self, value):
        recognisedTypes = ["str", "FASTA", "ZS_SeqIO.FASTA", "FastASeq", "ZS_SeqIO.FastASeq"]
        assert type(value).__name__ in recognisedTypes, \
            "Target value type not handled by Exonerate class"
        self._target = value
    
    @property
    def exhaustive(self):
        return self._exhaustive
    
    @exhaustive.setter
    def exhaustive(self, value):
        assert isinstance(value, bool), \
            "Exhaustive property is a boolean switch; must be True or False"
        self._exhaustive = value
    
    @property
    def score(self):
        return self._score
    
    @score.setter
    def score(self, value):
        assert isinstance(value, int), \
            "Score value must be an integer"
        assert value >= 0, \
            "Score value must be >= 0"
        
        self._score = value
    
    @property
    def showtargetgff(self):
        return self._showtargetgff
    
    @showtargetgff.setter
    def showtargetgff(self, value):
        assert isinstance(value, bool), \
            "showtargetgff property is a boolean switch; must be True or False"
        self._showtargetgff = value
    
    @staticmethod
    def convert_windows_to_wsl_path(windowsPath):
        '''
        Provides simple functionality to infer the WSL path from
        a windows path, provided as a string.
        
        Parameters:
            windowsPath -- a string indicating the full path to a file
                           or directory of interest; this MUST include
                           the root character e.g., 'D:\\' or 'C:\\'
        Returns:
            wslPath -- a string indicating the inferred full path to the
                       given file or directory using WSL formatting
        '''
        # Check that the path is something we can work with
        driveRegex = re.compile(r"^([A-Za-z]{1}):\\")
        assert driveRegex.match(windowsPath) != None, \
            f"'{windowsPath}' is not recognised as a full, root drive inclusive path"
        assert os.path.exists(windowsPath), \
            f"'{windowsPath}' is not recognised as an existing path"
        
        # If it is, convert it
        driveLetter = driveRegex.match(windowsPath).group(1)
        wslPath = "/{0}".format("/".join(
            [
                "mnt",
                driveLetter.lower(),
                *windowsPath.split("\\")[1:]
            ]
        ))
        
        return wslPath
    
    def _get_filename_for_query_or_target(self, qt):
        '''
        Hidden method for use when boiling down one of the three data types (FASTA, FastASeq, and string)
        into a string representing the query file name.
        
        If it's already a string for a file, this does nothing. Otherwise, it will make sure a file exists
        with the FASTA data contents for use by exonerate.
        
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
                       else qt.fileOrder[0][0] if (type(qt).__name__ == "ZS_SeqIO.FASTA" or type(qt).__name__ == "FASTA") and qt.fileOrder != [] \
                       else str(qt) if (type(qt).__name__ == "ZS_SeqIO.FASTA" or type(qt).__name__ == "FASTA") and qt.fileOrder == []  \
                       else qt.id
        tmpHash = hashlib.sha256(bytes(hashForTmp + str(time.time()) + str(random.randint(0, 100000)), 'utf-8') ).hexdigest()
        
        # If qt is a string but not a file, make it a FastASeq object
        if isinstance(qt, str) and not os.path.isfile(qt):
            qt = FastASeq("tmpID", qt)
        
        # If qt is a FastASeq, make it a FASTA object
        if type(qt).__name__ == "FastASeq" or type(qt).__name__ == "ZS_SeqIO.FastASeq":
            tmpFastaName = _tmp_file_name_gen("qt_fasta_tmp" + tmpHash[0:20], "fasta")
            with open(tmpFastaName, "w") as fileOut:
                fileOut.write(">{0}\n{1}\n".format(qt.id, qt.seq))
            qt = FASTA(tmpFastaName)
            os.unlink(tmpFastaName)
        
        # If qt is a FASTA, make it into a file
        isTemporary = False
        if type(qt).__name__ == "FASTA" or type(qt).__name__ == "ZS_SeqIO.FASTA":
            tmpQtName = _tmp_file_name_gen("exonerate_tmp" + tmpHash[0:20], "fasta")
            qt.write(tmpQtName)
            qt = tmpQtName # after this point, qt will be a string indicating a FASTA file name
            isTemporary = True # if we set this, qt was not originally a string
        
        return qt, isTemporary
    
    def _format_exonerate_cmd(self, queryFile, targetFile):
        '''
        Hidden helper function for getting a list amenable to subprocess.run()
        that sets all relevant cmd tags depending on this object's properties.
        It expects that the query and target files have already been subjected
        to ._get_filename_for_query_or_target(), and their return values are
        to be given to this function.
        
        Parameters:
            queryFile -- a string pointing to a FASTA file that exists
            targetFile -- a string pointing to a FASTA file that exists
        Returns:
            cmds -- a list amenable to subprocess.run()
        '''
        cmds = []
        if platform.system() == "Windows":
            cmds += ["wsl", "~", "-e"]
        
        cmds += [
            self.exonerateExe, "--model", self.model,
            "--score", str(self.score)
        ]
        
        if self.showtargetgff is True:
            cmds += ["--showtargetgff", "yes"]
        if self.exhaustive is True:
            cmds += ["--exhaustive", "yes"]
        
        if platform.system() == "Windows":
            cmds += [
                Exonerate.convert_windows_to_wsl_path(os.path.abspath(queryFile)),
                Exonerate.convert_windows_to_wsl_path(os.path.abspath(targetFile))
            ]
        else:
            cmds += [queryFile, targetFile]
        
        return cmds
    
    def run_exonerate(self):
        '''
        Performs the exonerate operation using the parameters already specified during/after
        creation of this object. This function pawns off the handling to one of two hidden
        subfunctions depending on whether we're on Windows or not.
        
        Returns:
            features -- a list containing ZS_GFF3IO.Feature objects corresponding to
                        all sequence matches reported by exonerate
        '''
        # Get file names for query and target after data type coercion
        q, qIsTemporary = self._get_filename_for_query_or_target(self.query)
        t, tIsTemporary = self._get_filename_for_query_or_target(self.target)
        
        # Run exonerate
        cmds = self._format_exonerate_cmd(q, t)
        exonerate = subprocess.run(cmds, capture_output=True)
        
        # Parse exonerate results into gene features
        features = Exonerate.parse_exonerate_gff_stdout(exonerate.stdout.decode())
        
        # Clean up and return
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
        return features
    
    @staticmethod
    def parse_exonerate_gff_stdout(exonerateGffStdout):
        '''
        This function is capable of reading the stdout produced by running exonerate
        with --showtargetgff yes and producing GFF3 Features that encapsulate the
        details of aligned models reported by exonerate. This stdout can be obtained
        via subprocess or just by reading the file into a string (hopefully it's not
        too big!)
        
        Parameters:
            exonerateGffStdout -- a string containing all of the stdout from running
                                  exonerate with "--showtargetgff yes".
        Returns
            features -- a list containing ZS_GFF3IO.Feature objects corresponding to
                        all sequence matches reported by exonerate
        '''
        featureDict = {"gene": {}, "mRNA": {}}
        geneIDDict = {} # just for naming purposes for this script
        for line in exonerateGffStdout.split("\n"):
            # Skip irrelevant lines
            try:
                contig, source, featureType, start, end, \
                        score, strand, frame, attributes \
                        = line.rstrip('\t\n').split('\t')
                start = int(start)
                end = int(end)
                assert strand in ['+', '-']
            except:
                continue  # If any of the above fail we know this isn't a GFF line
            
            # Parse attributes
            splitAttributes = []
            for a in attributes.split(" ; "):
                splitAttributes += a.split(" ", maxsplit=1)
            attributesDict = {splitAttributes[i]: splitAttributes[i+1] for i in range(0, len(splitAttributes)-(len(splitAttributes)%2), 2)}
            
            # Reformat detail lines to be GFF3-style
            if featureType == 'gene':
                geneID, geneIDDict = Exonerate._exonerate_geneid_produce(contig, attributesDict['sequence'], geneIDDict)        # This will carry over into CDS/exon lines below this
                newAttributesDict = {
                    "ID": geneID,
                    "Name": f"exonerate_{geneID}",
                    "Sequence": attributesDict['sequence'],
                    "identity": attributesDict['identity'],
                    "similarity": attributesDict['similarity']
                }
                exonCount = 1
                cdsCount = 1
            elif featureType == 'cds':
                featureType = 'CDS'
                newAttributesDict = {
                    "ID": f"{geneID}.mrna1.cds{cdsCount}", # geneID carries over from the gene line
                    "Parent": f"{geneID}.mrna1",
                }
                cdsCount += 1
            elif featureType == 'exon':
                newAttributesDict = {
                    "ID": f"{geneID}.mrna1.exon{exonCount}",
                    "Parent": f"{geneID}.mrna1"
                }
                exonCount += 1
            else:
                continue # Skip any other lines; these include 'similarity', 'intron', and 'splice5'/'splice3'
            
            # Create a feature from what we have and associate it appropriately
            if featureType == 'gene':
                # Create the gene feature
                geneFeature = Feature()
                geneFeature.add_attributes(newAttributesDict)
                geneFeature.add_attributes({
                    "contig": contig, "source": source, "type": featureType,
                    "start": int(start), "end": int(end), "coords": [int(start), int(end)],
                    "score": score, "strand": strand, "frame": frame
                })
                featureDict["gene"][geneID] = geneFeature
                
                # Create an mRNA feature and nestle it as a child under its parental gene
                "So caring... so nurturing..."
                mrnaFeature = Feature()
                mrnaFeature.add_attributes({
                    "ID": f"{geneID}.mrna1",
                    "Name": f"exonerate_{geneID}",
                    "Sequence": attributesDict['sequence'],
                    "identity": attributesDict['identity'],
                    "similarity": attributesDict['similarity']
                })
                mrnaFeature.add_attributes({
                    "contig": contig, "source": source, "type": "mRNA",
                    "start": int(start), "end": int(end), "coords": [int(start), int(end)],
                    "score": score, "strand": strand, "frame": frame
                })
                featureDict["mRNA"][f"{geneID}.mrna1"] = mrnaFeature
                geneFeature.add_child(mrnaFeature)
            else:
                # Create the CDS/exon feature
                subFeature = Feature()
                subFeature.add_attributes(newAttributesDict)
                subFeature.add_attributes({
                    "contig": contig, "source": source, "type": featureType,
                    "start": int(start), "end": int(end), "coords": [int(start), int(end)],
                    "score": score, "strand": strand, "frame": frame
                })
                featureDict["mRNA"][newAttributesDict["Parent"]].add_child(subFeature)
        
        return [feature for feature in featureDict["gene"].values()]
    
    @staticmethod
    def _exonerate_geneid_produce(contigID, sequenceID, idDict):
        '''
        This is legacy code from the exonerate_gene_find.py program. I don't want
        to touch it even if it is ugly.
        '''
        # Produce the basic ID prefix
        sequenceBit = sequenceID.split('|')
        sequenceBit.sort(key=len, reverse=True) # This should help to handle sequence IDs like eg|asdf|c1; we assume the longest part is the most informative which should be true with Trinity and GenBank/Swiss-Prot IDs
        # Specifically handle older Trinity-style IDs
        if len(sequenceBit) > 1:
            if sequenceBit[1].startswith('TR') and sequenceBit[1][2:].isdigit():
                sequenceBit[0] = sequenceBit[1] + '_' + sequenceBit[0]
        # Specifically handle ToxProt-style IDs [Note that, normally, the longest bit in a UniProt ID is what we want, but with toxprot the format differs e.g., with "toxprot_sp|P25660|VKT9_BUNFA" the longest section is ambiguous and might return the toxprot bit]
        if len(sequenceBit) == 3 and sequenceID.split('|')[0] == 'toxprot_sp':
            sequenceBit[0] = sequenceID.split('|')[2]
        # Format gene ID
        geneID = contigID + '.' + sequenceBit[0]
        if geneID not in idDict:
            idDict[geneID] = 1
        # Produce the final geneID and iterate idDict's contents
        outGeneID = geneID + '.' + str(idDict[geneID])
        idDict[geneID] += 1
        return outGeneID, idDict

if __name__ == "__main__":
    pass

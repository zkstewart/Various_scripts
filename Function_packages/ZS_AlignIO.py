#! python3
# ZS_AlignIO.py
# Contains various Classes to perform manipulations involving
# FASTA objects that are to be aligned or have been aligned.

import os, platform, sys, re, subprocess, shutil
from pathlib import Path
from copy import deepcopy

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import ZS_Utility
import ZS_SeqIO
from ZS_GFF3IO import Feature, GFF3

import parasail
if platform.system() != 'Windows':
    from skbio.alignment import StripedSmithWaterman

def _create_file_from_FASTA_object(FASTA_obj, useExistingFile=True):
    '''
    Hidden helper function for the ZS_AlignIO package. It performs actions specific to
    this module, so I won't move it to ZS_Utility or ZS_SeqIO.Conversion.
    
    Parameters:
        useExistingFile -- a boolean indicating whether an existing file with the
                            expected file name hash should be used as-is, or if we
                            should create a new file instead.
    Returns:
        fastaFileName -- a string indicating the location of an existing FASTA file
        fileWasCreated -- a boolean indicating whether a new file was created (True)
                          or not (False)
    '''
    assert hasattr(FASTA_obj, "isFASTA") and FASTA_obj.isFASTA is True, \
        "Value provided to _get_filename_for_FASTA_object is not actually a FASTA object"
    assert isinstance(useExistingFile, bool), \
        "useExistingFile value must be True or False"
    
    # Get hash string for file creation
    tmpHash = ZS_SeqIO.Conversion.get_hash_for_input_sequences(FASTA_obj, randomHash=False, maxLength=20)
    
    # Get our expected file name
    if useExistingFile is True:
        fastaFileName = f"alignIO_tmp_{tmpHash}.fasta"
    else:
        fastaFileName = ZS_Utility.tmp_file_name_gen("alignIO_tmp_" + tmpHash, "fasta")
    
    # If file already exists, just return that file name
    "If the file exists, that means it already exists, AND that useExistingFile is True"
    if os.path.isfile(fastaFileName):
        return fastaFileName, False # file was not created
    
    # Otherwise, create a new file
    "If the file does not exist, that means the file wasn't already there OR useExistingFile is False"
    FASTA_obj.write(fastaFileName)
    return fastaFileName, True # file was created

##############################

class MSA:
    '''
    A namespace for containing static methods that are useful for manipulating
    multiple sequence alignments. It is expected that most functions will receive
    a FASTA file or ZS_SeqIO.FASTA object as input and return a ZS_SeqIO.FASTA
    object as output.
    '''
    @staticmethod
    def locate_variants_from_msa(FASTA_obj, isNucleotide, asCodons, reportUntilStop=False):
        '''
        Parameters:
            FASTA_obj -- a ZS_SeqIO.FASTA object containing aligned sequences
            isNucleotide -- a boolean indicating whether the sequences are nucleotide
                            (True) or protein (False)
            asCodons -- a boolean indicating whether to consider codons as the unit of
                        comparison for nucleotide sequences; only relevant if isNucleotide
                        is True.
            reportUntilStop -- OPTIONAL; a boolean indicating whether to report variants
                               after a stop codon is encountered. Default == False i.e.,
                               variants will be reported for a sequence even after a stop
                               codon has been encountered.
        Returns:
            variantDict -- a dictionary with structure like:
                        {
                            positionNumber: {
                                "consensus": "consensusResidue",
                                "variants": {
                                    "variantResidue1": ["seqID1", "seqID2", ...],
                                    "variantResidue2": ["seqID3", "seqID4", ...],
                                    ...
                                }
                            }
                        }
        '''
        def translate_codon(codon):
            'Private function of this method to translate a codon into a protein residue'
            residue = ZS_SeqIO.FastASeq.dna_to_protein(codon.rstrip("-"))
            if residue == "": # this might happen if the codon is incomplete and doesn't translate
                residue = "-"
            return residue
        
        if asCodons:
            assert isNucleotide, "locate_variants_from_msa: asCodons should only be True if isNucleotide is True"
        
        variantDict = {}
        
        # Generate a consensus sequence
        consensusSeq = FASTA_obj.generate_consensus(asCodons=asCodons)

        # Iterate over positions in the alignment to locate variants
        seqsToSkip = set()
        stepSize = 3 if asCodons else 1
        for i in range(0, len(consensusSeq), stepSize):
            # Get the codon-normalised position index
            positionNumber = int(i / stepSize)
            
            # Obtain the consensus residue at this position for comparison
            consensusResidue = consensusSeq[i:i+stepSize]
            consensusResidue = consensusResidue if not (isNucleotide and asCodons) \
                                else translate_codon(consensusResidue)
            
            # Iterate over each sequence in the alignment and look at this position
            for FastASeq_obj in FASTA_obj.seqs:
                if FastASeq_obj.id in seqsToSkip:
                    continue
                
                residue = FastASeq_obj.gap_seq[i:i+stepSize]
                residue = residue if not (isNucleotide and asCodons) \
                            else translate_codon(residue)
                
                # Store data if a variant exists
                if residue != consensusResidue:
                    variantDict.setdefault(positionNumber, {
                        "consensus": consensusResidue,
                        "variants": {}
                    })
                    variantDict[positionNumber]["variants"].setdefault(residue, [])
                    variantDict[positionNumber]["variants"][residue].append(FastASeq_obj.id)
                
                # If we're not reporting after stop codons, check if we've hit one
                if residue == "*" and reportUntilStop:
                    seqsToSkip.add(FastASeq_obj.id)
        
        return variantDict
    
    def trim(FASTA_obj, pctTrim=0.70):
        '''
        Receives an aligned ZS_SeqIO.FASTA object and trims the sequences by removing columns
        from the start and end of the alignment where the proportion of sequences with a residue
        present is less than the pctTrim threshold, stopping when the threshold is met.
        
        Parameters:
            FASTA_obj -- a ZS_SeqIO.FASTA object containing aligned sequences.
            pctTrim -- OPTIONAL; a float from 0->1 exclusive indicating the minimum proportion
                       of sequences that must have a residue present in a single column
                       for that column to be used as the left and right borders of the MSA;
                       default == 0.70 i.e., 70% of sequences must have a residue present.
        Returns:
            trimmedFASTA_obj -- a NEW ZS_SeqIO.FASTA object containing the sequences after trimming.
        '''
        assert FASTA_obj.isAligned, "ERROR: MSA.trim() cannot trim input FASTA object since .isAligned is False"
        assert 0 < pctTrim < 1, "ERROR: MSA.trim() pctTrim must be a float between 0 and 1 exclusive"
        
        # Find the start point where pctTrim is satisfied
        for startIndex in range(len(FASTA_obj.seqs[0].gap_seq)): # all seqs should be same length
            columnResidues = [ FastASeq_obj.gap_seq[startIndex] for FastASeq_obj in FASTA_obj.seqs ]
            pctResidues = 1 - (columnResidues.count('-') / len(columnResidues))
            if pctResidues <= pctTrim:
                continue
            break
        
        # Find the end point similarly
        for endIndex in range(len(FASTA_obj.seqs[0].gap_seq)-1, 0, -1):
            columnResidues = [ FastASeq_obj.gap_seq[endIndex] for FastASeq_obj in FASTA_obj.seqs ]
            pctResidues = 1 - (columnResidues.count('-') / len(columnResidues))
            if pctResidues <= pctTrim:
                continue
            break
        
        # Check our values to ensure they're sensible
        if startIndex >= endIndex: # if we'd trim everything
            raise Exception("ERROR: MSA.trim() pctTrim is too high; would trim all sequences")
        
        # Trim our FASTA object
        trimmedFASTA_obj = FASTA_obj.slice_cols(startIndex, endIndex)
        
        return trimmedFASTA_obj
    
    def drop_gappy_seqs(FASTA_obj, allowedGappiness=0.50, inPlace=False):
        '''
        Receives an aligned ZS_SeqIO.FASTA object and drops sequences that have a proportion
        of gaps above the allowedGappiness threshold. The threshold is calculated relative to
        the consensus sequence of the alignment.
        
        Parameters:
            FASTA_obj -- a ZS_SeqIO.FASTA object containing aligned sequences.
            allowedGappiness -- OPTIONAL; a float from 0->1 (but less than 1) indicating the maximum allowed
                                proportion of residues that can be gaps before a sequence will be dropped; value
                                is calculated relative to the consensus sequence of the alignment; default == 0.50
                                i.e., up to 50% of residues can be gaps with any excess causing dropping.
            inPlace -- OPTIONAL; a boolean indicating whether to modify the input FASTA object
                       (True) or return a new FASTA object with the empty columns removed (False);
                       default == False.
        Returns:
            prunedFASTA_obj -- a new or modified ZS_SeqIO.FASTA object with sequences potentially dropped if they
                               did not meet the gappyPct threshold; whether it's new or modified depends on inPlace.
        '''
        assert FASTA_obj.isAligned, "ERROR: MSA.drop_gappy_seqs() cannot process input FASTA object since .isAligned is False"
        assert 0 <= allowedGappiness < 1, "ERROR: MSA.drop_gappy_seqs() gappyPct must be a float >=0 and <1"
        
        # Create a deepcopy of the FASTA object to work with
        prunedFASTA_obj = deepcopy(FASTA_obj)
        
        # Generate a consensus sequence
        consensusSeq = prunedFASTA_obj.generate_consensus()
        
        # Check each sequence against the consensus for gappiness
        toDelete = []
        for index, FastASeq_obj in enumerate(prunedFASTA_obj.seqs):
            nonGapResidues = sum([1 for i in range(len(FastASeq_obj.gap_seq)) if FastASeq_obj.gap_seq[i] != "-"])
            gappyPct = 1 - (nonGapResidues / len(FastASeq_obj.gap_seq))
            if gappyPct > allowedGappiness:
                toDelete.append(index)
        
        # Check our values to ensure they're sensible
        if len(toDelete) == len(prunedFASTA_obj): # if we'd drop everything
            raise Exception("ERROR: MSA.drop_gappy_seqs() allowedGappiness is too low; would drop all sequences")
        
        # Drop the sequences that didn't meet the threshold
        for index in toDelete[::-1]: # reverse order to avoid index shifting
            del prunedFASTA_obj.seqs[index]
        
        # If we dropped sequences, make sure we don't have any blank columns now
        if toDelete != []:
            prunedFASTA_obj = MSA.drop_empty_columns(prunedFASTA_obj, inPlace=True)
        
        return prunedFASTA_obj

    def drop_empty_columns(FASTA_obj, inPlace=False):
        '''
        Parameters:
            FASTA_obj -- a ZS_SeqIO.FASTA object containing aligned sequences.
            inPlace -- OPTIONAL; a boolean indicating whether to modify the input FASTA object
                       (True) or return a new FASTA object with the empty columns removed (False);
                       default == False.
        Returns:
            cleanedFASTA_obj -- a new or modified ZS_SeqIO.FASTA object with any empty columns removed from each
                                FastASeq's .gap_seq attribute; whether it's new or modified depends on inPlace.
        '''
        assert FASTA_obj.isAligned, "ERROR: MSA.drop_empty_columns() cannot process input FASTA object since .isAligned is False"
        
        # Create a deepcopy of the FASTA object to work with (if applicable)
        if inPlace == False:
            cleanedFASTA_obj = deepcopy(FASTA_obj)
        else:
            cleanedFASTA_obj = FASTA_obj
        
        # Find empty columns
        colsToDelete = []
        for index in range(len(cleanedFASTA_obj.seqs[0].gap_seq)):
            columnResidues = [ FastASeq_obj.gap_seq[index] for FastASeq_obj in cleanedFASTA_obj.seqs ]
            if set(columnResidues) == {"-"}:
                colsToDelete.append(index)
        
        # Fix the sequences
        for index in colsToDelete[::-1]: # reverse order to avoid index shifting
            for FastASeq_obj in cleanedFASTA_obj.seqs:
                FastASeq_obj.gap_seq = FastASeq_obj.gap_seq[:index] + FastASeq_obj.gap_seq[index+1:]
        
        return cleanedFASTA_obj
    
class MAFFT:
    '''
    A mostly fully features wrapper program for MAFFT. Rather than allowing use of the mafft.bat
    on Windows, this class instead assumes that MAFFT is available through WSL.
    '''
    ALG_CONVERT = {
            "auto": ["--auto"],
            "einsi": ["--genafpair"],
            "linsi": ["--localpair"],
            "ginsi": ["--globalpair"],
            "fftns1": ["--retree", "1"],
            "fftns2": ["--retree", "2"],
            "fftnsi": None
        }
    ACCEPTED_ALGS = ["auto", "einsi", "linsi", "ginsi", "fftns1", "fftns2", "fftnsi"]
    ACCEPTED_MOLECULES = ["auto", "nucleotide", "protein"]
    
    def __init__(self, mafftPath, algorithm="auto", thread=1, maxiterate=0, molecule="auto"):
        self.exe = mafftPath
        self.algorithm = algorithm
        self.thread = thread
        self.maxiterate = maxiterate
        self.molecule = molecule
    
    @property
    def exe(self):
        return self._exe
    
    @exe.setter
    def exe(self, value):
        convertedValue = ZS_Utility.convert_to_wsl_if_not_unix(value)
        assert ZS_Utility.wsl_isfile(convertedValue), \
            f"MAFFT executable not found at '{convertedValue}' after WSL compatibility conversion"
        self._exe = convertedValue
    
    @property
    def algorithm(self):
        return self._algorithm
    
    @algorithm.setter
    def algorithm(self, value):
        assert value in MAFFT.ACCEPTED_ALGS, \
            f"algorithm '{value}' not recognized; please choose from {MAFFT.ACCEPTED_ALGS}"
        self._algorithm = value
    
    @property
    def molecule(self):
        return self._molecule
    
    @molecule.setter
    def molecule(self, value):
        assert value in MAFFT.ACCEPTED_MOLECULES, \
            f"molecule '{value}' not recognized; please choose from {MAFFT.ACCEPTED_MOLECULES}"
        self._molecule = value
    
    @property
    def thread(self):
        return self._thread
    
    @thread.setter
    def thread(self, value):
        assert isinstance(value, int), "thread must be an integer"
        assert value > 0, "thread must be greater than 0"
        self._thread = value
    
    @property
    def maxiterate(self):
        return self._maxiterate
    
    @maxiterate.setter
    def maxiterate(self, value):
        assert isinstance(value, int), "maxiterate must be an integer"
        assert value >= 0, "maxiterate must be greater than or equal to 0"
        self._maxiterate = value
    
    def _align(self, fastaFileName, reorder=False, keeplength=False):
        '''
        Simple wrapper function to run MAFFT on a FASTA file and write output
        to the specified file. For a fully-features handler of ZS_SeqIO.FASTA objects
        and such, you should wrap something around this wrapper.
        
        Parameters:
            fastaFileName -- a string indicating the location of a FASTA file to align
            reorder -- OPTIONAL; a boolean indicating whether to pass along the --reorder
                       arg to MAFFT. Default == False.
            keeplength -- OPTIONAL; a boolean indicating whether to pass along the
                          --keeplength arg to MAFFT. Default == False.
        Returns:
            msa -- a string indicating the multiple sequence alignment result as derived
                   directly from MAFFT stdout
        '''
        assert os.path.isfile(fastaFileName), \
            f"ERROR: MAFFT._align() could not find FASTA file '{fastaFileName}'"
        
        # Construct the cmd for subprocess
        cmd = ZS_Utility.base_subprocess_cmd(self.exe)
        
        # Handle algorithm specification
        if self.algorithm != "fftnsi":
            cmd += MAFFT.ALG_CONVERT[self.algorithm]
        else:
            assert self.maxiterate >= 2, "maxiterate must be >= 2 for fftnsi"
        
        # Handle maxiterate specification
        if self.algorithm != "auto":
            cmd += ["--maxiterate", str(self.maxiterate)]
        
        # Handle molecule specification
        if self.molecule != "auto":
            if self.molecule == "nucleotide":
                cmd += ["--nuc"]
            else:
                cmd += ["--amino"]
        
        # Handle optional flags
        if reorder:
            cmd.append("--reorder")
        if keeplength:
            cmd.append("--keeplength")
        
        # Handle threads and files to align
        cmd += [
            "--quiet", # necessary to prevent unusual MAFFT errors; see https://mafft.cbrc.jp/alignment/software/stderr.html
            "--thread", str(self.thread),
            ZS_Utility.convert_to_wsl_if_not_unix(fastaFileName)
        ]
        
        if platform.system() != "Windows":
            cmd = " ".join(cmd)
        
        # Run the command
        run_mafft = subprocess.Popen(cmd, shell = True,
                                     stdout = subprocess.PIPE,
                                     stderr = subprocess.PIPE)
        mafftout, maffterr = run_mafft.communicate()
        
        # Check to see if there was an error
        if (not mafftout.decode("utf-8").startswith(">")):
            raise Exception(("ERROR: MAFFT._align() encountered an error; have a look " +
                            f'at the stdout ({mafftout.decode("utf-8")}) and stderr ' + 
                            f'({maffterr.decode("utf-8")}) to make sense of this.'))
        
        # Parse out the MSA result
        msa = mafftout.decode("utf-8").replace("\r", "").rstrip("\n ")
        return msa
    
    def align(self, fasta):
        '''
        Handles the execution of MAFFT alignment from start to finish. It will:
        
            1) Convert the fasta value into a file to align (if you've provided 
            a ZS_SeqIO.FASTA object)
            2) Perform the alignment with MAFFT
            3) Parse the MSA output into a ZS_SeqIO.FASTA object to return
        
        Before calling this, you might consider making sure all values are set
        appropriately e.g., algorithm, thread, and maxiterate.
        
        Parameters:
            fasta -- a string indicating a file path, OR an object of ZS_SeqIO.FASTA class
        Returns:
            msa -- a ZS_SeqIO.FASTA object containing the aligned sequences
        '''
        # Validate input value type
        if isinstance(fasta, str):
            assert os.path.isfile(fasta), f"ERROR: MAFFT.align() could not find the FASTA file '{fasta}'"
            fastaFileName, fastaIsTemporary = fasta, False
        elif hasattr(fasta, "isFASTA") and fasta.isFASTA is True:
            fastaFileName, fastaIsTemporary = ZS_SeqIO.Conversion.get_filename_for_input_sequences(fasta)
        else:
            raise Exception(f"ERROR: MAFFT.align() requires a FASTA file or FASTA object as input; did not understand '{fasta}'")
        
        # Run MAFFT alignment
        msaString = self._align(fastaFileName)
        
        # Parse output into a FASTA object
        FASTA_obj = ZS_SeqIO.FASTA(msaString, isAligned=True)
        
        # Clean up temporary file
        if fastaIsTemporary:
            os.unlink(fastaFileName)
        
        return FASTA_obj
    
    def _add(self, originalFile, toAddFile, addType, reorder=False, keeplength=False):
        '''
        Simple wrapper function to run MAFFT to add a FASTA file into an existing alignment
        FASTA. Intended to be wrapped by a surrounding function to handle file format conversions.
        
        Parameters:
            originalFile -- a string indicating the location of the original FASTA file
                            to have sequenced added into
            toAddFile -- a string indicating the location of the FASTA file to add into the
                         originalFile
            addType -- a string indicating the type of adding to perform; must be "add", 
                       "addfragments", or "addlong"
            reorder -- OPTIONAL; a boolean indicating whether to pass along the --reorder
                       arg to MAFFT. Default == False.
            keeplength -- OPTIONAL; a boolean indicating whether to pass along the
                          --keeplength arg to MAFFT. Default == False.
        Returns:
            msa -- a string indicating the multiple sequence alignment result as derived
                   directly from MAFFT stdout
        '''
        assert os.path.isfile(originalFile), \
            f"ERROR: MAFFT._add() could not find FASTA file '{originalFile}'"
        assert os.path.isfile(toAddFile), \
            f"ERROR: MAFFT._add() could not find FASTA file '{toAddFile}'"
        assert addType in ["add", "addfragments", "addlong"], \
            f"ERROR: MAFFT._add() addType '{addType}' not recognized"
        
        # Construct the cmd for subprocess
        cmd = ZS_Utility.base_subprocess_cmd(self.exe)
        
        # Handle algorithm specification
        if self.algorithm != "fftnsi":
            cmd += MAFFT.ALG_CONVERT[self.algorithm]
        elif self.algorithm == "einsi":
            raise NotImplementedError(("MAFFT._add() does not support E-INSi for " + 
                                       "adding sequences since MAFFT itself does not!"))
        else:
            assert self.maxiterate >= 2, "maxiterate must be >= 2 for fftnsi"
        
        # Handle maxiterate specification
        if self.algorithm != "auto":
            cmd += ["--maxiterate", str(self.maxiterate)]
        
        # Handle optional flags
        if reorder:
            cmd.append("--reorder")
        if keeplength:
            cmd.append("--keeplength")
        
        # Handle threads and files to align
        cmd += [
            "--quiet",
            "--thread", str(self.thread),
            f"--{addType}", ZS_Utility.convert_to_wsl_if_not_unix(toAddFile),
            ZS_Utility.convert_to_wsl_if_not_unix(originalFile)
        ]
        
        if platform.system() != "Windows":
            cmd = " ".join(cmd)
        
        # Run the command
        run_mafft = subprocess.Popen(cmd, shell = True,
                                     stdout = subprocess.PIPE,
                                     stderr = subprocess.PIPE)
        mafftout, maffterr = run_mafft.communicate()
        
        # Check to see if there was an error
        if (not mafftout.decode("utf-8").startswith(">")):
            raise Exception(("ERROR: MAFFT._add() encountered an error; have a look " +
                            f'at the stdout ({mafftout.decode("utf-8")}) and stderr ' + 
                            f'({maffterr.decode("utf-8")}) to make sense of this.'))
        
        # Parse out the MSA result
        msa = mafftout.decode("utf-8").replace("\r", "").rstrip("\n ")
        return msa
    
    def add(self, originalFasta, toAddFasta, addType, reorder=False, keeplength=False):
        '''
        Handles the execution of MAFFT alignment from start to finish, but with
        the intention of adding sequences into an existing alignment. It will:
        
            1) Convert the fasta values into files to align (if you've provided 
            a ZS_SeqIO.FASTA object(s))
            2) Perform the alignment with MAFFT
            3) Parse the MSA output into a ZS_SeqIO.FASTA object to return
        
        Before calling this, you might consider making sure all values are set
        appropriately e.g., algorithm, thread, and maxiterate.
        
        Parameters:
            originalFasta -- a string indicating a file path, OR an object of ZS_SeqIO.FASTA class
                             to have sequences added into
            toAddFasta -- a string indicating a file path, OR an object of ZS_SeqIO.FASTA class
                          to add into the originalFasta
            addType -- a string indicating the type of adding to perform; must be "add",
                       "addfragments", or "addlong"
        Returns:
            msa -- a ZS_SeqIO.FASTA object containing the aligned sequences
        '''
        # Validate input value types
        assert addType in ["add", "addfragments", "addlong"], \
            f"ERROR: MAFFT.add() addType '{addType}' not recognized"
        
        if isinstance(originalFasta, str):
            assert os.path.isfile(originalFasta), \
                f"ERROR: MAFFT.add() could not find the FASTA file '{originalFasta}'"
            originalFileName, originalIsTemporary = originalFasta, False
        elif hasattr(originalFasta, "isFASTA") and originalFasta.isFASTA is True:
            originalFileName, originalIsTemporary = ZS_SeqIO.Conversion.get_filename_for_input_sequences(originalFasta, asAligned=True)
        else:
            raise Exception(f"ERROR: MAFFT.run() requires a FASTA file or FASTA object as input; did not understand '{originalFasta}'")
        
        try:
            if isinstance(toAddFasta, str):
                assert os.path.isfile(toAddFasta), \
                    f"ERROR: MAFFT.add() could not find the FASTA file '{toAddFasta}'"
                toAddFileName, toAddIsTemporary = toAddFasta, False
            elif hasattr(toAddFasta, "isFASTA") and toAddFasta.isFASTA is True:
                toAddFileName, toAddIsTemporary = ZS_SeqIO.Conversion.get_filename_for_input_sequences(toAddFasta)
            else:
                raise Exception(f"ERROR: MAFFT.run() requires a FASTA file or FASTA object as input; did not understand '{toAddFasta}'")
        except Exception as e:
            # Clean up first file conversion if necessary
            if originalIsTemporary: 
                os.unlink(originalFileName)
            # Now allow error to propagate
            raise e
        
        # Run MAFFT add alignment
        try:
            msaString = self._add(originalFileName, toAddFileName, addType, reorder, keeplength)
        except Exception as e:
            # Clean up file conversions if necessary
            if originalIsTemporary: 
                os.unlink(originalFileName)
            if toAddIsTemporary:
                os.unlink(toAddFileName)
            # Now allow error to propagate
            raise e
        
        # Parse output into a FASTA object
        FASTA_obj = ZS_SeqIO.FASTA(msaString, isAligned=True)
        
        # Clean up temporary files
        if originalIsTemporary: 
            os.unlink(originalFileName)
        if toAddIsTemporary:
            os.unlink(toAddFileName)
        
        return FASTA_obj
    
    def align_as_protein(self, fasta, strands=None, frames=None):
        '''
        Handles the execution of MAFFT alignment from start to finish, but by first
        translating nucleotides into proteins, and then de-translating them. It will:
        
            1) Translate nucleotides into their corresponding protein sequence
            2) Convert the fasta value into a file to align (if you've provided 
            a ZS_SeqIO.FASTA object)
            3) Perform the alignment with MAFFT
            4) De-translate sequences back into nucleotides
            5) Parse the MSA output into a ZS_SeqIO.FASTA object to return
        
        In effect, it allows nucleotide sequences to be aligned as codons rather than
        as individual base pairs, which might be more biologically relevant when dealing
        with closely related sequences.
        
        Parameters:
            fasta -- a string indicating a file path, OR an object of ZS_SeqIO.FASTA class
            strands -- OPTIONAL; a list of integers indicating the strand(s) to use for
                       for optimal translation of each sequence. Must be the same length
                       as the number of sequences in the FASTA object. If not provided,
                       the program will attempt to find the best strand for each sequence.
            frames -- OPTIONAL; a list of integers indicating the frame(s) to use for
                      for optimal translation of each sequence. Must be the same length
                      as the number of sequences in the FASTA object. If not provided,
                      the program will attempt to find the best frame for each sequence.
        Returns:
            msa -- a new ZS_SeqIO.FASTA object containing the aligned sequences if you
                   specified a file, or a modified version of the original FASTA object
                   (if you provided such an object type)
        '''
        # Validate input value type
        if isinstance(fasta, str):
            assert os.path.isfile(fasta), f"ERROR: MAFFT.run() could not find the FASTA file '{fasta}'"
            FASTA_obj = ZS_SeqIO.FASTA(fasta, isAligned=False)
        elif hasattr(fasta, "isFASTA") and fasta.isFASTA is True:
            FASTA_obj = fasta
        else:
            raise Exception(("ERROR: MAFFT.align_as_protein() requires a FASTA file or FASTA " + 
                             f"object as input; did not understand '{fasta}'"))
        
        if strands != None:
            assert len(strands) == len(FASTA_obj), \
                "strands must be the same length as the number of sequences in the FASTA object"
        if frames != None:
            assert len(frames) == len(FASTA_obj), \
                "frames must be the same length as the number of sequences in the FASTA object"
        
        # Create a new FASTA object with our translations
        forAligningFASTA = ZS_SeqIO.FASTA(None)
        newStrands = []
        newFrames = []
        for i in range(len(FASTA_obj)):
            FastASeq_obj = FASTA_obj.seqs[i]
            
            # Run the translation
            strand=None if strands == None else strands[i] if strands[i] != None else None
            frame=None if strand == None else None if frames == None else frames[i] if frames[i] != None else None
            findBestFrame=True if (strand == None or frame == None) else False
            
            protein, _strand, _frame = FastASeq_obj.get_translation(
                findBestFrame=findBestFrame,
                strand=strand,
                frame=frame)
            
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
            newStrands.append(_strand)
            newFrames.append(_frame)
            forAligningFASTA.add(ZS_SeqIO.FastASeq(FastASeq_obj.id, seq=protein))
        
        # Run MAFFT alignment
        tempFasta, _ = ZS_SeqIO.Conversion.get_filename_for_input_sequences(forAligningFASTA)
        try:
            msaString = self._align(os.path.abspath(tempFasta))
        except Exception as e:
            os.unlink(tempFasta)
            raise e
        os.unlink(tempFasta)
        
        # Parse alignment back into a FASTA object
        alignedFASTA_obj = ZS_SeqIO.FASTA(msaString, isAligned=True)
        
        # Map back aligned proteins to their original nucleotide counterparts
        trimmedBits = []
        for i in range(len(alignedFASTA_obj)):
            _strand = newStrands[i]
            _frame = newFrames[i]
            nuc = FASTA_obj[i].seq if _strand == 1 else FASTA_obj[i].get_reverse_complement()
            alignedProt = alignedFASTA_obj[i].gap_seq
            
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
        
        return FASTA_obj
    
    def __repr__(self):
        return "<MAFFT(exe='{0}',algorithm='{1}',thread={2},maxiterate={3})>".format(
            self.exe, self.algorithm, self.thread, self.maxiterate
        )
    
    def __str__(self):
        return "MAFFT(exe='{0}',algorithm='{1}',thread={2},maxiterate={3})".format(
            self.exe, self.algorithm, self.thread, self.maxiterate
        )

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
    def ssw_parasail(queryString, targetString, molecule, openPenalty=10, extendPenalty=1):
        '''
        molecule -- a string equal to "nucleotide" or "protein"; determines which
                        substitution matrix to use
        '''
        assert molecule in ["nucleotide", "protein"], \
            "ssw_parasail must be told what type of molecule it is aligning!"
        
        alignment = parasail.sw_trace_striped_sat(queryString, targetString,
                                                  openPenalty, extendPenalty,
                                                  parasail.nuc44 if molecule == "nucleotide"
                                                  else parasail.blosum62
                                                  )
        queryAlign = alignment.traceback.query
        targetAlign = alignment.traceback.ref
        
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
            self._exonerateExe = ZS_Utility.convert_windows_to_wsl_path(value)
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
    def filter_exonerate_resultsDict(resultsDict, num_hits=0, identity=0.0, similarity=0.0):
        '''
        Utility function to receive a resultsDict as produced by run_exonerate() and
        filter results to only maintain relevant ones. It will also sort the results
        in descending order of their 1) similarity, 2) identity, and 3) sequence length.
        
        It will set a .coding_length attribute on the gene and mRNA subfeature for each
        feature, which is used for filtering here and can be useful for external, specialisd
        filtering.
        
        And, it will automatically remove any features which lack CDS prediction, if that's
        even possible for exonerate to do...?
        
        Parameters:
            resultsDict -- a dictionary associating ZS_GFF3IO.Feature objects to their
                           query sequence with structure like:
                           {
                               'querySeqID1': [feature, feature, ...],
                               'querySeqID2': [ ... ],
                               ...
                           }
            num_hits -- an integer >= 0 indicating how many of the "top hits" we want to retain;
                        a value of 0 means there is no maximum
            identity -- a float or integer in the range 0->100 (inclusive) setting the minimum
                        identity value we will retain
            similarity -- a float or integer in the range 0->100 (inclusive) setting the minimum
                          similarity value we will retain
        Returns:
            newResultsDict -- a new dictionary with the same structure as the input parameter,
                              but with only results which pass filtration remaining.
        '''
        # Validate input types
        assert isinstance(num_hits, int) and num_hits >= 0, \
            "num_hits must be an integer >= 0"
        assert (isinstance(identity, int) or isinstance(identity, float)) and 100 >= identity >= 0, \
            "identity must be a float or integer from 0 to 100"
        assert (isinstance(similarity, int) or isinstance(similarity, float)) and 100 >= similarity >= 0, \
            "similarity must be a float or integer from 0 to 100"
        
        newResultsDict = {}
        for queryID, features in resultsDict.items():
            # Filter
            for feature in features:
                if float(feature.identity) >= identity and float(feature.similarity) >= similarity and hasattr(feature.mRNA[0], "CDS"):
                    # Calculate coding length
                    cdsCoords = GFF3._get_feature_coords(feature.mRNA[0], "CDS")
                    cdsLength = 0
                    for start, end in cdsCoords[0]:
                        cdsLength += (end - start + 1)
                    feature.coding_length = cdsLength
                    feature.mRNA[0].coding_length = cdsLength
                    
                    # Store in dictionary
                    newResultsDict.setdefault(queryID, [])
                    newResultsDict[queryID].append(feature)
            
            # Sort
            if queryID in newResultsDict:
                newResultsDict[queryID].sort(key = lambda x: 
                    (-float(x.identity), -float(x.similarity), -x.coding_length)
                )
            
            # Retain top num_hit results
            if queryID in newResultsDict:
                if num_hits > 0:
                    newResultsDict[queryID] = newResultsDict[queryID][0:num_hits]
        return newResultsDict
    
    def _format_exonerate_cmd(self, queryFile, targetFile):
        '''
        Hidden helper function for getting a list amenable to subprocess.run()
        that sets all relevant cmd tags depending on this object's properties.
        It expects that the query and target files have already been subjected
        to Conversion.get_filename_for_input_sequences(), and their return values
        are to be given to this function.
        
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
                ZS_Utility.convert_windows_to_wsl_path(os.path.abspath(queryFile)),
                ZS_Utility.convert_windows_to_wsl_path(os.path.abspath(targetFile))
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
        q, qIsTemporary = ZS_SeqIO.Conversion.get_filename_for_input_sequences(self.query)
        t, tIsTemporary = ZS_SeqIO.Conversion.get_filename_for_input_sequences(self.target)
        
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
            resultsDict -- a dictionary associating ZS_GFF3IO.Feature objects to their
                           query sequence with structure like:
                           {
                               'querySeqID1': [feature, feature, ...],
                               'querySeqID2': [ ... ],
                               ...
                           }
        '''
        featureDict = {"gene": {}, "mRNA": {}}
        geneIDDict = {} # just for naming purposes for this script
        for line in exonerateGffStdout.split("\n"):
            # Skip irrelevant lines
            try:
                contig, source, featureType, start, end, \
                        score, strand, frame, attributes \
                        = line.rstrip('\n').split('\t')
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
        
        # Disentangle results so that the originating sequence can be easily found for each feature
        resultsDict = {}
        for feature in featureDict["gene"].values():
            resultsDict.setdefault(feature.Sequence, [])
            resultsDict[feature.Sequence].append(feature)
        
        return resultsDict
    
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

class GMAP:
    '''
    The GMAP Class is intended to be used alongside the ZS_SeqIO.FASTA Class to
    handle the alignment of FASTA sequences with OOP magics. In short, it provides
    easy access to GMAP to perform various forms of alignment, returning a 
    parsed result.
    
    The gmapDir parameter expects a directory to be passed that contains the gmap
    and gmap_build executable files.
    '''
    def __init__(self, gmapDir, query, target):
        self.gmapDir = gmapDir
        self.query = query
        self.target = target
        self.dbWasCreated = False
        
        # Set default GMAP parameters
        self.npaths = 5 # GMAP default
        self.batch = 5 # GMAP default is 2, but I always use this
        self.max_intronlength_middle = 500000 # GMAP default
        self.max_intronlength_ends = 500000 # GMAP default is 10000, but I always use this
        self.nthreads = 1
    
    @property
    def gmapDir(self):
        return self._gmapDir
    
    @gmapDir.setter
    def gmapDir(self, value):
        assert isinstance(value, str) and os.path.isdir(value), \
            f"'{value}' is not a string or does not point to a directory"
        assert os.path.isfile(os.path.join(value, "gmap")), \
            f"gmap executable not found at location'{value}'"
        assert os.path.isfile(os.path.join(value, "gmap_build")), \
            f"gmap_build executable not found at location'{value}'"
        
        if platform.system() == "Windows":
            print("GMAP Class on Windows assumes GMAP was built using WSL; carrying on...")
            self._gmapDir = ZS_Utility.convert_windows_to_wsl_path(value)
        else:
            self._gmapDir = value
    
    @property
    def query(self):
        return self._query
    
    @query.setter
    def query(self, value):
        recognisedTypes = ["str", "FASTA", "ZS_SeqIO.FASTA", "FastASeq", "ZS_SeqIO.FastASeq"]
        assert type(value).__name__ in recognisedTypes, \
            "Query value type not handled by GMAP class"
        
        # Derive an existing file string, or non-existing file FASTA
        tmpQuery, _ = ZS_SeqIO.Conversion._intermediate_conversion(value) # throw away tmpHash return
        
        # If it's an existing file string, store it as-is and mark it as non-temporary
        if isinstance(tmpQuery, str): # if this is True, os.path.isfile(tmpQuery) is implicitly also True
            self._query = tmpQuery
            self._queryFile = value
            self._queryIsTemporary = False
        
        # Otherwise, store the FASTA object, create a file for it, and mark it as temporary
        else:
            self._query = tmpQuery 
            self._queryFile, self._queryIsTemporary = _create_file_from_FASTA_object(tmpQuery, useExistingFile=True)
    
    @property
    def target(self):
        return self._target
    
    @target.setter
    def target(self, value):
        recognisedTypes = ["str", "FASTA", "ZS_SeqIO.FASTA", "FastASeq", "ZS_SeqIO.FastASeq"]
        assert type(value).__name__ in recognisedTypes, \
            "Query value type not handled by GMAP class"
        
        # Derive an existing file string, or non-existing file FASTA
        tmpTarget, _ = ZS_SeqIO.Conversion._intermediate_conversion(value) # throw away tmpHash return
        
        # If it's an existing file string, store it as-is and mark it as non-temporary
        if isinstance(tmpTarget, str): # if this is True, os.path.isfile(tmpTarget) is implicitly also True
            self._target = tmpTarget
            self._targetFile = value
            self._targetIsTemporary = False
        
        # Otherwise, store the FASTA object, create a file for it, and mark it as temporary
        else:
            self._target = tmpTarget 
            self._targetFile, self._targetIsTemporary = _create_file_from_FASTA_object(tmpTarget, useExistingFile=True)
    
    @property
    def queryFile(self):
        return self._queryFile
    
    @property
    def targetFile(self):
        return self._targetFile
    
    @property
    def dbWasCreated(self):
        return self._dbWasCreated
    
    @dbWasCreated.setter
    def dbWasCreated(self, value):
        assert isinstance(value, bool), \
            "dbWasCreated value must be a boolean"
        
        self._dbWasCreated = value
    
    @property
    def npaths(self):
        return self._npaths
    
    @npaths.setter
    def npaths(self, value):
        assert isinstance(value, int), \
            "npaths value must be an integer"
        assert 0 <= value, \
            "npaths value must be greater than or equal to zero"
        
        self._npaths = value
    
    @property
    def batch(self):
        return self._batch
    
    @batch.setter
    def batch(self, value):
        assert isinstance(value, int), \
            "batch value must be an integer"
        assert 0 <= value <= 5, \
            "batch value must be in the range of 0->5"
        
        self._batch = value
    
    @property
    def max_intronlength_middle(self):
        return self._max_intronlength_middle
    
    @max_intronlength_middle.setter
    def max_intronlength_middle(self, value):
        assert isinstance(value, int), \
            "max_intronlength_middle value must be an integer"
        assert 0 <= value, \
            "max_intronlength_middle value must be greater than or equal to zero"
        
        self._max_intronlength_middle = value
    
    @property
    def max_intronlength_ends(self):
        return self._max_intronlength_ends
    
    @max_intronlength_ends.setter
    def max_intronlength_ends(self, value):
        assert isinstance(value, int), \
            "max_intronlength_ends value must be an integer"
        assert 0 <= value, \
            "max_intronlength_ends value must be greater than or equal to zero"
        
        self._max_intronlength_ends = value
    
    @property
    def nthreads(self):
        return self._nthreads
    
    @nthreads.setter
    def nthreads(self, value):
        assert isinstance(value, int), \
            "nthreads value must be an integer"
        assert 0 <= value, \
            "nthreads value must be greater than or equal to zero"
        
        self._nthreads = value
    
    @property
    def targetDB(self):
        return f"{self.targetFile}.gmap"
    
    def clean(self, query=False, target=False, db=False):
        '''
        This Class creates several temporary files. Call this method when you're done
        using it to remove the left over junk.
        
        Parameters:
            query -- a boolean indicating whether the query file should be cleaned up
                     (if applicable); if the query file was not made by this script,
                     it will NOT be removed
            target -- a boolean indicating whether the target file should be cleaned up
                      (if applicable); if the target file was not made by this script,
                      it will NOT be removed
            db -- a boolean indicating whether the gmap_build directory should be cleaned up
                  (if applicable); if this directory was not made by this script,
                  it will NOT be removed
        '''
        if query is True and self._queryIsTemporary is True:
            os.unlink(self._queryFile)
        
        if target is True and self._targetIsTemporary is True:
            os.unlink(self._targetFile)
        
        if db is True and self.dbWasCreated is True:
            shutil.rmtree(self.targetDB)
    
    def _format_build_cmd(self):
        '''
        Hidden helper function for getting a list amenable to subprocess.run()
        that sets all relevant cmd tags depending on this object's properties in
        order to run gmap_build and index a target FASTA file.
        
        Returns:
            cmds -- a list amenable to subprocess.run()
        '''
        cmds = []
        if platform.system() == "Windows":
            cmds += ["wsl", "~", "-e"]
        
        cmds += [
            Path(self.gmapDir, "gmap_build").as_posix(),
            
            "-D", ZS_Utility.convert_windows_to_wsl_path(os.path.abspath(os.path.dirname(self.targetFile))) \
                if platform.system() == "Windows" else os.path.dirname(self.targetFile),
            
            "-d", f"{os.path.basename(self.targetDB)}",
            
            ZS_Utility.convert_windows_to_wsl_path(os.path.abspath(self.targetFile)) if platform.system() == "Windows" \
                else self.targetFile
        ]
        
        return cmds
    
    def _target_build_exists(self):
        '''
        Hidden helper function for determining whether gmap_build has already been run
        on the .target value.
        
        Returns:
            buildExists -- a boolean indicating whether a gmap_build result already exists
                           (True) or not (False)
        '''
        return os.path.isdir(self.targetDB)
    
    def gmap_build(self):
        '''
        Makes a GMAP database out of the .target value. Skips doing so if the output
        folder already exists
        '''
        # Check if calling this function is necessary
        if self._target_build_exists():
            print("gmap_build results already exist")
            return
        
        # Run GMAP build
        cmds = self._format_build_cmd()
        gmapBuild = subprocess.run(cmds, capture_output=True)
        
        # Validate that build was successful
        stderr = gmapBuild.stderr.decode()
        assert stderr.rstrip("\n").split("\n")[-1] == "Done", \
            f"Unexpected error occuring during gmap_build run; stderr == {stderr}"
        
        # Set flag so we can clean up at the end if necessary
        self.dbWasCreated = True
    
    def _format_search_cmd(self):
        '''
        Hidden helper function for getting a list amenable to subprocess.run()
        that sets all relevant cmd tags depending on this object's properties in
        order to run GMAP and get GFF3 formatted results.
        
        Returns:
            cmds -- a list amenable to subprocess.run() that will perform a GMAP search
        '''
        cmds = []
        if platform.system() == "Windows":
            cmds += ["wsl", "~", "-e"]
        
        cmds += [
            Path(self.gmapDir, "gmap").as_posix(),
            
            "-D", ZS_Utility.convert_windows_to_wsl_path(os.path.dirname(os.path.abspath(self.targetDB))) \
                if platform.system() == "Windows" else os.path.dirname(os.path.abspath(self.targetDB)),
            
            "-d", f"{os.path.basename(self.targetDB)}",
            
            "-f", "2",
            "-n", str(self.npaths),
            "-t", str(self.nthreads),
            "-B", str(self.batch),
            
            f"--max-intronlength-middle={self.max_intronlength_middle}",
            
            f"--max-intronlength-ends={self.max_intronlength_ends}",
            
            ZS_Utility.convert_windows_to_wsl_path(os.path.abspath(self.queryFile)) if platform.system() == "Windows" \
                else os.path.abspath(self.queryFile)
        ]
        
        return cmds
    
    def run_gmap(self):
        '''
        Runs GMAP using the parameters set in this object.
        
        Returns:
            gmapGFF3 -- a ZS_GFF3IO.GFF3 object containing all results reported by
                        GMAP.
        '''
        # Make sure target file is ready for GMAP
        if not self._target_build_exists():
            self.gmap_build()
        
        # Get temporary results file name
        tmpHash = ZS_SeqIO.Conversion.get_hash_for_input_sequences(self._queryFile + self._targetFile,
                                                                     randomHash=True, maxLength=20)
        tmpResultsName = ZS_Utility.tmp_file_name_gen(f"alignIO_tmp_" + tmpHash, "gff3")
        
        # Get cmds and run
        cmds = self._format_search_cmd()
        with open(tmpResultsName, "w") as fileOut:
            gmap = subprocess.run(cmds, stderr=subprocess.PIPE, stdout=fileOut)
        
        # Validate that search was successful
        stderr = gmap.stderr.decode()
        assert re.findall(r"Processed \d{1,10} queries in", stderr) != [], \
            f"Unexpected error occuring during gmap run; stderr == {stderr}"
        
        # Parse result & clean up
        gmapGFF3 = GFF3(tmpResultsName)
        os.unlink(tmpResultsName)
        
        return gmapGFF3

class Minimap2:
    '''
    The Minimap2 Class provides access to the minimap2 search functionality.
    
    Attributes:
        query (REQUIRED) -- a string pointing to a FASTA file to query against the target.
        target (REQUIRED) -- a string point to a FASTA file to search against/be the reference.
        preset (REQUIRED) -- a string value specifying the preset to use for the search.
        minimap2Exe (REQUIRED) -- a string pointing to the location where the mmseqs executable
                                  is found.
        threads (OPTIONAL) -- an integer value specifying the number of threads to run when
                              performing MMseqs2 search. Defaults to 1.
    '''
    PRESETS = [
        "sr", # short read
        "map-ont", "ava-ont", "lr:hq", # nanopore
        "asm5", "asm10", "asm20", # assembly-to-assembly
        "map-pb", "map-hifi", "ava-pb", "lr:pacbio", # pacbio
        "splice", "splice:sr", "splice:hq", # spliced alignment
        ]
    FORMATS = ["sam", "paf"]
    
    def __init__(self, query, target, preset, minimap2Exe, threads=1):
        # Set default attributes
        self.threads = threads
        self.cigar = False # set this now as preset may change it
        self.format = "sam" # set this now as preset may change it
        
        # Set input attributes
        self.query = query
        self.target = target
        self.preset = preset
        self.minimap2Exe = minimap2Exe
        
        # Set helper attributes
        self.isMinimap2 = True
    
    @property
    def query(self):
        return self._query
    
    @query.setter
    def query(self, value):
        if not os.path.isfile(value):
            raise FileNotFoundError(f"query value '{value}' is not a file")
        
        self._query = os.path.abspath(value)
    
    @property
    def target(self):
        return self._target
    
    @target.setter
    def target(self, value):
        if not os.path.isfile(value):
            raise FileNotFoundError(f"target value '{value}' is not a file")
        
        self._target = os.path.abspath(value)
    
    @property
    def minimap2Exe(self):
        return self._minimap2Exe
    
    @minimap2Exe.setter
    def minimap2Exe(self, value):
        if not os.path.isfile(value):
            raise FileNotFoundError(f"minimap2 executable not found at location '{value}'")
        
        self.mmseqsExe = os.path.abspath(value)
    
    @property
    def threads(self):
        return self._threads
    
    @threads.setter
    def threads(self, value):
        assert isinstance(value, int), "threads must be an integer"
        assert 0 < value, "threads must be a positive integer"

        self._threads = value
    
    @property
    def preset(self):
        return self._preset
    
    @preset.setter
    def preset(self, value):
        assert value in Minimap2.PRESETS, \
            f"preset must be a value in the list {Minimap2.PRESETS}"
        
        self._preset = value
        
        if value in ["asm5", "asm10", "asm20"]:
            self.format = "paf"
            self.cigar = True
        elif value in ["ava-pb", "ava-ont"]:
            self.format = "paf"
    
    @property
    def format(self):
        return self._format
    
    @format.setter
    def format(self, value):
        assert value in Minimap2.FORMATS, \
            f"format must be a value in the list {Minimap2.FORMATS}"
        
        self._format = value
    
    @property
    def cigar(self):
        return self._cigar
    
    @cigar.setter
    def cigar(self, value):
        assert isinstance(value, bool), \
            "cigar must be True or False"
        
        self._cigar = value
    
    def minimap2(self, outFile, force=False):
        '''
        Performs the minimap2 alignment operation.
        
        Parameters:
            outFile -- a string indicating the location to write the output file to.
            force -- a boolean indicating whether to overwrite the output file if it already
                     exists; default is False (do not overwrite)
        '''
        # Check if output file already exists
        if force == False and os.path.exists(outFile):
            raise FileExistsError(f"minimap2 would overwrite file '{outFile}' and force is False")
        
        # Format I/O string
        io = f"-{'a' if self.format == 'sam' else ''}{'c' if self.cigar else ''}x"
        
        # Format minimap2 command
        cmd = [
            self.mmseqsExe, io, self.preset, "-t", str(self.threads), "-o", outFile,
            self.target, self.query
        ]
        
        # Run minimap2
        if platform.system() != "Windows":
            run_minimap2 = subprocess.Popen(" ".join(cmd), shell = True,
                                            stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        else:
            run_minimap2 = subprocess.Popen(cmd, shell = True,
                                            stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        minimap2out, minimap2err = run_minimap2.communicate()
        if not "Real time:" in minimap2err.decode("utf-8"):
            raise Exception('minimap2 error text below\n' + minimap2err.decode("utf-8"))

if __name__ == "__main__":
    pass

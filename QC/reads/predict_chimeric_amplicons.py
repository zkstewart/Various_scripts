#! python 3
# predict_chimeric_transcripts.py
# A script to receive a multiple sequence alignment of reference sequences
# alongside a FASTQ file of amplicon reads which are on the same strand as
# the reference. This script will then use smith waterman alignment and
# a rough heuristic to determine whether the read may be chimeric. The
# input FASTQ will be split into chimeric and non-chimeric output FASTQ files.

import os, argparse, sys, gzip, subprocess, platform
from contextlib import contextmanager
from Bio import SeqIO
from Levenshtein import distance
from math import inf

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from Function_packages import ZS_SeqIO, ZS_Utility

class MinimalMUSCLE:
    '''
    Provides a minimal interface to the MUSCLE aligner, for the purpose of adding
    sequences into an existing alignment. Its minimisation is to try to keep its
    footprint small for speed purposes, so it won't perform any kind of validation
    since it is assumed that the input sequences are already validated via their
    direct implementation in this script.
    '''
    def __init__(self, musclePath, gapOpen=-1, gapExtend=2):
        self.exe = musclePath
        self.gapOpen = gapOpen
        self.gapExtend = gapExtend
    
    def add(self, originalMSA, toAdd):
        # Construct the cmd for subprocess
        cmd = [self.exe, "-profile", "-in1", originalMSA, "-in2", toAdd,
               "-gapopen", str(self.gapOpen), "-gapextend", str(self.gapExtend)]
        
        if platform.system() != "Windows":
            cmd = " ".join(cmd)
        
        # Run the command
        run_muscle = subprocess.Popen(cmd, shell = True,
                                stdout = subprocess.PIPE,
                                stderr = subprocess.PIPE)
        muscleout, muscleerr = run_muscle.communicate()
        
        # Check to see if there was an error
        if (not "Writing output" in muscleerr.decode("utf-8")):
            raise Exception(("ERROR: MinimalMUSCLE.add() encountered an error; have a look " +
                            f'at the stdout ({muscleout.decode("utf-8")}) and stderr ' + 
                            f'({muscleerr.decode("utf-8")}) to make sense of this.'))
        
        # Parse out the MSA result
        msa = muscleout.decode("utf-8").replace("\r", "")
        return msa

def validate_args(args):
    # Validate input data locations
    if not os.path.isfile(args.referenceFile):
        print(f'I am unable to locate the reference FASTA file ({args.referenceFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.fastqFile):
        print(f'I am unable to locate the reads FASTQ file ({args.fastqFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    
    # Validate MUSCLE location
    if not os.path.isfile(args.muscle):
        print(f'I am unable to locate the MUSCLE executable ({args.muscle})')
        print("Make sure you've specified the correct location, and try again.")
        quit()
    
    # Validate numeric arguments
    if args.gapOpen > 0:
        print("gapOpen value must be 0 or negative; positive values induce weird behaviour in MUSCLE")
        quit()
    
    # Handle file output
    args.chimerOut = os.path.join(args.outputDirectory, "chimeras.fastq")
    args.nonChimerOut = os.path.join(args.outputDirectory, "non_chimeras.fastq")
    
    if os.path.isfile(args.outputDirectory):
        print('The specified output directory already exists as a file! This program will not allowing overwriting.')
        print('You should specify a directory that does not exist, or an empty directory that does exist.')
        print("Program will exit now.")
        quit()
    elif os.path.isdir(args.outputDirectory):
        print(f"The specified output directory '{args.outputDirectory}' already exists.")
        if any([ os.path.exists(f) for f in [args.chimerOut, args.nonChimerOut] ]):
            print(f"I want to write '{args.chimerOut}' and '{args.nonChimerOut}' files, " +
                  "but one or both already exist.")
            print('This program will not allowing overwriting. Fix this issue then try again.')
            quit()
        else:
            print("I will write the output files there.")
    else:
        try:
            os.mkdir(args.outputDirectory)
            print(f"'{os.path.abspath(args.outputDirectory)}' was created as part of argument validation")
        except:
            print(f"'{os.path.abspath(args.outputDirectory)}' could not be created as part of argument validation")
            print("This may be because you've specified a directory nested inside another " +
                  "which does not already exist.")
            print("Specify a different (existing, perhaps) directory and try again.")
            quit()

@contextmanager
def open_gz_file(filename):
    if filename.endswith(".gz"):
        with gzip.open(filename, "rt") as f:
            yield f
    else:
        with open(filename) as f:
            yield f

def get_ref_haplotype_code(fastaObj, variantDict):
    '''
    Using the variantDict, determine for each reference sequence
    what the haplotype code
    
    Paramters:
        variantDict -- a dictionary with structure like:
                       {
                           pos1: ['nuc1', 'nuc2', 'nuc3'],
                           pos2: ['nuc1', 'nuc2', 'nuc3'],
                           ...
                       }
    Returns:
        refHaplotypes -- a list containing strings akin to:
                         [
                             'CGGGGAGCTT',
                             'T----AACCT',
                             'C----CAAC-'
                         ]
    '''
    refHaplotypes = []
    for i in range(len(fastaObj)):
        refHaplotypes.append("")
        for varList in variantDict.values():
            refHaplotypes[-1] += varList[i]
    return refHaplotypes

def validate_ref_is_distinguishable(refHaplotypes):
    '''
    Function to figure out if the reference sequences (in the 
    fasta object) have variants that distinguish them from each other.
    If there are identical sequences this script cannot work. Also
    this function could probably be replaced by checking the hash
    of each reference sequence for equality but, this is what I did.
    
    Raises an error an and exits program if this check fails.
    
    Parameters:
        refHaplotypes -- a list containing strings akin to:
                         [
                             'CGGGGAGCTT',
                             'T----AACCT',
                             'C----CAAC-'
                         ]
    '''
    if len(set(refHaplotypes)) < len(refHaplotypes):
        print("Some of your reference sequences are not unique, and hence indistiguishable")
        print("These are...")
        for i in range(len(refHaplotypes)):
            numOfThisRefHaplotype = refHaplotypes.count(refHaplotypes[i])
            if numOfThisRefHaplotype > 1:
                print(f"> Sequence #{i+1}")
        print("You should handle this before running this script!")
        print("Program will exit now.")
        quit()

def trim_amplicon_insertions(amplicon, originalMSA, muscleObj):
    '''
    Peforms MUSCLE profile alignment of an amplicon sequence 
    against a reference MSA. Cleans up the alignment to make it
    uniform in length against the original alignment, most notably
    by removing insertions in the amplicon relative to the full
    reference MSA.
    
    Parameters:
        amplicon -- a string of nucleotides corresponding to an amplicon read
        originalMSA -- a string indicating the file location of the original
                       reference alleles MSA file.
        muscleObj -- a MinimalMUSCLE object for aligning sequences
    Returns:
        queryAlign -- a string with gaps indicated where relevant for
                      the alignment of the amplicon against the consensus
    '''
    # Write the record to a temporary FASTA file
    tmpHash = ZS_SeqIO.Conversion.get_hash_for_input_sequences(amplicon)
    tmpFileName = ZS_Utility.tmp_file_name_gen("chimeras_muscle_" + tmpHash, "fasta")
    with open(tmpFileName, "w") as fileOut:
        fileOut.write(f">{tmpHash}\n{amplicon}\n")
    
    # Run MUSCLE to add the record to the reference MSA
    msa = muscleObj.add(originalMSA, tmpFileName)
    
    # Clean up the temporary file
    os.remove(tmpFileName)
    
    # Parse the MSA into a FASTA object
    alignedFASTA = ZS_SeqIO.FASTA(msa, isAligned=True)
    
    # Drop any positions that don't have a residue in the original MSA
    refSeqs = [ x.gap_seq for x in alignedFASTA.seqs if x.id != tmpHash ]
    recordSeq = alignedFASTA[tmpHash].gap_seq
    
    adjustedRefs = [ "" for x in refSeqs ]
    adjustedRecord = ""
    for posIndex in range(len(recordSeq)):
        refResidues = set([ x[posIndex] for x in refSeqs ])
        if refResidues != {"-"}:
            for i in range(len(refSeqs)):
                adjustedRefs[i] += refSeqs[i][posIndex]
            adjustedRecord += recordSeq[posIndex]
    
    return adjustedRecord

def align_amplicon_to_best_match(amplicon, alleleSeqs, originalMSA, muscleObj):
    '''
    Similar to align_record_to_reference(), this peforms MUSCLE profile alignment
    of an amplicon sequence, but it aligns instead against the best matching
    reference sequence. This may help to address the arbitrariness of some gap
    positions in the alignment. Aftewards, it will  clean up the alignment to make
    it uniform in length against the reference sequence.
    
    Parameters:
        amplicon -- a string of nucleotides corresponding to an amplicon read
        alleleSeqs -- a dictionary with structure like:
                      {
                          0: ['referenceAlleleSequence0', 'fileLocation0'],
                          1: ['referenceAlleleSequence1', 'fileLocation1'],
                          ...
                      }
        originalMSA -- a string indicating the file location of the original
                       reference alleles MSA file from which the alleleSeqs
                       object was created
        muscleObj -- a MinimalMUSCLE object for aligning sequences
    Returns:
        queryAlign -- a string with gaps indicated where relevant for
                      the alignment of the amplicon against the consensus
    '''
    # Write the record to a temporary FASTA file
    tmpHash = ZS_SeqIO.Conversion.get_hash_for_input_sequences(amplicon)
    tmpFileName = ZS_Utility.tmp_file_name_gen("chimeras_muscle_" + tmpHash, "fasta")
    with open(tmpFileName, "w") as fileOut:
        fileOut.write(f">{tmpHash}\n{amplicon}\n")
    
    # Find the best matching reference sequence
    bestMatch = [inf, None]
    for index, alleleSeqPair in alleleSeqs.items():
        alleleSeq, alleleFile = alleleSeqPair
        seqDistance = distance(amplicon, alleleSeq)
        if seqDistance < bestMatch[0]:
            bestMatch = [seqDistance, alleleFile]
    
    # Run MUSCLE to add the record to the best matching reference sequence MSA
    msa = muscleObj.add(bestMatch[1], tmpFileName)
    
    # Clean up the temporary file
    os.remove(tmpFileName)
    
    # Parse the MUSCLE MSA into a FASTA object
    alignedFASTA = ZS_SeqIO.FASTA(msa, isAligned=True)
    
    # Add the MUSCLE aligned sequence into the full MSA
    #fullFASTA = ZS_SeqIO.FASTA(originalMSA, isAligned=True)
    #fullFASTA.add(alignedFASTA[record.id], isAligned=True)
    
    # Drop any positions that don't have a residue in the original MSA
    refSeqs = [ x.gap_seq for x in alignedFASTA.seqs if x.id != tmpHash ]
    recordSeq = alignedFASTA[tmpHash].gap_seq
    
    adjustedRefs = [ "" for x in refSeqs ]
    adjustedRecord = ""
    for posIndex in range(len(recordSeq)):
        refResidues = set([ x[posIndex] for x in refSeqs ])
        if refResidues != {"-"}:
            for i in range(len(refSeqs)):
                adjustedRefs[i] += refSeqs[i][posIndex]
            adjustedRecord += recordSeq[posIndex]
    
    return adjustedRecord

def get_vote_data(refHaplotypes, readHaplotype, turnOffSmoothing=True):
    '''
    Performs some data smoothing of how the read's haplotype
    relates to the reference haplotypes. Ideally, the resulting
    'votes' lists allow one to determine where breakpoints occur
    between 'vote lines' which would indicate amplicon chimerism.
    
    Parameters:
        refHaplotypes -- a list containing strings akin to:
                         [
                             'CGGGGAGCTT',
                             'T----AACCT',
                             'C----CAAC-'
                         ]
        readHaplotype -- a string akin to 'C---GCAACT'
        turnOffSmoothing -- OPTIONAL; a boolean indicating whether
                        vote smoothing over a sliding window should
                        be turned off. It is believed that smoothing
                        will make the algorithm more tolerant to
                        rare substitution errors. Default is True,
                        i.e., smoothing will not be applied.
    Returns:
        voteLines -- a list of strings akin to:
                     [
                         '--X-XX-X-X',
                         'XX---XXXXX',
                         '-X----X-X-'
                     ]
        voteSmoothed -- a list of lists akin to:
                        [
                            [1, 0, 0, 0, 0, 0, 0, 0, 0, 1],
                            [0, 1, 1, 1, 0, 0, 0, 1, 1, 1],
                            [1, 1, 1, 1, 1, 1, 1, 1, 1, 0]
                        ]
    '''
    WINDOW_SIZE = 1
    
    # Create "vote line" list for finding where variants are shared
    voteLines = [ '' for i in range(len(refHaplotypes)) ]
    for i in range(len(readHaplotype)):
        nuc = readHaplotype[i]
        for x in range(len(refHaplotypes)):
            if refHaplotypes[x][i].upper() == nuc.upper():
                voteLines[x] += "X"
            else:
                voteLines[x] += '-'
    
    # If we're not smoothing, return a simple vote line
    if turnOffSmoothing:
        return voteLines, [ [1 if x == "X" else 0 for x in line] for line in voteLines ]
    
    # Otherwise, smooth over a moving window for optimal alignment site
    voteSmoothed = []
    
    for voteLine in voteLines:
        voteSmoothed.append([])
        # Harshly weight first SNP
        if voteLine[0] == "X":
            voteSmoothed[-1].append(1)
        else:
            voteSmoothed[-1].append(0)
        
        # Calculate for body SNPs
        for i in range(WINDOW_SIZE, len(voteLine) - WINDOW_SIZE):
            if voteLine[i] == "X": # don't smooth if this position is the same as a reference type
                voteSmoothed[-1].append(1)
            else:
                window = voteLine[i-WINDOW_SIZE:i+WINDOW_SIZE+1]
                
                if window.count("X") / len(window) > 0.5:
                    voteSmoothed[-1].append(1)
                else:
                    voteSmoothed[-1].append(0)
        
        # Harshly weight last SNP
        if voteLine[-1] == "X":
            voteSmoothed[-1].append(1)
        else:
            voteSmoothed[-1].append(0)
    
    return voteLines, voteSmoothed

def find_last_index(inputList, value):
    return len(inputList) - inputList[::-1].index(value) - 1

def simple_chimera_detection(voteSmoothed):
    '''
    Employs a relatively simple heuristic upon the smoothed vote list to
    see if an amplicon contains variants from a mixture of different reference alleles.
    
    Parameters:
        voteSmoothed -- a list of lists akin to:
                        [
                            [1, 0, 0, 0, 0, 0, 0, 0, 0, 1],
                            [0, 1, 1, 1, 0, 0, 0, 1, 1, 1],
                            [1, 1, 1, 1, 1, 1, 1, 1, 1, 0]
                        ]
    Returns:
        isMaybeChimera -- a boolean of True if the amplicon might be chimeric, else False
    '''
    # Locate the best vote line
    voteNums = [ sum(voteSmoothed[i]) for i in range(len(voteSmoothed)) ]
    bestVote = voteNums.index(max(voteNums)) # it's okay if this has ties
    
    # If our best vote has no evidence of switching, return it now
    if not 0 in voteSmoothed[bestVote]:
        return False
    
    # Otherwise, see if reference switching may have occurred at a 0 position
    for x in range(len(voteSmoothed[bestVote])):
        # At a 0 position ...
        if voteSmoothed[bestVote][x] == 0:
            # ... see if the other vote lines have a 1
            if any([ voteSmoothed[i][x] == 1 for i in range(len(voteSmoothed)) if i != bestVote ]):
                return True
    
    # If not, just return False
    return False

def complex_chimera_detection(voteSmoothed, disallowAmbiguity=False):
    '''
    Employs a more complex heuristic upon the smoothed vote list to
    see if there may be a reference switch occurring in the amplicon sequence
    at its left and/or right sides.
    
    This heuristic may fail to discover complex chimerism (e.g., where the middle
    section is chimeric, but the left and right are the same reference) but, I'm
    operating under an assumption that that is exceedingly rare.
    
    Parameters:
        voteSmoothed -- a list of lists akin to:
                        [
                            [1, 0, 0, 0, 0, 0, 0, 0, 0, 1],
                            [0, 1, 1, 1, 0, 0, 0, 1, 1, 1],
                            [1, 1, 1, 1, 1, 1, 1, 1, 1, 0]
                        ]
        disallowAmbiguity -- OPTIONAL; a boolean indicating whether to, in situations
                             where chimerism is ambiguous and a read has no clear breakpoint,
                             to consider it as a chimera. Default is False, i.e., ambiguous
                             amplicons will not be considered as chimeras.
    Returns:
        isMaybeChimera -- a boolean of True if the amplicon might be chimeric, else False
    '''
    # Locate the best vote line
    voteNums = [ sum(voteSmoothed[i]) for i in range(len(voteSmoothed)) ]
    bestVote = voteNums.index(max(voteNums)) # it's okay if this has ties
    
    # Get the best smoothed vote line
    bestVoteSmoothed = voteSmoothed[bestVote]
    otherVoteSmoothed = [ voteSmoothed[i] for i in range(len(voteSmoothed)) if i != bestVote ]
    
    # If our best vote has no evidence of switching, return it now
    if not 0 in voteSmoothed[bestVote]:
        return False
    
    # Otherwise, see if reference switching at an internal cut line improves the best vote line
    for cutIndex in range(1, len(bestVoteSmoothed) - 1):
        bestLeft, bestRight = bestVoteSmoothed[0: cutIndex], bestVoteSmoothed[cutIndex: ]
        for otherVote in otherVoteSmoothed:
            otherLeft, otherRight = otherVote[0: cutIndex], otherVote[cutIndex: ]
            
            if disallowAmbiguity:
                # If the other left is possibly better than the best left
                if (0 in bestLeft and bestLeft != otherLeft and sum(otherLeft) >= sum(bestLeft)):
                    # AND the other right is possibly worse than the best right
                    "This fulfills our assumption of chimerism"
                    if (0 in otherRight and bestRight != otherRight and sum(otherRight) <= sum(bestRight)):
                        return True
                
                # Otherwise, if the other right is possibly better than the best right
                elif (0 in bestRight and bestRight != otherRight and sum(otherRight) >= sum(bestRight)):
                    # AND the other left is possibly worse than the best left
                    "This fulfills our assumption of chimerism"
                    if (0 in otherLeft and bestLeft != otherLeft and sum(otherLeft) <= sum(bestLeft)):
                        return True
            else:
                if sum(otherLeft) > sum(bestLeft) or sum(otherRight) > sum(bestRight):
                    return True
    
    # If not, just return False
    return False

def get_read_haplotype(sequence, variantDict):
    '''
    Parameters:
        sequence -- a string of nucleotides
        variantDict -- a dictionary with structure like:
                       {
                           pos1: ['nuc1', 'nuc2', 'nuc3'],
                           pos2: ['nuc1', 'nuc2', 'nuc3'],
                           ...
                       }
    Returns:
        readHaplotype -- a string of nucleotides akin to 'C---GCAACT'
                         indicating the haplotype of the read at each variant
                         position
    '''
    return "".join([ sequence[x] for x in variantDict.keys() ])

def main():
    usage = """%(prog)s will predict chimeric amplicons by using MUSCLE to add them into
    a reference alleles multiple sequence alignment. Variant sites in the reference alleles
    are identified in the amplicon, and this information is used to determine if the amplicon
    is chimeric or not. Reads will be output split into two FASTQ files for those reads deemed
    to be chimeric and non-chimeric.
    
    The mode to run this program in is specified by the --mode argument. The 'simple' mode
    will utilise a coarse heuristic for chimera detection that is likely to have the best recall
    but lowest precision. The 'complex' mode will utilise a more advanced heuristic that is likely
    to have the lowest recall but greatest precision. The 'complex-balanced' mode will employ the
    complex heuristic, but will also consider reads with ambiguous chimerism to be chimeras.
    This mode is expected to have a balanced recall and precision.
    
    Some notes include:
    
    1) The reference FASTA must have more than one sequence, and be 
    a multiple sequence alignment. Try to keep it to just the amplicon
    region.
    2) The input FASTQ(.gz) must have all reads in the same orientation
    as the reference sequence. You can ensure that by first running
    complement_amplicons.py.
    3) MUSCLE parameters are user specifiable, but the defaults appear to be necessary
    for this script to work well.
    4) --turnOnSmoothing enables some heuristics that make the script more tolerant
    to rare but genuine sequencing errors. Turning this on will make the script more
    relaxed (lower recall, more precision).
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-fa", dest="referenceFile",
                   required=True,
                   help="Specify the location of the reference FASTA MSA file")
    p.add_argument("-fq", dest="fastqFile",
                   required=True,
                   help="Specify the location of the reads FASTQ(.gz) file")
    p.add_argument("-o", dest="outputDirectory",
                   required=True,
                   help="Output directory where chimer/non-chimer FASTQ files will be written")
    p.add_argument("-m", dest="muscle",
                   required=True,
                   help="Specify the location of the MUSCLE executable")
    p.add_argument("--mode", dest="mode",
                   required=True,
                   choices=["simple", "complex", "complex-balanced"],
                   help="Indicate which mode you want to run this program in")
    # Opts (behavioural)
    p.add_argument("--turnOnSmoothing", dest="turnOnSmoothing",
                   required=False,
                   action="store_true",
                   help="""Optionally, specify this flag to enable smoothing
                   of the vote data. This will make all modes less tolerant to rare
                   substitution errors.""",
                   default=False)
    # Opts (MUSCLE)
    p.add_argument("--gapOpen", dest="gapOpen",
                   required=False,
                   type=int,
                   help="""Optionally, specify the gap opening penalty MUSCLE
                   should use when aligning an amplicon against the reference alleles
                   MSA; default==-1, which appears to work best for this script.""",
                   default=-1)
    p.add_argument("--gapExtend", dest="gapExtend",
                   required=False,
                   type=int,
                   help="""Optionally, specify the gap extension penalty MUSCLE
                   should use when aligning an amplicon against the reference alleles
                   MSA; default==2, which appears to work best for this script.""",
                   default=2)
    
    args = p.parse_args()
    validate_args(args)
    
    # Load reference fasta as MSA
    fastaObj = ZS_SeqIO.FASTA(args.referenceFile, isAligned=True)
    
    # Locate variants in the reference
    variantDict = {}
    seqs = [ x.gap_seq for x in fastaObj.seqs ]
    
    for index, positions in enumerate(zip(*seqs)):
        # Skip if there are no variants in this position
        if len(set(positions)) == 1:
            continue
        
        # Record the position of the variant and what's found in each sequence
        else:
            variantDict[index] = list(positions)
    
    # Validate that this script can work
    refHaplotypes = get_ref_haplotype_code(fastaObj, variantDict)
    validate_ref_is_distinguishable(refHaplotypes)
    
    # Establish quick access to individual reference allele sequences and files
    alleleSeqs = {}
    for index, FastASeq_obj in enumerate(fastaObj.seqs):
        alleleSeqs[index] = [ FastASeq_obj.gap_seq.upper(), os.path.join(args.outputDirectory, f"ref_{index}.fasta") ]
        with open(alleleSeqs[index][1], "w") as fileOut:
            fileOut.write(f">{FastASeq_obj.id}\n{FastASeq_obj.gap_seq}\n")
    
    # Iterate over FASTQ reads to determine if they are chimeric
    numChimeras = 0
    numReads = 0
    cachedResults = {}
    muscleObj = MinimalMUSCLE(args.muscle, args.gapOpen, args.gapExtend)
    with open_gz_file(args.fastqFile) as fileIn, open(args.chimerOut, "w") as chimeraOut, open(args.nonChimerOut, "w") as nonChimerOut:
        records = SeqIO.parse(fileIn, "fastq")
        for record in records:
            numReads += 1
            
            # If this amplicon is identical a previously processed one, use the cached result
            if str(record.seq) in cachedResults:
                isMaybeChimera = cachedResults[str(record.seq)]
            
            # Otherwise:
            else:
                # Align record against reference
                queryAlign = align_amplicon_to_best_match(str(record.seq), alleleSeqs,
                                                          args.referenceFile, muscleObj)
                
                # Get the variant at all relevant positions
                readHaplotype = get_read_haplotype(queryAlign, variantDict)
                
                # Obtain 'vote' data structures
                voteLines, voteSmoothed = get_vote_data(refHaplotypes, readHaplotype,
                                                        turnOffSmoothing = not args.turnOnSmoothing)
                
                # Use a heuristic to find if the read may be chimeric
                if args.mode == "simple":
                    isMaybeChimera = simple_chimera_detection(voteSmoothed)
                elif args.mode == "complex":
                    isMaybeChimera = complex_chimera_detection(voteSmoothed,
                                                               disallowAmbiguity = False)
                elif args.mode == "complex-balanced":
                    isMaybeChimera = complex_chimera_detection(voteSmoothed,
                                                               disallowAmbiguity = True)
                
                # Cache the result
                cachedResults[str(record.seq)] = isMaybeChimera
            
            # If we have a chimera, count its existence and write to file
            if isMaybeChimera:
                chimeraOut.write(record.format("fastq"))
                numChimeras += 1
            # If we don't have a chimera, write to its respective file
            else:
                nonChimerOut.write(record.format("fastq"))
    
    # Clean up the reference allele files
    for alleleFile in [ alleleSeqPair[1] for alleleSeqPair in alleleSeqs.values() ]:
        if os.path.isfile(alleleFile):
            os.remove(alleleFile)
    
    # Print chimera statistics then end program
    print("# predict_chimeric_amplicons output statistics")
    print(f"# Processing of file: '{args.fastqFile}'")
    print(f"# > Of {numReads} amplicons, {numChimeras} were detected as potential chimeras")
    print(f"# > This equates to {(numChimeras / numReads) * 100}% of amplicons")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

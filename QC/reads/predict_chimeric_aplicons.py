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

def most_common_position(positions):
    '''
    Checks a tuple of nucleotides to find the most common one which is not
    a gap ('-')
    
    Parameters:
        positions -- a tuple with more than one nucleotide as a single character
    Returns:
        nucleotide -- the most common nucleotide at this position that is not a gap.
    '''
    # Get a set of nucleotides sans gap characters
    positionsSet = set(positions)
    try:
        positionsSet.remove("-")
    except:
        pass
    
    # If there is only one nucleotide, return it here
    if len(positionsSet) == 1:
        return positionsSet.pop()
    
    # Otherwise, find the most frequent nucleotide
    else:
        'Ties result in the first found being the representative nucleotide'
        mostFrequent = [0, None]
        for position in positionsSet:
            numOfThisPosition = positions.count(position)
            if numOfThisPosition > mostFrequent[0]:
                mostFrequent = [numOfThisPosition, position]
        return mostFrequent[1]

def align_record_to_reference(record, originalMSA, muscleObj):
    '''
    Peforms MAFFT L-insi --add of the record (which should be an amplicon
    sequence) against a reference MSA. Cleans up the alignment to make it
    uniform in length against the original alignment.
    
    Parameters:
        record -- a SeqIO.SeqRecord object with a nucleotide amplicon read
        originalMSA -- a string indicating the file location of the original
                       reference alleles MSA file.
        muscleObj -- a MinimalMUSCLE object for aligning sequences
    Returns:
        queryAlign -- a string with gaps indicated where relevant for
                      the alignment of the amplicon against the consensus
    '''
    # Write the record to a temporary FASTA file
    tmpFileName = ZS_Utility.tmp_file_name_gen("chimeras_muscle", "fasta")
    with open(tmpFileName, "w") as fileOut:
        SeqIO.write(record, fileOut, "fasta")
    
    # Run MUSCLE to add the record to the reference MSA
    msa = muscleObj.add(originalMSA, tmpFileName)
    
    # Clean up the temporary file
    os.remove(tmpFileName)
    
    # Parse the MSA into a FASTA object
    alignedFASTA = ZS_SeqIO.FASTA(msa, isAligned=True)
    
    # Drop any positions that don't have a residue in the original MSA
    refSeqs = [ x.gap_seq for x in alignedFASTA.seqs if x.id != record.id ]
    recordSeq = alignedFASTA[record.id].gap_seq
    
    adjustedRefs = [ "" for x in refSeqs ]
    adjustedRecord = ""
    for posIndex in range(len(recordSeq)):
        refResidues = set([ x[posIndex] for x in refSeqs ])
        if refResidues != {"-"}:
            for i in range(len(refSeqs)):
                adjustedRefs[i] += refSeqs[i][posIndex]
            adjustedRecord += recordSeq[posIndex]
    
    return adjustedRecord

def get_vote_data(refHaplotypes, readHaplotype, turnOffSmoothing=False):
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
                        rare substitution errors. Default is False,
                        i.e., smoothing will be applied.
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
            if refHaplotypes[x][i] == nuc:
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

def is_possible_chimera(voteSmoothed, voteLines, disallowAmbiguity=False):
    '''
    Employs a relatively simple heuristic upon the smoothed vote list to
    see if there may be a reference switch occurring in the amplicon sequence
    at its left and/or right sides.
    
    This heuristic would fail to discover complex chimerism (e.g., where the middle
    section is chimeric, but the left and right are the same reference) but, I'm
    operating under an assumption that that is exceedingly rare.
    
    Parameters:
        voteSmoothed -- a list of lists akin to:
                        [
                            [1, 0, 0, 0, 0, 0, 0, 0, 0, 1],
                            [0, 1, 1, 1, 0, 0, 0, 1, 1, 1],
                            [1, 1, 1, 1, 1, 1, 1, 1, 1, 0]
                        ]
        voteLines -- a list of strings akin to:
                     [
                         '--X-XX-X-X',
                         'XX---XXXXX',
                         '-X----X-X-'
                     ]
        disallowAmbiguity -- OPTIONAL; a boolean indicating whether to, in situations
                             where chimerism is ambiguous and a read has no clear breakpoint,
                             to consider it as a chimera. Default is False, i.e., ambiguous
                             amplicons will not be considered as chimeras.
    Returns:
        isMaybeChimera -- a boolean of True if the amplicon might be chimeric, else False
    '''
    voteNums = [ sum(x) for x in voteSmoothed ]
    
    # If there are no ties, take the best vote at face value
    if voteNums.count(max(voteNums)) == 1:
        bestVote = voteNums.index(max(voteNums))
    
    # Tie break equal voteNums by voteLines
    else:
        voteNums = [ voteNums[i] + voteLines[i].count("X") for i in range(len(voteNums)) ]
        
        # If this breaks the tie, use it
        if voteNums.count(max(voteNums)) == 1:
            bestVote = voteNums.index(max(voteNums))
        
        # Tie break by left and right "X" status
        else:
            for x in range(0, len(voteLines[0]) - 1):
                voteNums = [
                    voteNums[i] +
                        (1 if voteLines[i][x] == "X" else 0) +
                        (1 if voteLines[i][::-1][x] == "X" else 0)
                    for i in range(len(voteNums)) 
                ]
                
                if voteNums.count(max(voteNums)) == 1:
                    break
            
            # Take the best vote at face value now - ties be damned if the above didn't solve it
            bestVote = voteNums.index(max(voteNums))
    
    # Get the best smoothed vote line
    bestVoteSmoothed = voteSmoothed[bestVote]
    otherVoteSmoothed = [ voteSmoothed[i] for i in range(len(voteSmoothed)) if i != bestVote ]
    
    # If our best vote has no evidence of switching, return it now
    if not 0 in bestVoteSmoothed:
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

def main():
    usage = """%(prog)s will predict chimeric amplicons by using MUSCLE to add them into
    a reference alleles multiple sequence alignment. Variant sites in the reference alleles
    are identified in the amplicon, and this information is used to determine if the amplicon
    is chimeric or not. Reads will be output split into two FASTQ files for those reads deemed
    to be chimeric and non-chimeric. Some notes include:
    
    1) The reference FASTA must have more than one sequence, and be 
    a multiple sequence alignment. Try to keep it to just the amplicon
    region.
    2) The input FASTQ(.gz) must have all reads in the same orientation
    as the reference sequence. You can ensure that by first running
    complement_amplicons.py.
    3) MUSCLE parameters are user specifiable, but the defaults appear to be necessary
    for this script to work well.
    4) --turnOffSmoothing disables some heuristics that make the script more tolerant
    to rare but genuine sequencing errors. Turning this off will make the script more
    strict.
    5) --noAmbiguity will consider any read whose chimerism is ambiguous to be a chimera.
    In this context, ambiguous chimerism occurs when there is a breakpoint in the read where
    the read does not match its best-matching reference allele properly, and it matches another
    allele equally well. Hence, the read could optimally be a chimer or a non-chimer.
    This option will make the script more strict by considering these reads as chimeras.
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
    # Opts (behavioural)
    p.add_argument("--turnOffSmoothing", dest="turnOffSmoothing",
                   required=False,
                   action="store_true",
                   help="""Optionally, specify this flag to disable the smoothing
                   of the vote data. This will produce a stricter output which is
                   less tolerant to rare substitution errors.""",
                   default=False)
    p.add_argument("--noAmbiguity", dest="noAmbiguity",
                   required=False,
                   action="store_true",
                   help="""Optionally, specify this flag if you'd like any reads
                   whose chimerism is ambiguous (e.g., there is no clear breakpoint,
                   but the read could optimally be a chimer or a non-chimer) to be
                   considered as chimeras. This will produce a stricter output which
                   is less tolerant to rare indel errors.""",
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
    
    consensusReference = ""
    for index, positions in enumerate(zip(*seqs)):
        # Generate a consensus reference sans any gaps
        nuc = most_common_position(positions)
        consensusReference += nuc
        
        # Skip if there are no variants in this position
        if len(set(positions)) == 1:
            continue
        
        # Record the position of the variant and what's found in each sequence
        else:
            variantDict[index] = list(positions)
    
    # Validate that this script can work
    refHaplotypes = get_ref_haplotype_code(fastaObj, variantDict)
    validate_ref_is_distinguishable(refHaplotypes)
    
    # Iterate over FASTQ reads to determine if they are chimeric
    numChimeras = 0
    numReads = 0
    muscleObj = MinimalMUSCLE(args.muscle, args.gapOpen, args.gapExtend)
    with open_gz_file(args.fastqFile) as fileIn, open(args.chimerOut, "w") as chimeraOut, open(args.nonChimerOut, "w") as nonChimerOut:
        records = SeqIO.parse(fileIn, "fastq")
        for record in records:
            numReads += 1
            
            # Align record against reference
            queryAlign = align_record_to_reference(record, args.referenceFile,
                                                   muscleObj)
            
            # Get the variant at all relevant positions
            readHaplotype = "".join([ queryAlign[x] for x in variantDict.keys() ])
            
            # Obtain 'vote' data structures
            voteLines, voteSmoothed = get_vote_data(refHaplotypes, readHaplotype,
                                                    turnOffSmoothing = args.turnOffSmoothing)
            
            # Use a heuristic to find if the read may be chimeric
            ## Find out which reference the amplicon primarily belongs to
            isMaybeChimera = is_possible_chimera(voteSmoothed, voteLines,
                                                 disallowAmbiguity = args.noAmbiguity)
            
            # If we have a chimera, count its existence and write to file
            if isMaybeChimera:
                chimeraOut.write(record.format("fastq"))
                numChimeras += 1
            # If we don't have a chimera, write to its respective file
            else:
                nonChimerOut.write(record.format("fastq"))
    
    # Print chimera statistics then end program
    print("# predict_chimeric_amplicons output statistics")
    print(f"# > Of {numReads} amplicons, {numChimeras} were detected as potential chimeras")
    print(f"# > This equates to {(numChimeras / numReads) * 100}% of amplicons")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

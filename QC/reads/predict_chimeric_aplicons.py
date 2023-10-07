#! python 3
# predict_chimeric_transcripts.py
# A script to receive a multiple sequence alignment of reference sequences
# alongside a FASTQ file of amplicon reads which are on the same strand as
# the reference. This script will then use smith waterman alignment and
# a rough heuristic to determine whether the read may be chimeric. The
# input FASTQ will be split into chimeric and non-chimeric output FASTQ files.

import os, argparse, sys, gzip
from contextlib import contextmanager
from Bio import SeqIO

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from Function_packages import ZS_SeqIO, ZS_AlignIO

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

def align_record_to_consensus(record, consensusReference):
    '''
    Peforms striped smith waterman alignment of the record
    (which should be an amplicon sequence) against a consensus
    sequence. Cleans up the alignment to make it uniform in length
    against the original consensus.
    
    Parameters:
        record -- a SeqIO.SeqRecord object with a nucleotide amplicon read
        consensusReference -- a string of a nucleotide generated via
                              consensus generation from a MSA.
    Returns:
        queryAlign -- a string with gaps indicated where relevant for
                      the alignment of the amplicon against the consensus
        consensusAlign -- a string of how the consensus was aligned against
                          the amplicon sequence.
    '''
    alignResult = ZS_AlignIO.SSW.ssw_parasail(
        str(record.seq), consensusReference, "nucleotide"
    ) # TBD: Alignment needs to be against a static reference
    
    # Adjust the query alignment to have the same positioning as the consensus
    queryAlign, consensusAlign = "", ""
    for qNuc, tNuc in zip(alignResult.queryAlign, alignResult.targetAlign):
        if tNuc != "-":
            queryAlign += qNuc
            consensusAlign += tNuc
    
    # Add any gaps to the start and end for uniformity of length
    queryAlign = '-'*alignResult.targetStartIndex + queryAlign
    queryAlign += '-'*(len(consensusReference) - len(queryAlign))
    
    return queryAlign, consensusAlign

def get_vote_data(refHaplotypes, readHaplotype):
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
    # Create "vote line" list for finding where variants are shared
    voteLines = [ '' for i in range(len(refHaplotypes)) ]
    for i in range(len(readHaplotype)):
        nuc = readHaplotype[i]
        for x in range(len(refHaplotypes)):
            if refHaplotypes[x][i] == nuc:
                voteLines[x] += "X"
            else:
                voteLines[x] += '-'
    
    # Calculate a moving window for optimal alignment site
    'This smoothes the vote line'
    voteSmoothed = []
    windowSize = 1
    for voteLine in voteLines:
        voteSmoothed.append([])
        # Harshly weight first SNP
        if voteLine[0] == "X":
            voteSmoothed[-1].append(1)
        else:
            voteSmoothed[-1].append(0)
        
        # Calculate for body SNPs
        for i in range(windowSize, len(voteLine) - windowSize):
            window = voteLine[i-windowSize:i+windowSize+1]
            if window.count("X") / len(window) > 0.5:
                voteSmoothed[-1].append(1)
            else:
                voteSmoothed[-1].append(0)
            # voteSmoothed[-1].append(window.count("X") / len(window))
        
        # Harshly weight last SNP
        if voteLine[-1] == "X":
            voteSmoothed[-1].append(1)
        else:
            voteSmoothed[-1].append(0)
    
    return voteLines, voteSmoothed

def find_last_index(inputList, value):
    return len(inputList) - inputList[::-1].index(value) - 1

def is_possible_chimera(voteSmoothed, voteLines):
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
    
    # Get the best smooted vote line
    bestVoteSmoothed = voteSmoothed[bestVote]
    otherVoteSmoothed = [ voteSmoothed[i] for i in range(len(voteSmoothed)) if i != bestVote ]
    
    # If our best vote has no evidence of switching, return it now
    if not 0 in bestVoteSmoothed:
        return False
    
    # Otherwise, see if reference switching improves the best vote line
    bestVoteLeft = bestVoteSmoothed[0: bestVoteSmoothed.index(0)]
    leftMaybeChimera = any([
        len(otherVote[0: otherVote.index(0)]) > len(bestVoteLeft)
        for otherVote in otherVoteSmoothed
    ])
    
    bestVoteRight = bestVoteSmoothed[find_last_index(bestVoteSmoothed, 0)+1: ]
    rightMaybeChimera = any([
        len(otherVote[find_last_index(otherVote, 0)+1: ]) > len(bestVoteRight)
        for otherVote in otherVoteSmoothed
    ])
    
    return leftMaybeChimera or rightMaybeChimera

def change_point_detection(voteSmoothed):
    '''
    Unused - was initially experimented upon but it is not a good
    fit for this analysis.
    '''
    import ruptures as rpt
    import numpy as np

    signal = np.array(voteSmoothed).T
    algo = rpt.Binseg(model="rbf").fit(signal)
    result = algo.predict(pen=2)
    # rpt.display(signal, result)
    isChimera = any([ r for r in result if r < len(signal) ])
    return isChimera

def main():
    usage = """%(prog)s will predict chimeric amplicons by ...
    
    1) The reference FASTA must have more than one sequence, and be 
    a multiple sequence alignment. Try to keep it to just the amplicon
    region.
    2) The input FASTQ must have all reads in the same orientation
    as the reference sequence. You can ensure that by first running
    unify_complement_seqs.py.
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-fa", dest="referenceFile",
                   required=True,
                   help="Specify the location of the reference FASTA file")
    p.add_argument("-fq", dest="fastqFile",
                   required=True,
                   help="Specify the location of the reads FASTQ file")
    p.add_argument("-o", dest="outputDirectory",
                   required=True,
                   help="Output directory where QC files will be written")
    
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
    with open_gz_file(args.fastqFile) as fileIn, open(args.chimerOut, "w") as chimeraOut, open(args.nonChimerOut, "w") as nonChimerOut:
        records = SeqIO.parse(fileIn, "fastq")
        for record in records:
            numReads += 1
            # Align record against consensus
            queryAlign, consensusAlign = align_record_to_consensus(record, consensusReference)
            
            # Get the variant at all relevant positions
            readHaplotype = "".join([ queryAlign[x] for x in variantDict.keys() ])
            
            # Obtain 'vote' data structures
            voteLines, voteSmoothed = get_vote_data(refHaplotypes, readHaplotype)
            
            # Use a heuristic to find if the read may be chimeric
            ## Find out which reference the amplicon primarily belongs to
            isMaybeChimera = is_possible_chimera(voteSmoothed, voteLines)
            
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

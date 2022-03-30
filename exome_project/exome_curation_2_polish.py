#! python3
# exome_curation_autoPolish.py
# Follows up on exome_curation_prep.py to perform some additional
# polishing of the MSAs including trimming and removal of indel errors.

import sys, argparse, os
sys.path.append(os.path.dirname(os.path.dirname(__file__))) # 2 dirs up is where we find dependencies
from Function_packages import ZS_SeqIO
from exome_liftover import ssw_parasail

def validate_args(args):
    # Validate input data location
    if not os.path.isdir(args.alignmentsDir):
        print('I am unable to locate the directory where the alignments files are (' + args.alignmentsDir + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Handle file output
    if os.path.isdir(args.outputDir):
        if os.listdir(args.outputDir) != []:
            print(args.outputDir + ' already contains files. Specify a new location or move any existing files elsewhere.')
            quit()
    else:
        try:
            os.mkdir(args.outputDir)
            print("Created '{0}' directory as part of argument validation".format(args.outputDir))
        except:
            print("Wasn't able to create '{0}' directory; does '{1}' actually exist?".format(args.outputDir, os.path.dirname(args.outputDir)))

def solve_translation_frames(FASTA_obj):
    '''
    This function will make a best attempt guess at the starting frame for
    sequences within the FASTA object.
    
    NOTE: This function assumes the strand for the translation is ALWAYS
    the positive strand. If that's not the case, you're gonna have problems.
    Since we're handling MSAs, they should already all be in the same strand,
    so maybe you ought to reverse complement them all before running this
    function.
    
    Params:
        FASTA_obj -- a ZS_SeqIO.FASTA object
    Returns:
        solutionDict -- a dictionary with [index] = [seq, frame, hasStopCodon (bool)] structure
                        where index maps to the FASTA_obj.seqs list.
    '''
    resultsDict = {} # contains [index] = [[seq, frame, hasStopCodon], ... +2 more frames]
    
    # Obtain translations in the three strand=1 frames
    for i in range(len(FASTA_obj)):
        FastASeq_obj = FASTA_obj[i]
        
        results = [] # contains triples of [seq, frame, hasStopCodon (bool)]
        for frame in range(0, 3):
            seq, strand, frame = FastASeq_obj.get_translation(strand=1, frame=frame)
            results.append([seq, frame, "*" in seq[:-1]]) # don't count the last position since that can be normal
        resultsDict[i] = results
    
    # Loop back through and find easy-to-solve and hard-to-solve scenarios
    solutionDict = {} # contains [index] = [seq, frame, hasStopCodon]
    problemDict = {} # contains [index] = [[seq, frame, hasStopCodon], ... +2 more frames]
    for i in range(len(FASTA_obj)):
        results = resultsDict[i]
        
        # Handle easy-to-solve scenarios
        numWithoutStopCodons = sum([1 for r in results if r[2] == False])
        if numWithoutStopCodons == 1:
            solutionDict[i] = [r for r in results if r[2] == False][0]
            continue
        
        # Note hard-to-solve scenarios
        else:
            problemDict[i] = resultsDict[i]
    
    # Find the best solution to problem sequences
    for i, results in problemDict.items():
        scores = [[], [], []] # will store scores for frame 0, 1, and 2
        for j in range(0, 3):
            problemSeq, _, _ = results[j]
            for x, solution in solutionDict.items():
                solutionSeq, _, _ = solution
                _, _, _, score = ssw_parasail(problemSeq, solutionSeq)
                scores[j].append(score)
        
        # Calculate two indicative metrics
        maxIndividualScores = [[max(scores[frame]), frame] for frame in range(0, 3)] # Gives [[score, frame], ... +2 frames]
        totalScores = [[sum(scores[frame]), frame] for frame in range(0, 3)] # Gives [[score, frame], ... +2 frames]
        
        maxIndividualScoresFrame = max(maxIndividualScores, key = lambda x: x[0])[1]
        maxTotalScoresFrame = max(totalScores, key = lambda x: x[0])[1]
        
        # Find the best frame from two indicative metrics
        if maxIndividualScoresFrame == maxTotalScoresFrame:
            bestFrame = maxIndividualScoresFrame
        # Find the best frame by picking the most informative of the two metrics
        else:
            individualScoreDifference = max(maxIndividualScores, key = lambda x: x[0])[0] / min(maxIndividualScores, key = lambda x: x[0])[0]
            totalScoreDifference = max(totalScores, key = lambda x: x[0])[0] / min(totalScores, key = lambda x: x[0])[0]

            individualScoreDifference = round(individualScoreDifference, 2)
            totalScoreDifference = round(totalScoreDifference, 2)
            
            if totalScoreDifference > individualScoreDifference: # we're checking to see which metric separates the data most
                bestFrame = maxTotalScoresFrame
            elif individualScoreDifference > totalScoreDifference:
                bestFrame = maxIndividualScoresFrame
            # Just pick a best-guess frame and alert the user to this occurrence
            else:
                bestFrame = maxTotalScoresFrame # assume total score sum should be most indicative
                print("Failed to find a good frame for FASTA based on {0}".format(FASTA_obj.fileOrder[0][0]))
        
        solutionDict[i] = results[bestFrame]
    
    return solutionDict

def add_codon_seqs(FASTA_obj, solutionDict):
    '''
    This function will modify the FASTA object such that, for every sequence,
    a new partner sequence is added with sequence like:
    
        --12-3123--1---2-3
    
    Indicating codon positions 1, 2, and 3.
    
    Params:
        FASTA_obj -- a ZS_SeqIO.FASTA object
        solutionDict -- a dictionary with [index] = [seq, frame, hasStopCodon (bool)] structure
                        where index maps to the FASTA_obj.seqs list. Should have been computed
                        for this FASTA_obj by solve_translation_frames().
    '''
    assert len(FASTA_obj) == len(solutionDict), "solutionDict wasn't built from the provided FASTA_obj!"
    
    dummySeqs = []
    for i in range(len(FASTA_obj)):
        FastASeq_obj = FASTA_obj[i]
        solution = solutionDict[i]
        
        numbersSeq = ""
        currentCodon = solution[1] # this is the starting frame
        for codonIndex in range(0, len(FastASeq_obj.seq), 3):
            codon = FastASeq_obj.seq[codonIndex:codonIndex+3].upper()
            
            if codon == "NNN":
                numbersSeq += "---" # no codon number can be inferred for gap positions
            else:
                numbersSeq += "{0}{1}{2}".format(
                    ((currentCodon + 0) % 3) + 1, # +1 for 1-based numbering
                    ((currentCodon + 1) % 3) + 1,
                    ((currentCodon + 2) % 3) + 1,
                )
                currentCodon = (currentCodon + 3) % 3
        
        # Add real "gaps" back into numbersSeq
        gappedNumbersSeq = ""
        ongoingCount = 0
        for letter in FastASeq_obj.gap_seq:
            if letter == "-":
                gappedNumbersSeq += "-"
            else:
                gappedNumbersSeq += numbersSeq[ongoingCount]
                ongoingCount += 1
        
        # Add dummy sequences to FASTA object
        description = "codons_" + FastASeq_obj.description
        if FastASeq_obj.alt != None:
            alt = "codons_" + FastASeq_obj.alt
            dummyFastASeq_obj = ZS_SeqIO.FastASeq(description, alt=alt, gapSeq=gappedNumbersSeq)
        else:
            dummyFastASeq_obj = ZS_SeqIO.FastASeq(description, gapSeq=gappedNumbersSeq)
        
        dummySeqs.append(dummyFastASeq_obj)
    
    # Insert into list at appropriate positions
    ongoingCount = 0 # Tracks how many other things have already been inserted for our offset
    for i in range(len(FASTA_obj)):
        FASTA_obj.insert(i + 1 + ongoingCount, dummySeqs[i])
        ongoingCount += 1

if __name__ == "__main__":
    usage = """%(prog)s receives a directory full of aligned FASTA files as part of the
    Oz Mammals genome project. Its goal is to polish these alignments to make manual inspection
    less tedious. This includes trimming of ends (especially where it's likely to be intron
    sequence), and removal of indel errors.
    
    Note: This should be step 2 in the Oz Mammals project!
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-a", dest="alignmentsDir", required=True,
                help="Specify the directory where aligned FASTA files are located")
    p.add_argument("-o", dest="outputDir", required=True,
                help="Output directory location (default == \"2_prep\")",
                default="2_polish")
    # Opts
    p.add_argument("--codons", dest="showCodonPositions", required=False, action="store_true",
                help="Optionally, specify this flag if you want to produce dummy sequences containing codon numbers",
                default=False)

    args = p.parse_args()
    validate_args(args)
    
    # Locate all files
    files = [os.path.join(args.alignmentsDir, file) for file in os.listdir(args.alignmentsDir)]

    # Load FASTA files
    fastaObjs = []
    for file in files:
        f = ZS_SeqIO.FASTA(file, isAligned=True)
        fastaObjs.append(f)
    
    # Optional step: derive codon position sequences
    if args.showCodonPositions:
        for FASTA_obj in fastaObjs:
            solutionDict = solve_translation_frames(FASTA_obj)
            add_codon_seqs(FASTA_obj, solutionDict)
    
    # Write output files
    ## TBD!

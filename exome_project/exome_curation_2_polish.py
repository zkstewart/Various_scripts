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
    for file in args.gff3s:
        if not os.path.isfile(file):
            print('I am unable to locate the GFF3 file (' + file + ')')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    for dir in args.liftovers:
        if not os.path.isdir(dir):
            print('I am unable to locate the liftover directory (' + dir + ')')
            print('Make sure you\'ve typed the directory location correctly and try again.')
            quit()
    for file in args.genomes:
        if not os.path.isfile(file):
            print('I am unable to locate the genome FASTA file (' + file + ')')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    # Validate liftovers argument
    if len(args.gff3s) != len(args.liftovers):
        print("gff3s and liftovers arguments need to be paired, which means we should receive the same number of values.")
        print("You provided {0} gff3 and {1} liftovers arguments.".format(len(args.gff3s), len(args.liftovers)))
        print("Fix this issue and try again.")
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

def get_cds_coords(gff3File):
    '''
    Parses a GFF3 file and locates CDS start and stop coordinates from within
    each chromosome.
    
    Params:
        gff3File -- a string indicating the known location of a GFF3 file
    Returns:
        cdsCoords -- a dictionary with chromosome: [[start, end], ... [start, end]]
                      structure.
    '''
    cdsCoords = {}
    with open(gff3File, "r") as fileIn:
        for line in fileIn:
            # Handle header lines
            if line.startswith("#") or line == "\n" or line == "\r\n":
                continue
            
            # Extract relevant details
            l = line.rstrip("\r\n").split("\t")
            chrom = l[0]
            annotType = l[2]
            start = int(l[3])
            end = int(l[4])
            
            # Set up storage structure at the start of each chromosome
            if chrom not in cdsCoords:
                cdsCoords[chrom] = []

            # Store CDS details
            if annotType == "CDS":
                coords = [start, end]
                cdsCoords[chrom].append(coords)
            
    return cdsCoords

def trim_SeqIO_FASTA(alignFastaFile, FASTA_obj, cdsCoordsList, liftoverFilesList, genomesList, INTRON_PCT=0.33):
    '''
    trims a ZS_SeqIO.FASTA object to the coordinates range indicated by annotated CDS
    coordinates.
    
    Params:
        alignFastaFile -- a string indicating the FASTA file that generated FASTA_obj.
                          We need this because the file name prefix should match a file
                          in our liftoverFilesList, and that's how we know which sequence
                          to grab.
        FASTA_obj -- a ZS_SeqIO.FASTA instance, containing all the sequences
                     indicated by the idsList parameter.
        cdsCoordsList -- a list containing one or more cdsCoords dictionar(y/ies) which
                     itself has a chromosome: [[start, end], ... [start, end]] structure.
                     The number of dictionaries in this list must match the number of values
                     provided in the idsList parameter.
        liftoverFilesList -- a list containing one or more lists which provide the locations of
                             liftover-ed exon sequences from genomes. This should be paired with
                             the GFF3s that made the cdsCoordsList values.
        genomesList -- a list containing ZS_SeqIO.FASTA instances of the genome sequences paired
                       to the GFF3s that created the cdsCoordsList values.
        INTRON_PCT -- a float value indicating the threshold to use for flagging a sequence
                      as (probably) containing intron sequence; this dictates whether this
                      function will trim based on CDS coordinates, or based on some notion
                      of removing intron flanks.
    Returns:
        wasTrimmed -- a Boolean indicating whether trimming occurred successfully or not.
        reason -- a string explaining why we didn't trim (if relevant) or just None if we did trim,
    '''
    assert len(cdsCoordsList) == len(liftoverFilesList), "Can't trim FASTA if param lengths don't match"
    assert isinstance(INTRON_PCT, float) or isinstance(INTRON_PCT, int), "INTRON_PCT must be a float!"
    assert 0 <= INTRON_PCT <= 1, "INTRON_PCT must be a float between 0 and 1!"
        
    # Locate relevant liftover files
    fastaPrefix = os.path.basename(alignFastaFile).rsplit(".", maxsplit=1)[0]
    foundFiles = [None for _ in liftoverFilesList] # Keep files ordered same as cdsCoordsList
    for i in range(len(liftoverFilesList)):
        for file in liftoverFilesList[i]:
            if os.path.basename(file).startswith(fastaPrefix):
                foundFiles[i] = file
                break
    
    # Abort if we can't trim like this
    noneFound = all([value == None for value in foundFiles])
    if noneFound:
        return False, "No liftover files"
    
    # Parse liftover files
    seqs = [None if file == None else ZS_SeqIO.FASTA(file)[0] for file in foundFiles] # Keep seqs ordered same as cdsCoordsList

    # Get liftover sequence IDs and start:end coordinates
    coords = [None for _ in cdsCoordsList] # Keep coords ordered same as cdsCoordsList
    for i in range(len(seqs)):
        if seqs[i] == None:
            continue
        
        chrom = seqs[i].description.split("chr=")[1].split(" start=")[0]
        start, end = seqs[i].description.split("start=")[1].split(" end=")
        coords[i] = [chrom, int(start), int(end)]
    
    # See if any liftover coords match to CDS coords
    matches = [None for _ in cdsCoordsList] # Keep coords ordered same as cdsCoordsList
    intronicFlags = []
    for i in range(len(coords)):
        if coords[i] == None:
            continue
        chrom, start, end = coords[i]
        
        # Check through cdsCoords dictionary
        cdsCoords = cdsCoordsList[i]
        hits = []
        if chrom not in cdsCoords:
            continue
        else:
            for cdsStart, cdsEnd in cdsCoords[chrom]:
                if cdsStart <= end and start <= cdsEnd:
                    hit = [cdsStart, cdsEnd]
                    if hit not in hits:
                        hits.append(hit)
        
        # If there are multiple hits, pick the best hit
        if len(hits) > 1:
            bestHit = [None, 0]
            for hitStart, hitEnd in hits:
                overlappedPositions = sum([1 for x in range(start, end+1) if x >= hitStart and x <= hitEnd])
                if bestHit[0] == None:
                    bestHit = [[hitStart, hitEnd],  overlappedPositions]
                else:
                    if overlappedPositions > bestHit[1]:
                        bestHit = [[hitStart, hitEnd],  overlappedPositions]
            hits = [bestHit[0]]
        
        # Check to see if we're in what looks to be an intron
        if hits != []:
            "The below list comprehension will check if every sequence position is covered by a hit"
            overlappedPositions = sum([1 for x in range(start, end+1) for hitStart, hitEnd in hits if x >= hitStart and x <= hitEnd])
            seqLen = len(seqs[i].seq)
            pctNoOverlap = (seqLen - overlappedPositions) / seqLen
            if pctNoOverlap >= INTRON_PCT:
                intronicFlags.append(True)
            else:
                intronicFlags.append(False)
        
        # Store hit if possible
        if hits != []:
            matches[i] = hits[0]
    
    # Abort if we didn't find any matches
    noneFound = all([value == None for value in matches])
    if noneFound:
        return False, "No CDS hits"
    
    # If we did find matches, get our sequence IDs
    ids = [None if seqs[i] == None or matches[i] == None else seqs[i].description.split(" ")[1] for i in range(len(seqs))] # seq.description structure is ENSSH... Species_ID... chr=...

    # Locate relevant sequences from FASTA_obj by their ID
    fastaSeqs = [x for x in ids] # copy ids list
    for FastASeq_obj in FASTA_obj:
        if FastASeq_obj.description in fastaSeqs:
            index = fastaSeqs.index(FastASeq_obj.description) # keep things in order always!
            fastaSeqs[index] = FastASeq_obj
    
    # Enter into trimming phase
    if all(intronicFlags):
        print("WARNING: {0} likely contains intronic sequence; will attempt to trim...".format(os.path.basename(alignFastaFile)))
    return _intron_trimming(FASTA_obj, coords, matches, ids, fastaSeqs, genomesList)

def _intron_trimming(FASTA_obj, coords, matches, ids, fastaSeqs, genomesList, PROBLEM_THRESHOLD=0.25):
    '''
    Hidden function for use by trim_SeqIO_FASTA(). This will enact trimming behaviour to remove
    what is likely to be intronic sequence.
    
    Intron trimming assumes the exon is surrounded by one or more intron flanks of varying length.
    It's possible that we have two exons in the sequence, with an internal intron. However, we ...
    
    Params:
        PROBLEM_THRESHOLD -- a float value indicating what proportion of the genomic exon sequence
                             is allowed to NOT align against the MSA exon sequence.
    '''
    assert isinstance(PROBLEM_THRESHOLD, float)
    assert 0 <= PROBLEM_THRESHOLD <= 1
    
    # Extract the exon sequence(s) from the genome(s)
    exonSeqs = []
    for i in range(len(matches)):
        match = matches[i]
        if matches[i] == None:
            exonSeqs.append(None)
            continue
        genomeStart, genomeEnd = match
        chrom, _, _ = coords[i]
        
        seqString = genomesList[i][chrom].seq[genomeStart:genomeEnd+1]
        exonSeqs.append(seqString)
    
    # Align exonSeqs against fastaSeqs
    alignResults = []
    for i in range(len(exonSeqs)):
        exonSeqStr = exonSeqs[i]
        if exonSeqStr == None:
            alignResults.append(None)
            continue
        exonSeqStrRevComp = ZS_SeqIO.FastASeq.get_reverse_complement(None, exonSeqStr) # Don't want to bother figuring out genome strand
        fastaSeqStr = fastaSeqs[i].seq
        
        # Align forward and reverse complement
        fwd_fastaSeqAlign, fwd_exonAlign, _, fwd_score = ssw_parasail(exonSeqStr, fastaSeqStr)
        rev_fastaSeqAlign, rev_exonAlign, _, rev_score = ssw_parasail(exonSeqStrRevComp, fastaSeqStr)

        # Pick our orientation based on score maximisation
        if rev_score > fwd_score:
            fastaSeqAlign, exonAlign, score = rev_fastaSeqAlign, rev_exonAlign, rev_score
        else:
            fastaSeqAlign, exonAlign, score = fwd_fastaSeqAlign, fwd_exonAlign, fwd_score

        # Stop if things didn't work well
        unalignedPct = 1 - (len(exonAlign) / len(exonSeqStr))
        if unalignedPct >= PROBLEM_THRESHOLD:
            alignResults.append(None)
        # If things did work well, figure out where our boundaries should be
        else:
            fastaSeqAlign = fastaSeqAlign.replace("-", "")
            startIndex = fastaSeqStr.find(fastaSeqAlign)
            endIndex = startIndex + len(fastaSeqAlign)
            alignResults.append([startIndex, endIndex])
            
    # Return warning if we failed to find boundaries like this
    if all([r == None for r in alignResults]):
        return False, "Failed to get a good alignment from overlapping CDS"
    
    # Map ungapped boundaries to the aligned sequences
    trueBoundaries = []
    for i in range(len(alignResults)):
        result = alignResults[i]
        if result == None:
            # trueBoundaries.append(None) # don't need to do that here
            continue
        gappedFastaSeqStr = fastaSeqs[i].gap_seq
        
        trueStart, trueEnd = None, None
        sequenceOngoingCount = 0
        for x in range(len(gappedFastaSeqStr)):
            if x == len(gappedFastaSeqStr) - 1:
                if trueEnd == None:
                    trueEnd = x + 1 # +1 to offset range not being inclusive, or something
                break
            
            letter = gappedFastaSeqStr[x]
            if letter == "-":
                continue
            
            if sequenceOngoingCount == result[0]:
                trueStart = x
            elif sequenceOngoingCount == result[1]:
                trueEnd = x + 1 # +1 to offset range not being inclusive, or something

            sequenceOngoingCount += 1
        assert trueStart != None and trueEnd != None
        
        trueBoundaries.append([trueStart, trueEnd])
    
    # Get the most optimistic boundaries for one or more possible boundaries
    trueBoundary = [
        min([b[0] for b in trueBoundaries]),
        max([b[1] for b in trueBoundaries])
    ]
    
    # Trim any excess sequence
    startTrim = trueBoundary[0] # No need to change
    endTrim = len(gappedFastaSeqStr) - trueBoundary[1]
    
    FASTA_obj.trim_left(startTrim, asAligned=True)
    FASTA_obj.trim_right(endTrim, asAligned=True)
    
    return True, None # indicate that trimming worked successfully

def _CDS_coords_trimming(FASTA_obj, coords, matches, ids, fastaSeqs):
    '''
    Hidden function for use by trim_SeqIO_FASTA(). This will enact trimming behaviour based on CDS coords
    when the intron flag was not raised. As such, this function will simply trim the MSA to the best
    estimated coordinates for one or more of the genomic annotation's CDS regions.
    
    Currently deprecated, being kept around in case I want to scavenge its code for later.
    '''
    
    # Find optimistic boundaries for exon CDS
    '''
    i.e., if we have more than one match, pick the longest section of CDS
    based on all matches. We do the difference calculations to give us an offset
    from 0.
    
    For example, if the current start position of the sequence we're assessing
    is considered 0 in the MSA, a value of -4 means we'll start our MSA 4 positions to
    the left of the sequence and no further. If the end position is also considered 0,
    a value of 4 would mean we'll end our MSA 4 positions to the right of the sequence.
    So -ve == left, +ve == right.
    '''
    boundaries = [None for _ in ids]
    for i in range(len(matches)):
        match = matches[i]
        if match == None:
            continue
        chrom, startCoord, endCoord = coords[i]
        
        startDifference = match[0] - startCoord
        endDifference =  match[1] - endCoord
        
        boundaries[i] = [startDifference, endDifference]
    
    # Compare theoretical optimistic boundaries to actual gapped sequences to see where the most optimistic boundaries are
    '''
    Our theoretical boundaries were just computed. Because the gapped sequences might start at
    different positions, we need to relate their starts to the MSA's true coordinates to see
    where any trimming should occur, if at all.
    '''
    
    bestStart, bestEnd = None, None # more -ve is better for start, more +ve is better for end
    for i in range(len(boundaries)):
        fastaSeq = fastaSeqs[i]
        if fastaSeq == None:
            continue
        startOffset, endOffset = boundaries[i]
        
        # Find our new start, with offset adjusted to actual sequence positions excluding gaps
        firstLetter = fastaSeq.seq[0]
        msaStart = fastaSeq.gap_seq.find(firstLetter)
        sequenceOngoingCount = 0
        if startOffset < 0:
            iterRange = range(msaStart-1,-1,-1)
        else:
            iterRange = range(msaStart, len(fastaSeq.gap_seq))
        for x in iterRange:
            if sequenceOngoingCount == startOffset:
                newStart = x
                break
            
            letter = fastaSeq.gap_seq[x]
            if letter == "-":
                continue
            else:
                sequenceOngoingCount += 1
        #newStart = msaStart + startOffset
        
        lastLetter = fastaSeq.seq[-1]
        msaEnd = len(fastaSeq.gap_seq) - fastaSeq.gap_seq[::-1].find(lastLetter)
        sequenceOngoingCount = 0
        if endOffset < 0:
            iterRange = range(msaEnd-1,-1,-1)
        else:
            iterRange = range(msaEnd, len(fastaSeq.gap_seq))
        for x in iterRange:
            if sequenceOngoingCount == endOffset or x == len(fastaSeq.gap_seq) - 1:
                newEnd = x + 1 # counter act range not being inclusive, or something
                break
            
            letter = fastaSeq.gap_seq[x]
            if letter == "-":
                continue
            else:
                sequenceOngoingCount += 1
        #newEnd = msaEnd + endOffset
        
        if bestStart == None:
            bestStart = newStart
        else:
            bestStart = min([bestStart, newStart])
        
        if bestEnd == None:
            bestEnd = newEnd
        else:
            bestEnd = max([bestEnd, newEnd])
        
        gapSeqLength = len(fastaSeq.gap_seq) # Set this variable to use it below
    
    # Make start and end positions possible
    bestStart = 0 if bestStart < 0 else bestStart
    bestEnd = gapSeqLength if bestEnd > gapSeqLength else bestEnd
    
    # Trim any excess sequence
    startTrim = bestStart
    endTrim = gapSeqLength - bestEnd
    
    FASTA_obj.trim_left(startTrim, asAligned=True)
    FASTA_obj.trim_right(endTrim, asAligned=True)
    
    return True # indicate that trimming worked successfully

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
    p.add_argument("-g", dest="gff3s", required=True, nargs="+",
                help="Specify one or more GFF3s to provide CDS boundary information")
    p.add_argument("-lo", dest="liftovers", required=True, nargs="+",
                help="Specify one or more liftover exon FASTAs dirs paired to the gff3s")
    p.add_argument("-ge", dest="genomes", required=True, nargs="+",
                help="Specify one or more genomes to provide paired to the gff3s")
    p.add_argument("-o", dest="outputDir", required=True,
                help="Output directory location (default == \"2_prep\")",
                default="2_polish")
    # Opts
    p.add_argument("--codons", dest="showCodonPositions", required=False, action="store_true",
                help="Optionally, specify this flag if you want to produce dummy sequences containing codon numbers",
                default=False)

    args = p.parse_args()
    validate_args(args)
    
    # Locate all aligned FASTA files
    files = [os.path.join(args.alignmentsDir, file) for file in os.listdir(args.alignmentsDir)]

    # Load aligned FASTA files
    fastaObjs = []
    for file in files:
        f = ZS_SeqIO.FASTA(file, isAligned=True)
        fastaObjs.append(f)
    
    # Parse GFF3 files
    cdsCoordsList = []
    for file in args.gff3s:
        cdsCoordsList.append(get_cds_coords(file))
    
    # Locate liftover FASTA files
    liftoverFilesList = []
    for dir in args.liftovers:
        loFiles = [os.path.join(dir, file) for file in os.listdir(dir)]
        liftoverFilesList.append(loFiles)
    
    # Parse genomes
    genomesList = []
    for file in args.genomes:
        g = ZS_SeqIO.FASTA(file)
        genomesList.append(g)
        
    # Update CDS boundaries for each exon sequence (where possible)
    for i in range(len(files)):
        alignFastaFile = files[i] # i=5 for testing intron trim
        FASTA_obj = fastaObjs[i]
        trimmed, reason = trim_SeqIO_FASTA(alignFastaFile, FASTA_obj, cdsCoordsList, liftoverFilesList, genomesList)
        if trimmed == False:
            if reason == "No liftover files":
                print("Unable to trim {0}; no representative liftovers".format(alignFastaFile))
                ## TBD: Special trimming procedure without liftover files
            else:
                print("Unable to trim {0}; reason = {1}".format(alignFastaFile, reason))
                stophere # For development in IPython

    
    # # Optional step: derive codon position sequences
    # if args.showCodonPositions:
    #     for FASTA_obj in fastaObjs:
    #         solutionDict = solve_translation_frames(FASTA_obj)
    #         add_codon_seqs(FASTA_obj, solutionDict)
    
    # Write output files
    ## TBD!
    print("Program completed successfully!")

#! python3
# exome_curation_autoPolish.py
# Follows up on exome_curation_prep.py to perform some additional
# polishing of the MSAs including trimming and removal of indel errors.

import sys, argparse, os, subprocess, hashlib, time, random
import numpy as np
sys.path.append(os.path.dirname(os.path.dirname(__file__))) # 2 dirs up is where we find dependencies
from Function_packages import ZS_SeqIO
from exome_liftover import ssw_parasail
from copy import deepcopy

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
    for file in args.transcriptomes:
        if not os.path.isfile(file):
            print('I am unable to locate the transcriptome file (' + file + ')')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    # Validate liftovers argument
    if len(args.gff3s) != len(args.liftovers):
        print("gff3s and liftovers arguments need to be paired, which means we should receive the same number of values.")
        print("You provided {0} gff3 and {1} liftovers arguments.".format(len(args.gff3s), len(args.liftovers)))
        print("Fix this issue and try again.")
        quit()
    # Validate INTRON_CHAR argument
    if len(args.INTRON_CHAR) != 1:
        print("INTRON_CHAR must be exactly one character in length")
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
        if FastASeq_obj.seq == "": # skip empty/dummy seqs
            continue
        
        results = [] # contains triples of [seq, frame, hasStopCodon (bool)]
        for frame in range(0, 3):
            seq, strand, frame = FastASeq_obj.get_translation(strand=1, frame=frame)
            results.append([seq, frame, "*" in seq[:-1]]) # don't count the last position since that can be normal
        resultsDict[i] = results
    
    # Loop back through and find easy-to-solve and hard-to-solve scenarios
    solutionDict = {} # contains [index] = [seq, frame, hasStopCodon]
    problemDict = {} # contains [index] = [[seq, frame, hasStopCodon], ... +2 more frames]
    for i in range(len(FASTA_obj)):
        if i not in resultsDict:
            continue
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
            strand = l[6]
            
            # Set up storage structure at the start of each chromosome
            if chrom not in cdsCoords:
                cdsCoords[chrom] = []

            # Store CDS details
            if annotType == "CDS":
                coords = [start, end, strand]
                cdsCoords[chrom].append(coords)
            
    return cdsCoords

def trim_SeqIO_FASTA_evidenced(alignFastaFile, FASTA_obj, cdsCoordsList, liftoverFilesList, genomesList, INTRON_PCT=0.33):
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
        pctTrimmed -- a float value indicating what proportion of the alignment was trimmed.
        reason -- a string explaining why we didn't trim (if relevant) or just None if we did trim.
        intronFlag -- a Boolean indicating whether the intronic sequence flag was raised (i.e.,
                      if this sequence is very likely to contain a significant portion of intron sequence)
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
        return False, 0.0, "No liftover files", False
    
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
        return False, 0.0, "No CDS hits", False
    
    # If we did find matches, get our sequence IDs
    ids = [None if seqs[i] == None or matches[i] == None else seqs[i].description.split(" ")[1] for i in range(len(seqs))] # seq.description structure is ENSSH... Species_ID... chr=...

    # Locate relevant sequences from FASTA_obj by their ID
    fastaSeqs = [x for x in ids] # copy ids list
    for FastASeq_obj in FASTA_obj:
        if FastASeq_obj.description in fastaSeqs:
            index = fastaSeqs.index(FastASeq_obj.description) # keep things in order always!
            fastaSeqs[index] = FastASeq_obj
    
    # Enter into trimming phase
    wasTrimmed, pctTrimmed, reason = _intron_trimming(FASTA_obj, coords, matches, fastaSeqs, genomesList)
    return wasTrimmed, pctTrimmed, reason, all(intronicFlags)

def _intron_trimming(FASTA_obj, coords, matches, fastaSeqs, genomesList, PROBLEM_THRESHOLD=0.66):
    '''
    Hidden function for use by trim_SeqIO_FASTA_evidenced(). This will enact trimming behaviour to remove
    what is likely to be intronic sequence.
    
    Intron trimming assumes the exon is surrounded by one or more intron flanks of varying length.
    It's possible that we have two exons in the sequence, with an internal intron. However, we ...
    
    Params:
        PROBLEM_THRESHOLD -- a float value indicating what proportion of the genomic exon sequence
                             is allowed to NOT align against the MSA exon sequence.
    Returns:
        wasTrimmed -- a Boolean indicating whether trimming occurred successfully or not.
        pctTrimmed -- a float value indicating what proportion of the alignment was trimmed.
        reason -- a string explaining why we didn't trim (if relevant) or just None if we did trim.
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
        exonUnalignedPct = 1 - (len(exonAlign) / len(exonSeqStr))
        msaUnalignedPct = 1 - (len(fastaSeqAlign) / len(fastaSeqStr))
        if exonUnalignedPct >= PROBLEM_THRESHOLD and msaUnalignedPct >= PROBLEM_THRESHOLD:
            alignResults.append(None)
        # If things did work well, figure out where our boundaries should be
        else:
            fastaSeqAlign = fastaSeqAlign.replace("-", "")
            startIndex = fastaSeqStr.find(fastaSeqAlign)
            endIndex = startIndex + len(fastaSeqAlign)
            alignResults.append([startIndex, endIndex])
            
    # Return warning if we failed to find boundaries like this
    if all([r == None for r in alignResults]):
        return False, 0.0, "Failed to get a good alignment from overlapping CDS"
    
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
            elif sequenceOngoingCount + 1 == result[1]: # also +1 here for reasons
                trueEnd = x + 1 # +1 to offset range not being inclusive, or something
                break

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
    
    # Calculate how much was trimmed for tracking purposes
    pctTrimmed = (startTrim + endTrim) / len(gappedFastaSeqStr)
    
    return True, pctTrimmed, None # indicate that trimming worked successfully

def trim_SeqIO_FASTA_denovo(FASTA_obj, EXCLUSION_PCT=0.90, PROBLEM_THRESHOLD=0.66):
    '''
    Trims a ZS_SeqIO.FASTA object to the best guess boundaries of the CDS region.
    It does this without genomic evidence (hence "de novo") by assessment of how
    to get the longest uninterrupted (i.e., no stop codons) region from a MSA.
    
    Params:
        FASTA_obj -- a ZS_SeqIO.FASTA instance
        EXCLUSION_PCT -- a float value indicating the proportion of sequences that
                         must have their stop codons excluded within the selected region.
                         E.g., if 0.90, then 90% of the sequences must NOT contain a stop
                         codon within the selected region.
        PROBLEM_THRESHOLD -- a float value indicating what proportion of the MSA is allowed 
                             to be trimmed using this function. Higher proportions will
                             NOT be trimmed to prevent major malfunctions.
    Returns:
        wasTrimmed -- a Boolean indicating whether trimming occurred successfully or not.
        pctTrimmed -- a float value indicating what proportion of the alignment was trimmed.
        reason -- a string explaining why we didn't trim (if relevant) or just None if we did trim.
    '''
    assert isinstance(EXCLUSION_PCT, float) or isinstance(EXCLUSION_PCT, int)
    assert 0 <= EXCLUSION_PCT <= 1, "EXCLUSION_PCT must be between 0 or 1 (inclusive of 0 and 1)"
    
    assert isinstance(PROBLEM_THRESHOLD, float) or isinstance(PROBLEM_THRESHOLD, int)
    assert 0 <= PROBLEM_THRESHOLD <= 1, "PROBLEM_THRESHOLD must be between 0 or 1 (inclusive of 0 and 1)"
    
    EXCLUSION_PCT = int(EXCLUSION_PCT*100)
    
    # Get sequence translations
    solutionDict = solve_translation_frames(FASTA_obj)
    
    # Locate longest segment boundaries for each sequence that excludes stop codons
    boundaries = []
    for i in range(len(FASTA_obj)):
        if i not in solutionDict:
            continue
        
        FastASeq_obj = FASTA_obj[i]
        translationSeq, frame, hasStopCodon = solutionDict[i]
        
        # If there's no stop codons, just put the maximal sequence length in
        if not hasStopCodon:
            boundaries.append([0, len(FastASeq_obj.gap_seq)])
        # If there are stop codons, find the borders for the longest section excluding stop codons
        else:
            longestSection = max(translationSeq.split("*"), key=len)
            proteinStart = longestSection.find(longestSection)
            proteinEnd = proteinStart + len(longestSection)
            
            # Translate .seq positions to .gap_seq positions
            ongoingPositionCount = 0
            for x in range(len(FastASeq_obj.gap_seq)):
                letter = FastASeq_obj.gap_seq[x]
                if letter == "-":
                    continue
                
                if ongoingPositionCount == proteinStart*3:
                    startIndex = x
                if ongoingPositionCount == proteinEnd*3:
                    endIndex = x # proteinEnd marks the first position of the stop codon, setting this to endIndex will exclude it
                
                ongoingPositionCount += 1
            
            boundaries.append([startIndex, endIndex])
    
    # Find true boundaries which maximise sequence length according to EXCLUSION_PCT threshold
    '''
    It's not precise, but the below code intends to trim our MSA to only the section that best
    accounts for ~{EXCLUSION_PCT}% of the start and stop sites in boundaries. As such, if ~{EXCLUSION_PCT}%
    == 90%, then 90% of our boundaries should have their start site LESS THAN OR EQUAL TO the
    selected start site, and 90% of the boundaries should have their end site GREATER THAN OR
    EQUAL TO the selected end site. It's messy but it should do the job!
    '''
    trueStartIndex = np.percentile([x[0] for x in boundaries], EXCLUSION_PCT)
    trueEndIndex = np.percentile([x[1] for x in boundaries], 100-EXCLUSION_PCT) # Need to get percentile in reverse, kinda
    
    # Check to see how much we will trim
    startTrim = int(trueStartIndex) # No need to change
    endTrim = len(FASTA_obj[0].gap_seq) - int(trueEndIndex) # any FastASeq will do, they should all be the same length
    pctTrimmed = (startTrim + endTrim) / len(FASTA_obj[0].gap_seq)
    if pctTrimmed > PROBLEM_THRESHOLD:
        return False, 0.0, "De novo trimming would trim {0}%; exceeds {1}% threshold".format(round(pctTrimmed*100, 2), PROBLEM_THRESHOLD*100)
    
    # Trim to boundaries
    FASTA_obj.trim_left(startTrim, asAligned=True)
    FASTA_obj.trim_right(endTrim, asAligned=True)
        
    return True, pctTrimmed, None # indicate that trimming worked successfully

def trim_noninformative_flanks(FASTA_obj, INFO_DROP_PCT=0.98):
    '''
    Trims a ZS_SeqIO.FASTA object to remove non-informative flanks i.e.,
    sequence that is ALMOST ENTIRELY only gap ("-") or unknown ("N") across
    entire columns on the left and right sides of the MSA.
    
    Params:
        FASTA_obj -- a ZS_SeqIO.FASTA instance
        INFO_DROP_PCT -- a float value indicating what percentage of the MSA members
                         must be informative to avoid being trimmed. If left at the
                         default 0.98, then 2% or more of the sequences must be informative
                         to avoid having the head or tail position trimmed.
    Returns:
        pctTrimmed -- a float value indicating what proportion of the alignment was trimmed.
    '''
    assert isinstance(INFO_DROP_PCT, float) or isinstance(INFO_DROP_PCT, int)
    assert 0 <= INFO_DROP_PCT <= 1, "INFO_DROP_PCT must be between 0 or 1 (inclusive of 0 and 1)"
    
    startTrim = 0
    for i in range(len(FASTA_obj[0].gap_seq)):
        isNonInformative = [FastASeq_obj.gap_seq[i].lower() in ["n", "-"] for FastASeq_obj in FASTA_obj if FastASeq_obj.id != "Codons"]
        nonInfoPct = sum(isNonInformative) / len(isNonInformative)
        
        if nonInfoPct > INFO_DROP_PCT:
            startTrim += 1
        else:
            break
    
    endTrim = 0
    for i in range(len(FASTA_obj[0].gap_seq)-1, -1, -1):
        isNonInformative = [FastASeq_obj.gap_seq[i].lower() in ["n", "-"] for FastASeq_obj in FASTA_obj if FastASeq_obj.id != "Codons"]
        nonInfoPct = sum(isNonInformative) / len(isNonInformative)
        
        if nonInfoPct > INFO_DROP_PCT:
            endTrim += 1
        else:
            break
    
    # Calculate statistics
    pctTrimmed = (startTrim + endTrim) / len(FASTA_obj[0].gap_seq)
    
    # Trim to boundaries
    FASTA_obj.trim_left(startTrim, asAligned=True)
    FASTA_obj.trim_right(endTrim, asAligned=True)
    
    return pctTrimmed

def _tmp_file_name_gen(prefix, suffix):
        '''
        Hidden function for use by this script.
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

#### RE-DO OF INTRON LOCATING
def get_intron_locations(alignFastaFile, FASTA_obj, cdsCoordsList, liftoverFilesList, transcriptomesFilesList, INTRON_CHAR="4", PROBLEM_THRESHOLD=0.90):
    '''
    Receives a ZS_SeqIO.FASTA object representing a MSA. Using genomic annotation(s) for one or 
    more of the sequences present in the MSA, this function will attempt to locate MSA columnar
    positions that are likely to be intronic. If this isn't possible, this function will
    instead launch into an evidence-less approach to derive probable intron regions based on
    sequence conservation patterns.
    
    ## TBD-rethought params and returns
    
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
        transcriptomesFilesList -- a list containing one or more strings indicating the locations
                                   of transcriptome FASTA files.
        INTRON_CHAR -- a string containing a single character for how the intron position will
                       be represented within the dummy sequence output.
        PROBLEM_THRESHOLD -- a float value indicating what proportion of the MSA is allowed to
                             be predicted as being "intronic" before we suspect that the genomic
                             annotations might be flawed.
    Returns:
        If successful:
            dummyString -- a string value indicating the predicted intron positions as INTRON_CHAR
                           and CDS positions as "-"
            pctIntron -- a float value indicating what proportion of the alignment was predicted
                         to be intronic.
            mode -- a string indicating whether "evidenced" prediction or "evidenceless" prediction
                    was used.
    '''
    assert len(cdsCoordsList) == len(liftoverFilesList), "Can't trim FASTA if param lengths don't match"
    
    # Locate relevant liftover files
    fastaPrefix = os.path.basename(alignFastaFile).rsplit(".", maxsplit=1)[0]
    foundFiles = [None for _ in liftoverFilesList] # Keep files ordered same as cdsCoordsList
    for i in range(len(liftoverFilesList)):
        for file in liftoverFilesList[i]:
            if os.path.basename(file).startswith(fastaPrefix):
                foundFiles[i] = file
                break
    
    # Split operation into evidenced mode and evidence-less mode depending on whether we have genomic evidence
    noneFound = all([value == None for value in foundFiles])
    if noneFound:
        return _evidenceless_intron_locations(FASTA_obj, transcriptomesFilesList, INTRON_CHAR)
    
    # Parse liftover files
    seqs = [None if file == None else ZS_SeqIO.FASTA(file)[0] for file in foundFiles] # Keep seqs ordered same as cdsCoordsList

    # Get liftover sequence IDs and start:end coordinates
    liftoverCoords = [None for _ in cdsCoordsList] # Keep coords ordered same as cdsCoordsList
    for i in range(len(seqs)):
        if seqs[i] == None:
            continue
        
        chrom = seqs[i].description.split("chr=")[1].split(" start=")[0]
        start, end = seqs[i].description.split("start=")[1].split(" end=")
        liftoverCoords[i] = [chrom, int(start), int(end)]
    
    # See if any liftover coords match to GFF3 CDS coords
    gff3Coords = [None for _ in liftoverCoords] # Keep coords ordered same as cdsCoordsList
    for i in range(len(liftoverCoords)):
        if liftoverCoords[i] == None:
            continue
        chrom, start, end = liftoverCoords[i]
        
        # Check through cdsCoords dictionary
        cdsCoords = cdsCoordsList[i]
        hits = []
        if chrom not in cdsCoords:
            continue
        else:
            for cdsStart, cdsEnd, strand in cdsCoords[chrom]:
                if cdsStart <= end and start <= cdsEnd:
                    hit = [cdsStart, cdsEnd, strand]
                    if hit not in hits:
                        hits.append(hit)
        
        # Store hits if possible
        if hits != []:
            gff3Coords[i] = hits
    
    # Once again, split operation into evidenced mode or evidenceless depending on whether we have genomic evidence
    noneFound = all([value == None for value in gff3Coords])
    if noneFound:
        return _evidenceless_intron_locations(FASTA_obj, transcriptomesFilesList, INTRON_CHAR)
    
    # If we did find matches, get our sequence IDs
    ids = [None if seqs[i] == None or gff3Coords[i] == None else seqs[i].description.split(" ")[1] for i in range(len(seqs))] # seq.description structure is ENSSH... Species_ID... chr=...
    
    # Locate relevant sequences from FASTA_obj by their ID
    fastaSeqs = [x for x in ids] # copy ids list
    for FastASeq_obj in FASTA_obj:
        if FastASeq_obj.description in fastaSeqs:
            index = fastaSeqs.index(FastASeq_obj.description) # keep things in order always!
            fastaSeqs[index] = FastASeq_obj
    
    # Enter into evidenced mode of operation
    evidencedResult = evidencedString, evidencedPctIntron, evidencedMode = _evidenced_intron_locations(liftoverCoords, gff3Coords, fastaSeqs, INTRON_CHAR)
    if evidencedPctIntron > PROBLEM_THRESHOLD:
        evidencelessResult = _evidenceless_intron_locations(FASTA_obj, transcriptomesFilesList, INTRON_CHAR)
        if evidencelessResult != False:
            evidencelessString, evidencelessPctIntron, evidencelessMode = evidencelessResult
            if evidencelessPctIntron < evidencedPctIntron:
                return evidencelessResult # If evidenceless succeeds and predicts less intron, we'll go with that
    return evidencedResult # _evidenced function will ALWAYS return a value, but _evidenceless won't necessarily do that

def _evidenced_intron_locations(liftoverCoords, gff3Coords, fastaSeqs, INTRON_CHAR):
    '''
    Hidden function for use by get_intron_locations(). This will enact intron prediction behaviour
    based on the genome annotations.
    
    Mental note: gff3Coords translates to the GFF3 CDS matches!
                 liftoverCoords translates to the genomic position of the HMMER extracted sequence!
    
    Returns:
        If successful:
            dummyString -- a string value indicating the predicted intron positions as INTRON_CHAR
                           and CDS positions as "-"
            pctIntron -- a float value indicating what proportion of the alignment was predicted
                         to be intronic.
            mode -- a string with the value "evidenced" to indicate that trimming occurred using
                    this function, and not with _evidenceless_intron_locations().
        If unsuccessful:
            False -- just a False boolean indicating that this process failed.
    '''

    # Map match coordinates to the MSA sequence positions
    msaCoords = []
    for i in range(len(gff3Coords)):
        if gff3Coords[i] == None:
            continue
        
        # Get relevant details for this match mapping iteration
        match = gff3Coords[i]
        _, coordStart, coordEnd = liftoverCoords[i] # _ == chromosome ID, irrelevant
        
        # Adjust match coordinates to be bounded by the MSA coordinates
        '''
        This is so when we do the below mapping and run things like
            adjustedMatchStart = matchStart - coordStart
        ..., we're always going to find the sequence index within the MSA
        and msaStart and msaEnd (below) should always be set to a non-None value.
        In theory at least.
        '''
        for j in range(len(match)):
            matchStart, matchEnd, strand = match[j]
            match[j] = [max(matchStart, coordStart), min(matchEnd, coordEnd), strand]
        
        # Iterate through gap_seq and perform the mapping of coordinates
        gappedFastaSeqStr = fastaSeqs[i].gap_seq
        fastaSeqStr = fastaSeqs[i].seq
        msaStart, msaEnd = None, None
        for matchStart, matchEnd, strand in match:
            sequenceOngoingCount = 0
            
            adjustedMatchStart = matchStart - coordStart
            adjustedMatchEnd = matchEnd - coordStart
            
            if strand == "-":
                '''
                If we're on the genomic -ve strand, our coordinates need to be inverted to
                accommodate this. The below process does this. Don't ask me to explain it,
                I just validated it by comparison to the old trimming approach that used
                ssw alignment to derive all those coordinate details.
                '''
                _adjustedMatchStart = len(fastaSeqStr) - adjustedMatchEnd
                _adjustedMatchEnd = len(fastaSeqStr) - adjustedMatchStart - 1
                adjustedMatchStart, adjustedMatchEnd = _adjustedMatchStart, _adjustedMatchEnd
            
            for x in range(len(gappedFastaSeqStr)):
                letter = gappedFastaSeqStr[x]
                if letter == "-":
                    continue
                
                if sequenceOngoingCount == adjustedMatchStart:
                    if msaStart == None or x < msaStart:
                        msaStart = x
                if sequenceOngoingCount == adjustedMatchEnd: # also +1 here for reasons ## Turned off +1 to ongoingCount, problems??
                    if msaEnd == None or x > msaEnd:
                        msaEnd = x + 1 # +1 to offset range not being inclusive, or something
                    break
                
                sequenceOngoingCount += 1
            assert msaStart != None and msaEnd != None
        
        msaCoords.append([msaStart, msaEnd])
    
    # Merge MSA match coordinates into a flat list of contiguous coordinates
    msaCoords = _merge_coords_list(msaCoords)
    
    # Create a string with ${INTRON_CHAR}'s for introns, and gaps for CDS positions
    dummyString = ""
    for i in range(len(gappedFastaSeqStr)): # doesn't matter which seq we use, gappedFastaSeqStr is already defined above and is not None
        isCds = any([True for matchStart,matchEnd in msaCoords if i<=matchEnd and i>=matchStart])
        if isCds:
            dummyString += "-"
        else:
            dummyString += INTRON_CHAR
    
    # Calculate how much would be marked as intron
    pctIntron = dummyString.count(INTRON_CHAR) / len(dummyString)
    
    # Return result string otherwise and pct for logging purposes
    return dummyString, pctIntron, "evidenced"

def _evidenceless_intron_locations(FASTA_obj, transcriptomesFilesList, INTRON_CHAR, REPRESENTATIVE_PERCENTILE=90, BLAST_EVALUE_CUTOFF=1e-5, IDENTITY_PCT_CUTOFF=75.00):
    '''
    Hidden function for use by get_intron_locations(). This will enact intron prediction behaviour
    based on the genome annotations.

    Params:
        REPRESENTATIVE_PERCENTILE -- an integer value between 0 and 100 (inclusive) to indicate what
                                     portion of sequences, sorted by length, we should consider to be
                                     representative of the MSA as a whole.
        BLAST_EVALUE_CUTOFF -- a float/int value representing the E-value to use as a cut-off threshold
                               when retrieving BLAST results for potential exon regions.
        IDENTITY_PCT_CUTOFF -- a float/int value representing the BLAST identity percentage to use as a
                               cut-off threshold for retrieving BLAST results for potential exon regions.
    Returns:
        If successful:
            dummyString -- a string value indicating the predicted intron positions as INTRON_CHAR
                           and CDS positions as "-"
            pctIntron -- a float value indicating what proportion of the alignment was predicted
                         to be intronic.
            mode -- a string with the value "evidenced" to indicate that trimming occurred using
                    this function, and not with _evidenceless_intron_locations().
        If unsuccessful:
            False -- just a False boolean indicating that this process failed.
    '''
    assert isinstance(REPRESENTATIVE_PERCENTILE, int)
    assert 0 <= REPRESENTATIVE_PERCENTILE <= 100, "PROBLEM_THRESHOLD must be between 0 or 100 (inclusive of 0 and 1)"
    
    assert isinstance(BLAST_EVALUE_CUTOFF, float) or isinstance(BLAST_EVALUE_CUTOFF, int)
    assert 0 <= BLAST_EVALUE_CUTOFF, "BLAST_EVALUE_CUTOFF must be a positive number"
    
    assert isinstance(IDENTITY_PCT_CUTOFF, float) or isinstance(IDENTITY_PCT_CUTOFF, int)
    assert 0 <= IDENTITY_PCT_CUTOFF, "IDENTITY_PCT_CUTOFF must be a positive number; set to 0 if you don't want an identity percentage cut-off"
    
    # Make sure transcriptome files are ready for BLAST
    for transcriptomeFile in transcriptomesFilesList:
        if not blastdb_exists(transcriptomeFile):
            makeblastdb(transcriptomeFile, "nucl")
    
    # Assess the length distribution of sequences in the MSA
    lengths = []
    for FastASeq_obj in FASTA_obj:
        seqLen = len(FastASeq_obj.seq.upper().replace("-","").replace("N","")) # we want the informative length, not absolute length
        lengths.append(seqLen)
    lengthCutoff = np.percentile([l for l in lengths if l > 0], REPRESENTATIVE_PERCENTILE) # we don't want to bias the distribution with 0's
    
    # Get best representative sequences for BLAST
    representatives = []
    for i in range(len(FASTA_obj)):
        if lengths[i] < lengthCutoff:
            continue
        
        representatives.append(FASTA_obj[i])
    
    # Reduce representatives to a more manageable number with CD-HIT
    dummy_FASTA_obj = deepcopy(FASTA_obj)
    dummy_FASTA_obj.seqs = representatives # create a FASTA and sidestep the init requirement that it be a file... yes, I'm aware of the problem there
    params = [0.8, 5, 1, 0.9, 0.9, 4000, 1] # [0.9, 5, 1, 0.0, 0.0, 400, 1] is default, we want to change aS and aL specifically, but also the similarity param
    clustered_FASTA_obj = get_cdhit_results(r"D:\Bioinformatics\Protein_analysis\cdhit-master\cdhit-master", dummy_FASTA_obj, params)
    
    # Get BLAST results
    blastDict = {}
    for transcriptomeFile in transcriptomesFilesList:
        _tmpBlastDict = get_blast_results(clustered_FASTA_obj, transcriptomeFile, "blastn", BLAST_EVALUE_CUTOFF, 2)
        for key, value in _tmpBlastDict.items():
            if key not in blastDict:
                blastDict[key] = value
            else:
                blastDict[key] += value
    if blastDict == {}:
        return False # if we have no results, we can't use this method of intron prediction
    
    # Predict CDS from BLAST results
    seqCoords = []
    for key, value in blastDict.items():
        for targetID, identityPct, qstart, qend, tstart, tend, evalue in value:
            identityPct, qstart, qend, tstart, tend, evalue = float(identityPct), int(qstart), int(qend), int(tstart), int(tend), float(evalue)
            if identityPct < IDENTITY_PCT_CUTOFF or evalue > BLAST_EVALUE_CUTOFF:
                continue
            # qstart, qend = min(qstart, qend), max(qstart, qend) # shouldn't be necessary, BLAST always orders query coords and does strand stuff to target...?
            seqCoords.append([key, qstart, qend])
    
    # Map CDS coordinates to the MSA sequence positions
    '''
    This works similarly to what's seen in _evidenced_intron_locations(), but we
    don't have to do as complex a set of operations here.
    '''
    msaCoords = []
    for i in range(len(seqCoords)):
        queryID, qstart, qend = seqCoords[i]
        
        # Get the query sequence
        for FastASeq_obj in FASTA_obj:
            if FastASeq_obj.id == queryID:
                break
        assert FastASeq_obj.id == queryID
        gappedFastaSeqStr = FastASeq_obj.gap_seq
        
        # Iterate through gap_seq and perform the mapping of coordinates
        msaStart, msaEnd = None, None
        sequenceOngoingCount = 0
        
        for x in range(len(gappedFastaSeqStr)):
            letter = gappedFastaSeqStr[x]
            if letter == "-":
                continue
            
            if sequenceOngoingCount == qstart - 1: # -1 here since qstart is 1-based
                msaStart = x
            elif sequenceOngoingCount == qend - 1: # -1 here since qend is 1-based
                msaEnd = x + 1 # +1 to offset range not being inclusive, or something
                break
            
            sequenceOngoingCount += 1
        assert msaStart != None and msaEnd != None
        
        msaCoords.append([msaStart, msaEnd])
    
    # Merge MSA match coordinates into a flat list of contiguous coordinates
    msaCoords = _merge_coords_list(msaCoords)
    
    # Create a string with ${INTRON_CHAR}'s for introns, and gaps for CDS positions
    dummyString = ""
    for i in range(len(gappedFastaSeqStr)): # doesn't matter which seq we use, gappedFastaSeqStr is already defined above and is not None
        isCds = any([True for matchStart,matchEnd in msaCoords if i<=matchEnd and i>=matchStart])
        if isCds:
            dummyString += "-"
        else:
            dummyString += INTRON_CHAR
    
    # Calculate how much would be marked as intron
    pctIntron = dummyString.count(INTRON_CHAR) / len(dummyString)
    
    # Return result string and pct for logging purposes
    return dummyString, pctIntron, "evidenceless"

def _merge_coords_list(coords):
    '''
    Helper function pulled out because it's used in _evidenced... and _evidenceless...
    functions.
    
    Because there may be multiple genomes/transcriptomes being used, we might have redundant
    coordinate ranges. This redundancy can be a benefit since we can merge the coordinates
    optimistically to get the largest possible non-intronic regions.
    '''
    mergeIndex = 0
    while mergeIndex < len(coords):
        iterate = True
        
        for i in range(mergeIndex + 1, len(coords)):
            match1Start, match1End = coords[mergeIndex]
            match2Start, match2End = coords[i]
            if match1Start <= match2End and match2Start <= match1End: # i.e., if they overlap
                newMatchStart = min(match1Start, match2Start)
                newMatchEnd = max(match1End, match2End)
                coords[mergeIndex] = [newMatchStart, newMatchEnd]
                del coords[i]
                iterate = False
                break
        
        if iterate:
            mergeIndex += 1
    
    return coords

## BLAST-specific functions
# It would be good for me to make a ZS_BlastIO, wouldn't it?
def blastdb_exists(fastaFile):
    '''
    Relies on a simple assumption that a .nsq file's presence
    indicates that a BLAST database was successfully created
    from the FASTA file.
    '''
    return os.path.isfile("{0}.nsq".format(fastaFile))

def makeblastdb(dbFastaFile, dbType):
    assert dbType.lower() in ['nucl', 'nucleotide', 'prot', 'protein'] or dbType.lower() in ['blastp', 'blastn', 'tblastn', 'tblastx']
    
    # Handle BLAST mode dbType values
    if dbType.lower() in ['blastp', 'blastn', 'tblastn', 'tblastx']:
        if dbType.lower() == "blastp":
            dbType = "prot"
        elif dbType == "blastn":
            dbType = "nucl"
        elif dbType == "tblastn":
            dbType = "nucl"
        elif dbType == "tblastx": # I think?
            dbType = "prot"
    
    cmd = 'makeblastdb -in "{}" -dbtype {} -out "{}"'.format(dbFastaFile, dbType, dbFastaFile)
    run_makedb = subprocess.Popen(cmd, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
    makedbout, makedberr = run_makedb.communicate()
    if makedberr.decode("utf-8") != '':
            raise Exception('Makeblastdb error text below\n' + makedberr.decode("utf-8")) 

def blast(queryFasta, dbFastaFile, blastType, evalue, threads, outFile):
    blastType = blastType.lower()
    assert blastType in ['blastp', 'blastn', 'tblastn', 'tblastx']
    
    cmd = '{} -query "{}" -db "{}" -num_threads {} -evalue {} -out "{}" -outfmt 6'.format(blastType, queryFasta, dbFastaFile, threads, evalue, outFile)
    run_blast = subprocess.Popen(cmd, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
    blastout, blasterr = run_blast.communicate()
    if blasterr.decode("utf-8") != '':
            raise Exception('BLAST error text below\n' + blasterr.decode("utf-8"))

def parse_blast_hit_coords(resultFile, evalueCutoff):
    # Set up
    blastDict = {}
    # Main loop
    with open(resultFile, 'r') as fileIn:
        for line in fileIn:
            # Extract details
            sl = line.split('\t')
            qid = sl[0]
            tid = sl[1]
            identityPct = sl[2]
            qstart = sl[6]
            qend = sl[7]
            tstart = sl[8]
            tend = sl[9]
            evalue = sl[10]
            # Skip if evalue isn't significant
            if float(evalue) > float(evalueCutoff):
                continue
            # Store result
            if qid not in blastDict:
                blastDict[qid] = [[tid, identityPct, qstart, qend, tstart, tend, evalue]]
            else:
                blastDict[qid].append([tid, identityPct, qstart, qend, tstart, tend, evalue])
    # Sort individual entries in blastDict
    for value in blastDict.values():
        value.sort(key = lambda x: float(x[6])) # sort by evalue
    # Return dict
    return blastDict

def get_blast_results(query, target, blastType, evalue, threads):
    '''
    This function pipelines the process of obtaining BLAST results. Intermediate files are
    deleted automatically, and hence this function will only result in the return of the
    blastDict object.
    
    Returns:
        blastDict -- a dict with structure:
            query_id: [[target_id, identity_pct, query_start, query_end, target_start, target_end, evalue], ...]
    '''
    assert blastType.lower() in ['blastp', 'blastn', 'tblastn', 'tblastx']
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
    
    # Get a hash for temporary file creation
    queryForHash = query if isinstance(query, str) \
                         else query.fileOrder[0][0] if type(query).__name__ == "ZS_SeqIO.FASTA" or type(query).__name__ == "FASTA" \
                         else query.id
    targetForHash = target if isinstance(target, str) \
                         else target.fileOrder[0][0] if type(target).__name__ == "ZS_SeqIO.FASTA" or type(target).__name__ == "FASTA" \
                         else target.id
    tmpHash = hashlib.sha256(bytes(queryForHash + targetForHash + str(time.time()) + str(random.randint(0, 100000)), 'utf-8') ).hexdigest()
    
    # If we've received a ZS_SeqIO.FastASeq query, make it into a SeqIO.FASTA object
    if type(query).__name__ == "FastASeq" or type(query).__name__ == "ZS_SeqIO.FastASeq":
        tmpFastaName = _tmp_file_name_gen("query_fasta_tmp" + tmpHash[0:20], "fasta")
        with open(tmpFastaName, "w") as fileOut:
            fileOut.write(">{0}\n{1}\n".format(query.id, query.seq))
        query = ZS_SeqIO.FASTA(tmpFastaName)
        os.unlink(tmpFastaName)
    
    # If we've received a ZS_SeqIO.FastASeq target, make it into a SeqIO.FASTA object
    if type(target).__name__ == "FastASeq" or type(target).__name__ == "ZS_SeqIO.FastASeq":
        tmpFastaName = _tmp_file_name_gen("target_fasta_tmp" + tmpHash[0:20], "fasta")
        with open(tmpFastaName, "w") as fileOut:
            fileOut.write(">{0}\n{1}\n".format(target.id, target.seq))
        target = ZS_SeqIO.FASTA(tmpFastaName)
        os.unlink(tmpFastaName)
    
    # If we've received a ZS_SeqIO.FASTA query, make it into a file
    '''
    If I adapt this code into a BlastIO function, it should do magics with the 
    FASTA object to set alt IDs to unique, non-messy values, then associate the
    results back by index or something. I'm not doing that here.
    '''
    tmpQueryName = None
    if type(query).__name__ == "FASTA" or type(query).__name__ == "ZS_SeqIO.FASTA":
        tmpQueryName = _tmp_file_name_gen("blast_query_tmp" + tmpHash[0:20], "fasta")
        query.write(tmpQueryName)
        
    # If we've received a ZS_SeqIO.FASTA target, make it into a file
    tmpTargetName = None
    if type(target).__name__ == "FASTA" or type(target).__name__ == "ZS_SeqIO.FASTA":
        tmpTargetName = _tmp_file_name_gen("blast_target_tmp" + tmpHash[0:20], "fasta")
        target.write(tmpTargetName)
    
    # Make sure target file is ready for BLAST
    if not blastdb_exists(target):
        makeblastdb(target, blastType.lower())
    
    # Run BLAST
    tmpResultName = _tmp_file_name_gen("blast_result_tmp" + tmpHash[0:20], "outfmt6")
    blast(
        tmpQueryName if tmpQueryName != None else query,
        tmpTargetName if tmpTargetName != None else target,
        blastType.lower(),
        evalue, threads, tmpResultName
    )
    
    # Parse BLAST results
    blastDict = parse_blast_hit_coords(tmpResultName, evalue)
    
    # Clean up & return
    if tmpQueryName != None:
        os.unlink(tmpQueryName)
    if tmpTargetName != None:
        os.unlink(tmpTargetName)
    os.unlink(tmpResultName)
    return blastDict

## CDHIT-specific functions
# It would be good for me to make a ZS_ClustIO, wouldn't it?
def validate_cdhit_params(params):
    assert isinstance(params, list), "params should be a list"
    assert len(params) == 7, "params should be a list containing 7 values"
    c, n, G, aS, aL, M, T = params
    assert isinstance(c, float) and 0<=c<=1, "c param should be a float between 0 and 1 (inclusive)"
    assert isinstance(n, int) and 2<=n<=5, "n param should be an int between 2 and 5 (inclusive)"
    assert isinstance(G, int) and (G == 0 or G == 1), "G param should be an int equal to 0 or 1"
    assert isinstance(aS, float) and 0<=aS<=1, "aS param should be a float between 0 and 1 (inclusive)"
    assert isinstance(aL, float) and 0<=aL<=1, "aL param should be a float between 0 and 1 (inclusive)"
    assert isinstance(M, int) and M >= 100, "M param should be an int of reasonable value (at least greater than 100 Mbytes)"
    assert isinstance(T, int), "T param should be an int"

def cdhit(cdhitDir, outputDir, inputFasta, outputFasta, params):
    '''
    If no params value is provided, defaults will be set as per the below information.
    
    Params = [-c, -n, -G, -aS, -aL, -M, -T]
        or...
    Params = [
        sequence identity threshold (default == 0.9),
        word length (default == 5),
        global sequence identity (0 == local, 1 == global; default == 1),
        alignment coverage for shorter sequence (default == 0.0),
        alignment coverage for longer sequence (default == 0.0),
        max memory in megabytes (default == 400),
        threads (0 uses all available; default == 1)
    ]
    '''
    assert os.path.isfile(os.path.join(cdhitDir, 'cd-hit')) or os.path.isfile(os.path.join(cdhitDir, 'cd-hit.exe')), "cd-hit executable not found at {}".format(cdhitDir)
    assert os.path.isdir(outputDir), "output directory does not exist; I want you to do that for me"
    assert os.path.isfile(inputFasta), "input fasta file does not exist"
    assert os.path.basename(outputFasta) == outputFasta, "output fasta file needs to be just the file name; its location is specified in the outputDir method parameter"
    assert not os.path.isfile(outputFasta), "output fasta file already exists; I don't want to overwrite it"
    if params == None:
        params = [0.9, 5, 1, 0.0, 0.0, 400, 1]
    else:
        validate_cdhit_params(params)
    
    cmd = "{} -i {} -o {} -c {} -n {} -G {} -aS {} -aL {} -M {} -T {}".format(os.path.join(cdhitDir, 'cd-hit'), inputFasta, os.path.join(outputDir, outputFasta), *params)
    run_cdhit = subprocess.Popen(cmd, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
    cdout, cderr = run_cdhit.communicate()
    if cderr.decode("utf-8") != '':
        raise Exception('CD-HIT Error text below' + str(cderr.decode("utf-8")))

def get_cdhit_results(cdhitDir, input, params):
    '''
    If I do this as a ClustIO class or something, params will do away and we'll refer to
    self attributes.
    
    Returns:
        FASTA_obj -- a ZS_SeqIO.FASTA object of the clustered CD-HIT results
    '''
    assert type(input).__name__ == "str" \
        or type(input).__name__ == "FASTA" \
        or type(input).__name__ == "ZS_SeqIO.FASTA" \
        or type(input).__name__ == "FastASeq" \
        or type(input).__name__ == "ZS_SeqIO.FastASeq"
    validate_cdhit_params(params)
    
    # Get a hash for temporary file creation
    inputForHash = input if isinstance(input, str) \
                         else input.fileOrder[0][0] if type(input).__name__ == "ZS_SeqIO.FASTA" or type(input).__name__ == "FASTA" \
                         else input.id
    tmpHash = hashlib.sha256(bytes(inputForHash + str(time.time()) + str(random.randint(0, 100000)), 'utf-8') ).hexdigest()
    
    # If we've received a ZS_SeqIO.FastASeq input, make it into a SeqIO.FASTA object
    if type(input).__name__ == "FastASeq" or type(input).__name__ == "ZS_SeqIO.FastASeq":
        tmpFastaName = _tmp_file_name_gen("input_fasta_tmp" + tmpHash[0:20], "fasta")
        with open(tmpFastaName, "w") as fileOut:
            fileOut.write(">{0}\n{1}\n".format(input.id, input.seq))
        input = ZS_SeqIO.FASTA(tmpFastaName)
        os.unlink(tmpFastaName)
    
    # If we've received a ZS_SeqIO.FASTA input, make it into a file
    '''
    If I adapt this code into a ClustIO function, it should do magics with the 
    FASTA object to set alt IDs to unique, non-messy values, then associate the
    results back by index or something. I'm not doing that here.
    '''
    tmpInputName = None
    if type(input).__name__ == "FASTA" or type(input).__name__ == "ZS_SeqIO.FASTA":
        tmpInputName = _tmp_file_name_gen("cdhit_input_tmp" + tmpHash[0:20], "fasta")
        input.write(tmpInputName)
    
    # Run CD-HIT
    tmpResultName = _tmp_file_name_gen("cdhit_result_tmp" + tmpHash[0:20], "fasta")
    cdhit(
        cdhitDir,
        ".", # "" for working dir being current one
        tmpInputName if tmpInputName != None else input,
        tmpResultName,
        params
    )
    
    # Parse CD-HIT results
    result_FASTA_obj = ZS_SeqIO.FASTA(tmpResultName)

    # Clean up & return
    if tmpInputName != None:
        os.unlink(tmpInputName)
    os.unlink(tmpResultName)
    os.unlink(tmpResultName + ".clstr")
    
    return result_FASTA_obj

#### END RE-DO

if __name__ == "__main__":
    usage = """%(prog)s receives a directory full of aligned FASTA files as part of the
    Oz Mammals genome project. Its goal is to polish these alignments to make manual inspection
    less tedious. This includes prediction of intron regions and inclusion of a >Codons dummy
    sequence at the start of the FASTA with these positions indicated. It additionally trims
    the borders of the MSA to remove non-informative flank regions.
    
    Some programs are assumed to be available in the system PATH. These include the BLAST program
    suite, as well as cd-hit.
    
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
    p.add_argument("-o", dest="outputDir", required=True,
                help="Output directory location (default == \"2_prep\")",
                default="2_polish")
    p.add_argument("-t", dest="transcriptomes", required=True, nargs="+",
                help="Specify one or more transcriptomes for intron prediction in absence of liftover data")
    # Opts
    p.add_argument("--INTRON_CHAR", dest="INTRON_CHAR", required=False,
                help="Optionally, specify what character should be used to denote intron positions (default==\"4\")",
                default="4")
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
    
    # Intron prediction
    log = [[
        "fileName", "predictionMode", "pctIntron",
        "initialLength", "trimmedLength", "pctTrimmed"
    ]]
    for i in range(len(files)):
        # Get details for this MSA
        alignFastaFile = files[i] # i=5 for testing intron trim
        FASTA_obj = fastaObjs[i]
        
        # Perform intron prediction
        result = get_intron_locations(alignFastaFile, FASTA_obj, cdsCoordsList, liftoverFilesList, args.transcriptomes, INTRON_CHAR=args.INTRON_CHAR)
        if result != False:
            dummyString, pctIntron, mode = result
            
            # Insert dummy sequence into FASTA_obj
            dummyFastASeq_obj = ZS_SeqIO.FastASeq("Codons", gapSeq=dummyString)
            FASTA_obj.insert(0, dummyFastASeq_obj)
            if mode == "evidenceless":
                print("evidenceless; i={0}... pctIntron={1}".format(i, result[1]))
        else:
            print("failed to work; i={0}...".format(i))
            stophere
        
        # Perform trimming to remove non-informative sites
        initialLength = len(FASTA_obj[0].gap_seq) # hold onto for later statistics
        pctTrimmed = trim_noninformative_flanks(FASTA_obj)
        endLength = len(FASTA_obj[0].gap_seq)
        
        # Statistics and logging
        log.append([
            os.path.basename(alignFastaFile),
            mode if result != False else ".",
            str(round(pctIntron*100, 2)) if result != False else ".",
            str(initialLength),
            str(endLength),
            str(round(pctTrimmed*100, 2)),
        ])
        
        # Write output FASTA file
        outputFileName = os.path.join(args.outputDir, os.path.basename(alignFastaFile))
        FASTA_obj.write(outputFileName, withDescription=True, asAligned=True)
        
    # Write output logging file
    logFileName = _tmp_file_name_gen("2_polish_log", "txt")
    with open(logFileName, "w") as fileOut:
        fileOut.write("\n".join(["\t".join(l) for l in log]))
    
    print("Program completed successfully!")

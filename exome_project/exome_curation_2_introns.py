#! python3
# exome_curation_polish.py
# Follows up on exome_curation_prep.py to perform some additional
# polishing of the MSAs including prediction of intron regions.

import sys, argparse, os, re, math, random
import numpy as np
sys.path.append(os.path.dirname(os.path.dirname(__file__))) # 2 dirs up is where we find dependencies
from Function_packages import ZS_SeqIO, ZS_BlastIO
from exome_liftover import ssw_parasail

def validate_args(args):
    # Validate input data location
    if not os.path.isdir(args.alignmentsDir):
        print('I am unable to locate the directory where the alignments files are (' + args.alignmentsDir + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.transcriptomeFile):
        print('I am unable to locate the transcriptome file (' + args.transcriptomeFile + ')')
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

## Intron locating based on genomic evidence
def get_intron_locations_genomic(alignFastaFile, FASTA_obj, cdsCoordsList, liftoverFilesList, INTRON_CHAR="4"):
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
        INTRON_CHAR -- a string containing a single character for how the intron position will
                       be represented within the dummy sequence output.
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
    
    # If we have no genomic exons, we can't predict introns; stop here
    noneFound = all([value == None for value in foundFiles])
    if noneFound:
        return False
    
    # Parse liftover files
    seqs = [None if file == None else ZS_SeqIO.FASTA(file)[0] for file in foundFiles] # Keep seqs ordered same as cdsCoordsList
    
    # Get liftover sequence IDs and start:end coordinates
    '''
    The coordinates system of 1_prep work properly through mass testing. We can trust
    these and hence don't need to reobtain them from the genome (saves time and memory)
    '''
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
    
    # If we don't have genomic evidence, we can't predict introns; stop here
    noneFound = all([value == None for value in gff3Coords])
    if noneFound:
        return False
    
    # If we did find matches, get our sequence IDs
    ids = [None if seqs[i] == None or gff3Coords[i] == None else seqs[i].description.split(" ")[1] for i in range(len(seqs))] # seq.description structure is ENSSH... Species_ID... chr=...
    
    # Locate relevant sequences from FASTA_obj by their ID
    fastaSeqs = [x for x in ids] # copy ids list
    for FastASeq_obj in FASTA_obj:
        if FastASeq_obj.description in fastaSeqs:
            index = fastaSeqs.index(FastASeq_obj.description) # keep things in order always!
            fastaSeqs[index] = FastASeq_obj
    
    # Run intron locating method
    '''This pattern lets me plug and play different versions for testing, and
    reduces mental burden by keeping this method to just the locating of relevant
    evidence. The downstream function does the heavy lifting and logical work.'''
    evidencedString, evidencedPctIntron = evidencedString, evidencedPctIntron = _evidenced_intron_locations(liftoverCoords, gff3Coords, fastaSeqs, INTRON_CHAR) # always returns a result
    return evidencedString, evidencedPctIntron, "evidenced"

def _evidenced_intron_locations(liftoverCoords, gff3Coords, fastaSeqs, INTRON_CHAR):
    '''
    Hidden function for use by get_intron_locations_genomic(). This will enact intron prediction behaviour
    based on the genome annotations.
    
    It will also TRIM the MSA to the region covered by at least one genomic exon,
    since we can't know what, outside of this region, is intronic or not.
    
    Mental note: gff3Coords translates to the GFF3 CDS matches!
                 liftoverCoords translates to the genomic position of the HMMER extracted sequence!
    
    Returns:
        If successful:
            dummyString -- a string value indicating the predicted intron positions as INTRON_CHAR
                           and CDS positions as "-"
            pctIntron -- a float value indicating what proportion of the alignment was predicted
                         to be intronic.
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
    
    # Find out how much to trim from the left and right
    '''
    The goal here is to trim the MSA to just the regions which have a genomic exon
    represented. It's NOT to trim to the outer borders of the predicted CDS region,
    since we actually want those regions to be noted as intronic if they exist.
    '''
    startTrim = 0
    for i in range(len(FASTA_obj[0].gap_seq)):
        isGap = [FastASeq_obj.gap_seq[i] == "-" for FastASeq_obj in fastaSeqs if FastASeq_obj != None]
        
        if all(isGap):
            startTrim += 1
        else:
            break
    endTrim = 0
    for i in range(len(FASTA_obj[0].gap_seq)-1, -1, -1):
        isGap = [FastASeq_obj.gap_seq[i] == "-" for FastASeq_obj in fastaSeqs if FastASeq_obj != None]
        
        if all(isGap):
            endTrim += 1
        else:
            break
    
    # Trim the MSA, and the dummy string
    FASTA_obj.trim_left(startTrim, asAligned=True)
    FASTA_obj.trim_right(endTrim, asAligned=True)
    if endTrim == 0:
        dummyString = dummyString[startTrim:]
    else:
        dummyString = dummyString[startTrim:-endTrim]
    
    # Calculate how much would be marked as intron
    pctIntron = dummyString.count(INTRON_CHAR) / len(dummyString)
    
    # Return result string pct for logging purposes
    return dummyString, pctIntron

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

## Intron locating based on parsimony
def solve_translation_frames(FASTA_obj, transcriptomeFile):
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
        transcriptomeFile -- a string indicating the location of a transcriptome
                             file which we can BLAST against when solving hard
                             scenarios.
    Returns:
        solutionDict -- a dictionary with structure like:
            {
                sequence_id: [seq, frame, hasStopCodon (bool)]
            }
    '''
    # Obtain translations in the three strand=1 frames
    resultsDict = {} # same structure as solutionDict, but +2 more frames
    for FastASeq_obj in FASTA_obj:
        if FastASeq_obj.id == "Codons": # skip the >Codons line
            continue
        elif FastASeq_obj.seq == "": # skip empty/dummy seqs
            continue
        
        results = [] # contains triples of [seq, frame, hasStopCodon (bool)]
        for frame in range(0, 3):
            seq, strand, frame = FastASeq_obj.get_translation(strand=1, frame=frame)
            results.append([seq, frame, "*" in seq[:-1]]) # don't count the last position since that can be normal
        resultsDict[FastASeq_obj.id] = results
    
    # Loop back through our resultsDict to find easy-to-solve and hard-to-solve scenarios
    solutionDict, problemDict = _simple_scenario_handler(FASTA_obj, resultsDict)
    
    # If we didn't find many easy-to-solve scenarios with simple handling, get more advanced
    DESIRABLE_NUMBER=5 # arbitrary, should work for the oz mammals exomes project
    if len(solutionDict) <= DESIRABLE_NUMBER:
        solutionDict, problemDict = _advanced_scenario_handler(FASTA_obj, resultsDict, transcriptomeFile)
    
    # Enforce consistency in solutionDict values
    '''
    Sometimes, the true translation can have a stop codon towards the start/end of the sequence
    but otherwise show the same signs as an "easy-to-solve" sequence. It's a bit of a conundrum.
    
    Also, note that memoryDict will behave like a dynamic programming table with structure like:
        {
            thisIndex: {
                frame1: {
                    x1: [...],
                    x2: [...],
                    ...
                },
                frame2: {
                    ...
                },
                frame3: {
                    ...
                }
            },
            otherKey: ...
        }
    '''
    memoryDict = {} # this will speed up program execution to remember solved comparisons
    if len(solutionDict) > 2: # won't work unless we have a few here
        for seqID in solutionDict.keys():
            scores = _calculate_solution_scores(resultsDict, solutionDict, memoryDict, seqID) # will store scores for frame 0, 1, and 2
            
            # Calculate metrics to find the best frame
            bestFrame = _use_scores_metrics_to_get_best_frame(scores, resultsDict[seqID])
            
            # Update the solutionDict value if applicable
            if bestFrame != solutionDict[seqID][1]:
                solutionDict[seqID] = resultsDict[seqID][bestFrame] # we do this since the sequence has changed
                _reset_memory_dict_index(memoryDict, seqID) # and hence any scores associated to it are irrelevant
    
    # Find the best solution to problem sequences
    memoryDict = {} # reset memory dict since solution->solution comparisons don't happen below
    if solutionDict != {}: # prevent error when using an empty solutionDict
        for seqID in problemDict.keys():
            scores = _calculate_solution_scores(resultsDict, solutionDict, memoryDict, seqID) # will store scores for frame 0, 1, and 2
            
            # Calculate metrics to find the best frame
            bestFrame = _use_scores_metrics_to_get_best_frame(scores, resultsDict[seqID])
            
            # Add the problem sequence into our solutionDict
            solutionDict[seqID] = resultsDict[seqID][bestFrame]
    
    return solutionDict

def _simple_scenario_handler(FASTA_obj, resultsDict):
    # Loop back through and find easy-to-solve and hard-to-solve scenarios
    solutionDict = {}
    problemDict = {} # same structure as solutionDict, but +2 more frames
    for FastASeq_obj in FASTA_obj:
        if FastASeq_obj.id not in resultsDict:
            continue
        results = resultsDict[FastASeq_obj.id]
        
        # Skip blank sequences
        if all([r[0] == "" for r in results]) or all([r[0].upper().replace("X","") == "" for r in results]) or all([r[0].upper().replace("N","") == "" for r in results]):
            continue
        
        # Handle easy-to-solve scenarios
        numWithoutStopCodons = sum([1 for r in results if r[2] == False])
        if numWithoutStopCodons == 1:
            solutionDict[FastASeq_obj.id] = [r for r in results if r[2] == False][0]
            continue

        # Note hard-to-solve scenarios
        else:
            problemDict[FastASeq_obj.id] = resultsDict[FastASeq_obj.id]
    
    return solutionDict, problemDict

def _advanced_scenario_handler(FASTA_obj, resultsDict, transcriptomeFile, DESIRABLE_NUMBER=50):
    '''
    Function to find solutions in a more advanced way via BLAST. It will randomly iterate through
    our FASTA object and exert effort to find up to DESIRABLE_NUMBER sequences for our solutionDict.
    The randomness should ensure that we get full representation of the sequences in the MSA.
    For that to prove true, DESIRABLE_NUMBER needs to be set to an appropriate value.
    
    Params:
        DESIRABLE_NUMBER -- an integer value indicating the number of solutions we want to exert
                            BLAST effots to find a solution for.
    '''
    SHORT_SEQ_LEN = 10 # short lengths aren't going to resolve well in BLAST
    
    # Loop back through and find easy-to-solve and hard-to-solve scenarios
    solutionDict = {}
    problemDict = {} # same structure as solutionDict, but +2 more frames
    for FastASeq_obj in sorted(FASTA_obj.seqs,key=lambda _: random.random()):
        if FastASeq_obj.id not in resultsDict:
            continue
        results = resultsDict[FastASeq_obj.id]
        
        # Skip blank sequences
        if all([r[0] == "" for r in results]) or all([r[0].upper().replace("X","") == "" for r in results]) or all([r[0].upper().replace("N","") == "" for r in results]):
            continue
        
        # Handle easy-to-solve scenarios
        numWithoutStopCodons = sum([1 for r in results if r[2] == False])
        if numWithoutStopCodons == 1:
            solutionDict[FastASeq_obj.id] = [r for r in results if r[2] == False][0] # r[2] == hasStopCodons bool value
            continue
        
        # Handle hard-to-solve scenarios
        if len(solutionDict) < DESIRABLE_NUMBER:
            blastResults = [] # holds the best E-value, or +inf if no result found
            for seq, _, _ in results:
                # Skip blank sequences and _almost blank_ sequences
                if len(seq) < SHORT_SEQ_LEN:
                    blastResults.append(math.inf)
                    continue
                # Set up our BLAST handler
                blaster = ZS_BlastIO.BLAST(ZS_SeqIO.FastASeq("eg", seq=seq), transcriptomeFile, "blastp")
                blaster.set_threads(4)
                # Run BLAST
                blastDict, _ = blaster.get_blast_results() # throw away the tmpBlastName since it will be None
                blastResults.append(blastDict["eg"][0][6] if "eg" in blastDict else math.inf)
            if min(blastResults) != math.inf:
                bestFrame = blastResults.index(min(blastResults))
                solutionDict[FastASeq_obj.id] = results[bestFrame]
                continue

        # Note unsolveable (/not solved if DESIRABLE_NUMBER is met) scenarios
        problemDict[FastASeq_obj.id] = resultsDict[FastASeq_obj.id]
    
    return solutionDict, problemDict

def _reset_memory_dict_index(memoryDict, thisSeqID):
    memoryDict[thisSeqID] = {} # forget the parent index
    indicesToForget = []
    for key, frameDict in memoryDict.items():
        for frame, scoreDict in frameDict.items():
            for otherSeqID, score in scoreDict.items():
                if otherSeqID == thisSeqID:
                    indicesToForget.append([key, frame, otherSeqID])
    for key, frame, otherSeqID in indicesToForget:
        del memoryDict[key][frame][otherSeqID] # forget the results for other parents to this index

def _calculate_solution_scores(resultsDict, solutionDict, memoryDict, thisSeqID):
    scores = [[], [], []] # will store scores for frame 0, 1, and 2
    results = resultsDict[thisSeqID]
    for j in range(0, 3):
        targetSeq, _, _ = results[j]
        if targetSeq == "":
            scores[j].append(-math.inf)
            continue
        for otherSeqID, solution in solutionDict.items():
            # Skip self-comparison
            if otherSeqID == thisSeqID: # skip self comparison if relevant
                continue
            # Speed up operation with memoryDict
            try:
                solutionFrame = solution[1]
                try:
                    score = memoryDict[thisSeqID][j][otherSeqID][solutionFrame]
                except:
                    score = memoryDict[otherSeqID][solutionFrame][thisSeqID][j]
            except:
            # Otherwise, do it anew
                solutionSeq, solutionFrame, _ = solution
                _, _, _, score = ssw_parasail(targetSeq, solutionSeq)
                # Store it into memoryDict
                memoryDict.setdefault(thisSeqID, {})
                memoryDict[thisSeqID].setdefault(j, {})
                memoryDict[thisSeqID][j].setdefault(otherSeqID, {})
                memoryDict[thisSeqID][j][otherSeqID][solutionFrame] = score
            scores[j].append(score)
    
    return scores

def _use_scores_metrics_to_get_best_frame(scores, results):
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
        
        # If the two metrics are similar, try to base it off of the number of stop codons
        individualCodonsCount = results[maxIndividualScoresFrame][0].count("*")
        totalCodonsCount = results[maxTotalScoresFrame][0].count("*")
        bestFrame = None
        if (individualScoreDifference / 2) < totalScoreDifference and (totalScoreDifference / 2) < individualScoreDifference:
            if individualCodonsCount < totalCodonsCount:
                bestFrame = maxIndividualScoresFrame
            elif totalCodonsCount < individualCodonsCount:
                bestFrame = maxTotalScoresFrame
        # If the two metrics are dissimilar, or codons count doesn't separarate them, find the best one
        if bestFrame == None and totalScoreDifference > individualScoreDifference: # we're checking to see which metric separates the data most
            bestFrame = maxTotalScoresFrame
        elif bestFrame == None and individualScoreDifference > totalScoreDifference:
            bestFrame = maxIndividualScoresFrame
        # Just pick a best-guess frame [but silently now]
        elif bestFrame == None:
            bestFrame = maxTotalScoresFrame # assume total score sum should be most indicative
    return bestFrame

def trim_intron_locations_denovo(FASTA_obj, transcriptomeFile, EXCLUSION_PCT=0.90):
    '''
    Trims a ZS_SeqIO.FASTA object to the best guess boundaries of the CDS region.
    It does this without genomic evidence (hence "de novo") by assessment of how
    to get the longest uninterrupted (i.e., no stop codons) region from a MSA.
    
    This behaves differently to get_intron_locations_genomic() in that we are
    actually TRIMMING to the predicted CDS region, rather than just noting where
    the intron sequence is with the intron characters ("4"s usually)
    
    Params:
        FASTA_obj -- a ZS_SeqIO.FASTA instance
        transcriptomeFile -- a string indicating the location of a transcriptome file.
                             This will be used to identify the best translation frame for
                             sequences where this isn't readily apparent.
        EXCLUSION_PCT -- a float value indicating the proportion of sequences that
                         must have their stop codons excluded within the predicted CDS region.
                         E.g., if 0.90, then 90% of the sequences must NOT contain a stop
                         codon within the region.
    Returns:
        dummyString -- a string value with "-" for every character position; intended to
                       mimic get_intron_locations_genomic()'s return.
        pctIntron -- a float value indicating what proportion of the alignment was predicted
                     to be intronic and was hence trimmed
        mode -- a string indicating that denovo prediction occurred
    '''
    assert isinstance(EXCLUSION_PCT, float) or isinstance(EXCLUSION_PCT, int)
    assert 0 <= EXCLUSION_PCT <= 1, "EXCLUSION_PCT must be between 0 or 1 (inclusive of 0 and 1)"
    
    EXCLUSION_PCT = int(EXCLUSION_PCT*100)
    
    # Get sequence translations
    solutionDict = solve_translation_frames(FASTA_obj, transcriptomeFile)
    
    # Locate longest segment boundaries for each sequence that excludes stop codons
    boundaries = _get_segment_boundaries(FASTA_obj, solutionDict)
    
    # Find true boundaries which maximise sequence length according to EXCLUSION_PCT threshold
    '''
    It's not precise, but the below code intends to trim our MSA to only the section that best
    accounts for ~{EXCLUSION_PCT}% of the start and stop sites in boundaries. As such, if ${EXCLUSION_PCT}%
    == 90%, then 90% of our boundaries should have their start site LESS THAN OR EQUAL TO the
    selected start site, and 90% of the boundaries should have their end site GREATER THAN OR
    EQUAL TO the selected end site. It's messy but it should do the job!
    '''
    trueStartIndex = np.percentile([x[0] for x in boundaries], EXCLUSION_PCT)
    trueEndIndex = np.percentile([x[1] for x in boundaries], 100-EXCLUSION_PCT) # Need to get percentile in reverse, kinda
    
    # Correct mismatched circumstances
    '''
    A mismatched circumstance is when the startIndex is greater than the endIndex. This
    can occur when a MSA primarily consists of fragments on either side of the MSA
    (start frags and end frags) without any consistent "middle" to the MSA. To correct
    this, we basically invert how we use EXCLUSION_PCT to be much, much more lenient
    than we'd otherwise be.
    '''
    if trueStartIndex >= trueEndIndex:
        trueStartIndex = np.percentile([x[0] for x in boundaries], 100-EXCLUSION_PCT)
        trueEndIndex = np.percentile([x[1] for x in boundaries], EXCLUSION_PCT)
        assert trueStartIndex < trueEndIndex, "Mismatched circumstances still isn't handled, Zac, fix this pls"
    
    # Find out how much to trim from the left and right
    startTrim = int(trueStartIndex) # No need to change
    endTrim = len(FASTA_obj[0].gap_seq) - int(trueEndIndex) # any FastASeq will do, they should all be the same length
    pctIntron = (startTrim + endTrim) / len(FASTA_obj[0].gap_seq) # calculate before we trim
    
    # Trim the MSA, and the dummy string
    FASTA_obj.trim_left(startTrim, asAligned=True)
    FASTA_obj.trim_right(endTrim, asAligned=True)
    
    # Create a dummy string with only gaps to mimic genomic intron location prediction
    dummyString = "-"*len(FASTA_obj[0].gap_seq)
        
    # Return result string pct for logging purposes, as well as method of operation
    return dummyString, pctIntron, "denovo"

def _get_segment_boundaries(FASTA_obj, solutionDict, NATURAL_PCT=0.25):
    '''
    Hidden method of trim_intron_locations_denovo() for getting the bounded region
    in which no stop codons exist for each sequence in the FASTA_obj, depending
    on the results of solve_translation_frames() [i.e., solutionDict].
    
    This is pulled out into a separate function to 1) reduce the mental burden
    of interpretting trim_intron_locations_denovo(), but moreso so 2) I can reuse
    the function in 3_polish.
    
    Params:
        NATURAL_PCT -- a float ratio between 0 and 1 (inclusive) indicating the
                       percentage of sequences that must naturally reach the end
                       i.e., have their best ORF reach up to the end of the MSA,
                       in order for us to extend any truncated sequences to the end
                       as well. Basically, it's just a way of controlling truncated
                       sequences so they can "back up" any sequences which naturally
                       reach the end, or have them provide evidence against those
                       sequences when there aren't enough of them.
    '''
    boundaries = []
    boundaries_with_extension = []
    '''
    When determining boundaries, we generally want to have truncated sequences
    (i.e., those whose ORF does not end in the MSA) support a full length MSA.
    We do this with the boundaries_with_extension list, and determine whether we
    want the extensions to be supported depending on our NATURAL_PCT cut-off.
    
    We use our trunacted_ lists to see 1) how many sequences have nucleotides at
    the start and/or end of the sequence, and 2) of those sequences, how many
    have an ORF which extends beyond those positions? This information tells
    us if we should have truncated sequences which do NOT have nucleotides at
    the start and/or end of the sequence to have their boundaries assumed to
    be at the boundaries of the MSA. This matters because it will, using the
    np.percentile cut-off system seen here and in exome_curation_3_polish.py,
    lead to us "soft-trimming" less of the sequence than we might otherwise.
    '''
    truncated_start = []
    truncated_end = []
    for FastASeq_obj in FASTA_obj:
        if FastASeq_obj.id not in solutionDict:
            continue
        translationSeq, frame, hasStopCodon = solutionDict[FastASeq_obj.id]
        
        # Find the borders for the longest section excluding stop codons
        longestSection = max(translationSeq.split("*"), key=len)
        proteinStart = translationSeq.find(longestSection)
        proteinEnd = proteinStart + len(longestSection)
        
        # Translate .seq positions to .gap_seq positions
        ongoingPositionCount = 0
        startIndex, endIndex = None, None
        for x in range(len(FastASeq_obj.gap_seq)):
            letter = FastASeq_obj.gap_seq[x]
            if letter == "-":
                continue
            
            if ongoingPositionCount == proteinStart*3:
                startIndex = x
            if ongoingPositionCount+1 >= proteinEnd*3: # +1 to handle sequences that end early and have all gap at end
                "The >= is important above for various reasons I've come to semi-understand"
                endIndex = x+1 # proteinEnd marks the first position of the stop codon, setting this to endIndex will exclude it
            
            ongoingPositionCount += 1
        assert startIndex != None and endIndex != None
        
        # Store the unextended boundary details now
        boundaries.append([startIndex, endIndex])
        
        # Figure out if our MSA is touching the starts and ends, and if it is, is it truncated?
        if FastASeq_obj.gap_seq[0] != "-":
            if startIndex == 0:
                truncated_start.append(True)
            else:
                truncated_start.append(False)
        
        if FastASeq_obj.gap_seq[-1] != "-":
            if endIndex == len(FastASeq_obj.gap_seq):
                truncated_end.append(True)
            else:
                truncated_end.append(False)
        
        # Adjust startIndex when the longestSection starts at the beginning of the translationSeq
        if proteinStart == 0 or not hasStopCodon:
            startIndex = 0
        
        # Adjust endIndex when translationSeq does not end in a stop codon
        "This means it's probably been truncated, or it's just an exon MSA stopped at the correct position"
        if proteinEnd == len(translationSeq) or not hasStopCodon:
            endIndex = len(FastASeq_obj.gap_seq)
        
        # Store the extended boundary details now
        boundaries_with_extension.append([startIndex, endIndex])
    
    # Adjust our boundaries to include or exclude extensions appropriately
    if len(truncated_start) > 0:
        if (sum(truncated_start) / len(truncated_start)) > NATURAL_PCT:
            for j in range(len(boundaries)):
                boundaries[j][0] = boundaries_with_extension[j][0]
    
    if len(truncated_end) > 0:
        if (sum(truncated_end) / len(truncated_end)) > NATURAL_PCT:
            for j in range(len(boundaries)):
                boundaries[j][1] = boundaries_with_extension[j][1]

    return boundaries

def check_if_genome_annotation_is_good(dummyString, FASTA_obj, INTRON_CHAR="4"):
    '''
    Extra validation function to check to see if get_intron_locations_genomic()
    is being overzealous with the intron prediction. If it's cutting out _too much_
    sequence, we'll assume that the genome annotation isn't giving us good information
    and we'll disregard it.
    
    The logic goes like this: we want to see how much sequence is gappy in the
    predicted intron region versus the predicted exon region. If our exon region
    is much gappier than our intron region, then we're gonna think something has
    gone wrong and maybe we should fall back to a denovo method which focuses on
    trying to capture the maximum amount of useful information, rather than trusting
    an annotation that isn't guaranteed to be correct.
    
    Params:
        dummyString -- a string representation of where the intron is predicted
                       to be as derived from get_intron_locations_genomic()
        FASTA_obj -- a ZS_SeqIO.FASTA object.
    Returns:
        isGood -- a boolean indicating whether the results from get_intron_locations_genomic()
                  are good (True) or if we should disregard them (False).
    '''
    intronGappiness = []
    exonGappiness = []
    
    for i in range(len(dummyString)):
        intronLetter = dummyString[i]
        gappiness = [1 if FastASeq_obj.gap_seq[i] != "-" else 0 for FastASeq_obj in FASTA_obj]
        gappinessPct = sum(gappiness) / len(gappiness)
        if intronLetter == INTRON_CHAR:
            intronGappiness.append(gappinessPct)
        else:
            exonGappiness.append(gappinessPct)
    
    if intronGappiness == []:
        return True # if no intron is predicted, then there should be no problems
    
    intronGappinessPct = sum(intronGappiness) / len(intronGappiness)
    
    try:
        exonGappinessPct = sum(exonGappiness) / len(exonGappiness)
    except:
        print("check_if_genome_annotation_is_good had an issue with {0}".format(FASTA_obj.fileOrder[0][0]))
        return False # if this fails, then exonGappiness == [] which should never happen...?
    
    if (intronGappinessPct - exonGappinessPct) > exonGappinessPct:
        '''
        This is a simple heuristic. Basically, if the pct's end up being approximately equal,
        then this should never prove true. But if there's a large difference between them,
        then we don't want to rely on the genome's intron prediction for fear of it
        ruining things downstream.
        '''
        return False # not isGood
    else:
        return True # isGood

if __name__ == "__main__":
    usage = """%(prog)s receives a directory full of aligned FASTA files as part of the
    Oz Mammals genome project. Its goal is to annotate where probable intronic sequences
    are in the alignment. When a lack of genomic evidence is available, it will instead
    instead trim the MSA to only the region that most probably contains CDS, excluding
    regions which have an abundance of stop codons indicating that they're non-coding.
    
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
    p.add_argument("-t", dest="transcriptomeFile", required=True,
                help="Specify the location of a single representative (protein) transcriptome file")
    p.add_argument("-o", dest="outputDir", required=True,
                help="Output directory location (default == \"2_prep\")",
                default="2_polish")
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
    
    # Main processing loop
    logFileName = _tmp_file_name_gen("2_introns_log", "txt")
    with open(logFileName, "w") as logFileOut:
        # Begin writing output logging file
        logFileOut.write("{0}\n".format(
            "\t".join([
                    "fileName", "predictionMode", "pctIntron", "initialLength", "trimmedLength", "pctTrimmed"
                ])
            )
        )
        # Intron prediction
        for i in range(len(files)):
            # Get details for this MSA
            alignFastaFile = files[i] # i=5 for testing intron trim
            FASTA_obj = fastaObjs[i]
            initialLength = len(FASTA_obj[0].gap_seq) # hold onto for later statistics
            
            # Perform intron prediction
            result = get_intron_locations_genomic(alignFastaFile, FASTA_obj, cdsCoordsList, liftoverFilesList, INTRON_CHAR=args.INTRON_CHAR)
            if result != False:
                isGood = check_if_genome_annotation_is_good(result[0], FASTA_obj, INTRON_CHAR=args.INTRON_CHAR)
            else:
                isGood = False # prevent poential undefined variable error for first loop
            
            if result == False or not isGood:
                result = trim_intron_locations_denovo(FASTA_obj, args.transcriptomeFile)
            
            dummyString, pctIntron, mode = result
            
            # Insert dummy sequence into FASTA_obj
            dummyFastASeq_obj = ZS_SeqIO.FastASeq("Codons", gapSeq=dummyString)
            FASTA_obj.insert(0, dummyFastASeq_obj)
            
            # Write output FASTA file
            outputFileName = os.path.join(args.outputDir, os.path.basename(alignFastaFile))
            FASTA_obj.write(outputFileName, withDescription=True, asAligned=True)
            
            # Write statistics and logging information
            trimmedLength = len(FASTA_obj[0].gap_seq)
            pctTrimmed = (initialLength - trimmedLength) / initialLength
            logFileOut.write("{0}\n".format(
                "\t".join([
                        os.path.basename(alignFastaFile),
                        mode if result != False else ".",
                        str(round(pctIntron*100, 2)) if result != False else ".",
                        str(initialLength),
                        str(trimmedLength),
                        str(round(pctTrimmed*100, 2))
                    ])
                )
            )
    
    print("Program completed successfully!")

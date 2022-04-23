#! python3
# exome_curation_polish.py
# Follows up on exome_curation_prep.py to perform some additional
# polishing of the MSAs including prediction of intron regions.

import sys, argparse, os, re
import numpy as np
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
def get_intron_locations_genomic(alignFastaFile, FASTA_obj, cdsCoordsList, liftoverFilesList, genomesList, INTRON_CHAR="4"):
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
        genomesList -- a list containing ZS_SeqIO.FASTA instances of the genome sequences paired
                       to the GFF3s that created the cdsCoordsList values.
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
    I don't think my coordinates are reliable, I've probably made a mistake. I'm going
    to directly retrieve them from the genome to prevent any mishaps.
    '''
    liftoverCoords = [None for _ in cdsCoordsList] # Keep coords ordered same as cdsCoordsList
    for i in range(len(seqs)):
        if seqs[i] == None:
            continue
        chrom = seqs[i].description.split("chr=")[1].split(" start=")[0]
        
        # Find the start and end coordinates from the genome sequence
        start, end = None, None
        search = re.search(seqs[i].seq, genomesList[i][chrom].seq, re.IGNORECASE)
        if search == None:
            search = re.search(seqs[i].get_reverse_complement(), genomesList[i][chrom].seq, re.IGNORECASE)
        assert search != None, "Exon sequence doesn't match; Zac has messed something up"
        start, end = search.span()
        start += 1 # we want to keep things 0-based
        
        liftoverCoords[i] = [chrom, start, end]
        
        # Testing section
        'This is to see if my coordinate prediction IS actually faulty'
        seqStart, seqEnd = seqs[i].description.split("start=")[1].split(" end=")
        if int(seqStart) != start or end != int(seqEnd):
            print("These coords don't match; file={0}, start={1} vs {2}, end={3} vs {4}".format(alignFastaFile, start, seqStart, end, seqEnd))
    
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

def get_intron_locations_denovo(FASTA_obj, INTRON_CHAR="4", EXCLUSION_PCT=0.90):
    '''
    Trims a ZS_SeqIO.FASTA object to the best guess boundaries of the CDS region.
    It does this without genomic evidence (hence "de novo") by assessment of how
    to get the longest uninterrupted (i.e., no stop codons) region from a MSA.
    
    Params:
        FASTA_obj -- a ZS_SeqIO.FASTA instance
        INTRON_CHAR -- a string containing a single character for how the intron position will
                       be represented within the dummy sequence output.
        EXCLUSION_PCT -- a float value indicating the proportion of sequences that
                         must have their stop codons excluded within the predicted CDS region.
                         E.g., if 0.90, then 90% of the sequences must NOT contain a stop
                         codon within the region.
    Returns:
        dummyString -- a string value indicating the predicted intron positions as INTRON_CHAR
                       and CDS positions as "-"
        pctIntron -- a float value indicating what proportion of the alignment was predicted
                     to be intronic.
        mode -- a string indicating that denovo prediction occurred
    '''
    assert isinstance(EXCLUSION_PCT, float) or isinstance(EXCLUSION_PCT, int)
    assert 0 <= EXCLUSION_PCT <= 1, "EXCLUSION_PCT must be between 0 or 1 (inclusive of 0 and 1)"
    
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
    accounts for ~{EXCLUSION_PCT}% of the start and stop sites in boundaries. As such, if ${EXCLUSION_PCT}%
    == 90%, then 90% of our boundaries should have their start site LESS THAN OR EQUAL TO the
    selected start site, and 90% of the boundaries should have their end site GREATER THAN OR
    EQUAL TO the selected end site. It's messy but it should do the job!
    '''
    trueStartIndex = np.percentile([x[0] for x in boundaries], EXCLUSION_PCT)
    trueEndIndex = np.percentile([x[1] for x in boundaries], 100-EXCLUSION_PCT) # Need to get percentile in reverse, kinda
    
    # Create a string with ${INTRON_CHAR}'s for introns, and gaps for CDS positions
    dummyString = "".join(["-" if i<trueEndIndex and i>=trueStartIndex else "4" for i in range(len(FASTA_obj[0].gap_seq))])
    
    # Calculate how much would be marked as intron
    pctIntron = dummyString.count(INTRON_CHAR) / len(dummyString)
    
    # Return result string pct for logging purposes, as well as method of operation
    return dummyString, pctIntron, "denovo"

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
    p.add_argument("-ge", dest="genomes", required=True, nargs="+",
                help="Specify one or more genomes to provide paired to the gff3s")
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
    
    # Parse genomes
    genomesList = []
    for file in args.genomes:
        g = ZS_SeqIO.FASTA(file)
        genomesList.append(g)
    
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
        initialLength = len(FASTA_obj[0].gap_seq) # hold onto for later statistics
        
        # Perform intron prediction
        result = get_intron_locations_genomic(alignFastaFile, FASTA_obj, cdsCoordsList, liftoverFilesList, genomesList, INTRON_CHAR=args.INTRON_CHAR)
        if result == False:
            result = get_intron_locations_denovo(FASTA_obj, INTRON_CHAR=args.INTRON_CHAR)
        dummyString, pctIntron, mode = result
        
        # Insert dummy sequence into FASTA_obj
        dummyFastASeq_obj = ZS_SeqIO.FastASeq("Codons", gapSeq=dummyString)
        FASTA_obj.insert(0, dummyFastASeq_obj)
        
        # Statistics and logging
        trimmedLength = len(FASTA_obj[0].gap_seq)
        log.append([
            os.path.basename(alignFastaFile),
            mode if result != False else ".",
            str(round(pctIntron*100, 2)) if result != False else ".",
            initialLength,
            trimmedLength,
            (initialLength - trimmedLength) / initialLength # gives % trimmed
        ])
        
        # Write output FASTA file
        outputFileName = os.path.join(args.outputDir, os.path.basename(alignFastaFile))
        FASTA_obj.write(outputFileName, withDescription=True, asAligned=True)
    
    # Write output logging file
    logFileName = _tmp_file_name_gen("2_polish_log", "txt")
    with open(logFileName, "w") as fileOut:
        fileOut.write("\n".join(["\t".join(l) for l in log]))
    
    print("Program completed successfully!")

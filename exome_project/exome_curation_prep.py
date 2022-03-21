#! python3
# exome_curation_prep.py
# Program to enable manual curation of exome sequencing
# to occur for the Oz Mammals Genomics initiative as part
# of Matthew Phillips and Andrew Baker (et. al.'s) group.

import sys, argparse, os, math
sys.path.append(os.path.dirname(os.path.dirname(__file__))) # 2 dirs up is where we find dependencies
from Function_packages import ZS_SeqIO
from exome_liftover import ssw_parasail

def validate_args(args):
    # Validate input data location
    if not os.path.isdir(args.alignmentsDir):
        print('I am unable to locate the directory where the alignments files are (' + args.alignmentsDir + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.metadataFile):
        print('I am unable to locate the metadata file (' + args.metadataFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate numeric inputs
    if args.chunkSize < 1:
        print("chunkSize should be at least 1")
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

def get_dasyurid_metadata_dict(metadataFile):
    metadataDict = {}
    with open(metadataFile, "r") as fileIn:
        for line in fileIn:
            l = line.rstrip("\r\n ").split("\t")
            metadataDict[l[0]] = l[1]
    return metadataDict

def get_missing_species(fastaFiles, metadataDict):
    '''
    Function is just for debugging and getting an updated metadata file
    '''
    missingIDs = set()
    for fastaFile in fastaFiles:
        FASTA_obj = ZS_SeqIO.FASTA(fastaFile)
        for FastASeq_obj in FASTA_obj:
            description = FastASeq_obj.description
            metadataID = description.split(" ")[1].rsplit("_", maxsplit=1)[0] # Get just the middle part sans _S### suffix
            if metadataID in metadataDict:
                continue
            else:
                missingIDs.add(metadataID)
    return list(missingIDs)

def set_alts(FASTA_obj, metadataDict):
    '''
    This function will receive a single FASTA object and perform
    Oz Mammals-specific processing to get an alt ID appropriately
    set for each FastASeq object.
    '''
    for FastASeq_obj in FASTA_obj:
        description = FastASeq_obj.description
        metadataID = description.split(" ")[1].rsplit("_", maxsplit=1)[0] # Get just the middle part sans _S### suffix
        altID = metadataDict[metadataID]
        FastASeq_obj.alt = altID

def add_missing_seqs(FASTA_obj, sequenceIDs):
    '''
    This function will receive a single FASTA object and perform
    Oz Mammals-specific processing to add missing sequences into
    the FASTA. The sequenceIDs object should be a list derived
    from the metadataDict's keys.
    '''
    # Figure out which IDs are missing
    altIDs = [FastASeq_obj.alt for FastASeq_obj in FASTA_obj]
    missingIDs = [id for id in sequenceIDs if id not in altIDs]
    
    # Add dummy sequences to FASTA object
    mockSequence = "-" * len(FASTA_obj[0].gap_seq) # Mock up a fully-gapped sequence
    for mID in missingIDs:
        dummyFastASeq_obj = ZS_SeqIO.FastASeq(mID, alt=mID, gapSeq = mockSequence)
        FASTA_obj.insert(0, dummyFastASeq_obj) # Just insert at index = 0, we'll sort things later

def get_chunking_points(numberToChunk, chunkSize):
    '''
    This is a general purpose function to take in a number of "things"
    that you want to chunk, and find out how to chunk them evenly.
    
    Params:
        numberToChunk -- an integer value, possibly derived from a list length as example.
        chunkSize -- an integer value for the desired number of things per chunk.
    '''
    assert isinstance(numberToChunk, int)
    assert isinstance(chunkSize, int)
    if numberToChunk <= chunkSize:
        raise Exception("Chunking only valid if chunkSize is smaller than numberToChunk")
    
    numChunks = int(numberToChunk / chunkSize)
    rawNum = numberToChunk / numChunks # This line is more relevant in the multithreading code I took this from, but it's okay to just leave it.
    numRoundedUp = round((rawNum % 1) * numChunks, 0) # By taking the decimal place and multiplying it by the num of chunks, we can figure out how many chunks need to be rounded up
    
    chunkPoints = []
    ongoingCount = 0
    for i in range(numChunks):
        if i+1 <= numRoundedUp: # ngl I don't remember why this is needed; I'm borrowing this code from something I wrote a while back
            chunkPoints.append(math.ceil(rawNum) + ongoingCount) # Round up the rawNum, and also add our ongoingCount which corresponds to the number of things already put into a chunk
            ongoingCount += math.ceil(rawNum)
        else:
            chunkPoints.append(math.floor(rawNum) + ongoingCount)
            ongoingCount += math.floor(rawNum)
        if ongoingCount >= numberToChunk: # Without this check, if we have more chunks than things to chunk, we can end up with "extra" numbers in the list (e.g., [1, 2, 3, 4, 5, 6, 6, 6, 6, 6]).
            break  # This doesn't actually affect program function, but for aesthetic reasons and for clarity of how this function works, I prevent this from occurring.
    
    return chunkPoints

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

def get_mock_genename_sequence(FASTA_obj):
    '''
    Very specific to this Oz Mammals project, this function helps
    to format the gene name dummy sequence that Matt would like
    in the FASTA file
    '''
    geneName = os.path.basename(FASTA_obj.fileOrder[0][0]).split("-mx.fa")[0]
    return geneName.ljust(len(FASTA_obj[0].gap_seq), '-')

if __name__ == "__main__":
    usage = """%(prog)s receives a directory full of aligned FASTA files as part of the
    Oz Mammals genome project. Its goal is to transform these alignments into a format that
    is friendly for manual inspection and modification.
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-a", dest="alignmentsDir", required=True,
                help="Specify the directory where aligned FASTA files are located")
    p.add_argument("-m", dest="metadataFile", required=True,
                help="Specify the metadata file location")
    p.add_argument("-o", dest="outputDir", required=True,
                help="Output directory location")
    # Opts
    p.add_argument("-c", dest="chunkSize", type=int, required=False,
                help="Optionally, specify how many exons should be concatenated in a file",
                default=50)
    p.add_argument("--codons", dest="showCodonPositions", required=False, action="store_true",
                help="Optionally, specify this flag if you want to produce dummy sequences containing codon numbers",
                default=False)

    args = p.parse_args()
    validate_args(args)
    
    # Locate all files
    files = [os.path.join(args.alignmentsDir, file) for file in os.listdir(args.alignmentsDir)]
    
    # Parse metadata file
    metadataDict = get_dasyurid_metadata_dict(args.metadataFile)
    
    ## Debugging: find missing species IDs
    #missingIDs = get_missing_species(files, metadataDict)
    #pyperclip.copy("\n".join(missingIDs))
    
    # Load FASTA files
    fastaObjs = []
    for file in files:
        f = ZS_SeqIO.FASTA(file, isAligned=True)
        fastaObjs.append(f)
    
    # Set FASTA alt IDs
    for FASTA_obj in fastaObjs:
        set_alts(FASTA_obj, metadataDict)
    
    # Optional step: derive codon position sequences
    if args.showCodonPositions:
        for FASTA_obj in fastaObjs:
            solutionDict = solve_translation_frames(FASTA_obj)
            add_codon_seqs(FASTA_obj, solutionDict)
    
    # Add dummy missing sequences
    sequenceIDs = list(metadataDict.values())
    for FASTA_obj in fastaObjs:
        add_missing_seqs(FASTA_obj, sequenceIDs)
    
    # Sort FASTA objects to have consistent internal ordering
    for FASTA_obj in fastaObjs:
        FASTA_obj.seqs.sort(key = lambda x: sequenceIDs.index(x.alt))
    
    # Sort FASTA object list by conservation proportion & GC content
    gcList = []
    for FASTA_obj in fastaObjs:
        gcList.append(FASTA_obj.gc_content())
    
    conserveList = []
    for FASTA_obj in fastaObjs:
        FASTA_obj.generate_consensus()
        conserveList.append(FASTA_obj.conserved_proportion())
    
    sortList = [[i, round(gcList[i], 2), conserveList[i]] for i in range(len(conserveList))]
    sortList.sort(key = lambda x: (-x[1], -x[2]))
    fastaObjs = [fastaObjs[index] for index, _, _ in sortList]
    
    # Figure out how to concatenate FASTA files in chunks of ~50 (or whatever args.chunkSize is)
    chunkPoints = get_chunking_points(len(fastaObjs), args.chunkSize)
    
    # Perform the chunking
    concatFastaObjs = []
    prevChunkStart = 0 # Init as 0 for first loop iteration
    for i in range(0, len(chunkPoints)):
        baseFASTA = fastaObjs[prevChunkStart] # We use the first FASTA in this chunk as our base and concat to it
        geneNameSeq = get_mock_genename_sequence(baseFASTA)
        
        for x in range(prevChunkStart+1, chunkPoints[i]): # +1 since we're already using the first file as our base
            concatFASTA = fastaObjs[x]
            geneNameSeq += get_mock_genename_sequence(concatFASTA)
            for seqIndex in range(len(concatFASTA)):
                baseFASTA.seqs[seqIndex].extend(concatFASTA[seqIndex].gap_seq)
        
        dummyFastASeq_obj = ZS_SeqIO.FastASeq("GeneName", alt="GeneName", gapSeq = geneNameSeq)
        baseFASTA.insert(0, dummyFastASeq_obj)
        
        concatFastaObjs.append(baseFASTA) # Store our modified FASTA with all the concatenation performed
        prevChunkStart = chunkPoints[i]

    # Write output files
    for i in range(0, len(chunkPoints)):
        outputFileName = os.path.join(args.outputDir, "exons_chunk_{0}.fa".format(i+1))
        FASTA_obj = concatFastaObjs[i]
        FASTA_obj.write(outputFileName, withAlt=True, asAligned=True, withConsensus=False)
    

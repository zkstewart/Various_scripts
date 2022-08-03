#! python3
# exome_curation_3_polish.py
# Script which receives an MSA with a >Codons line as its first entry
# and attempts to fix indel errors present in a subset of the alignments.
# If they all have indels, god help you since I can't.

import sys, argparse, os, platform, statistics
import concurrent.futures
import numpy as np
from collections import Counter
from copy import deepcopy
from itertools import repeat

sys.path.append(os.path.dirname(os.path.dirname(__file__))) # 2 dirs up is where we find dependencies
from Function_packages import ZS_SeqIO, ZS_AlignIO, ZS_ORF, ZS_Indel

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
    if not os.path.isdir(args.mafftDir):
        print('I am unable to locate the directory where the MAFFT executables are (' + args.alignmentsDir + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    else:
        if platform.system() == "Windows":
            if not os.path.isfile(os.path.join(args.mafftDir, "mafft.bat")):
                raise Exception("{0} does not exist".format(os.path.join(args.mafftDir, "mafft.bat")))
        else:
            if not os.path.isfile(os.path.join(args.mafftDir, "mafft")) and not os.path.isfile(os.path.join(args.mafftDir, "mafft.exe")):
                raise Exception("mafft or mafft.exe does not exist at {0}".format(args.mafftDir))
    # Validate numeric inputs
    if args.threads < 1:
        print("threads argument must be a positive integer greater than 0")
        quit()
    # Handle file output
    if os.path.isdir(args.outputDir):
        if os.listdir(args.outputDir) != []:
            print(args.outputDir + ' already contains files. I am assuming that we are resuming a run.')
    else:
        try:
            os.mkdir(args.outputDir)
            print("Created '{0}' directory as part of argument validation".format(args.outputDir))
        except:
            print("Wasn't able to create '{0}' directory; does '{1}' actually exist?".format(args.outputDir, os.path.dirname(args.outputDir)))

def polish_MSA_denovo(FASTA_obj, transcriptomeFile, mafftDir):
    '''
    Polishes a ZS_SeqIO.FASTA object to remove probable indel errors from sequences.
    It does this without genomic evidence (hence "de novo") by assessment of how
    a subset of sequences containing stop codons differ to the majority which lack
    that stop codon. When a simple error can be identified, it will rectify it.
    
    It will also perform MSA realignment after codon fixing using codon-based alignment
    with MAFFT.
    
    Params:
        FASTA_obj -- a ZS_SeqIO.FASTA instance
        transcriptomeFile -- a string indicating the location of a transcriptome
                             file which we can BLAST against when solving hard
                             scenarios.
        mafftDir -- a string indicating the location of the MAFFT executable files.
    Returns:
        result_FASTA_obj -- a new ZS_SeqIO.FASTA instance with indels polished and
                            the MSA realigned by codons.
        exonFrames -- a list of dictionaries containing the predicted translation frame
                      for all the sequences this could be derived for. Mostly useful
                      for codon numbering. Key = seqID, value = frame (int). List indexes
                      correspond to exon regions ordered as per _get_exon_coords().
    '''
    assert FASTA_obj[0].id == "Codons", "FASTA lacks a Codons line at its start!"
    mafftAligner = ZS_AlignIO.MAFFT(mafftDir) # set up here for use later
        
    # Get the coordinate spans of exons
    exonCoords = _get_exon_coords(FASTA_obj)
    
    # Get copies of all exon regions from within the FASTA object
    exon_FASTA_objs = []
    for start, end in exonCoords:
        exon_FASTA_obj = deepcopy(FASTA_obj)
        
        # Find out how much to trim from the left and right
        startTrim = start # No need to change
        endTrim = len(exon_FASTA_obj[0].gap_seq) - int(end) # any FastASeq will do, they should all be the same length
        
        # Trim the exon MSA
        exon_FASTA_obj.trim_left(startTrim, asAligned=True)
        exon_FASTA_obj.trim_right(endTrim, asAligned=True)
        exon_FASTA_objs.append(exon_FASTA_obj)
    
    # For each exon region, find and polish indels!
    codonsProblemLeft = []
    codonsProblemRight = []
    exonFrames = []
    for x in range(len(exon_FASTA_objs)):
        exon_FASTA_obj = exon_FASTA_objs[x]
        
        # Remove the codons line
        "It gets in the way now when doing translations, we can put it in again later"
        codons_FastASeq_obj = exon_FASTA_obj[0]
        exon_FASTA_obj.seqs = exon_FASTA_obj.seqs[1:]
        
        # Get sequence translations
        solutionDict = ZS_ORF.ORF.solve_translation_frames(exon_FASTA_obj, transcriptomeFile)
        
        # If solutionDict is empty, abort any handling of this exon
        '''
        If we find an exon where we simply can't get a good CDS, we can't work with it in
        any way. So, we need to bypass a few steps in the process without upsetting anything
        along the way.
        '''
        if solutionDict == {}:
            exonLength = len(exon_FASTA_obj[0].gap_seq)
            startTrim = int(exonLength/2)
            endTrim = exonLength - startTrim
            exon_FASTA_obj.trim_left(startTrim, asAligned=True)
            exon_FASTA_obj.trim_right(endTrim, asAligned=True)
            codonsProblemLeft.append(startTrim)
            codonsProblemRight.append(endTrim)
            continue
        
        # Get the boundaries for this exon region that excludes stop codons
        "We do this here since it wasn't applied to genomic evidenced intron predictions"
        boundaries = ZS_ORF.MSA_ORF.get_segment_boundaries(exon_FASTA_obj, solutionDict)
        
        # Find true boundaries which maximise sequence length according to EXCLUSION_PCT threshold
        "Refer to exome_curation_2_introns for detailed comment if interested"
        EXCLUSION_PCT = 90 # Hard-code this, it isn't changing in 2_introns at the present moment
        trueStartIndex = np.percentile([x[0] for x in boundaries], EXCLUSION_PCT)
        trueEndIndex = np.percentile([x[1] for x in boundaries], 100-EXCLUSION_PCT)
        
        # Correct mismatched circumstances
        "Refer to commenting in exome_curation_2_introns.py; this code is copied from there"
        if trueStartIndex >= trueEndIndex:
            trueStartIndex = np.percentile([x[0] for x in boundaries], 100-EXCLUSION_PCT)
            trueEndIndex = np.percentile([x[1] for x in boundaries], EXCLUSION_PCT)
            assert trueStartIndex < trueEndIndex, "Mismatched circumstances still isn't handled, Zac, fix this pls"
        
        # Trim "problem areas" to focus only on the best CDS region
        '''
        "Problem areas" are those that persist beyond 2_introns operations, which is expected
        to be exclusively found in genomic predicted sequences. The genome annotations might not
        give us excellent boundaries for the CDS region, but they might not be "introns" per se.
        We'll just note these as a special character (HARD-CODED FOR NOW) and ignore them otherwise.
        '''
        startTrim = int(trueStartIndex) # No need to change
        endTrim = len(exon_FASTA_obj[0].gap_seq) - int(trueEndIndex) # any FastASeq will do, they should all be the same length
        exon_FASTA_obj.trim_left(startTrim, asAligned=True)
        exon_FASTA_obj.trim_right(endTrim, asAligned=True)
        
        # Take note of how much we trimmed WITHIN the exon region
        codonsProblemLeft.append(startTrim)
        codonsProblemRight.append(endTrim)
        
        # Predict our solutionDict again if needed
        "If we've changed our sequence region, we might need to update our translations"
        if startTrim != 0 or endTrim != 0:
            solutionDict = ZS_ORF.ORF.solve_translation_frames(exon_FASTA_obj, transcriptomeFile)
        
        # Abort if we found no solutions after trimming
        if solutionDict == None:
            continue
        
        # Loop through solutionDict and polish sequences that need it
        polishedSequences = False
        for problemSeqID, value in solutionDict.items():
            protSeq, frame, hasStopCodon = value
            # Decide if we can skip polishing this sequence
            '''
            This is a "newer" addition to the code. The goal is to detect cryptic indels
            that would normally coast through here otherwise.
            '''
            if not hasStopCodon:
                problemFrames = ZS_Indel.FastASeqAlignmentFrames(exon_FASTA_obj[problemSeqID], solutionDict)
                indelRegions = ZS_Indel.IndelPredictor.predict_indel_regions(problemFrames)
                if indelRegions == []:
                    continue
            
            # Find any fixes suggested by matching against solutionDict values
            fix = _get_fixes(exon_FASTA_obj, problemSeqID, solutionDict)
            if fix == []:
                continue
            
            # If we found a fix to make, get the edited sequence
            polishedSequences = True # if we get to here, we'll want to recompute the solutionDict again
            editedSeq = ZS_Indel.IndelPredictor.use_fix_on_sequence(fix, exon_FASTA_obj[problemSeqID].seq)
            
            # Then, update the sequence in our exon_FASTA_obj
            exon_FASTA_obj[problemSeqID].seq = editedSeq
            exon_FASTA_obj[problemSeqID].gap_seq = None # make sure we know that our gap_seq is no longer valid

        # Recompute solutionDict if necessary
        if polishedSequences == True:
            solutionDict = ZS_ORF.ORF.solve_translation_frames(exon_FASTA_obj, transcriptomeFile)
        
        # Get the frames to translate the sequences into from solutionDict
        frames, thisExonFrames = _get_frames_from_solutionDict(exon_FASTA_obj, solutionDict)
        exonFrames.append(thisExonFrames)
        
        # Align edited exon region
        mafftAligner.run_nucleotide_as_protein(exon_FASTA_obj, strand=1, frame=frames) # always search on strand=1 for an ORF
        
        # Perform extra polishing after alignment
        """
        Once we've performed alignment, we can use the information from the MSA to identify
        cases where the alignment is of poor quality. In these cases, it's likely that we have
        cryptic indels that evaded previous detection due to the lack of a stop codon. This
        extra step will hopefully address most of these situations.
        """
        polishedSequences = _outlier_fix_up(exon_FASTA_obj, solutionDict)
        if polishedSequences == True:
            solutionDict = ZS_ORF.ORF.solve_translation_frames(exon_FASTA_obj, transcriptomeFile)
            frames, thisExonFrames = _get_frames_from_solutionDict(exon_FASTA_obj, solutionDict)
            mafftAligner.run_nucleotide_as_protein(exon_FASTA_obj, strand=1, frame=frames) # re-align again
    
    # Merge the edited exon regions back into the overall sequence
    for x in range(len(exonCoords)-1, -1, -1): # iterate backwards through coords since we know they're ordered ascendingly
        # Get values for this iteration
        exonStart, exonEnd = exonCoords[x]
        problemLeft = codonsProblemLeft[x]
        problemRight = codonsProblemRight[x]
        exon_FASTA_obj = exon_FASTA_objs[x]
        
        # First, get our Codons sequence for this exon region and insert it
        codons_FastASeq_obj = ZS_SeqIO.FastASeq("Codons", gapSeq="-"*len(exon_FASTA_obj[0].gap_seq))
        exon_FASTA_obj.insert(0, codons_FastASeq_obj)
        
        # Then, get our trimmed slices of the original sequence
        "+ and - the problem left/right bits since they modify our coordinates of the exon region"
        left_FASTA_obj = FASTA_obj.slice_cols(0, exonStart + problemLeft)
        left_FASTA_obj[0].gap_seq = left_FASTA_obj[0].gap_seq[0: exonStart] + "5"*problemLeft # replace problemLeft length with 5's
        right_FASTA_obj = FASTA_obj.slice_cols(exonEnd - problemRight, len(FASTA_obj[0].gap_seq))
        right_FASTA_obj[0].gap_seq = "5"*problemRight + right_FASTA_obj[0].gap_seq[problemRight: ] # replace problemRight length with 5's
        
        # Write the exon and right FASTA objects to temporary files
        "We only do this because of how the FASTA.concat() method is implemented"
        exonTmpFileName = _tmp_file_name_gen("exon", "fasta")
        rightTmpFileName = _tmp_file_name_gen("right", "fasta")
        exon_FASTA_obj.write(exonTmpFileName, asAligned=True)
        right_FASTA_obj.write(rightTmpFileName, asAligned=True)
        
        # Merge the left, right, and exon FASTAs
        left_FASTA_obj.concat(exonTmpFileName)
        left_FASTA_obj.concat(rightTmpFileName)
        
        # Overwrite FASTA_obj with the new exon-edited version & clean up
        FASTA_obj = left_FASTA_obj
        os.unlink(exonTmpFileName)
        os.unlink(rightTmpFileName)
        
    return FASTA_obj, exonFrames # this has the same variable name, but its value is DIFFERENT than the input

def _get_exon_coords(FASTA_obj):
    '''
    Hidden method for getting exon coordinates based on the Codons line.
    
    NOTE: This will give 0-based numbers appropriate for range i.e.,
    end not inclusive
    '''
    exonCoords = []
    start, end = None, None
    for i in range(len(FASTA_obj[0].gap_seq)):
        letter = FASTA_obj[0].gap_seq[i]
        if start == None and (letter == "-" or letter in ["1", "2", "3"]):
            start = i
        elif start != None and (letter != "-" and letter not in ["1", "2", "3"]):
            end = i # since we've gone past the border of the CDS region, taking 'i' will be range acceptable
        if start != None and end != None:
            exonCoords.append([start, end])
            start, end = None, None # reset, so this gives us the chance to find more exons
    if start != None: # this means we found a start but never got to the stop
        exonCoords.append([start, i+1]) # +1 to make end non-inclusive
    return exonCoords

def _get_frames_from_solutionDict(exon_FASTA_obj, solutionDict):
    '''
    Hidden helper function to look through sequences in our exon FASTA object
    and identify the translation frame for these sequences from the solutionDict.
    
    Returns:
        frames -- a list containing frames (as 0-based integers) for sequences
                  ordered as in exon_FASTA_obj.
        thisExonFrames -- a dict containing the same information as in the frames
                          list but instead it's indexed by the sequence ID.
    '''
    frames = []
    thisExonFrames = {}
    for FastASeq_obj in exon_FASTA_obj:
        seqID = FastASeq_obj.id
        if seqID not in solutionDict:
            frames.append(None)
        else:
            frames.append(solutionDict[seqID][1])
            thisExonFrames[seqID] = solutionDict[seqID][1]
    return frames, thisExonFrames

def _get_fixes(FASTA_obj, problemSeqID, solutionDict):
    # Check matches to see what fix they suggest
    nuclSeq = FASTA_obj[problemSeqID].seq
    matchesFiltered = ZS_Indel.IndelPredictor.obtain_matches_via_ssw(FASTA_obj, FASTA_obj[problemSeqID], solutionDict, FILTER=True)
    matchesAll = ZS_Indel.IndelPredictor.obtain_matches_via_ssw(FASTA_obj, FASTA_obj[problemSeqID], solutionDict, FILTER=False)
    
    fixFiltered = ZS_Indel.IndelPredictor.get_fix_from_ssw_matches(matchesFiltered, nuclSeq)
    fixAll = ZS_Indel.IndelPredictor.get_fix_from_ssw_matches(matchesAll, nuclSeq)
    
    # If no fixes are suggested, return a blank list
    if fixFiltered == []:
        "Only check filtered fixes since we want to know the high-quality matches suggest a fix"
        return []
    # Edit sequences using the fixes
    fixFilteredSeq = ZS_Indel.IndelPredictor.use_fix_on_sequence(fixFiltered, nuclSeq)
    fixFilteredFastASeq_obj = ZS_SeqIO.FastASeq(id=problemSeqID, seq=fixFilteredSeq)
    
    fixAllSeq = ZS_Indel.IndelPredictor.use_fix_on_sequence(fixAll, nuclSeq)
    fixAllFastASeq_obj = ZS_SeqIO.FastASeq(id=problemSeqID, seq=fixAllSeq)
    
    # Ensure that fixing with all doesn't ruin things
    fixFilteredTranslation = fixFilteredFastASeq_obj.get_translation(findBestFrame=True)[0]
    fixAllTranslation = fixAllFastASeq_obj.get_translation(findBestFrame=True)[0]
    
    allIsBad = True if "*" in fixAllTranslation and not "*" in fixFilteredTranslation else False

    # Ensure that fix is sensible
    "Code originally here was cut"
    filteredIsBad = False
    
    # Return the all fix if it worked out, otherwise return the filtered fix
    return fixFiltered if allIsBad else fixAll if not filteredIsBad else []

def _outlier_fix_up(exon_FASTA_obj, solutionDict):
    '''
    Function is under development and testing and uncertain to be part of the final module.
    '''
    # Get the position frequencies for all MSA columns
    frequencies = {}
    for x in range(len(exon_FASTA_obj[0].gap_seq)): # doesn't matter which one we use
        letters = [FastASeq_obj.gap_seq[x] for FastASeq_obj in exon_FASTA_obj.seqs]
        counts = Counter(letters)
        _frequency = {}
        for _letter, _count in counts.items():
            _frequency[_letter] = _count / len(letters)
        frequencies[x] = _frequency
    
    # Loop through sequences and score them
    scores = {}
    for FastASeq_obj in exon_FASTA_obj:
        if FastASeq_obj.id == "Codons":
            continue
        
        exonSeq = FastASeq_obj.gap_seq
        if len(exonSeq.upper().replace("-", "").replace("N","")) < 3: # no point doing anything to short sequences
            continue
        
        # Find the sequence's details for the region it's present in
        startPosition = exonSeq.find(FastASeq_obj.seq[0])
        endPosition = [j for j in range(len(exonSeq)) if exonSeq[j] != "-"][-1]
        
        # Score the region & store it
        thisScore = []
        for x in range(startPosition, endPosition+1): # only score the region of sequence that's present
            letter = exonSeq[x]
            if letter == "N": # don't penalise N's, we want to find non-matching nucleotides
                continue
            thisScore.append(frequencies[x][letter])
        scores[FastASeq_obj.id] = sum(thisScore) / len(thisScore)
    
    # Find outliers from score distribution
    medianScore = statistics.median(scores.values())
    maxScore = max(scores.values())
    maxMedianDiff = maxScore - medianScore
    lowerCutoff = medianScore - maxMedianDiff
    possibleOutliers = [seqID for seqID, score in scores.items() if score < lowerCutoff]
    
    # Try to find indel errors in the outliers
    polishedSequences = False
    for problemSeqID in possibleOutliers:
        if problemSeqID not in solutionDict: # can't help this
            continue
        
        # Find any fixes suggested by matching against solutionDict values
        fix = _get_fixes(exon_FASTA_obj, problemSeqID, solutionDict)
        if fix == []:
            continue
        
        # If we found a fix to make, get the edited sequence
        polishedSequences = True # if we get to here, we'll want to recompute the solutionDict again
        editedSeq = ZS_Indel.IndelPredictor.use_fix_on_sequence(fix, exon_FASTA_obj[problemSeqID].seq)
        
        # Then, update the sequence in our exon_FASTA_obj
        exon_FASTA_obj[problemSeqID].seq = editedSeq
        exon_FASTA_obj[problemSeqID].gap_seq = None
    return polishedSequences

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

def add_codon_numbers(FASTA_obj, exonFrames):
    '''
    This function will modify the FASTA object's Codons line to have codon numbering.
    
    Developer's note: I used to have a bug in here with relation to the frame.
    Conceptually, if we are obtaining frame 2 (0-based), we're only trimming 1
    position from the start. And if we're obtaining frame 1, we're trimming 2
    positions from the start. Frame 0 involves no trimming.
    
    Returns no values since it modifies FASTA_obj directly.
    
    Params:
        FASTA_obj -- a ZS_SeqIO.FASTA object.
        exonFrames -- a list of dictionaries containing the predicted translation frame
                      for all the sequences this could be derived for via the polishing
                      function.
    '''
    
    # Get the exon coordinates for numbering
    exonCoords = _get_exon_coords(FASTA_obj)
    
    # Loop through exon coordinates to figure out which frame the region is in
    exonBestFrames = []
    for j in range(len(exonCoords)):
        exonStart, exonEnd = exonCoords[j]
        # Get a slice of the exon region
        exon_FASTA_obj = FASTA_obj.slice_cols(exonStart, exonEnd)
        
        # Check what the starting frame is predicted to be for each sequence
        startingFrames = []
        for exon_FastASeq_obj in exon_FASTA_obj:
            if exon_FastASeq_obj.id == "Codons":
                continue
            elif len(exon_FastASeq_obj.gap_seq.replace("-","")) < 3: # skip empty sequences / ones without a single codon
                continue
            
            # Get the translation frame
            if exon_FastASeq_obj.id in exonFrames[j]:
                frame = exonFrames[j][exon_FastASeq_obj.id]
            else:
                _, _, frame = exon_FastASeq_obj.get_translation(True)
            
            # Conceptually map the "frame" to how we want it displayed
            '''
            This is where I've struggled with bugs at various places in this project.
            When we say a translation starts at frame 2 (1-based), it's not correct
            to number the first location of the sequence as position 2 of the codon!
            In actual fact, the codon position would be the third in a codon. Why?
            Because if the translation STARTS at frame 2 of the sequence, we essentially
            need to TRIM ONE position from the start to get there. The amount that we trim
            ends up being the codon position.
            '''
            frame = 0 if frame == 0 else 1 if frame == 2 else 2
            
            # Map the translation start to MSA coordinates
            '''
            Not every sequence will start at the very border of the exon region.
            In truncated sequences, we need to know where it is starting.
            '''
            for x in range(len(exon_FastASeq_obj.gap_seq)):
                letter = exon_FastASeq_obj.gap_seq[x]
                if letter != "-":
                    msaStart = x
                    break
            
            # Backtrack the frame to the start of the MSA if necessary
            '''
            For truncated sequences, we can derive the intended frame for the first
            MSA position based on where its start position and translation frame is.
            '''
            for x in range(msaStart-1, -1, -1):
                if frame == 0:
                    frame = 2
                else:
                    frame -= 1
            
            startingFrames.append(frame)
        
        # Figure out which frame is most supported by consensus
        mostCommonStartingFrame = Counter(startingFrames).most_common(1)[0][0]
        exonBestFrames.append(mostCommonStartingFrame)
    
    # For each exon region, add the frame numbering in
    for i in range(len(exonCoords)):
        exonStart, exonEnd = exonCoords[i]
        frame = exonBestFrames[i]
        
        ongoingFrameCount = frame + 1 # we'll write it 1-based to the >Codons line
        for x in range(exonStart, exonEnd):
            FASTA_obj[0].gap_seq = FASTA_obj[0].gap_seq[:x] + str(ongoingFrameCount) + FASTA_obj[0].gap_seq[x+1:]
            ongoingFrameCount = 1 if ongoingFrameCount == 3 else ongoingFrameCount + 1

def polishing_handler(args, alignFastaFile, FASTA_obj):
    # Get details for this MSA
    outputFileName = os.path.join(args.outputDir, os.path.basename(alignFastaFile))
    
    # Skip if we've already processed this file
    if os.path.isfile(outputFileName):
        return
    
    # Perform polishing procedure
    FASTA_obj, exonFrames = polish_MSA_denovo(FASTA_obj, args.transcriptomeFile, args.mafftDir)
    
    # Number codons
    add_codon_numbers(FASTA_obj, exonFrames)
    
    # Write output FASTA file
    FASTA_obj.write(outputFileName, withDescription=True, asAligned=True)

def main():
    usage = """%(prog)s receives a directory full of aligned FASTA files as part of the
    Oz Mammals genome project. Its goal is to remove indel errors from a MSA that has
    been previously subjected to exome_curation_2_introns.py. Specifically, it uses the
    information of which regions are coding to identify probable indel errors present in
    a minority of the sequences.
    
    It will also realign the MSA using a codon-aware method (i.e., translation to protein,
    alignment, then untranslating the alignment) and add the codon numbering into the
    >Codons line.
    
    Note: This should be step 3 in the Oz Mammals project!
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-a", dest="alignmentsDir", required=True,
                help="Specify the directory where aligned FASTA files are located")
    p.add_argument("-t", dest="transcriptomeFile", required=True,
                help="Specify the location of a single representative (protein) transcriptome file")
    p.add_argument("-m", dest="mafftDir", required=True,
                help="Specify the directory where MAFFT executables are located")
    p.add_argument("-o", dest="outputDir", required=True,
                help="Output directory location (default == \"3_polish\")",
                default="3_polish")
    # Opts
    p.add_argument("--threads", dest="threads", required=False, type=int, default=1,
                help="Optionally specify how many threads to run with")
    args = p.parse_args()
    validate_args(args)
    
    # Locate all aligned FASTA files
    files = [os.path.join(args.alignmentsDir, file) for file in os.listdir(args.alignmentsDir)]
    
    # Load aligned FASTA files
    fastaObjs = []
    for file in files:
        f = ZS_SeqIO.FASTA(file, isAligned=True)
        f.make_uppercase()
        fastaObjs.append(f)
    
    # Polishing
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as executor:
        executor.map(polishing_handler, repeat(args, len(files)), files, fastaObjs)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

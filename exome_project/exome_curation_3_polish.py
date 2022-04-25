#! python3
# exome_curation_3_polish.py
# Script which receives an MSA with a >Codons line as its first entry
# and attempts to fix indel errors present in a subset of the alignments.
# If they all have indels, god help you since I can't.

import sys, argparse, os, re, platform
import numpy as np
sys.path.append(os.path.dirname(os.path.dirname(__file__))) # 2 dirs up is where we find dependencies
from Function_packages import ZS_SeqIO, ZS_AlignIO
from exome_liftover import ssw_parasail
from copy import deepcopy
from exome_curation_2_introns import solve_translation_frames, _get_segment_boundaries

def validate_args(args):
    # Validate input data location
    if not os.path.isdir(args.alignmentsDir):
        print('I am unable to locate the directory where the alignments files are (' + args.alignmentsDir + ')')
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

def polish_MSA_denovo(FASTA_obj, mafftDir):
    '''
    Polishes a ZS_SeqIO.FASTA object to remove probable indel errors from sequences.
    It does this without genomic evidence (hence "de novo") by assessment of how
    a subset of sequences containing stop codons differ to the majority which lack
    that stop codon. When a simple error can be identified, it will rectify it.
    
    It will also perform MSA realignment after codon fixing using codon-based alignment
    with MAFFT.
    
    Params:
        FASTA_obj -- a ZS_SeqIO.FASTA instance
        mafftDir -- a string indicating the location of the MAFFT executable files.
    Returns:
        result_FASTA_obj -- a new ZS_SeqIO.FASTA instance with indels polished and
                            the MSA realigned by codons.
    '''
    
    assert FASTA_obj[0].id == "Codons", "FASTA lacks a Codons line at its start!"
    mafftAligner = ZS_AlignIO.MAFFT(mafftDir) # set up here for use later
    
    # Get the coordinate spans of exons
    'NOTE: This will give 0-based numbers appropriate for range i.e., end not inclusive'
    exonCoords = []
    start, end = None, None
    for i in range(len(FASTA_obj[0].gap_seq)):
        letter = FASTA_obj[0].gap_seq[i]
        if start == None and letter == "-": # any other letter at this point will be the intron letter
            start = i
        elif start != None and letter != "-":
            end = i # since we've gone past the border of the CDS region, taking 'i' will be range acceptable
        if start != None and end != None:
            exonCoords.append([start, end])
            start, end = None, None # reset, so this gives us the chance to find more exons
    if start != None: # this means we found a start but never got to the stop
        exonCoords.append([start, i+1]) # +1 to make end non-inclusive
    
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
    for x in range(len(exon_FASTA_objs)):
        exon_FASTA_obj = exon_FASTA_objs[x]
        
        # Remove the codons line
        "It gets in the way now when doing translations, we can put it in again later"
        codons_FastASeq_obj = exon_FASTA_obj[0]
        exon_FASTA_obj.seqs = exon_FASTA_obj.seqs[1:]
        
        # Get sequence translations
        solutionDict = solve_translation_frames(exon_FASTA_obj)
        
        # Get the boundaries for this exon region that excludes stop codons
        "We do this here since it wasn't applied to genomic evidenced intron predictions"
        boundaries = _get_segment_boundaries(exon_FASTA_obj, solutionDict)
        
        # Find true boundaries which maximise sequence length according to EXCLUSION_PCT threshold
        "Refer to exome_curation_2_introns for detailed comment if interested"
        EXCLUSION_PCT = 90 # Hard-code this, it isn't changing in 2_introns at the present moment
        trueStartIndex = np.percentile([x[0] for x in boundaries], EXCLUSION_PCT)
        trueEndIndex = np.percentile([x[1] for x in boundaries], 100-EXCLUSION_PCT)
        
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

        # Predict our solutionDict again if needed
        "If we've changed our sequence region, we might need to update our translations"
        if startTrim != 0 or endTrim != 0:
             solutionDict = solve_translation_frames(exon_FASTA_obj)
        
        # Loop through solutionDict and polish sequences that need it
        for problemIndex, value in solutionDict.items():
            protSeq, frame, hasStopCodon = value
            nuclSeq = exon_FASTA_obj[problemIndex].seq
            if not hasStopCodon: # no changes needed here
                continue
            
            # Loop through all solutionDict values and find the best match to this problem one
            matches = []
            for targetIndex, _value in solutionDict.items():
                if targetIndex == problemIndex:
                    continue
                _, _, targetHasStopCodon = _value
                if targetHasStopCodon:
                    continue
                else:
                    targetNuclSeq = exon_FASTA_obj[targetIndex].seq
                    targetAlign, problemAlign, startIndex, score = ssw_parasail(nuclSeq, targetNuclSeq) # this function has poorly ordered outputs
                    matches.append([problemAlign, targetAlign, startIndex, score])
            matches.sort(key = lambda x: -x[3]) # order by score
            
            # Check all matches to see what fix they suggest
            SCORE_CUTOFF = np.percentile([x[3] for x in matches], 100-EXCLUSION_PCT) # borrow EXCLUSION_PCT which is 90
            fix = _get_suggested_fix_from_ssw_matches(matches, nuclSeq, SCORE_CUTOFF)
            if fix == []:
                continue
            
            # If we found a fix to make, get the edited sequence
            editedSeq = _enact_fix_to_seq(fix, exon_FASTA_obj[problemIndex].seq)
            
            # Then, update the sequence in our exon_FASTA_obj
            exon_FASTA_obj[problemIndex].seq = editedSeq
            exon_FASTA_obj[problemIndex].gap_seq = None # make sure we know that our gap_seq is no longer valid
        
        # Align edited exon region
        mafftAligner.run_nucleotide_as_protein(exon_FASTA_obj, findBestFrame=True)
        
        # Briefly fix up any MAFFT weirdness
        "I think my codon system kinda fucks up the last codon by missing 1-2 gaps at times"
        maxLen = max([len(FastASeq_obj.gap_seq) for FastASeq_obj in exon_FASTA_obj])
        for FastASeq_obj in exon_FASTA_obj:
            if len(FastASeq_obj.gap_seq) != maxLen:
                FastASeq_obj.gap_seq += "-"*(maxLen - len(FastASeq_obj.gap_seq))
        
        # Take note of how much we trimmed WITHIN the exon region
        codonsProblemLeft.append(startTrim)
        codonsProblemRight.append(endTrim)
        
    # Merge the edited exon regions back into the overall sequence
    for x in range(len(exonCoords)-1, -1, -1): # iterate backwards through coords since we know they're ordered ascendingly
        # Get values for this iteration
        exonStart, exonEnd = exonCoords[x]
        problemLeft = codonsProblemLeft[x]
        problemRight = codonsProblemRight[x]
        exon_FASTA_obj = exon_FASTA_objs[x]
        
        # First, get our Codons sequence for this exon region and insert it
        codons_FastASeq_obj = ZS_SeqIO.FastASeq("Codons", gapSeq="5"*problemLeft + "-"*len(exon_FASTA_obj[0].gap_seq) + "5"*problemRight)
        exon_FASTA_obj.insert(0, codons_FastASeq_obj)
        
        # Then, get our trimmed slices of the original sequence
        left_FASTA_obj = FASTA_obj.slice_cols(0, exonStart)
        right_FASTA_obj = FASTA_obj.slice_cols(exonEnd, len(FASTA_obj[0].gap_seq))
        
        # Write the exon and right FASTA objects to temporary files
        "We only do this because of how the FASTA.concat() method is implemented"
        exonTmpFileName = _tmp_file_name_gen("exon", "fasta")
        rightTmpFileName = _tmp_file_name_gen("right", "fasta")
        exon_FASTA_obj.write(exonTmpFileName, asAligned=True)
        right_FASTA_obj.write(rightTmpFileName, asAligned=True)
        
        # Merge the left, right, and exon FASTAs
        left_FASTA_obj.concat(exonTmpFileName)
        left_FASTA_obj.concat(rightTmpFileName)
        
        # Overwrite FASTA_obj with the new exon-edited version
        FASTA_obj = left_FASTA_obj
        
    return FASTA_obj # this has the same variable name, but its value is DIFFERENT than the input

def _get_suggested_fix_from_ssw_matches(matches, nuclSeq, score_cutoff, GOOD_ALIGN_PCT=0.80, GOOD_GAPS_NUM=2):
    '''
    Hidden function for use by polish_MSA_denovo(). It's been pulled aside here since it's quite
    complex and I don't want it being a mental burden when interpretting the parent function.

    Importantly, the fixes will never overlap.
    
    Parameters:
        GOOD_ALIGN_PCT -- arbitary, hard-coded magic number. I don't think we should change this.
        GOOD_GAPS_NUM -- arbitary, hard-coded magic number. 2 sounds right to me, yknow?
    Returns:
        fixes -- a list with format of: [[start, end, nLength], ...]. It should
                 used to make changes to the sequence, eventually.
    '''
    fixes = [] # format of 
    for problemAlign, targetAlign, startIndex, score in matches:
        # Limit ourselves to only good matches for indel fixing
        ## 1) Skip anything with a score that doesn't meet cut-off
        if score < score_cutoff:
            continue
        
        ## 2) Check if the alignment is good based on % overlap
        pctOverlap = len(problemAlign) / len(nuclSeq)
        if pctOverlap < GOOD_ALIGN_PCT:
            continue
        
        ## 3) Check if it's only 1 or 2 (max?) gap opens in either sequence
        problemGaps = list(re.finditer(r"-+", problemAlign))
        targetGaps = list(re.finditer(r"-+", targetAlign))
        numGaps = len(problemGaps) + len(targetGaps)
        if numGaps > GOOD_GAPS_NUM or numGaps == 0: # if we have no gaps, the stop codon has to be substitution related
            continue # we're not going to fix substitution errors, only indels
        
        # If we found a good match (i.e., we get here), see if there's anything to fix
        fix = []
        for gap in problemGaps:
            '''
            Here, the rationale is that gaps in the problem sequence mean it's MISSING
            something. Under that assumption, it's easy for us to just add N's into this
            gap region that maintain the reading frame, and everything should just work.
            '''
            gapLen = len(gap.group())
            gapFrame = gapLen % 3 # this will give 0, 1, or 2
            # If gapFrame isn't 0 (gapLen not divisible by 3), add N's to address the situation
            if gapFrame != 0:
                #gapPatch = "N"*(gapFrame) # gives a length of N's to insert
                start = gap.span()[0] + startIndex # + startIndex to give a consistent position in the sequence
                end = gap.span()[0] + gapFrame + startIndex # we're replacing only the length of "-"s that we want to be "n"s
                fix.append([start, end, gapFrame])
                # ongoingCount = 0
                # for j in range(gap.span()[0], gap.span()[1]): # iterate through gap positions
                #     if ongoingCount == gapFrame: # exit condition
                #         break
                #     newProblemAlign = newProblemAlign[0:j] + "N" + newProblemAlign[j+1:]
                #     ongoingCount += 1
        for gap in targetGaps:
            '''
            Here, the rationale is that gaps in the target sequence mean the problem sequence
            has ADDED something incorrectly. Under that assumption, it's non-trivial finding
            which position(s) have the problem especially if the gap region is longer than 3 base
            pairs. The safe solution is to REPLACE this entire region in the problem sequence
            with N's that fix the reading frame, and everything should still work.
            '''
            gapLen = len(gap.group())
            gapFrame = gapLen % 3 # this will give 0, 1, or 2
            # If gapFrame isn't 0 (gapLen not divisible by 3), use N's to address the situation
            if gapFrame != 0:
                #gapPatch = "N"*(gapLen-gapFrame) # gives a length of N's with appropriate shortening for good reading frame
                start = gap.span()[0] + startIndex # we're just going to replace the entire gap region
                end = gap.span()[1] + startIndex
                fix.append([start, end, gapLen-gapFrame])
                #newProblemAlign = newProblemAlign[0:gap.span()[0]] + gapPatch + newProblemAlign[gap.span()[1]:]
        fixes.append(fix)
    
    # Find the most consistently supported fix
    fixesCounts = {}
    ongoingCount = 0
    for fix in fixes:
        if ongoingCount == 0:
            fixesCounts[ongoingCount] = [fix, 1]
            ongoingCount += 1
        else:
            found = False
            for i in range(0, ongoingCount):
                if fixesCounts[i][0] == fix:
                    fixesCounts[i][1] += 1
                    found = True
                    break
            if found == False:
                fixesCounts[ongoingCount] = [fix, 1]
    
    bestFix = [None, 0]
    for value in fixesCounts.values():
        if value[1] > bestFix[1]:
            bestFix = value
    
    # Return the most supported fix
    return bestFix[0]

def _enact_fix_to_seq(fix, seq):
    '''
    Hidden function for use by polish_MSA_denovo(). Its goal is to take a fix identified by
    _get_suggested_fix_from_ssw_matches() and make the changes to the given sequence.
    
    NOTE: This is NOT the .gap_seq value, it is the .seq value!
    
    Return:
        editedGapSeq -- the fixed gapSeq.
    '''
    
    # Sort fixes in descending order
    fix.sort(key = lambda x: -x[1]) # doesn't matter if we sort by start or end since it's non-overlapping
    
    # Iterate through fixes and make all suggested changes
    for start, end, nLength in fix:
        seq = seq[0:start] + "N"*nLength + seq[end:]
    return seq

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

if __name__ == "__main__":
    usage = """%(prog)s receives a directory full of aligned FASTA files as part of the
    Oz Mammals genome project. Its goal is to remove indel errors from a MSA that has
    been previously subjected to exome_curation_2_introns.py. Specifically, it uses the
    information of which regions are coding to identify probable indel errors present in
    a minority of the sequences.
    
    Note: This should be step 3 in the Oz Mammals project!
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-a", dest="alignmentsDir", required=True,
                help="Specify the directory where aligned FASTA files are located")
    p.add_argument("-m", dest="mafftDir", required=True,
                help="Specify the directory where MAFFT executables are located")
    p.add_argument("-o", dest="outputDir", required=True,
                help="Output directory location (default == \"3_polish\")",
                default="3_polish")
    args = p.parse_args()
    validate_args(args)
    
    # Locate all aligned FASTA files
    files = [os.path.join(args.alignmentsDir, file) for file in os.listdir(args.alignmentsDir)]
    
    # Load aligned FASTA files
    fastaObjs = []
    for file in files:
        f = ZS_SeqIO.FASTA(file, isAligned=True)
        fastaObjs.append(f)
    
    # Polishing
    for i in range(len(files)):
        # Get details for this MSA
        alignFastaFile = files[i] # i=8 for first bug with fix.sort, fix is None somehow
        FASTA_obj = fastaObjs[i]
        
        # Perform polishing procedure
        FASTA_obj = polish_MSA_denovo(FASTA_obj, args.mafftDir)
        
        # Number codons
        ## TBD
        
        # Write output FASTA file
        outputFileName = os.path.join(args.outputDir, os.path.basename(alignFastaFile))
        FASTA_obj.write(outputFileName, withDescription=True, asAligned=True)
    
    print("Program completed successfully!")
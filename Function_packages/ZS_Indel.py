#! python3
# ZS_Indel.py
# Contains Class(es) to locate spurious
# indel errors in sequences.

import os, sys, re
import numpy as np
from copy import deepcopy

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from ZS_AlignIO import SSW
from ZS_ORF import _merge_coords_list, plummetdet

def plateaudet(lineChart, plummets):
    '''
    Another function to do things with a line chart (1D array) and try to find
    plateaus. Specifically, this works differently by taking the results of
    plummetdet() and using these as boundaries for potential plateau regions.
    This function then just finds the maximum within any potential region and
    returns the index of those positions for later plateau extension.
    
    Parameters:
        lineChart -- a list or 1D np.array() containing numeric values to search for
                     plateau regions in.
        plummets -- a list containing the indices where plummetdet() found plummets
                    within the lineChart value.
    Returns:
        maxindices -- a np.array() containing the indices for potential plateau maximums.
    '''
    lineChart = np.array(lineChart) # make sure it's always the same format
    
    # Get potential plateau regions
    potentialRegions = []
    for i in range(len(plummets)):
        thisValue = plummets[i]
        
        # Handle first
        if i == 0:
            potentialRegions.append([0, max(1, thisValue-1)])
        
        # Handle intermediate points
        else:
            start = min(plummets[i-1]+1, thisValue-1) # do this to handle adjacent plummets
            end = max(plummets[i-1]+1, thisValue-1)
            potentialRegions.append([start, end])
    
    # Handle last
    if plummets != []:
        potentialRegions.append([thisValue+1, len(lineChart)])
    
    # Get the maximum index for any potential regions
    maxindices = []
    for start, end in potentialRegions:
        maxindices.append(
            [
                np.where(np.array(lineChart[start:end]) == \
                    np.max(lineChart[start:end]))[0][0] + start,
                np.max(lineChart[start:end])
            ]
        )
    
    return np.array(maxindices)

class FastASeqAlignmentFrames:
    '''
    This Class represents the three-frame translation of a nucleotide FastASeq object.
    Specifically, it allows for mathematical operations to be performed based on the
    amount of alignment coverage it receives
    
    Note that we specifically handle only THREE frames. The sequence is assumed to be
    in the correct 5'->3' orientation already.
    
    Parameters:
        solutionDict -- a dictionary with structure like:
            {
                sequence_id: [seq, frame, hasStopCodon (bool), target_sequence_id]
            }
    '''
    def __init__(self, FastASeq_obj, solutionDict):
        assert type(FastASeq_obj).__name__ == "FastASeq" or type(FastASeq_obj).__name__ == "ZS_SeqIO.FastASeq"
        self.FastASeq = deepcopy(FastASeq_obj)
        self.FastASeq.seq = self.FastASeq.seq.strip("nN") # remove leading/trailing ambiguous characters
        self.solutionDict = solutionDict
        
        self.framing()
        self.frame_position_numbering()
    
    def framing(self):
        '''
        Performs the three frame translation.
        '''
        self.frame_1, _, _ = self.FastASeq.get_translation(strand=1, frame=0) # drop the strand and frame returns as '_'
        self.frame_2, _, _ = self.FastASeq.get_translation(strand=1, frame=1)
        self.frame_3, _, _ = self.FastASeq.get_translation(strand=1, frame=2)
    
    def frame_position_numbering(self):
        '''
        Numbers each position in each frame according to the amount of alignment coverage it receives.
        '''
        self.numbers_1 = np.zeros(len(self.FastASeq.seq))
        self.numbers_2 = np.zeros(len(self.FastASeq.seq))
        self.numbers_3 = np.zeros(len(self.FastASeq.seq))
        sequence = self.FastASeq.seq
        
        # Prevent bugs with empty frames
        if self.frame_1.upper().replace("X", "") == "" or \
            self.frame_2.upper().replace("X", "") == "" or \
            self.frame_3.upper().replace("X", "") == "":
                return
        
        # Run SSW against all solution sequences and tally numbers
        for x in range(0, 3):
            numbers = [self.numbers_1, self.numbers_2, self.numbers_3][x] # reuse code by selecting our numbers_# values here
            frame = [self.frame_1, self.frame_2, self.frame_3][x] # as above for frame_# values
            
            for value in self.solutionDict.values():
                solutionSeq, _, _ = value # drop translationFrame and hasStopCodon from value
                
                # Prevent bugs from empty solution sequences
                if solutionSeq.upper().replace("X", "") == "":
                    continue
                
                sswResult = SSW.ssw_parasail(frame, solutionSeq) # query, target
                
                # Find our nucleotide positions for the alignment
                queryStartIndex = sswResult.queryStartIndex * 3
                endIndex = queryStartIndex + (len(sswResult.queryAlign) * 3)
                
                # Iterate through our nucleotide sequence and note positions that are in the alignment
                ongoingCount = 0
                for nucleotide in sequence:                    
                    if nucleotide != "-":
                        if ongoingCount in range(queryStartIndex, endIndex + 1):
                            numbers[ongoingCount] += 1
                        ongoingCount += 1
    
    def max(self, position):
        '''
        Get the best ORF length value for the selected position by checking the underlying
        three frames.
        
        Parameters:
            position -- a zero-based integer for the position to check ORF lengths for.
        Returns:
            value -- an integer value for the best ORF length represented by this position.
        '''
        return max(self.numbers_1[position], self.numbers_2[position], self.numbers_3[position])

## Hidden methods for IndelPredictor
def _plateau_extens(lineChart, plateaus, coverages, plummets, delta, allowedDecrease=0.3):
    '''
    Parameters:
        delta -- the float value used for peakdet.
        plummets -- a list of indices corresponding to positions that should mark the
                    boundaries of any plateau extension
        allowedDecrease -- a float value indicating how much a min-max normalised
                            value can decrease before we no longer want to extend
                            the plateau. Works in tandem with plummet to limit plateau
                            extension to only sensible boundaries.
    '''
    for i in range(len(plateaus)):
        cutoff = coverages[i] - allowedDecrease
        
        # Look back
        prevCov = lineChart[plateaus[i][0]]
        newStart = plateaus[i][0]
        for x in range(plateaus[i][0]-1, -1, -1):
            indexCov = lineChart[x]
            # Plummet check
            if indexCov in plummets:
                break
            # Increasing check
            if indexCov > prevCov + delta:  # This means we're leading up to another peak, and should stop extending this plateau
                break
            # Decreasing cut-off
            if indexCov >= cutoff:
                newStart = x
                prevCov = indexCov
                continue
            else:
                break
        plateaus[i][0] = newStart
        
        # Look forward
        prevCov = lineChart[plateaus[i][1]]
        newEnd = plateaus[i][1]
        for x in range(plateaus[i][1]+1, len(lineChart)):
            indexCov = lineChart[x]
            # Plummet check
            if indexCov in plummets:
                break
            # Increasing check
            if indexCov > prevCov + delta:  # This means we're leading up to another peak, and should stop extending this plateau
                break
            # Decreasing cut-off
            if indexCov >= cutoff:
                newEnd = x
                prevCov = indexCov
                continue
            else:
                break
        plateaus[i][1] = newEnd
    return plateaus

def _filter_matches(matches, DESIRABLE_NUMBER=20):
    '''
    This function attempts to reduce the noise associated with fix finding from matches.
    It's important to configure DESIRABLE_NUMBER to be something that makes sense. It used
    to be 20 by default, but it seemed like this could introduce problems. I've upped the
    default, and also made it so filtering is an optional component of obtain_matches_via_ssw().
    '''
    # Find a good cut-off to use
    if len(matches) <= DESIRABLE_NUMBER:
         SCORE_CUTOFF = np.percentile([m.score for m in matches], 10) # drop the worst 10% only
    else:
        # Find a cut-off that ensures at least ~{DESIRABLE_NUMBER} matches
        for x in [90, 70, 50, 30, 20, 15]: # 15 will always be our fallback condition
            SCORE_CUTOFF = np.percentile([m.score for m in matches], x)
            numMatchesWithCutoff = sum([1 for m in matches if m.score >= SCORE_CUTOFF])
            if numMatchesWithCutoff >= DESIRABLE_NUMBER:
                break
    
    # Filter matches & return
    return [match for match in matches if match.score >= SCORE_CUTOFF]

class Match:
    '''
    Simple object to act as a container for match results
    '''
    def __init__(self, queryAlign, targetAlign, score, queryStartIndex, targetStartIndex, targetSeqID):
        self.queryAlign = queryAlign
        self.targetAlign = targetAlign
        self.score = score
        self.queryStartIndex = queryStartIndex
        self.targetStartIndex = targetStartIndex
        self.targetSeqID = targetSeqID
    
    def __repr__(self):
        return "<Match object; queryAlign={0}, targetAlign={1}, score={2}".format(
            self.queryAlign, self.targetAlign, self.score
        )

class IndelPredictor:
    '''
    This exists as a container for static methods that are useful for indel prediction.
    '''
    @staticmethod
    def predict_indel_regions(FastASeqAlignmentFrames_obj, delta=0.01, allowedDecrease=0.1):
        '''
        This method attempts to predict regions likely to contain indels within a sequence
        using a peak detection algorithm.
        
        In short, each position in the sequence is numbered based on the amount of times it
        aligned against another sequence. The maximum value for each position across the
        three frames is taken, and that gives us a kind of line chart we can perform peak
        detection from.
        
        Parameters:
            FastASeqAlignmentFrames_obj -- a FastASeqAlignmentFrames object.
            delta -- a float value indicating what change in the min-max normalised values must
                     occur before a new peak is detected
            allowedDecrease -- a float value indicating what change in the min-max normalised
                               values is permitted for extending a plateau. If this value is
                               less than delta, plateaus will be allowed to merge.
        Returns:
            indelRegions -- a list containing lists with structure like:
                              [
                                  [regionStart_1, regionEnd_1],
                                  [regionStart_2, regionEnd_2],
                                  ...
                              ]
        '''
        # Magic number declaration
        MIN_LENGTH_FOR_TRIM = 90 # this should do it...?
        TRIM_LENGTH = 20
        MIN_PLATEAU_LENGTH = 10
        BIAS_CUTOFF = 0.05 # ehhhh, it should work
        LOTS_OF_NS_PCT = 0.33 # ehh, idk?
        
        # Get our line chart
        seqLength = len(FastASeqAlignmentFrames_obj.FastASeq.seq)
        
        lineChart = np.zeros(seqLength)
        for position in range(seqLength):
            lineChart[position] = FastASeqAlignmentFrames_obj.max(position)
        
        # Early exit for empty line charts
        if np.sum(lineChart) == 0.0:
            return []
        
        # Min-max normalise line chart values
        lineChart = [(value - np.min(lineChart)) / (np.max(lineChart) - np.min(lineChart)) for value in lineChart]
        
        # De-noise data
        SMOOTH_RANGE = 10
        lineChart = [min(lineChart[max(0, i-SMOOTH_RANGE):min(i+SMOOTH_RANGE, len(lineChart))]) for i in range(len(lineChart))]
        
        # Trim any potentially biased data
        '''
        The goal here is to identify when our head/tail regions of a line chart have
        low values because there's N's in the nucleotide sequence resulting in poor
        alignment. These can be flagged as plateaus unfairly and result in problems
        where other behaviours are contingent on some level of confidence that an indel
        truly exists.
        '''
        ## Trim tail
        if len(lineChart) >= MIN_LENGTH_FOR_TRIM:
            tailNsCount = FastASeqAlignmentFrames_obj.FastASeq.seq[-TRIM_LENGTH:].lower().count("n")
            tailNsPct = tailNsCount / TRIM_LENGTH
            tailMedianValue = np.median(lineChart[-TRIM_LENGTH:])
            
            if tailNsPct > LOTS_OF_NS_PCT and tailMedianValue < BIAS_CUTOFF:
                trimAmount = sum([1 for value in lineChart[-TRIM_LENGTH:] if value < BIAS_CUTOFF])
                lineChart = lineChart[:-trimAmount]
        ## Trim head
        if len(lineChart) >= MIN_LENGTH_FOR_TRIM:
            headNsCount = FastASeqAlignmentFrames_obj.FastASeq.seq[:TRIM_LENGTH].lower().count("n")
            headNsPct = headNsCount / TRIM_LENGTH
            headMedianValue = np.median(lineChart[:TRIM_LENGTH])
            
            if headNsPct > LOTS_OF_NS_PCT and headMedianValue < BIAS_CUTOFF:
                trimAmount = sum([1 for value in lineChart[:TRIM_LENGTH] if value < BIAS_CUTOFF])
                lineChart = lineChart[:-trimAmount]
        
        # Perform plateau detection
        plummets = plummetdet(lineChart)
        maxindices = plateaudet(lineChart, plummets)
        
        # Get plateau regions
        plateaus = []
        coverages = []
        for maximum in maxindices:
            index = int(maximum[0])
            coverage = maximum[1]
            # Look forward
            '''The peakdet values are always at the start of the plateau,
            so we don't need to look back, we just need to look forward to find
            where the plateau ends'''
            plat = None
            for i in range(index, len(lineChart)):
                if lineChart[i] == coverage:
                    continue
                else:
                    plat = [index,i-1] # This is 0-indexed, and we -1 since we want the previous i value
                    break
            if plat == None: # This acts as a check for plateaus that run to the end of the sequence
                plat = plat = [index,i] # We don't -1 here since i will be equal to the last position of the sequence (in 0-based notation)
            plateaus.append(plat)
            coverages.append(coverage)
        
        # Extend plateaus
        plateaus = _plateau_extens(lineChart, plateaus, coverages, plummets, delta, allowedDecrease)
        
        # Merge overlapping plateaus
        plateaus = _merge_coords_list(plateaus)
        
        # Drop any short plateaus
        "The ordering here is intentional; drop comes AFTER merging"
        dropIndices = []
        for i in range(len(plateaus)):
            start, end = plateaus[i]
            if end - start < MIN_PLATEAU_LENGTH:
                dropIndices.append(i)
        for index in dropIndices[::-1]:
            del plateaus[index]
        
        # Find indel regions inbetween plateaus
        plateaus.sort()
        indelRegions = []
        for i in range(1, len(plateaus)):
            start = plateaus[i-1][1]+1
            end = plateaus[i][0]-1
            start, end = min(start, end), max(start, end) # accommodates single base positions
            indelRegions.append([start, end])
        
        # Return
        return indelRegions
    
    @staticmethod
    def obtain_matches_via_ssw(FASTA_obj, query_FastASeq_obj, solutionDict, FILTER=True):
        '''
        Function to take a sequence and, looking through the solutionDict,
        find the best SSW alignment for the sequence to all the solved
        sequences in solutionDict.
        
        The query_FastASeq_obj does NOT have to be part of FASTA_obj, but if it is, we'll
        skip self-comparison via comparing the .id attribute of the query_FastASeq_obj
        to whatever is in the solutionDict.
        
        Parameters:
            FASTA_obj -- a ZS_SeqIO.FASTA object
            query_FastASeq_obj -- a FastASeq object to use as a query and find matches for
            solutionDict -- a dictionary with structure like:
                {
                    sequence_id: [seq, frame, hasStopCodon (bool), target_sequence_id]
                }
            FILTER -- a boolean to indicate whether matches should be automatically filtered
        Returns:
            matches -- a list with structure like:
                [
                    [aligned_problem_seq, aligned_target_seq, start_index, score],
                    ...
                ]
        '''
        # assert query_FastASeq_obj.id in FASTA_obj.ids, \
        #     "{0} doesn't exist in the provided FASTA object".format(query_FastASeq_obj.id)
        
        queryNuclSeq = query_FastASeq_obj.seq
        
        matches = []
        for targetSeqID, value in solutionDict.items():
            if targetSeqID == query_FastASeq_obj.id: # skip potential self comparison
                continue
            _, _, targetHasStopCodon = value
            if targetHasStopCodon:
                continue
            else:
                targetNuclSeq = FASTA_obj[targetSeqID].seq
                sswResult = SSW.ssw_parasail(queryNuclSeq, targetNuclSeq)
                matches.append(
                    Match(
                        sswResult.queryAlign, sswResult.targetAlign,
                        sswResult.score, sswResult.queryStartIndex,
                        sswResult.targetStartIndex, targetSeqID
                    )
                )
        
        matches.sort(key = lambda x: -x.score) # order by score
        return matches if FILTER is False else _filter_matches(matches)
    
    @staticmethod
    def get_fix_from_ssw_matches(matches, querySeq, GOOD_ALIGN_PCT=0.60, PERFECT_VOTE_PCT=0.95, GOOD_GAPS_NUM=3, FIX_PCT=0.05):
        '''
        Following on from obtain_matches_via_ssw() for example, this function will interpret
        matches into a list of fixes that can occur to rectify probable indel errors.
        
        Parameters:
            matches -- a list with structure like:
                [
                    [aligned_problem_seq, aligned_target_seq, start_index, score],
                    ...
                ]
            querySeq -- a string of the sequence used as query to obtain matches to using e.g.,
                        obtain_matches_via_ssw().
            GOOD_ALIGN_PCT -- magic number. Used to determine whether a match is good enough
                              to look at for any fixes.
            PERFECT_VOTE_PCT -- magic number. Used to weigh against finding fixes if a sequence
                                is an almost "perfect" match and no fixes are suggested.
            GOOD_GAPS_NUM -- magic number. 2 USED TO sound right to me, but now I'm ride
                             or die with 3. Basically, a fix needs to have this many or fewer
                             suggested changes to be considered.
            FIX_PCT -- magic number. It lets us control weird scenarios where, out
                       of 50 matches, only 1 suggests a fix. If set to 0.05, 5%
                       of the matches would need to propose some kind of fix for us to continue
                       to find a fix. PERFECT_VOTE_PCT acts in opposition to this.
        Returns:
            fixes -- a list with format of: [[gapStart, resume, nLength], ...]. It should be
                    used to make changes to the sequence, eventually. Note that the gapStart
                    and resume variables are intentionally named to help conceptualisation.
                    
                    The gapStart value is where the gap region BEGINS and hence should be
                    EXCLUDED.
                    
                    The resume value is where we should RESUME the sequence from after 
                    any deletions or insertions i.e., it should be INCLUDED.
        '''
        fixes = []
        perfectVotes = 0
        for match in matches:
            # Limit ourselves to only good matches for indel fixing
            
            ## 1) Check if the alignment is good based on % overlap
            pctOverlap = len(match.queryAlign.replace("-","")) / len(querySeq)
            if pctOverlap < GOOD_ALIGN_PCT:
                continue
            
            ## 2) Check if it's only 1 or 2 (max?) gap opens in either sequence
            problemGapsHits = list(re.finditer(r"-+", match.queryAlign))
            problemGaps = [] # we're going to drop any gaps divisible by 3 because they shouldn't matter
            for gap in problemGapsHits:
                if len(gap.group()) % 3 != 0:
                    problemGaps.append(gap)
            targetGapsHits = list(re.finditer(r"-+", match.targetAlign))
            targetGaps = []
            for gap in targetGapsHits:
                if len(gap.group()) % 3 != 0:
                    targetGaps.append(gap)
            numGaps = len(problemGaps) + len(targetGaps)
            if numGaps > GOOD_GAPS_NUM or numGaps == 0: # if we have no gaps, the stop codon has to be substitution related
                if pctOverlap >= PERFECT_VOTE_PCT and numGaps == 0:
                    perfectVotes += 1
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
                    gapStart = gap.span()[0] + match.queryStartIndex # + startIndex to give a consistent position in the sequence
                    resume = gap.span()[0] + match.queryStartIndex
                    fix.append([gapStart, resume, gapFrame])
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
                    gapStart = gap.span()[0] + match.queryStartIndex # we're just going to replace the entire gap region
                    resume = gap.span()[1] + match.queryStartIndex # trust me, we ONLY USE queryStartIndex
                    fix.append([gapStart, resume, gapLen-gapFrame])
            fixes.append(fix)
        
        # Abort operation if we've found no valid fixes
        if fixes == []:
            return fixes
        
        # Also abort operation if most matches don't propose any fixes
        if (len(fixes) / len(matches)) < FIX_PCT:
            return []
        
        # Also abort operation if we've got lots of perfect matches
        if perfectVotes > len(fixes):
            return []
        
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
                    ongoingCount += 1
        
        bestFix = [None, 0]
        for value in fixesCounts.values():
            if value[1] > bestFix[1]:
                bestFix = value
        
        # Return the most supported fix
        return bestFix[0]
    
    @staticmethod
    def use_fix_on_sequence(fix, seq):
        '''
        The goal of this function is to take a fix identified by
        get_fix_from_ssw_matches() and make the changes to the given sequence.
        
        Changes are made using ambiguous characters i.e., "n"s. Note that these will be
        lowercase. It's probably a good idea to make the rest of the sequence uppercase, so
        it's obvious where these fixes have occurred.
        
        NOTE: Sequence should NOT be the .gap_seq value from a FastASeq object,
        it should be the .seq value!
        
        Parameters:
            fix -- a list with format of:
                [
                    [gapStart, resume, nLength],
                    ...
                ]
            seq -- a string of a nucleotide sequence.
        Return:
            editedSeq -- the fixed nucleotide sequence.
        '''
        
        # Sort fixes in descending order
        fix.sort(key = lambda x: -x[1]) # doesn't matter if we sort by start or end since it's non-overlapping
        
        # Iterate through fixes and make all suggested changes
        for start, end, nLength in fix:
            '''
            Start will be the first index where the gap (-) occurs, so running up to but 
            NOT INCLUDING start gives us the desirable behaviour. End is range() index
            i.e., non inclusive, so we want to INCLUDE it here.
            '''
            seq = seq[0:start] + "n"*nLength + seq[end:]
        return seq

if __name__ == "__main__":
    pass

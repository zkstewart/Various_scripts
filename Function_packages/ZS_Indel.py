#! python3
# ZS_Indel.py
# Contains Class(es) to locate spurious
# indel errors in sequences.

import os, sys, re
import numpy as np

sys.path.append(os.path.dirname(__file__))
from ZS_AlignIO import SSW
from ZS_ORF import peakdet, _merge_coords_list, plummetdet

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
        self.FastASeq = FastASeq_obj
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
        
        # Run SSW against all solution sequences and tally numbers
        for x in range(0, 3):
            numbers = [self.numbers_1, self.numbers_2, self.numbers_3][x] # reuse code by selecting our numbers_# values here
            frame = [self.frame_1, self.frame_2, self.frame_3][x] # as above for frame_# values
            
            for value in self.solutionDict.values():
                solutionSeq, _, _ = value # drop translationFrame and hasStopCodon from value
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

def _filter_matches(matches):
    # Find a good cut-off to use
    if len(matches) <= 20:
         SCORE_CUTOFF = np.percentile([x[3] for x in matches], 10) # drop the worst 10% only
    else:
        # Find a cut-off that ensures at least ~20 matches [this is an arbitrary cut-off]
        DESIRABLE_NUMBER = 20
        for x in [90, 70, 50, 30, 20, 15]: # 15 will always be out fallback condition
            SCORE_CUTOFF = np.percentile([x[3] for x in matches], x)
            numMatchesWithCutoff = sum([1 for m in matches if m[3] >= SCORE_CUTOFF])
            if numMatchesWithCutoff >= DESIRABLE_NUMBER:
                break
    
    # Filter matches & return
    return [match for match in matches if match[3] >= SCORE_CUTOFF]

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
        seqLength = len(FastASeqAlignmentFrames_obj.FastASeq.seq)
        
        lineChart = np.zeros(seqLength)
        for position in range(seqLength):
            lineChart[position] = FastASeqAlignmentFrames_obj.max(position)
        
        # Min-max normalise line chart values
        lineChart = [(value - np.min(lineChart)) / (np.max(lineChart) - np.min(lineChart)) for value in lineChart]
        
        # Perform peak detection & plummet detection
        plummets = plummetdet(lineChart)
        maxindices, minindices = peakdet(lineChart, delta)
        
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
        SHORT_CUTOFF = 30
        dropIndices = []
        for i in range(len(plateaus)):
            start, end = plateaus[i]
            if end - start < SHORT_CUTOFF:
                dropIndices.append(i)
        for index in dropIndices[::-1]:
            del plateaus[index]
        
        # Find sections inbetween plateaus
        plateaus.sort()
        indelRegions = []
        for i in range(1, len(plateaus)):
            start = plateaus[i-1][1]+1
            end = plateaus[i][0]-1
            indelRegions.append([start, end])
        
        # Return
        return indelRegions
    
    @staticmethod
    def obtain_matches_via_ssw(FASTA_obj, querySeqID, solutionDict, FILTER=True):
        '''
        Function to take a sequence and, looking through the solutionDict,
        find the best SSW alignment for the sequence to all the solved
        sequences in solutionDict.
        
        Parameters:
            FASTA_obj -- a ZS_SeqIO.FASTA object
            querySeqID -- a string identifying a sequence within FASTA_obj to have
                          matches found for
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
        assert querySeqID in FASTA_obj.ids, \
            "{0} doesn't exist in the provided FASTA object".format(querySeqID)
        
        matches = []
        for targetSeqID, value in solutionDict.items():
            if targetSeqID == querySeqID:
                continue
            _, _, targetHasStopCodon = value
            if targetHasStopCodon:
                continue
            else:
                queryNuclSeq = FASTA_obj[querySeqID].seq
                targetNuclSeq = FASTA_obj[targetSeqID].seq
                sswResult = SSW.ssw_parasail(queryNuclSeq, targetNuclSeq) # this function has poorly ordered outputs
                matches.append(
                    [
                        sswResult.queryAlign, sswResult.targetAlign,
                        sswResult.score, sswResult.queryStartIndex,
                        sswResult.targetStartIndex, targetSeqID
                    ]
                )
        
        matches.sort(key = lambda x: -x[3]) # order by score
        return matches if FILTER is False else _filter_matches(matches)
    
    @staticmethod
    def get_fix_from_ssw_matches(matches, querySeq, GOOD_ALIGN_PCT=0.60, GOOD_GAPS_NUM=2):
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
            GOOD_ALIGN_PCT -- arbitary, hard-coded magic number. I don't think we should change this.
            GOOD_GAPS_NUM -- arbitary, hard-coded magic number. 2 sounds right to me.
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
        for problemAlign, targetAlign, startIndex, _, _ in matches:
            # Limit ourselves to only good matches for indel fixing
            
            ## 1) Check if the alignment is good based on % overlap
            pctOverlap = len(problemAlign.replace("-","")) / len(querySeq)
            if pctOverlap < GOOD_ALIGN_PCT:
                continue
            
            ## 2) Check if it's only 1 or 2 (max?) gap opens in either sequence
            problemGapsHits = list(re.finditer(r"-+", problemAlign))
            problemGaps = [] # we're going to drop any gaps divisible by 3 because they shouldn't matter
            for gap in problemGapsHits:
                if len(gap.group()) % 3 != 0:
                    problemGaps.append(gap)
            targetGapsHits = list(re.finditer(r"-+", targetAlign))
            targetGaps = []
            for gap in targetGapsHits:
                if len(gap.group()) % 3 != 0:
                    targetGaps.append(gap)
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
                    gapStart = gap.span()[0] + startIndex # + startIndex to give a consistent position in the sequence
                    resume = gap.span()[0] + startIndex
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
                    gapStart = gap.span()[0] + startIndex # we're just going to replace the entire gap region
                    resume = gap.span()[1] + startIndex
                    fix.append([gapStart, resume, gapLen-gapFrame])
            fixes.append(fix)
        
        # Abort operation if we've found no valid fixes
        if fixes == []:
            return fixes
        
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

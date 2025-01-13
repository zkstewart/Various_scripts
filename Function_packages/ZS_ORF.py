#! python3
# ZS_ORF.py
# Contains Class(es) to manipulate ORFs e.g., through
# the prediction of ORFs for each sequence or from a
# MSA.

import os, sys, inspect, random, math, re
import numpy as np
from copy import deepcopy
from Bio import SeqIO

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from ZS_BlastIO import BLAST
from ZS_SeqIO import FASTA, FastASeq
from ZS_AlignIO import SSW

def peakdet(v, delta, x = None):
    import sys
    from numpy import NaN, Inf, arange, isscalar, asarray, array
    """
    Converted from MATLAB script at http://billauer.co.il/peakdet.html
    
    Returns two arrays
    
    function [maxtab, mintab]=peakdet(v, delta, x)
    %PEAKDET Detect peaks in a vector
    %                [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
    %                maxima and minima ("peaks") in the vector V.
    %                MAXTAB and MINTAB consists of two columns. Column 1
    %                contains indices in V, and column 2 the found values.
    %          
    %                With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
    %                in MAXTAB and MINTAB are replaced with the corresponding
    %                X-values.
    %
    %                A point is considered a maximum peak if it has the maximal
    %                value, and was preceded (to the left) by a value lower by
    %                DELTA.
    
    % Eli Billauer, 3.4.05 (Explicitly not copyrighted).
    % This function is released to the public domain; Any use is allowed.
    
    """
    maxtab = []
    mintab = []
        
    if x is None:
        x = arange(len(v))
    
    v = asarray(v)
    
    if len(v) != len(x):
        sys.exit('Input vectors v and x must have same length')
    
    if not isscalar(delta):
        sys.exit('Input argument delta must be a scalar')
    
    if delta <= 0:
        sys.exit('Input argument delta must be positive')
    
    mn, mx = Inf, -Inf
    mnpos, mxpos = NaN, NaN
    
    lookformax = True
    
    for i in arange(len(v)):
        this = v[i]
        if this > mx:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]
        
        if lookformax:
            if this < mx-delta:
                maxtab.append((mxpos, mx))
                mn = this
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn+delta:
                mintab.append((mnpos, mn))
                mx = this
                mxpos = x[i]
                lookformax = True
    
    return array(maxtab), array(mintab)

def maxdet(v):
    """
    Complement to peakdet which will address problems where no peak was found due to
    how delta affects peakdet. Will just return the first position in any maximums.
    """
    maxindices = np.where(v == np.max(v))[0]
    outindices = [(maxindices[0], 1.0)] # make sure we always at least return the first maximum
    return np.array(
        outindices + [(maxindices[i], 1.0) for i in range(1, len(maxindices)) if (maxindices[i] > maxindices[i-1]+1)]
    )

def plateau_start_det(v, minimum=0.9):
    """
    Alternative to peakdet which will address problems where peaks are not found for reasons
    that have not been fully understood. This function will return the start of regions that should
    later become plateaus.
    
    Parameters:
        v -- a list or np.array() containing float values in the range of 0.0->1.0
        minimum -- a float value indicating the minimum value that a plateau should have
    """
    # Find where the values are greater than the minimum
    startindices = np.where(v >= np.float64(minimum))[0]
    startindices = [
        [si, v[si]]
        for si in startindices
        if si == 0 or v[si-1] < minimum
    ]
    
    # Push values back to where they first become 1 if possible
    outindices = []
    for si, sv in startindices:
        if sv == 1.0:
            outindices.append([si, sv])
        else:
            _si, _sv = si, sv
            for i in range(si, len(v)):
                if v[i] < minimum:
                    break
                if v[i] > _sv:
                    _si, _sv = i, v[i]
                    if v[i] == 1.0:
                        break
            outindices.append([_si, _sv])
    
    return np.array(outindices)

def outliers_z_score(ys, threshold=3):
    """
    Credit to http://colingorrie.github.io/outlier-detection.html
    """
    mean_y = np.mean(ys)
    stdev_y = np.std(ys)
    z_scores = [(y - mean_y) / stdev_y for y in ys]
    return np.where(np.abs(z_scores) > threshold)

def plummetdet(v):
    """
    Designed to behave kinda similarly to peakdet but instead, we want to find unusually large
    declines in min-max normalised values.
    
    Parameters:
        v -- a list or np.array() containing float values in the range of 0.0->1.0
    Returns:
        plummets -- a list of integers indicating the index of any positions that plummet
    """
    declines = []
    for i in range(1, len(v)):
        this = v[i]
        prev = v[i-1]
        declines.append(prev - this)
    
    declines = np.array(declines)
    
    # Calculate plummets using z-score
    if np.all(declines == 0): # avoid a runtime warning
        return []
    else:
        plummets = outliers_z_score(declines)
        return list(plummets[0])

class FastASeqFrames:
    '''
    This Class represents the three-frame translation of a nucleotide FastASeq object.
    Specifically, it allows for mathematical operations to be performed based on the
    length of any ORFs found across its three frames.
    
    Note that we specifically handle only THREE frames. The sequence is assumed to be
    in the correct 5'->3' orientation already.
    '''
    def __init__(self, FastASeq_obj):
        assert type(FastASeq_obj).__name__ == "FastASeq" or type(FastASeq_obj).__name__ == "ZS_SeqIO.FastASeq"
        
        self.framing(FastASeq_obj) # sets self.seq, self.frame_1, self.frame_2, self.frame_3
        self.frame_position_numbering()
    
    def framing(self, FastASeq_obj):
        '''
        Performs the three frame translation;
        drops the strand and frame returns as '_' internally.
        '''
        startRegex = re.compile(r"^[nN-]+")
        endRegex = re.compile(r"[nN-]+$")
        
        # Identify where the sequence is concealed by Ns or gaps
        thisSeq = FastASeq_obj.gap_seq
        
        startConcealed = startRegex.search(thisSeq)
        endConcealed = endRegex.search(thisSeq)
        
        startConcealed = len(startConcealed.group()) if startConcealed != None else 0
        endConcealed = len(endConcealed.group()) if endConcealed != None else 0
        
        # Generate a new sequence with the concealed regions replaced by Ns
        self.seq = "N"*startConcealed + thisSeq[0+startConcealed:len(thisSeq)-endConcealed] + "N"*endConcealed
        thisSeq = self.seq.replace("-", "") # get rid of gaps so we can use it for translation next
        
        # Translate the sequence
        self.frame_1 = FastASeq.dna_to_protein(thisSeq)
        self.frame_2 = FastASeq.dna_to_protein(thisSeq[1:])
        self.frame_3 = FastASeq.dna_to_protein(thisSeq[2:])
        
        # Note where the sequence is concealed
        self.concealed = np.array( [1]*startConcealed + [0]*(len(self.seq) - startConcealed - endConcealed) + [1]*endConcealed )
    
    def frame_position_numbering(self):
        '''
        Numbers each position in each frame according to the length of ORF it is participating in.
        '''
        self.numbers_1 = np.zeros(len(self.seq))
        self.numbers_2 = np.zeros(len(self.seq))
        self.numbers_3 = np.zeros(len(self.seq))
        
        for x in range(0, 3):
            numbers = [self.numbers_1, self.numbers_2, self.numbers_3][x] # reuse code by selecting our numbers_# values here
            frame = [self.frame_1, self.frame_2, self.frame_3][x] # as above for frame_# values
            
            # Fill in first numbers values cut off by this frame
            splitFrame = frame.split("*")
            numbers[0:x] = len(splitFrame[0])*3 # *3 for nucleotide length
            
            # Fill in the rest of the numbers values based on the ORF lengths
            ongoingCount = x
            for y in range(len(splitFrame)):
                frameSeq = splitFrame[y]
                orfLength = len(frameSeq)*3 # *3 for nucleotide length
                
                # Iterate through the ORF sequence
                for amino in frameSeq:
                    aminoCount = 0
                    for nucleotide in self.seq[ongoingCount:]:
                        if aminoCount == 3 or ongoingCount == len(numbers):
                            break
                         
                        numbers[ongoingCount] = orfLength
                        if nucleotide != "-":
                            aminoCount += 1
                        ongoingCount += 1
                
                # Handle a stop codon at the end if we didn't reach len(numbers)
                if ongoingCount != len(numbers):
                    stopCodonCount = 0
                    for nucleotide in self.seq[ongoingCount:]:
                        if stopCodonCount == 3 or ongoingCount == len(numbers):
                            break
                        
                        numbers[ongoingCount] = orfLength
                        if nucleotide != "-":
                            stopCodonCount += 1
                        ongoingCount += 1
        
        self.numbers_max = np.max((self.numbers_1, self.numbers_2, self.numbers_3))
    
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

def _merge_coords_list(coords):
    '''
    Helper function for merging overlapping coords lists
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

class MSA_ORF:
    '''
    TBD; unsure if this class will actually be a thing so no point documenting it
    heavily just yet.
    
    This Class will ONLY receive a ZS_SeqIO.FASTA object. It MUST contain nucleotide
    sequences - what's the point of looking for ORFs in a protein file?
    
    It will NOT modify the input FASTA. A deepcopy is made of the provided object.
    
    The FASTA object must already be aligned! This is a MSA-based approach to ORF
    finding.
    '''
    def __init__(self, FASTA_obj):
        # Validation
        assert type(FASTA_obj).__name__ == "FASTA" or type(FASTA_obj).__name__ == "ZS_SeqIO.FASTA"
        
        if not FASTA_obj.isAligned:
            raise Exception("FASTA object isn't flagged as being aligned; cant use it with MSA_ORF")
        
        for FastASeq_obj in FASTA_obj:
            if FastASeq_obj.gap_seq == None:
                raise Exception(inspect.cleandoc("""
                                Sequence with ID {0} lacks a gap seq value; 
                                can't use it with MSA_ORF""".format(FastASeq_obj.id)))
        
        prevLen = None
        for FastASeq_obj in FASTA_obj:
            if prevLen == None:
                pass
            else:
                if prevLen != len(FastASeq_obj.gap_seq):
                    raise Exception("Not all FASTA sequences are identical in length; can't use it with MSA_ORF")
            prevLen = len(FastASeq_obj.gap_seq)
        
        # Store relevant values
        self.FASTA = deepcopy(FASTA_obj)
        self.peaksCoordinates = None
        self.minimalStopsCoordinates = None
    
    def _plateau_extens(self, plateaus, coverages, plummets, delta, allowedDecrease=0.3):
        '''
        Function to help mitigate the effect of internal stop codons messing with an
        ORF region because one or two sequences had one which affected the ORF
        length calculation.
        
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
            prevCov = self.lineChart[plateaus[i][0]]
            newStart = plateaus[i][0]
            for x in range(plateaus[i][0]-1, -1, -1):
                indexCov = self.lineChart[x]
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
            prevCov = self.lineChart[plateaus[i][1]]
            newEnd = plateaus[i][1]
            for x in range(plateaus[i][1]+1, len(self.lineChart)):
                indexCov = self.lineChart[x]
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
    
    def denovo_prediction_peaks(self, delta=0.005, allowedDecrease=0.2, PEAK_MINIMUM=0.90, CONCEALED_COUNT_RATIO=0.5):
        '''
        This method attempts to predict multiple ORFs from within the MSA using a peak 
        detection algorithm.
        
        In short, each position in each sequence is numbered based on the maximum length of ORF
        that position participates in across the 3 frames. Those values are summed per position
        and that gives us a kind of line chart we can perform peak detection from.
        
        Parameters:
            delta -- a float value indicating what change in the min-max normalised values must
                     occur before a new peak is detected
            allowedDecrease -- a float value indicating what change in the min-max normalised
                               values is permitted for extending a plateau. If this value is
                               less than delta, plateaus will be allowed to merge.
            PEAK_MINIMUM -- a float value in the range of 0->1 which will limit what regions
                            of the min-max array can be detected as a peak. This will prevent
                            low-quality peaks from being found.
            CONCEALED_COUNT_RATIO -- a float value in the range of 0->1 which will determine whether
                                     a sequence that is concealed at its start or end (i.e., has gaps or
                                     N's exclusively) will have the benefit of the doubt applied when it
                                     comes to contributing to the line chart. Specifically, if the ratio
                                     of non-concealed sequences has a full length ORF at a position, then
                                     the concealed sequence will be treated as if it has a full length ORF
                                     at that position.
        Returns:
            orfCoordinates -- a list containing lists with structure like:
                              [
                                  [orfStart_1, orfEnd_1],
                                  [orfStart_2, orfEnd_2],
                                  ...
                              ]
        '''
        # Init values
        msaLength = len(self.FASTA[0].gap_seq) # all seqs should be same length due to __init__ validations
        frames = [FastASeqFrames(FastASeq_obj) for FastASeq_obj in self.FASTA]
        
        # Generate data structure for peak detection
        self.lineChart = np.zeros(msaLength)
        for position in range(msaLength):
            # Determine whether we should treat concealed sequences as full length
            isFullLength = [
                1 if 
                  frame.concealed[position] == 0 # concealed sequences are ambiguous if they're full length
                  and frame.max(position) == frame.numbers_max
                else 0
                for frame in frames
                if 0 in frame.concealed # skip sequences that are all concealed
            ]
            isFullLengthRatio = sum(isFullLength) / len(isFullLength)
            
            # Generate line chart values based on whether we should treat concealed sequences as full length
            if isFullLengthRatio >= CONCEALED_COUNT_RATIO:
                self.lineChart[position] = sum([frame.max(position) for frame in frames])
            else:
                self.lineChart[position] = sum([frame.max(position) for frame in frames if frame.concealed[position] == 0])
        
        # End detection if no peaks are possible
        "i.e., if the line chart is just a flat line"
        if np.min(self.lineChart) == np.max(self.lineChart):
            self.peaksCoordinates = [[0, msaLength]]
            return self.peaksCoordinates
        
        # End detection if trough is really minor
        "i.e., if the difference between the min and max value is small and being contributed by a minority of sequences"
        if (np.min(self.lineChart) + np.max(self.lineChart)*0.01) >= np.max(self.lineChart):
            self.peaksCoordinates = [[0, msaLength]]
            return self.peaksCoordinates
        
        # Min-max normalise line chart values
        self.lineChart = [(value - np.min(self.lineChart)) / (np.max(self.lineChart) - np.min(self.lineChart)) for value in self.lineChart]
        
        # Perform peak detection & plummet detection
        plummets = plummetdet(self.lineChart) # plummets are actually any outlier sudden increase or decrease
        maxindices = plateau_start_det(self.lineChart, PEAK_MINIMUM)
        # maxindices, minindices = peakdet(self.lineChart, delta)
        # if len(maxindices) == 0:
        #     maxindices = maxdet(self.lineChart)
        # else:
        #     maxindices = maxindices[np.where(np.abs(maxindices[:,1]) > PEAK_MINIMUM)] # remove crappy peaks
        
        # Get plateau regions
        plateaus = []
        coverages = []
        for maximum in maxindices:
            index = int(maximum[0])
            coverage = maximum[1]
            # Look forward
            '''The peakdet/plateau_start_det values are always at the start of the plateau,
            so we don't need to look back, we just need to look forward to find
            where the plateau ends'''
            plat = None
            for i in range(index, len(self.lineChart)):
                if self.lineChart[i] == coverage:
                    continue
                else:
                    plat = [index,i-1] # This is 0-indexed, and we -1 since we want the previous i value
                    break
            if plat == None: # This acts as a check for plateaus that run to the end of the sequence
                plat = [index,i] # We don't -1 here since i will be equal to the last position of the sequence (in 0-based notation)
            plateaus.append(plat)
            coverages.append(coverage)
        
        # Extend plateaus to deal with rare stop codons
        "Think of a stop codon in just one or two sequences as a small barrier to our finding an ideal ORF region"
        plateaus = self._plateau_extens(plateaus, coverages, plummets, delta, allowedDecrease)
        
        # Drop any short plateaus
        SHORT_CUTOFF = 30
        dropIndices = []
        for i in range(len(plateaus)):
            start, end = plateaus[i]
            if end - start < SHORT_CUTOFF:
                dropIndices.append(i)
        for index in dropIndices[::-1]:
            del plateaus[index]
        
        # Merge overlapping plateaus
        plateaus = _merge_coords_list(plateaus)
        
        # Return
        orfCoordinates = plateaus # just for conceptual purposes and matching our method string
        self.peaksCoordinates = orfCoordinates
        return orfCoordinates
    
    def denovo_prediction_minimalStops(self, transcriptomeFile, EXCLUSION_PCT=0.90):
        '''
        Predicts the best ORF regions across all sequences using combined evidence.
        It does this without genomic evidence by assessment of how to get the
        longest uninterrupted (i.e., no stop codons) region from a MSA.
        
        Due to how this works, it will find only one ORF per MSA.
        
        Params:
            transcriptomeFile -- a string indicating the location of a transcriptome file.
                                This will be used to identify the best translation frame for
                                sequences where this isn't readily apparent.
            EXCLUSION_PCT -- a float value indicating the proportion of sequences that
                            must have their stop codons excluded within the predicted CDS region.
                            E.g., if 0.90, then 90% of the sequences must NOT contain a stop
                            codon within the region.
        Returns:
            orfCoordinates -- a list containing lists with structure like:
                              [
                                  [orfStart_1, orfEnd_1],
                                  [orfStart_2, orfEnd_2],
                                  ...
                              ]
        '''
        assert isinstance(EXCLUSION_PCT, float) or isinstance(EXCLUSION_PCT, int)
        assert 0 <= EXCLUSION_PCT <= 1, "EXCLUSION_PCT must be between 0 or 1 (inclusive of 0 and 1)"
        
        EXCLUSION_PCT = int(EXCLUSION_PCT*100)
        
        # Get sequence translations
        solutionDict = ORF.solve_translation_frames(self.FASTA, transcriptomeFile)
        
        # Locate longest segment boundaries for each sequence that excludes stop codons
        boundaries = MSA_ORF.get_segment_boundaries(self.FASTA, solutionDict)
        
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
        
        # Return
        orfCoordinates = [[trueStartIndex, trueEndIndex]] # just for conceptual purposes and matching our method string
        self.minimalStopsCoordinates = orfCoordinates
        return orfCoordinates

    @staticmethod
    def get_segment_boundaries(FASTA_obj, solutionDict, NATURAL_PCT=0.25):
        '''
        Method for getting the bounded region in which no stop codons exist
        for each sequence in the FASTA_obj, depending on the results of 
        solve_translation_frames() [i.e., solutionDict].
        
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
                
                if ongoingPositionCount+1 >= proteinEnd*3: # +1 to handle sequences that end early and have all gap at end
                    "The >= is important above for various reasons I've come to semi-understand"
                    endIndex = x+1 # proteinEnd marks the first position of the stop codon, setting this to endIndex will exclude it
                
                if letter == "-":
                    "Must skip gaps AFTER checking for the end of the sequence in case"
                    "the sequence ends in a gap"
                    continue
                
                if ongoingPositionCount == proteinStart*3:
                    "OK to check start position after skipping gaps since the start"
                    "position will never be a gap"
                    startIndex = x
                
                ongoingPositionCount += 1
            
            if endIndex == None:
                assert letter != "-" # check that the last letter isn't a gap
                endIndex = x+1 # +1 to include the last position when using range() type logic
            
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

class ORF:
    '''
    Provides static methods for the prediction of ORFs for each sequence in a FASTA object.
    '''
    @staticmethod
    def _simple_scenario_handler(FASTA_obj, resultsDict):
        '''
        Helper function for solve_translation_frames()
        '''
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
    
    @staticmethod
    def _advanced_scenario_handler(FASTA_obj, resultsDict, transcriptomeFile, DESIRABLE_NUMBER=50):
        '''
        Helper function for solve_translation_frames()
        
        Finds solutions in a more advanced way via BLAST. It will randomly iterate through
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
                    if len(seq.upper().replace("X","")) < SHORT_SEQ_LEN:
                        blastResults.append(math.inf)
                        continue
                    # Set up our BLAST handler
                    blaster = BLAST(FastASeq("eg", seq=seq), transcriptomeFile, "blastp")
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
    
    @staticmethod
    def _reset_memory_dict_index(memoryDict, thisSeqID):
        '''
        Helper function for solve_translation_frames()
        '''
        memoryDict[thisSeqID] = {} # forget the parent index
        indicesToForget = []
        for key, frameDict in memoryDict.items():
            for frame, scoreDict in frameDict.items():
                for otherSeqID, score in scoreDict.items():
                    if otherSeqID == thisSeqID:
                        indicesToForget.append([key, frame, otherSeqID])
        for key, frame, otherSeqID in indicesToForget:
            del memoryDict[key][frame][otherSeqID] # forget the results for other parents to this index

    @staticmethod
    def _calculate_solution_scores(resultsDict, solutionDict, memoryDict, thisSeqID):
        '''
        Helper function for solve_translation_frames()
        '''
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
                # Otherwise, do it anew
                except:
                    solutionSeq, solutionFrame, _ = solution
                    sswResult = SSW.ssw_parasail(solutionSeq, targetSeq, "protein")
                    score = sswResult.score
                    # Store it into memoryDict
                    memoryDict.setdefault(thisSeqID, {})
                    memoryDict[thisSeqID].setdefault(j, {})
                    memoryDict[thisSeqID][j].setdefault(otherSeqID, {})
                    memoryDict[thisSeqID][j][otherSeqID][solutionFrame] = score
                scores[j].append(score)
        
        return scores

    @staticmethod
    def _use_scores_metrics_to_get_best_frame(scores, results):
        '''
        Helper function for solve_translation_frames()
        '''
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
    
    @staticmethod
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
            transcriptomeFile -- a string indicating the location of a transcriptome
                                file which we can BLAST against when solving hard
                                scenarios.
        Returns:
            solutionDict -- a dictionary with structure like:
                {
                    sequence_id: [seq, frame, hasStopCodon (bool), target_sequence_id]
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
        solutionDict, problemDict = ORF._simple_scenario_handler(FASTA_obj, resultsDict)
        
        # If we didn't find many easy-to-solve scenarios with simple handling, get more advanced
        DESIRABLE_NUMBER=5 # arbitrary, should work for the oz mammals exomes project
        if len(solutionDict) <= DESIRABLE_NUMBER:
            solutionDict, problemDict = ORF._advanced_scenario_handler(FASTA_obj, resultsDict, transcriptomeFile)
        
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
                scores = ORF._calculate_solution_scores(resultsDict, solutionDict, memoryDict, seqID) # will store scores for frame 0, 1, and 2
                
                # Calculate metrics to find the best frame
                bestFrame = ORF._use_scores_metrics_to_get_best_frame(scores, resultsDict[seqID])
                
                # Update the solutionDict value if applicable
                if bestFrame != solutionDict[seqID][1]:
                    solutionDict[seqID] = resultsDict[seqID][bestFrame] # we do this since the sequence has changed
                    ORF._reset_memory_dict_index(memoryDict, seqID) # and hence any scores associated to it are irrelevant
        
        # Find the best solution to problem sequences
        memoryDict = {} # reset memory dict since solution->solution comparisons don't happen below
        if solutionDict != {}: # prevent error when using an empty solutionDict
            for seqID in problemDict.keys():
                scores = ORF._calculate_solution_scores(resultsDict, solutionDict, memoryDict, seqID) # will store scores for frame 0, 1, and 2
                
                # Calculate metrics to find the best frame
                bestFrame = ORF._use_scores_metrics_to_get_best_frame(scores, resultsDict[seqID])
                
                # Add the problem sequence into our solutionDict
                solutionDict[seqID] = resultsDict[seqID][bestFrame]
        
        ORF.solutionDict = solutionDict
        return solutionDict

class ORF_Find:
    '''
    Copy-pasted reimplementation of the biopython_orf_find.py script's logic with
    some minor modifications to remove unneeded functionality.
    '''
    def __init__(self, fastaFile):
        '''
        Params:
            fastaFile -- a string indicating the location of a FASTA file.
        '''
        self.fastaFile = fastaFile
        self.startCodon = re.compile(r'^.*?(M.*)')
        self.xRegex = re.compile(r'X+')
        
        # Set default values
        self._minProLen = 30
        self._maxProLen = 0
        self._hitsToPull = 3
        self._altCodonStringency = 49
        self._noCodonStringency = 99
        self._sequenceType = "prot"
        self._unresolvedCodon = 0
        self._translationTable = 1
    
    @property
    def minProLen(self):
        return self._minProLen
    
    @minProLen.setter
    def minProLen(self, value):
        assert isinstance(value, int)
        assert value > 0, "minProLen must be greater than 0"
        self._minProLen = value
    
    @property
    def maxProLen(self):
        return self._maxProLen
    
    @maxProLen.setter
    def maxProLen(self, value):
        assert isinstance(value, int)
        assert value >= 0, "maxProLen must be greater than or equal to 0"
        self._maxProLen = value
    
    @property
    def hitsToPull(self):
        return self._hitsToPull
    
    @hitsToPull.setter
    def hitsToPull(self, value):
        assert isinstance(value, int)
        assert value > 0, "hitsToPull must be greater than 0"
        self._hitsToPull = value
    
    @property
    def altCodonStringency(self):
        return self._altCodonStringency
    
    @altCodonStringency.setter
    def altCodonStringency(self, value):
        assert isinstance(value, int)
        assert value >= 0, "altCodonStringency must be greater than or equal to 0"
        self._altCodonStringency = value
    
    @property
    def noCodonStringency(self):
        return self._noCodonStringency
    
    @noCodonStringency.setter
    def noCodonStringency(self, value):
        assert isinstance(value, int)
        assert value >= 0, "noCodonStringency must be greater than or equal to 0"
        self._noCodonStringency = value
    
    @property
    def unresolvedCodon(self):
        return self._unresolvedCodon
    
    @unresolvedCodon.setter
    def unresolvedCodon(self, value):
        assert isinstance(value, int)
        assert value >= 0, "unresolvedCodon must be greater than or equal to 0"
        self._unresolvedCodon = value
    
    @property
    def translationTable(self):
        return self._translationTable
    
    @translationTable.setter
    def translationTable(self, value):
        assert isinstance(value, int)
        assert value > 0, "translationTable must be greater than 0"
        self._translationTable = value
    
    def process(self):
        records = SeqIO.parse(open(self.fastaFile, 'r'), 'fasta')
        
        # Iterate over transcript records
        for record in records:
            # Declare output holding values that should reset for each transcript/record
            tempMProt = []
            tempMNucl = []
            tempAltProt = []
            tempAltNucl = []
            tempNoneProt = []
            tempNoneNucl = []
            
            # Iterate over strands
            for strand, nuc in [(+1, record.seq), (-1, record.seq.reverse_complement())]:
                # Iterate over reading frames
                for frame in range(3):
                    length = 3 * ((len(record)-frame) // 3)
                    frameNuc = str(nuc[frame:frame+length])
                    frameProt = str(nuc[frame:frame+length].translate(table=self.translationTable))
                    
                    # Split protein/nucleotide into corresponding ORFs
                    ongoingLength = 0
                    splitNucleotide = []
                    splitProtein = []
                    frameProt = frameProt.split('*')
                    for i in range(len(frameProt)):
                        if len(frameProt) == 1 or i + 1 == len(frameProt):
                            splitProtein.append(frameProt[i])
                            splitNucleotide.append(frameNuc[ongoingLength:ongoingLength+len(frameProt[i])*3])
                            ongoingLength += len(frameProt[i])*3
                        else:
                            splitProtein.append(frameProt[i] + '*')
                            splitNucleotide.append(frameNuc[ongoingLength:ongoingLength+len(frameProt[i] + '*')*3])       
                            ongoingLength += (len(frameProt[i]) + 1)*3
                    
                    # Fix unresolved regions
                    resolvedProt = []
                    resolvedNuc = []
                    indicesForDel = []
                    for i in range(len(splitProtein)):
                        if 'X' in splitProtein[i]:
                            posProt = []
                            for x in re.finditer(self.xRegex, splitProtein[i]):
                                if x.end() - x.start() > self.unresolvedCodon:
                                    posProt += [x.start(), x.end()]
                            if posProt == []:
                                continue
                            indicesForDel.insert(0, i)
                            
                            # Pull out resolved regions
                            resolvedProt.append(splitProtein[i][:posProt[0]])
                            resolvedNuc.append(splitNucleotide[i][:posProt[0]*3])
                            for x in range(1, len(posProt)-1, 2):
                                start = posProt[x]
                                end = posProt[x+1]
                                resolvedProt.append(splitProtein[i][start:end])
                                resolvedNuc.append(splitNucleotide[i][start*3:end*3])
                            resolvedProt.append(splitProtein[i][posProt[-1]:])
                            resolvedNuc.append(splitNucleotide[i][posProt[-1]*3:])
                    
                    # Delete old entries and add resolved entries
                    for index in indicesForDel:
                        del splitProtein[index]
                        del splitNucleotide[index]
                    splitProtein += resolvedProt
                    splitNucleotide += resolvedNuc
                    
                    # Enter the main processing loop with our resolved regions
                    for i in range(len(splitProtein)):                              # Note that I have done a 'for i in range...' loop rather than a 'for value in splitProtein' loop which would have been simpler for a reason explained below on the 'elif i + 1 ==' line
                        # Declare blank values needed for each potential ORF region so we can tell which things were 'found'
                        mPro = None
                        altPro = None
                        nonePro = None
                        codonIndex = None
                        noneCodonContingency = None
                        
                        # Process sequences to determine whether we're ignoring this, or adding an asterisk for length counts
                        if len(splitProtein[i]) < self.minProLen:            # Disregard sequences that won't meet the size requirement
                            continue
                        acceptedPro = str(splitProtein[i])
                        
                        # Alternative start coding      
                        nucSeqOfProt = splitNucleotide[i]               # Don't need to do it, but old version of script extensively uses this value and cbf changing it
                        codons = re.findall('..?.?', nucSeqOfProt)          # Pulls out a list of codons from the nucleotide
                        for codon in codons:                    # Cycle through this list of codons to find the first alternative start of the normal class (GTG and TTG) and the rare class (CTG)
                            if codon == 'GTG' or codon == 'TTG':
                                codonIndex = codons.index(codon)    # This will save the position of the first GTG or TTG encountered. Note that by breaking after this,  we stop looking for CTG as it is irrelevant after this
                                break
                            elif codon == 'CTG':
                                if noneCodonContingency == None:    # noneCodonContingency is set to None at the end of each loop. Thus, this line of code will 'capture' the position of the first CTG in a sequence if a GTG or TTG was not encountered first
                                    noneCodonContingency = codons.index(codon)
                        
                        # Get the three ORF versions from each region inbetween stop codons
                        if 'M' in str(acceptedPro):                 # Obtains a traditional methionine initiated ORF starting from the first methionine if there is one in the sequence
                            mPro = self.startCodon.search(str(acceptedPro)).groups()[0]  # Note that startCodon was declared at the start of this file     
                        
                        if codonIndex != None:                  # Gets the start position of the protein if we found a likely alternative start (aka a 'GTG' or 'TTG')
                            altPro = acceptedPro[codonIndex:]
                        elif noneCodonContingency != None:              # This will match an alternative start to 'CTG' only if 'TTG' or 'GTG' are not present
                            altPro = acceptedPro[codonIndex:]
                        
                        nonePro = acceptedPro
                        
                        # Store if sequence lengths are within the desired range
                        if mPro != None:
                            if (len(mPro) >= self.minProLen) and (self.maxProLen == 0 or len(mPro) <= self.maxProLen):
                                tempMProt.append(mPro)
                                newStartPosition = acceptedPro.find(mPro)
                                tempMNucl.append(str(nucSeqOfProt[newStartPosition*3:]))
                        if altPro != None:
                            if (len(altPro) >= self.minProLen) and (self.maxProLen == 0 or len(altPro) <= self.maxProLen):
                                tempAltProt.append(altPro)
                                newStartPosition = acceptedPro.find(altPro[1:]) - 1
                                tempAltNucl.append(str(nucSeqOfProt[newStartPosition*3:]))
                        if nonePro != None:
                            if (len(nonePro) >= self.minProLen) and (self.maxProLen == 0 or len(nonePro) <= self.maxProLen):
                                tempNoneProt.append(nonePro)
                                tempNoneNucl.append(str(nucSeqOfProt))
            
            # Sort our top hits from each inter-stop codon fragment by size and category (i.e. mPro or altPro?) and select the top X hits
            if len(tempMProt + tempAltProt + tempNoneProt) >= 1:
                # Append '-' entries to lists which have less entries than we want to pull to allow the below 'for' loops to run without exceptions
                ## Prot list
                for i in range(0, self.hitsToPull-len(tempMProt)):
                    tempMProt.append('-')
                for i in range(0, self.hitsToPull-len(tempAltProt)):
                    tempAltProt.append('-')
                for i in range(0, self.hitsToPull-len(tempNoneProt)):
                    tempNoneProt.append('-')
                ## Nucl list
                for i in range(0, self.hitsToPull-len(tempMNucl)):
                    tempMNucl.append('-')
                for i in range(0, self.hitsToPull-len(tempAltNucl)):
                    tempAltNucl.append('-')
                for i in range(0, self.hitsToPull-len(tempNoneNucl)):
                    tempNoneNucl.append('-')
                
                # Sort the lists by size (largest on the bottom to allow the .pop() method to remove a hit when accepted)
                ## Prot
                tempSortedMProt = sorted(tempMProt, key=len)
                tempSortedAltProt = sorted(tempAltProt, key=len)
                tempSortedNoneProt = sorted(tempNoneProt, key=len)
                ## Nucl
                tempSortedMNucl = sorted(tempMNucl, key=len)
                tempSortedAltNucl = sorted(tempAltNucl, key=len)
                tempSortedNoneNucl = sorted(tempNoneNucl, key=len)
                
                # Run a final size comparison to choose the best ORF(s).
                tempOverallProt = []
                tempOverallNucl = []
                for i in range(0, self.hitsToPull):
                    # Accept a none protein if it's longer+cutoff than the alternatives
                    if len(tempSortedNoneProt[-1]) > len(tempSortedAltProt[-1]) + self.noCodonStringency and len(tempSortedNoneProt[-1]) > len(tempSortedMProt[-1]) + self.noCodonStringency:     # Again, we add the stringency values to help with determining priority of ORF ordering. Since this script will often be returning either 1, 3, or 5 potential ORFs, it is important that we order these in the most logical way
                        tempOverallProt.append(tempSortedNoneProt[-1])
                        tempSortedNoneProt.pop()
                        
                        tempOverallNucl.append(tempSortedNoneNucl[-1])
                        tempSortedNoneNucl.pop()
                    # Accept an alt protein if it's longer+cutoff than the methionine protein
                    elif len(tempSortedAltProt[-1]) > len(tempSortedMProt[-1]) + self.altCodonStringency:
                        tempOverallProt.append(tempSortedAltProt[-1])
                        tempSortedAltProt.pop()
                        
                        tempOverallNucl.append(tempSortedAltNucl[-1])
                        tempSortedAltNucl.pop()
                    # Accept an M protein if it isn't a blank
                    elif tempSortedMProt[-1] != "-":
                        tempOverallProt.append(tempSortedMProt[-1])                                                             # By using this as the 'else' position, sequences with methionine starts will be selected in the majority of situations as the default stringency settings ensure that it is rare an alternative start is used instead of a methionine start
                        tempSortedMProt.pop()
                        
                        tempOverallNucl.append(tempSortedMNucl[-1])                                                             # By using this as the 'else' position, sequences with methionine starts will be selected in the majority of situations as the default stringency settings ensure that it is rare an alternative start is used instead of a methionine start
                        tempSortedMNucl.pop()
                    # Default to a none protein if nothing else is available and it isn't blank
                    elif tempSortedNoneProt[-1] != "-":
                        tempOverallProt.append(tempSortedNoneProt[-1])
                        tempSortedNoneProt.pop()
                        
                        tempOverallNucl.append(tempSortedNoneNucl[-1])
                        tempSortedNoneNucl.pop()
                
                # Yield results
                for i in range(0, self.hitsToPull):
                    if len(tempOverallProt) <= i:
                        break
                    if tempOverallProt[i] != "-":
                        yield record.id, tempOverallProt[i], tempOverallNucl[i]

if __name__ == "__main__":
    pass

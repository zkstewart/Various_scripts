#! python3
# exome_liftover.py
# Program to enable discovery of exon sequences
# from genome sequences on the basis of exome
# sequencing alignments

import sys, argparse, os, math, statistics, platform, hashlib, time, random
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) # 2 dirs up is where we find dependencies
from Function_packages import ZS_SeqIO, ZS_HmmIO, ZS_Utility
import parasail # need to import this always to handle modules importing our ssw_parasail function
if platform.system() != 'Windows':
    from skbio.alignment import StripedSmithWaterman

def validate_args(args):
    # Validate input data location
    if not os.path.isdir(args.alignmentsDir):
        print('I am unable to locate the directory where the alignments files are (' + args.alignmentsDir + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    for exe in ["hmmpress", "hmmbuild", "hmmsearch"]:
        if not os.path.isfile(os.path.join(args.hmmerDir, exe)) and not os.path.isfile(os.path.join(args.hmmerDir, exe + ".exe")):
            print("{0} does not exist at {1}".format(exe, args.hmmerDir))
            print('Make sure you\'ve typed the location correctly and try again.')
            quit()
    if not os.path.isfile(args.genomeFile):
        print('I am unable to locate the genome FASTA file (' + args.genomeFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if args.identifier == "":
        print('identifier cannot be empty!')
        quit()
    # Validate numeric arguments
    if args.Evalue < 0:
        print('Evalue should be greater than 0')
        quit()
    if args.threads < 1:
        print('Threads should be greater than or equal to 1')
        quit()
    # Handle file output
    if os.path.isdir(args.outputDir):
        print("""Output directory already exists; note that this program will NOT
              overwrite existing files and instead will opt to resume any steps not
              completed. If files are corrupted or truncated expect unhandled errors""")
    else:
        try:
            os.mkdir(args.outputDir)
            print("Created '{0}' directory as part of argument validation".format(args.outputDir))
        except:
            print("Wasn't able to create '{0}' directory; does '{1}' actually exist?".format(args.outputDir, os.path.dirname(args.outputDir)))

def get_prediction_from_domdict(domDict):
    '''
    This method is intended to perform half the work of the liftover operation.
    Its goal is to take a domDict with potentially many domain predictions, and return
    the single best domain if such a scenario exists. If multiple competing options 
    exist, no result will be returned to indicate the inability for us to realistically
    predict which exon is the target of our liftover.
    
    Downstream, it still remains to be seen if this exon is a good match for liftover
    to occur.
    
    Params:
        domDict -- a dictionary with chromosome:[predictions] structure as per hmmer_parse()
                   or nhmmer_parse().
    '''
    # Flatten dictionary structure to a list
    flatList = [[key] + v for key, value in domDict.items() for v in value] # Put the chromosome ID at start of list
    
    # Adjust negative stranded predictions
    '''nhmmer seems to write coordinates in an adjusted fashion when strand is negative
    so that where we'd normally expect it to be start=[2] and end=[3], which throws off
    some of our operations.'''
    for i in range(len(flatList)):
        if flatList[i][7] == "-":
            start, end = flatList[i][3], flatList[i][2]
            flatList[i][2], flatList[i][3] = start, end
    
    # Sort list for priority of coordinate ordering prior to merging
    flatList.sort(key = lambda x: (x[0], x[2], x[3]))
    
    # Merge segments separated by gaps if relevant
    GAP_LEN_CUTOFF=500 # This is derived through checking the median values of gap lengths
    HMM_DIFF_ALLOWANCE_RATIO = 0.3 # This is an arbitrary magic number
    
    newFlatList = []
    skipIndices = [] # This list will stop us from making truncated chains after making the longest chain possible
    for i in range(len(flatList)):
        # Skip values that have been merged
        if i in skipIndices:
            continue
        
        entry1 = flatList[i] # Get the "start" of a potentially gappy bit
        possibleGaps = [flatList[i]] # This is a list to hold all potential fragments
        for x in range(len(flatList)):
            # Skip self-comparison
            if i == x:
                continue
            entry2 = flatList[x]
            # Skip bits on different chromosomes
            if entry1[0] != entry2[0]:
                continue
            # Skip bits on different strands
            if entry1[7] != entry2[7]:
                continue
            # Skip bits that aren't in the gap length allowance
            "Now, we wan't to use the last entry of possibleGaps to allow chaining of multiple bits"
            if (possibleGaps[-1][3] + GAP_LEN_CUTOFF) < entry2[2]:
                continue
            gapDist = entry2[2] - possibleGaps[-1][3] # Our flatList is sorted so this always works            
            # Skip bits where the gap in the sequence doesn't roughly match the gap in the HMM
            hmmGapDist = entry2[5] - possibleGaps[-1][6]
            if hmmGapDist not in range(int(gapDist - (gapDist * HMM_DIFF_ALLOWANCE_RATIO)), int(gapDist + (gapDist * HMM_DIFF_ALLOWANCE_RATIO))):
                continue
            # If we get through all this, it's probably a HMM match with a gap in it
            possibleGaps.append(entry2)
            skipIndices.append(x) # Since we're merging this, it shouldn't show up in the output list
        
        # Merge any relevant possibleGaps values
        newValue = entry1 # By default, we just keep this sequence as-is
        if len(possibleGaps) > 1:
            newStart = possibleGaps[0][2]
            newEnd = possibleGaps[-1][3]
            newHmmStart = possibleGaps[0][5]
            newHmmEnd = possibleGaps[-1][6]
            strand = entry1[7]
            
            evalues = [g[4] for g in possibleGaps]
            if 0 in evalues:
                newEvalue = 0
            else:
                newExponent = int(sum([math.log10(e) for e in evalues]))
                newEvalue = eval("1e{0}".format(newExponent)) # This is wack but idk maths so...
            newValue = [entry1[0], entry1[1], newStart, newEnd, newEvalue, newHmmStart, newHmmEnd, strand]
        
        # Update our newFlatList for downstream usage
        newFlatList.append(newValue)
    
    # Sort list for E-value priority prior to heuristic application
    newFlatList.sort(key = lambda x: (x[4], x[2], x[3]))
    
    # Perform prediction operation with simple heuristics
    ## Heuristic 1: single result
    if len(newFlatList) == 1:
        return newFlatList[0] # Simple, a single result is always our best result
    else:
        bestEvalue = newFlatList[0][4]
        secondBestEvalue = newFlatList[1][4]
        
        ## Heuristic 2: best result has E=0, second best does not!
        if bestEvalue == 0 and secondBestEvalue != 0:
            return newFlatList[0]
        
        ## Heuristic 3: best and second best results are both E=0
        "If they're both 0, then we just decide it's too hard to pick a 'best'"
        if bestEvalue == 0 and secondBestEvalue == 0:
            return None
        
        ## Heuristic 4: many results, but a much better E-value as a scaling value
        BETTER_FACTOR = 1.5 # 1.5 times better is always significant amount
        if abs(math.log10(bestEvalue)) >= abs(math.log10(secondBestEvalue)) * BETTER_FACTOR:
            return newFlatList[0]
        
        ## Heuristic 5: many results, but a large E-value difference as a flat value
        BETTER_FLAT_VALUE = 20 # an E-value difference of 1e-20 is always a significant amount
        if abs(math.log10(bestEvalue)) >= abs(math.log10(secondBestEvalue)) + BETTER_FLAT_VALUE:
            return newFlatList[0]
        
        ## Heuristic 6: many results, but HMM alignment length is better as a scaling value
        LENGTH_BETTER_FACTOR = 1.5 # 1.5 times better is always significant amount
        length1 = newFlatList[0][6] - newFlatList[0][5]
        length2 = newFlatList[1][6] - newFlatList[1][5]
        if length1 >= length2 * LENGTH_BETTER_FACTOR:
            return newFlatList[0]
        
        ## If all above heuristics fail to check out, we can't know for sure which exon is "the" exon
        else:
            return None

def get_hmm_length(hmmFile):
    '''
    This function is a bit rough, but it should work to find out how many
    residues are coded in a HMM to get the "length" of the HMM. I don't know
    how best to parse its coding scheme esp. when it pertains to amino acid
    containing models, so I'm using some basic heuristics which seem to apply
    for the files I'm looking at right now.
    '''
    currentCount = 1
    with open(hmmFile, "r") as fileIn:
        for line in fileIn:
            sl = line.split()
            if sl == []:
                continue
            elif sl[0] != str(currentCount):
                continue
            elif "." not in str(sl[1]):
                continue
            else:
                currentCount += 1
    return currentCount - 1 # Since we add one every iteration, we'll end up adding 1 too many by the end

def check_if_prediction_is_good(bestPrediction, hmmerObj, hmmFastaFile, 
                                genomeFastaFile, genome_FASTA_obj):
    '''
    This function will apply some general heuristics to see if the exon we've
    predicted seems to match the HMM well. If there's something minor wrong 
    with it, we'll fix it up here. Otherwise, we'll discard it.
    
    Params:
        bestPrediction -- a list created by get_prediction_from_domdict()
        hmmerObj -- a ZS_HmmIO.HMMER object configured with the HMM
        hmmFastaFile -- a string indicating the location of the FASTA file that
                        the HMM is based on.
        genomeFastaFile -- a string indicating the location of the FASTA file of
                           the genome sequences.
        genome_FASTA_obj -- a ZS_SeqIO.FASTA object containing the genome sequences.
    Returns:
        isGood -- a boolean indicating whether the prediction is good or not.
    '''    
    # Heuristic 1: Check if the exon mostly aligns end-to-end with the HMM
    '''
    From testing, this heuristic can't be applied because some of the alignments
    are horrible. I could do this if I applied my MSA trimming protocol first, but
    it makes things more complicated for not enough benefit. The E-value based
    heuristic below should deal with problematic MSAs like ENSSHAP00000000088-mx
    '''
    # ALLOWED_MISSING_RATIO = 0.15 # don't allow more than 15% of the sequence to not match the HMM
    # hmmLength = get_hmm_length(hmmer.HMM.hmmFile)
    # predictionHmmStart, predictionHmmEnd = bestPrediction[5:7]
    # overlapLength = hmmLength - ((hmmLength - predictionHmmEnd) - (predictionHmmStart - 1))
    # if overlapLength < hmmLength - (hmmLength * ALLOWED_MISSING_RATIO):
    #     return False
    
    # Heuristic 2: Check if the exon has a comparable E-value to sequences taken from the HMM itself
    ## Run HMMER again against the HMM's own sequences
    tmpHash = hashlib.sha256(bytes(str(hmmFastaFile) + str(time.time()) + str(random.randint(0, 100000)), 'utf-8') ).hexdigest()
    tmpTbloutName = ZS_Utility.tmp_file_name_gen("goodPred" + tmpHash[0:20], "tblout")
    
    hmmerObj.run(genomeFastaFile, tmpTbloutName, isNucleotide=True)
    domDict = hmmerObj.domDict
    
    os.unlink(tmpTbloutName) # Clean up temporary files now since we've saved .domDict
    
    ## Get the statistical distribution of E-values
    evalues = [ value[3] for chrom in domDict.keys() for value in domDict[chrom] ]
    evalues.sort()
    assert evalues != [], "HMMER search against own HMM fails to find results, alignment file is out of sync probably"
    
    medianEvalue = statistics.median(evalues)
    stdev = statistics.stdev(evalues)
    
    ## Calculate the minimum allowed floor
    STDEV_RANGE_FACTOR = 2 # we'll allow the exponent to increase or decrease by 2x the st.dev.'s exponent
    stdevExponent = math.log10(stdev)
    if medianEvalue != 0:
        medianExponent = math.log10(medianEvalue)
    else:
        medianExponent = -250 # zero has no exponent, so we can just set a very good E-value as our base
    worstAllowedEvalue = eval("1e{0}".format(int(medianExponent - (stdevExponent * STDEV_RANGE_FACTOR))))
    
    # If we fail because our E-value is worse than the floor for heuristic passing...
    if bestPrediction[4] > worstAllowedEvalue:
        # ... try to rescue it first to see if it's a gap sequence messing with things
        
        ## Rescue step 1: Find out how much sequence might be left out due to gap regions
        hmmLength = get_hmm_length(hmmerObj.hmmFile)
        predictionHmmStart, predictionHmmEnd = bestPrediction[5:7]
        
        potentialNewStart = bestPrediction[2] - predictionHmmStart + 1
        startDifference = bestPrediction[2] - potentialNewStart
        
        potentialNewEnd = bestPrediction[3] + (hmmLength - predictionHmmEnd)
        endDifference = potentialNewEnd - bestPrediction[3]
        
        ## Rescue step 2: Find the new sequence section we want to work with
        '''
        We only extend the head OR the tail, not both. It gets too complicated to try to rescue
        both ends, and it's unlikely we're going to run into these situations anyway.
        '''
        for FastASeq_obj in genome_FASTA_obj:
            if FastASeq_obj.id == bestPrediction[0]:
                break
        assert FastASeq_obj.id == bestPrediction[0], "Get in and fix things pls Zac"
        
        if endDifference > startDifference:
            newSequence = FastASeq_obj.seq[bestPrediction[2]-1: potentialNewEnd] # -1 to bring it into 0-based indexing
        else:
            newSequence = FastASeq_obj.seq[potentialNewStart-1: bestPrediction[3]] # No +1 needed here
        
        ## Rescue step 3: Check if the new sequence section contains gap characters
        LONG_GAP_LENGTH = 20 # shorter gaps shouldn't really register... I think
        if endDifference > startDifference:
            additionalSequence = FastASeq_obj.seq[bestPrediction[3]: potentialNewEnd]
        else:
            additionalSequence = FastASeq_obj.seq[potentialNewStart: bestPrediction[2]]
        if "n" * LONG_GAP_LENGTH not in additionalSequence.lower():
            '''
            If there isn't a sufficiently long gap to explain the failure to capture this bit,
            it's probably not a gaps fault and is because the exon is dodgy.
            '''
            return False
        
        ## Rescue step 4: If the new sequence DOES contain gaps, align it to the exon and check for good alignment
        ## 4.1: Create a consensus sequence for the exon
        FASTA_obj = ZS_SeqIO.FASTA(hmmFastaFile, isAligned=True)
        consensus = FASTA_obj.generate_consensus().replace("-", "")
        ## 4.2: Get the additionalSequence bit we want to align (sans N's)
        querySequence = max(additionalSequence.lower().split("n" * LONG_GAP_LENGTH), key = lambda x: len(x)).lstrip("n").upper()
        if querySequence == "":
            '''It can't be rescued if this happens, since it means the extended
            sequence bit is purely N's'''
            return False
        ## 4.3: Align it
        if platform.system() == 'Windows':
            queryAlign, targetAlign, startIndex, score = ssw_parasail(consensus, querySequence)
        else:
            queryAlign, targetAlign, startIndex, score = ssw_skbio(consensus, querySequence)
        if queryAlign == None:
            '''If this has happened, then we're probably using skbio and the
            alignment just didn't work. Haven't looked too deep into this.'''
            return False
        ## 4.4: Check if the query aligns well, fail it if not
        ALLOWED_NONALIGNING_RATIO = 0.1 # can only miss 10% of the extra sequence
        querySequenceLen = len(querySequence)
        if len(queryAlign.replace("-", "")) < querySequenceLen - (querySequenceLen*ALLOWED_NONALIGNING_RATIO):
            return False
        ## 4.5: Check if the aligned region is where we expect it to be, fail if not
        if endDifference > startDifference: # i.e., if it's a tail extension
            if startIndex < predictionHmmEnd: # i.e., if we're starting earlier than the point the original sequence ends
                return False
        else:
            if startIndex > predictionHmmStart: # i.e., if we're starting later than the point the original sequence starts
                return False
        ## If we get to here, it looks like we've got something good, so let's update the sequence details
        if endDifference > startDifference:
            bestPrediction[3] = potentialNewEnd
            bestPrediction[6] = startIndex + len(queryAlign)
        else:
            bestPrediction[2] = potentialNewStart
            bestPrediction[5] = startIndex
    
    ## If all the above heuristics pass, we can conclude that this exon is "good"
    return True

def ssw_parasail(targetString, queryString):
    '''
    Special implementation of striped Smith Waterman alignment for exon liftover
    project.
    '''
    # Perform SSW with parasail implementation
    profile = parasail.profile_create_sat(targetString, parasail.blosum62)
    alignment = parasail.sw_trace_striped_profile_sat(profile, queryString, 10, 1)
    targetAlign = alignment.traceback.query
    queryAlign = alignment.traceback.ref
    # Figure out where we're starting in the target with this alignment
    startIndex = targetString.find(targetAlign.replace('-', ''))
    
    return [queryAlign, targetAlign, startIndex, alignment.score]

def ssw_skbio(targetString, queryString):
    # Perform SSW with scikit.bio implementation
    query = StripedSmithWaterman(targetString)
    alignment = query(queryString)
    targetAlign = alignment.aligned_query_sequence
    queryAlign = alignment.aligned_target_sequence
    # Figure out where we're starting in the target with this alignment
    startIndex = targetString.find(targetAlign.replace('-', ''))
    
    return [queryAlign, targetAlign, startIndex, alignment.optimal_alignment_score]

def write_prediction_to_fasta(prediction, genome_FASTA_obj, identifier, outputFileName):
    '''
    Function to take the bestPrediction list output of get_prediction_from_domdict()
    and write it to a FASTA file.
    
    It's tuned for use in exome_liftover.py specifically, so keep that in mind.
    
    Params:
        prediction -- a list containing [chrom, id, genomeStart, genomeEnd, evalue, hmmStart, hmmEnd, strand]
        genome_FASTA_obj -- a ZS_SeqIO.FASTA object of the genome FASTA file.
        identifier -- a string uniquely identifying the species for which exons are being predicted.
        outputFileName -- a string indicating the file name and location to write the FASTA to
    '''
    # Extract relevant details
    chrom, id, start, end, _, _, _, strand = prediction
    seq = genome_FASTA_obj[chrom].seq[start-1: end] # -1 for 0-based indexing
    if strand == "-":
        seq = ZS_SeqIO.FastASeq.get_reverse_complement(None, seq) # None fills the role of self
    
    # Produce an alt ID
    altID = "{0} {1} chr={2} start={3} end={4}".format(id.rsplit("-")[0], identifier, chrom, start, end) # kill the -mx suffix in {0}
    
    # Write result to FASTA
    FastASeq_obj = ZS_SeqIO.FastASeq(id, seq=seq, alt=altID)
    with open(outputFileName, "w") as fileOut:
        fileOut.write(FastASeq_obj.get_str(withAlt=True))

def main():
    usage = """%(prog)s receives a directory full of aligned FASTA files as part of the
    Oz Mammals genome project. Its goal is to transform these alignments into HMMs that can
    then be queried against a genome of interest to locate the relevant exon sequence from
    said genome.
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-a", dest="alignmentsDir", required=True,
                help="Specify the directory where aligned FASTA files are located")
    p.add_argument("-hmm", dest="hmmerDir", required=True,
                help="Specify the directory where HMMER executables are located")
    p.add_argument("-g", dest="genomeFile", required=True,
                help="Specify the location of the genome FASTA file to find exons in")
    p.add_argument("-o", dest="outputDir", required=True,
                help="Output directory location (working and final files go here)")
    p.add_argument("-id", dest="identifier", required=True,
                help="Specify the species identifier to embed in sequence ID outputs")
    # Opts
    p.add_argument("-e", dest="Evalue", required=False, type=float,
                help="Optionally, specify the E-value cut-off for HMMER results (default==1e-20)",
                default=1e-20)
    p.add_argument("-t", dest="threads", required=False, type=int,
                help="Optionally, specify the number of threads to run HMMER with (default==1)",
                default=1)
    
    args = p.parse_args()
    validate_args(args)
    
    # Locate all files
    files = [os.path.join(args.alignmentsDir, file) for file in os.listdir(args.alignmentsDir)]
    
    # Create HMMs from aligned files
    hmmsDir = os.path.join(args.outputDir, "hmms")
    hmmsList = []
    os.makedirs(hmmsDir, exist_ok=True)
    for f in files:
        hmmName = os.path.basename(f).rsplit(".", maxsplit=1)[0]
        hmmFile = os.path.join(hmmsDir, hmmName + ".hmm")
        
        # Create HMM if it doesn't exist
        if not os.path.isfile(hmmName):
            hmm = ZS_HmmIO.HMM(args.hmmerDir, isNucleotide=True)
            hmm.create(f, hmmFile, hmmName=hmmName, isNucleotide=True)
        
        # Store HMM file name for later use
        hmmsList.append(hmmFile)
    
    # Load the genome FASTA for later use
    genome_FASTA_obj = ZS_SeqIO.FASTA(args.genomeFile)
    
    # Use our HMMs to query the genome for possible exon hits
    tbloutsDir = os.path.join(args.outputDir, "tblouts")
    fastasDir = os.path.join(args.outputDir, "fastas")
    os.makedirs(tbloutsDir, exist_ok=True)
    os.makedirs(fastasDir, exist_ok=True)
    for i in range(len(hmmsList)):
        hmmFile = hmmsList[i]
        
        # Derive exon FASTA name and skip if it already exists
        exonFastaFile = os.path.join(fastasDir, os.path.basename(hmmFile).rsplit(".", maxsplit=1)[0]) + ".fasta"
        if os.path.isfile(exonFastaFile):
            continue
        
        # If tblout doesn't exist, run HMMER
        tbloutName = os.path.join(
            tbloutsDir,
            "{0}.tblout".format(os.path.basename(hmmFile).rsplit(".", maxsplit=1)[0])
        )
        if not os.path.isfile(tbloutName):
            hmmer = ZS_HmmIO.HMMER(args.hmmerDir, hmmFile,
                                   threads=args.threads, evalue=args.Evalue)
            hmmer.run(args.genomeFile, tbloutName, isNucleotide=True)
            domDict = hmmer.domDict # This is what we want out of HMMER
        # If it does exist, simply load it in
        else:
            hmmer = ZS_HmmIO.HMMER(args.hmmerDir, hmmFile, # we use this later
                                   threads=args.threads, evalue=args.Evalue)
            domDict = ZS_HmmIO.nhmmer_parse(tbloutName, args.Evalue, extendedDetails=True)
        
        # If domDict is empty, skip this exon since we've failed to find it
        if domDict == {}:
            continue
        
        # Parse domDict and find the best exon
        bestPrediction = get_prediction_from_domdict(domDict)
        if bestPrediction == None:
            continue
        
        # If we could find a single "best" exon, check if it's any good at all
        fastaFile = files[i]
        isGood = check_if_prediction_is_good(bestPrediction, hmmer, fastaFile,
                                             args.genomeFile, genome_FASTA_obj)
        if not isGood:
            continue
        
        # If it's good, write it to file
        write_prediction_to_fasta(bestPrediction, genome_FASTA_obj, args.identifier, exonFastaFile)

    print("Program completed successfully!")

if __name__ == "__main__":
    main()

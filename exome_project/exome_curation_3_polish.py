#! python3
# exome_curation_3_polish.py
# Script which receives an MSA with a >Codons line as its first entry
# and attempts to fix indel errors present in a subset of the alignments.
# If they all have indels, god help you since I can't.

import sys, argparse, os
import numpy as np
sys.path.append(os.path.dirname(os.path.dirname(__file__))) # 2 dirs up is where we find dependencies
from Function_packages import ZS_SeqIO
from exome_liftover import ssw_parasail
from copy import deepcopy
from exome_curation_2_introns import solve_translation_frames, _get_segment_boundaries

def validate_args(args):
    # Validate input data location
    if not os.path.isdir(args.alignmentsDir):
        print('I am unable to locate the directory where the alignments files are (' + args.alignmentsDir + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
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

def polish_MSA_denovo(FASTA_obj):
    '''
    Polishes a ZS_SeqIO.FASTA object to remove probable indel errors from sequences.
    It does this without genomic evidence (hence "de novo") by assessment of how
    a subset of sequences containing stop codons differ to the majority which lack
    that stop codon. When a simple error can be identifier, it will rectify it.
    
    Params:
        FASTA_obj -- a ZS_SeqIO.FASTA instance
    Returns:
        TBD -- unsure what this needs to return at preset
    '''
    
    assert FASTA_obj[0].id == "Codons", "FASTA lacks a Codons line at its start!"
    
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
    for exon_FASTA_obj in exon_FASTA_objs:
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
        
        codonsProblemLeft = "5"*startTrim # This will be left-appended to the codons sequence for this subregion
        codonsProblemRight = "5"*endTrim # And this will be right-appended to the ...
        
        # Predict our solutionDict again if needed
        "If we've changed our sequence region, we might need to update our translations"
        if startTrim != 0 or endTrim != 0:
             solutionDict = solve_translation_frames(exon_FASTA_obj)
        
        # 
        
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
                    problemAlign, targetAlign, startIndex, score = ssw_parasail(nuclSeq, targetNuclSeq)
                    matches.append([problemAlign, targetAlign, startIndex, score])
            matches.sort(key = lambda x: -x[3]) # order by score
            bestMatch = matches[0]
            
            # Try to see if there's any easily fixable indel errors
            ## TBD...
        
        
        # # Identify our best sequences
        # seqs = []
        # for i in range(len(exon_FASTA_obj)):
        #     FastASeq_obj = exon_FASTA_obj[i]
        #     if i not in solutionDict:
        #         continue
        #     else:
        #         seq, frame, hasStopCodon = solutionDict[i]
        #         if not hasStopCodon:
        #             seqs.append(seq)
        
        # lengthCutoff = np.percentile([len(s) for s in seqs], 100-EXCLUSION_PCT)
        # MAXIMUM_GOOD_SEQS = 10 # hard-coded, we just want to minimise wasted time when there's lots of good sequences
        # bestSeqs = []
        # for i in range(len(seqs)):
        #     if len(bestSeqs) >= MAXIMUM_GOOD_SEQS:
        #         break
            
        #     if len(seqs[i]) >= lengthCutoff:
        #         bestSeqs.append(seqs[i])
        # # 
        # '''
        # It's important to consider something at this point. solve_translation_frames() is not
        # guaranteed to provide a good answer. It will only provide a good answer if at least one
        # sequence does NOT contain stop codons in it. If that doesn't prove true, it simply
        # won't work.
        # '''
        # ongoingCount=0
        # for x in range(min(solutionDict.keys()), max(solutionDict.keys())+1):
        #     if x not in solutionDict:
        #         continue
        #     else:
        #         b = boundaries[ongoingCount]
        #         s = solutionDict[x]
        #         ongoingCount +=1
        

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
        alignFastaFile = files[i]
        FASTA_obj = fastaObjs[i]
        
        # Perform polishing procedure
        
        # Write output FASTA file
        outputFileName = os.path.join(args.outputDir, os.path.basename(alignFastaFile))
        FASTA_obj.write(outputFileName, withDescription=True, asAligned=True)
    
    print("Program completed successfully!")

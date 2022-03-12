#! python3
# exome_liftover.py
# Program to enable discovery of exon sequences
# from genome sequences on the basis of exome
# sequencing alignments

import sys, argparse, os, math
sys.path.append(os.path.dirname(os.path.dirname(__file__))) # 2 dirs up is where we find dependencies
from Function_packages import ZS_SeqIO, ZS_HmmIO

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

def predict_exon_from_domDict(domDict):
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
    flatList = [[key] + v for key, value in domDict.items() for v in value] # Put the chromosome ID at start of list
    flatList.sort(key = lambda x: (x[4], x[2], x[3]))
    
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
            
            evalues = [g[4] for g in possibleGaps]
            newExponent = int(sum([math.log10(e) for e in evalues]))
            newEvalue = eval("1e{0}".format(newExponent)) # This is wack but idk maths so...
            
            newValue = [entry1[0], entry1[1], newStart, newEnd, newEvalue, newHmmStart, newHmmEnd] # Drops strand here, who cares
        
        # Update our newFlatList for downstream usage
        newFlatList.append(newValue)
    newFlatList.sort(key = lambda x: (x[4], x[2], x[3]))
    
    # Perform prediction operation with simple heuristics
    ## Heuristic 1: single result
    if len(newFlatList) == 1:
        return newFlatList[0] # Simple, a single result is always our best result
    else:
        bestEvalue = newFlatList[0][4]
        secondBestEvalue = newFlatList[1][4]
        
        ## Heuristic 2: many results, but a much better E-value as a scaling value
        BETTER_FACTOR = 1.5 # 1.5 times better is always significant amount
        if abs(math.log10(bestEvalue)) >= abs(math.log10(secondBestEvalue)) * BETTER_FACTOR:
            return newFlatList[0]
        
        ## Heuristic 3: many results, but a large E-value difference as a flat value
        BETTER_FLAT_VALUE = 20 # an E-value difference of 1e-20 is always a significant amount
        if abs(math.log10(bestEvalue)) >= abs(math.log10(secondBestEvalue)) + BETTER_FLAT_VALUE:
            return newFlatList[0]
        
        ## Heuristic 4: many results, but HMM alignment length is better as a scaling value
        LENGTH_BETTER_FACTOR = 1.5 # 1.5 times better is always significant amount
        length1 = newFlatList[0][6] - newFlatList[0][5]
        length2 = newFlatList[1][6] - newFlatList[1][5]
        if length1 >= length2 * LENGTH_BETTER_FACTOR:
            return newFlatList[0]
        
        ## If all above heuristics fail to check out, we can't know for sure which exon is "the" exon
        else:
            return None

if __name__ == "__main__":
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
        hmmName = os.path.join(hmmsDir, os.path.basename(f).rsplit(".", maxsplit=1)[0] + ".hmm")
        hmm = ZS_HmmIO.HMM(args.hmmerDir, isNucleotide=True)
        
        # If HMM exists, load it in
        if os.path.isfile(hmmName):
            hmm.load_HMM_file(hmmName)
        # Create HMM if it doesn't exist
        else:
            hmm.load_FASTA_from_file(f)
            hmm.create_HMM(hmmName, hmmBuildExtraArgs="--dna")
        hmmsList.append(hmm)
    
    # Use our HMMs to query the genome for possible exon hits
    tbloutsDir = os.path.join(args.outputDir, "tblouts")
    os.makedirs(tbloutsDir, exist_ok=True)
    exonPredictions = [] # This will store our successful hits
    for hmm in hmmsList:
        # Derive domtblout name
        tbloutName = os.path.join(
            tbloutsDir,
            "{0}.tblout".format(os.path.basename(hmm.hmmFile).rsplit(".", maxsplit=1)[0])
        )
        
        # If domtblout doesn't exist, run HMMER
        if not os.path.isfile(tbloutName):
            hmmer = ZS_HmmIO.HMMER(hmm)
            hmmer.load_FASTA_from_file(args.genomeFile)
            hmmer.set_output_name(tbloutName)
            hmmer.set_Evalue(args.Evalue)
            hmmer.set_threads(args.threads)
            hmmer.run_search()
            domDict = hmmer.domDict # This is what we want out of HMMER
        # If it does exist, simply load it in
        else:
            domDict = ZS_HmmIO.nhmmer_parse(tbloutName, args.Evalue, extendedDetails=True)

        # If domDict is empty, skip this exon since we've failed to find it
        if domDict == {}:
            continue
        
        # Locate the best exon prediction from the domDict
        for tbloutName in tbloutFiles:
            domDict = ZS_HmmIO.hmmer_parse(tbloutName, Evalue, extendedDetails=True) ## Testing
            if domDict == {}:
                continue
            bestPrediction = predict_exon_from_domDict(domDict)
            if bestPrediction == None:
                continue
            # If we could find a single "best" exon, check if it's any good at all
            ## TBD...
        
    print("Program completed successfully!")

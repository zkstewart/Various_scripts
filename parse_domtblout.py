#! python3
# parse_domtblout

# Python script to handle domain overlaps from HMMER domtblout results.
# This involves a two-step system. In the first step, identical domains are
# specifically handled by trimming, splitting, or merging these models based
# on the percentage overlap shared between them. The second system uses a
# process I refer to as "seed overlap resolving" which involes selecting the
# most significant domain prediction (by E-value) and trimming domains that
# only slightly overlap it, and removing any that overlap it beyond a user-specified
# percentage. After this, it moves to the next most significant prediction and
# performs the same procedure. The result is a domain list which prioritises
# the "best" domains and removes overlaps.

import argparse, os

# Define functions for later use

## Validate arguments
def validate_args(args):
        # Validate input file locations
        if not os.path.isfile(args.inputHmmer):
                print('I am unable to locate the HMMER domtblout file (' + args.inputHmmer + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Handle numerical arguments
        if args.evalue < 0.0:
                print('E-value cannot be < 0. Try again.')
                quit()
        if not 0 <= args.ovlCutoff <= 100.0:
                print('Percentage overlap cutoff must be given as a value <= 0 and >= 100. Try again.')
                quit()
        args.ovlCutoff = args.ovlCutoff / 100  
        # Handle file overwrites
        if os.path.isfile(args.outputFileName):
                print(args.outputFileName + ' already exists. Delete/move/rename this file and run the program again.')
                quit()
        return args

## Domain overlap handling
def findMiddle(input_list):             # https://stackoverflow.com/questions/38130895/find-middle-of-a-list
    middle = float(len(input_list))/2
    if middle % 2 != 0:
        return [input_list[int(middle - .5)]]
    else:
        return [input_list[int(middle)], input_list[int(middle-1)]]

def ovl_resolver(ovlCutoff, inputList):
        import copy
        seqHits = copy.deepcopy(inputList)
        seqHits.sort(key = lambda x: (x[3], x[1], x[2]))        # Format is [domainID, start, stop, E-value]
        while True:
                # Flow control: if no possibility for overlap, break
                if len(seqHits) == 1:
                        break
                # Flow control: if no overlaps remain, break
                overlapping = 'n'
                for y in range(len(seqHits)-1):
                        for z in range(y+1, len(seqHits)):
                                if seqHits[y][2] >= seqHits[z][1] and seqHits[z][2] >= seqHits[y][1]:           # i.e., if there is overlap
                                        overlapping = 'y'
                                        break
                if overlapping == 'n':
                        break
                # Handle overlaps through looping structure
                seqHits = seed_looping_structure(seqHits, ovlCutoff)
        seqHits.sort(key = lambda x: (x[1], x[2]))
        return seqHits

def seed_looping_structure(seqHits, ovlCutoff):
        # Set up
        import copy
        origSeqHits = copy.deepcopy(seqHits)    # This lets us compare our new domain against the original to make sure we haven't excessively cut and trimmed it
        origCutoff = 0.60       # Arbitrary; this means that, if the trimmed model is less than "a bit above" the length of the original, we drop it entirely
        # Main loop
        for y in range(len(seqHits)-1):
                z = y + 1
                while True:
                        # Exit condition
                        if z >= len(seqHits):   # len(seqHits)-1 would correspond to the final entry, this means we've gone at least one step further beyond
                                break
                        # Trim 1-bp overlaps [note the consistent design choice: the most N-proximal domain gets the extra AA position for no particular reason, we just need to handle this]
                        if seqHits[y][2] == seqHits[z][1]:
                                seqHits[z][1] =  seqHits[z][1] + 1
                        if seqHits[z][2] == seqHits[y][1]:
                                seqHits[y][1] =  seqHits[y][1] + 1
                        # If there is overlap, resolve this
                        if seqHits[y][2] > seqHits[z][1] and seqHits[z][2] > seqHits[y][1]:
                                # Get details of sequence overlap
                                sharedPos = set(range(max(seqHits[y][1], seqHits[z][1]), min(seqHits[y][2], seqHits[z][2]) + 1))
                                ovlLen = len(sharedPos)
                                seq1Perc = ovlLen / (seqHits[y][2] - seqHits[y][1] + 1)
                                seq2Perc = ovlLen / (seqHits[z][2] - seqHits[z][1] + 1)
                                bestEval = min(seqHits[y][3], seqHits[z][3])
                                # Handle slight mutual overlaps by trimming based on best E-value
                                if seq1Perc < ovlCutoff and seq2Perc < ovlCutoff:
                                        ## Identical E-values [mutual trimming]
                                        if seqHits[y][3] == seqHits[z][3]:
                                                posList = list(sharedPos)
                                                posList.sort()
                                                midPoint = findMiddle(posList)
                                                if seqHits[y][1] < seqHits[z][1]:
                                                        seqHits[y][2] = midPoint[0]
                                                        seqHits[z][1] = midPoint[0] + 1
                                                else:
                                                        seqHits[z][2] = midPoint[0]
                                                        seqHits[y][1] = midPoint[0] + 1
                                        ## Different E-values [trim lower E-value]
                                        elif bestEval == seqHits[y][3]:
                                                if seqHits[y][1] < seqHits[z][1]:
                                                        seqHits[z][1] = seqHits[y][2] + 1
                                                else:
                                                        seqHits[z][2] = seqHits[y][1] - 1
                                        else:
                                                if seqHits[y][1] < seqHits[z][1]:
                                                        seqHits[y][2] = seqHits[z][1] - 1
                                                else:
                                                        seqHits[y][1] = seqHits[z][2] + 1
                                        # If we've trimmed one of these domains too much, drop it
                                        assert origSeqHits[y][0] == seqHits[y][0]       # Make sure that things are working correctly
                                        assert origSeqHits[z][0] == seqHits[z][0]
                                        changed = False
                                        if (seqHits[z][2] - seqHits[z][1] + 1) / (origSeqHits[z][2] - origSeqHits[z][1] + 1) < origCutoff:      # Need to handle z first lest we upset the ordering
                                                del seqHits[z]
                                                del origSeqHits[z]
                                                changed = True
                                        if (seqHits[y][2] - seqHits[y][1] + 1) / (origSeqHits[y][2] - origSeqHits[y][1] + 1) < origCutoff:
                                                del seqHits[y]
                                                del origSeqHits[y]
                                                changed = True
                                        if changed == False:
                                                z += 1  # We've made the current pair compatible, now we can just move onto the next pairing
                                # Handle larger overlaps by deleting based on E-value
                                else:
                                        ## Identical E-values [delete the most C-proximal]
                                        if seqHits[y][3] == seqHits[z][3]:
                                                if seqHits[y][1] < seqHits[z][1]:
                                                        del seqHits[z]
                                                        del origSeqHits[z]      # Keep these lists equivalent
                                                else:
                                                        del seqHits[y]
                                                        del origSeqHits[y]
                                        ## Different E-values [delete the lowest E-value]
                                        elif bestEval == seqHits[y][3]:
                                                del seqHits[z]
                                                del origSeqHits[z]
                                        else:
                                                del seqHits[y]
                                                del origSeqHits[y]
                                                # We make no changes to our z value since we deleted a sequence
                        # If there is no overlap, continue the loop
                        else:
                                z += 1
        return seqHits

def split_middle(sharedPos, modelGroup, y):
        splitPos = list(sharedPos)
        splitPos.sort()
        middle = findMiddle(splitPos)
        if len(middle) == 1:
                modelGroup[y] = [*modelGroup[y][0:2], middle[0], modelGroup[y][3]]
                modelGroup[y+1] = [modelGroup[y+1][0], middle[0]+1, *modelGroup[y+1][2:]]
        else:
                modelGroup[y] = [*modelGroup[y][0:2], middle[1], modelGroup[y][3]]
                modelGroup[y+1] = [modelGroup[y+1][0], middle[0], *modelGroup[y+1][2:]]
        return modelGroup

def join_models(modelGroup, y):
        firstPos = modelGroup[y][1]
        lastPos = max(modelGroup[y][2], modelGroup[y+1][2])
        highestEval = max(modelGroup[y][3], modelGroup[y+1][3])                 # We want to associate the worst E-value to the joined model for the purpose of later handling
        modelGroup[y] = [modelGroup[y][0], firstPos, lastPos, highestEval]      # This is mostly because joining these models isn't entirely "natural" but it is more accurate than the significant overlap scenario that caused the joing to occur
        del modelGroup[y+1]
        return modelGroup

def single_database_domain_overlap_loop(domDict):
        # Setup
        finalDict = {}
        extensCutoff = 20       # This is arbitrary; seems to work well, don't see any reason why this should be variable by the user
        for key, value in domDict.items():
                # Figure out which domain models are associated with this gene model
                uniqueModels = []
                for val in value:
                        uniqueModels.append(val[0])
                uniqueModels = list(set(uniqueModels))
                # Collapse overlaps of identical domains
                collapsedIdentical = []
                for model in uniqueModels:
                        modelGroup = []
                        for val in value:
                                if val[0] == model:
                                        modelGroup.append(val)
                        modelGroup.sort(key = lambda x: (x[1], x[2]))   # Technically this should not be needed - the HMMER domtblout file is pre-sorted - but it's useful to put here _just in case_, and to make it clear that this script operates on the basis of this sorting
                        # Begin collapsing process
                        overlapping = 'y'
                        while True:
                                if len(modelGroup) == 1 or overlapping == 'n':
                                        break
                                for y in range(len(modelGroup)-1):
                                        if modelGroup[y+1][1] > modelGroup[y][2] and y != len(modelGroup)-2:    # i.e., if the start of seq2 > end of seq1, there is no overlap; we also want to skip this if it's the last pair we're inspecting since that will allow us to reach the final "else" condition and exit out of the loop
                                                continue
                                        elif modelGroup[y+1][1] == modelGroup[y][2]:                            # i.e., if the start of seq2 == end of seq1, there is 1 bp of overlap to handle
                                                modelGroup[y+1][1] =  modelGroup[y+1][1] + 1                    # Consistent design choice: the most N-proximal domain gets the extra AA position for no particular reason, we just need to handle this
                                                continue
                                        elif modelGroup[y+1][1] < modelGroup[y][2]:                             # i.e., if the start of seq2 < end of seq1, there is more than 1bp of overlap of handle
                                                # Calculate overlap proportion
                                                seq1Len = modelGroup[y][2] - modelGroup[y][1] + 1
                                                seq2Len = modelGroup[y+1][2] - modelGroup[y+1][1] + 1
                                                sharedPos = set(range(max(modelGroup[y][1], modelGroup[y+1][1]), min(modelGroup[y][2], modelGroup[y+1][2]) + 1))        # +1 to offset Python counting up-to but not including the last value in a range
                                                ovlLen = len(sharedPos)
                                                r1Perc = ovlLen / (seq1Len + 1)
                                                r2Perc = ovlLen / (seq2Len + 1)
                                                highest = max(r1Perc, r2Perc)
                                                lowest = min(r1Perc, r2Perc)
                                                # Determine the length of the sequence extension of the most-overlapped sequence
                                                if highest == 0.50:
                                                        longest = max(seq1Len, seq2Len)
                                                        if longest == seq1Len:
                                                                extension = seq2Len - ovlLen
                                                        else:
                                                                extension = seq1Len - ovlLen
                                                elif highest == r1Perc:
                                                        extension = seq1Len - ovlLen
                                                else:
                                                        extension = seq2Len - ovlLen
                                                ## Handle the various scenarios indicated by the highest/lowest values
                                                # Scenario 1: (TRIM BASED ON E-VALUE) small overlap of both sequences
                                                if highest <= 0.20:
                                                        if modelGroup[y][3] < modelGroup[y+1][3]:
                                                                # Trim y+1
                                                                modelGroup[y+1] = [modelGroup[y+1][0], modelGroup[y][2]+1, *modelGroup[y+1][2:]]
                                                        elif modelGroup[y+1][3] < modelGroup[y][3]:
                                                                # Trim y
                                                                modelGroup[y] = [*modelGroup[y][0:2], modelGroup[y+1][1]-1, modelGroup[y][3]]
                                                        else:
                                                                # If the two E-value are identical, we just split down the middle!
                                                                modelGroup = split_middle(sharedPos, modelGroup, y)
                                                        continue
                                                # Scenario 2: (SPLIT MIDDLE) intermediate overlap of one sequence with a significant length of sequence extension beyond the overlap region
                                                elif extension > extensCutoff and lowest <= 0.80:
                                                        modelGroup = split_middle(sharedPos, modelGroup, y)
                                                        continue
                                                # Scenario 3: (JOIN) intermediate or large overlap of one sequence with a short length of sequence extension beyond the overlap region
                                                else:
                                                        modelGroup = join_models(modelGroup, y)
                                                        break
                                        else:   # We need the y != check above since we need to set an exit condition when no more overlaps are present. The if/elif will always trigger depending on whether there is/is not an overlap UNLESS it's the second last entry and there is no overlap. In this case we finally reach this else clause, and we trigger an exit.
                                                overlapping = 'n'
                                                break
                        # Add corrected individual models to collapsedIdentical list
                        collapsedIdentical += modelGroup
                # Process collapsedIdentical list to get our list of domains annotated against the sequence from each individual database
                if len(collapsedIdentical) == 1:
                        if key not in finalDict:
                                finalDict[key] = collapsedIdentical
                        else:
                                finalDict[key].append(collapsedIdentical)
                else:
                        collapsedIdentical = ovl_resolver(args.ovlCutoff, collapsedIdentical)           # We've merged, joined, and trimmed identical domain models above. Now, we're looking at different domains from the same database.
                        if key not in finalDict:                                                        # We employ a similar approach here, but it's focused on E-values rather than on overlap proportions.
                                finalDict[key] = collapsedIdentical
                        else:
                                finalDict[key].append(collapsedIdentical)
        return finalDict

## Ensure accuracy
def dom_dict_check(finalDict):
        for key, value in finalDict.items():
                seqHits = value
                for y in range(len(seqHits)-1):
                        z = y + 1
                        while True:
                                # Exit condition
                                if z >= len(seqHits):   # len(seqHits)-1 would correspond to the final entry, this means we've gone at least one step further beyond
                                        break
                                # Checks
                                assert seqHits[y][2] != seqHits[z][1]
                                assert seqHits[z][2] != seqHits[y][1]
                                assert not (seqHits[y][2] > seqHits[z][1] and seqHits[z][2] > seqHits[y][1])
                                z += 1

## Domtblout parsing
def hmmer_parse(domtbloutFile, evalueCutoff):
        domDict = {}                            # We need to use a dictionary for later sorting since hmmsearch does not produce output that is ordered in the way we want to work with. hmmscan does, but it is SIGNIFICANTLY slower.
        with open(domtbloutFile, 'r') as fileIn:
                for line in fileIn:
                        # Skip unnecessary lines
                        if line.startswith('#') or line == '' or line == '\n' or line == '\r\n':
                                continue
                        # Parse line and skip if evalue is not significant
                        sl = line.rstrip('\r\n').split()
                        evalue = float(sl[12])
                        if evalue > float(evalueCutoff):
                                continue
                        # Get relevant details
                        pid = sl[0]
                        did = sl[3]
                        dstart = int(sl[17])
                        dend = int(sl[18])
                        # Add into domain dictionary
                        if pid not in domDict:
                                domDict[pid] = [[did, dstart, dend, evalue]]
                        else:
                                domDict[pid].append([did, dstart, dend, evalue])
        return domDict

## Output function
def output_func(inputDict, outputFileName):
        with open(outputFileName, 'w') as fileOut:
                for key, value in inputDict.items():
                        fileOut.write(key + '\t' + '\t'.join(list(map(str, value))) + '\n')

#### USER INPUT SECTION
usage = """%(prog)s reads .domtblout file and returns non-overlapped (or
specified percent of overlapping) domains below given e-value.
"""

# Reqs
p = argparse.ArgumentParser(description=usage)
p.add_argument("-i", "-input", dest="inputHmmer",
                   help="Input domtblout HMMER3 domtblout file.")
p.add_argument("-e", "-evalue", dest="evalue", type=float,
                   help="E-value significance cut-off for domain predictions (default == 1e-3)", default=1e-3)
p.add_argument("-p", "-percOvl", dest="ovlCutoff", type=float,
                   help="Percentage overlap cutoff (below == trimming to prevent overlap, above = deletion of lower E-value hit, default == 25.0).", default=25.0)
p.add_argument("-o", "-output", dest="outputFileName",
                   help="Output file name.")

args = p.parse_args()
args = validate_args(args)

# Parse hmmer domblout file
domDict = hmmer_parse(args.inputHmmer, args.evalue)

# Delve into parsed hmmer dictionary and sort out overlapping domain hits from different databases
finalDict = single_database_domain_overlap_loop(domDict)

# Check that the program worked correctly
dom_dict_check(finalDict)

# Generate output
output_func(finalDict, args.outputFileName)
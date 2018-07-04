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
def validate_args(args, dom_prefixes):
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
        # Handle domain prefix selection
        if args.databaseSelect != False:
                dom_prefixes = [args.databaseSelect]
                args.hmmdbScript = False        # These options are incompatible, make sure it's turned off here
        # Handle file overwrites
        if os.path.isfile(args.outputFileName):
                print(args.outputFileName + ' already exists. Delete/move/rename this file and run the program again.')
                quit()
        return args, dom_prefixes

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

def single_database_domain_overlap_loop(domDict, ovlCutoff):
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
                        collapsedIdentical = ovl_resolver(ovlCutoff, collapsedIdentical)           # We've merged, joined, and trimmed identical domain models above. Now, we're looking at different domains from the same database.
                        if key not in finalDict:                                                        # We employ a similar approach here, but it's focused on E-values rather than on overlap proportions.
                                finalDict[key] = collapsedIdentical
                        else:
                                finalDict[key].append(collapsedIdentical)
        return finalDict

def hmm_db_download_domain_overlap_loop(domDict, dom_prefixes, ovlCutoff, databaseSelect):
        # Setup
        finalDict = {}
        extensCutoff = 20       # This is arbitrary; seems to work well, don't see any reason why this should be variable by the user
        for key, value in domDict.items():
                # Compare models from within each domain database and handle overlaps
                for prefix in dom_prefixes:
                        prefixHits = []
                        for val in value:
                                if prefix == 'SUPERFAMILY':
                                        if val[0].isdigit():    # Remember, as mentioned above, SUPERFAMILY models are just digits. No other database has the same model naming scheme so we can detect these with this check
                                                prefixHits.append(val)
                                else:
                                        if val[0].startswith(prefix):
                                                prefixHits.append(val)
                        if prefixHits == []:
                                continue
                        # Figure out which domain models are associated with this gene model and this particular domain database
                        uniqueModels = []
                        for val in prefixHits:
                                uniqueModels.append(val[0])
                        uniqueModels = list(set(uniqueModels))
                        # Collapse overlaps of identical domains
                        collapsedIdentical = []
                        for model in uniqueModels:
                                modelGroup = []
                                for val in prefixHits:
                                        if val[0] == model:
                                                modelGroup.append(val)
                                modelGroup.sort(key = lambda x: (x[1], x[2]))                           # Technically this should not be needed - the HMMER domtblout file is pre-sorted - but it's useful to put here _just in case_, and to make it clear that this script operates on the basis of this sorting
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
                                        if databaseSelect == False:
                                                finalDict[key] = [collapsedIdentical]
                                        else:                                           # If databaseSelect != False, we are only getting 1 database's output, and thus we want a single list item and not one for each database
                                                finalDict[key] = collapsedIdentical
                                else:
                                        finalDict[key].append(collapsedIdentical)
                        else:
                                collapsedIdentical = ovl_resolver(ovlCutoff, collapsedIdentical)      # We've merged, joined, and trimmed identical domain models above. Now, we're looking at different domains from the same database.
                                if key not in finalDict:                                                        # We employ a similar approach here, but it's focused on E-values rather than on overlap proportions.
                                        if databaseSelect == False:
                                                finalDict[key] = [collapsedIdentical]
                                        else:
                                                finalDict[key] = collapsedIdentical
                                else:
                                        finalDict[key].append(collapsedIdentical)
        return finalDict

## Ensure accuracy
def dom_dict_check(finalDict, hmmdbDict):
        for key, value in finalDict.items():
                seqHits = value
                if hmmdbDict == True:
                        for val in seqHits:
                                for y in range(len(val)-1):
                                        z = y + 1
                                        while True:
                                                # Exit condition
                                                if z >= len(val):   # len(val)-1 would correspond to the final entry, this means we've gone at least one step further beyond
                                                        break
                                                # Checks
                                                assert val[y][2] != val[z][1]
                                                assert val[z][2] != val[y][1]
                                                assert not (val[y][2] > val[z][1] and val[z][2] > val[y][1])
                                                z += 1
                else:
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

## Output functions
def output_func(inputDict, outputFileName):
        with open(outputFileName, 'w') as fileOut:
                for key, value in inputDict.items():
                        fileOut.write(key + '\t' + '\t'.join(list(map(str, value))) + '\n')

def hmmdb_output_func(inputDict, outputFileName, ovlCutoff):
        # Set up
        dom_prefixes = ('cd', 'COG', 'KOG', 'LOAD', 'MTH', 'pfam', 'PHA', 'PRK', 'PTZ', 'sd', 'smart', 'TIGR', 'PLN', 'CHL', 'cath', 'SUPERFAMILY')
        # Main loop
        with open(outputFileName, 'w') as fileOut:
                for key, value in inputDict.items():
                        # Format hit receptacle list
                        hitReceptacle = ['-']*len(dom_prefixes)
                        for i in range(len(dom_prefixes)):
                                if dom_prefixes[i] != 'SUPERFAMILY':
                                        for val in value:
                                                if val[0][0].startswith(dom_prefixes[i]):
                                                        hitReceptacle[i] = '; '.join(list(map(str, val)))
                                else:
                                        for val in value:
                                                if val[0][0].isdigit():
                                                        hitReceptacle[i] = '; '.join(list(map(str, val)))
                        # Create a single column entry summarising all the different databases
                        seqHits = []
                        for val in value:
                                seqHits += val
                        seqHits.sort(key = lambda x: (x[1], x[2], x[3]))
                        if len(seqHits) != 1:
                                seqHits = ovl_resolver(ovlCutoff, seqHits)
                        hitReceptacle.insert(0, '; '.join(list(map(str, seqHits))))
                        # Write to file
                        fileOut.write(key + '\t' + '\t'.join(hitReceptacle) + '\n')

# Set up values needed for rest of script
dom_prefixes = ['cd', 'COG', 'KOG', 'LOAD', 'MTH', 'pfam', 'PHA', 'PRK', 'PTZ', 'sd', 'smart', 'TIGR', 'PLN', 'CHL', 'cath', 'SUPERFAMILY']    # These encompass the databases currently part of NCBI's CDD, and cath which I add to this resource. SUPERFAMILY is also included, but it is purely numbers so no prefix is applicable; if it lacks any of these prefixes, it's a SUPERFAMILY domain.

#### USER INPUT SECTION
usage = """%(prog)s reads a HMMER domtblout file and returns non-overlapping
domain predictions which pass user-specified E-value cut-off. Program can
operate in two modes. Without providing -hmm or -d argument, the input file
will be handled normally. Providing either of these two arguments will treat
the domtblout file as if it was produced by the hmm_db_download.py script
and can provide an alternative output format (-hmm; format is like that of
an annotation table; first column is combined database result with subsequent
columns being from individual databases) or produce output from a single
database (-d; extracts the specified databases' results to a normally-formatted
output file).
"""

# Reqs
p = argparse.ArgumentParser(description=usage)
p.add_argument("-i", "-input", dest="inputHmmer",
                   help="Input domtblout HMMER3 domtblout file.")
p.add_argument("-e", "-evalue", dest="evalue", type=float,
                   help="E-value significance cut-off for domain predictions (default == 1e-3)", default=1e-3)
p.add_argument("-p", "-percOvl", dest="ovlCutoff", type=float,
                   help="Percentage overlap cutoff (below == trimming to prevent overlap, above = deletion of lower E-value hit, default == 25.0).", default=25.0)
p.add_argument("-hmm", "-hmmdbScript", dest="hmmdbScript", action='store_true',
                   help="Optionally specify whether the HMM database was formatted by the hmm_db_download.py script and you want annotation table formatted results.", default=False)
p.add_argument("-d", "-databaseSelect", dest="databaseSelect", choices=dom_prefixes,
                   help="If the HMM database was formatted by the hmm_db_download.py script but you want to extract the results of just one database, choose here. Note that specifying this argument overrides -hmm", default=False)
p.add_argument("-o", "-output", dest="outputFileName",
                   help="Output file name.")

args = p.parse_args()
args, dom_prefixes = validate_args(args, dom_prefixes)

# Parse hmmer domblout file
domDict = hmmer_parse(args.inputHmmer, args.evalue)

# Delve into parsed hmmer dictionary and sort out overlapping domain hits from different databases
if args.hmmdbScript == False and args.databaseSelect == False:
        finalDict = single_database_domain_overlap_loop(domDict, args.ovlCutoff)
else:
        finalDict = hmm_db_download_domain_overlap_loop(domDict, dom_prefixes, args.ovlCutoff, args.databaseSelect)

# Check that the program worked correctly
dom_dict_check(finalDict, args.hmmdbScript)

# Generate output
if args.hmmdbScript == False:
        output_func(finalDict, args.outputFileName)
else:
        hmmdb_output_func(finalDict, args.outputFileName, args.ovlCutoff)

# All done!
print('Program completed successfully!')

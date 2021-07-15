#! python3
# get_cafe_family_ids_expand_contract.py
# Script to pull out the IDs of families that
# expanded/contracted, using the single best representative
# for the family as an output ID

import os, argparse
from Bio import SeqIO

# Define functions
def validate_args(args):
    # Validate that all arguments have been provided
    for key, value in vars(args).items():
        if value == None:
            print(key + ' argument was not specified. Fix this and try again.')
            quit()
    # Validate input file locations
    if not os.path.isfile(args.tabFileName):
        print('I am unable to locate the Base_change.tab file (' + args.tabFileName + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isdir(args.sequencesDir):
        print('The specified directory does not exist (' + args.sequencesDir + ')')
        print('Make sure you\'ve typed the location correctly and try again.')
        quit()
    # Validate output file location
    if os.path.isfile("{0}.expanded.txt".format(args.outputPrefix)):
        print('File already exists at output location (' + "{0}.expanded.txt".format(args.outputPrefix) + ')')
        print('Make sure you specify a unique file name and try again.')
        quit()
    if os.path.isfile("{0}.contracted.txt".format(args.outputPrefix)):
        print('File already exists at output location (' + "{0}.contracted.txt".format(args.outputPrefix) + ')')
        print('Make sure you specify a unique file name and try again.')
        quit()

def parse_basechange_tab_for_change(tabFileName, columnID):
    families = []
    found = False
    with open(tabFileName, "r") as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n").split("\t")
            # Handle header line
            if sl[0] == "FamilyID":
                for colIndex in range(1, len(sl)):
                    if sl[colIndex] == columnID: 
                        found = True # Marker for a good columnID
                        break # Break here so colIndex will be stored for later use
                if found == False:
                    print("Failed to find ID {0} in Base_change.tab file".format(columnID))
                    print("Make sure you specify the right input; program cannot continue; exiting now.")
                    quit()
                continue
            # Handle other lines
            num = int(sl[colIndex])
            direction = "expanded" if num > 0 else "contracted"
            if num != 0:
                groupID = sl[0]
                # Update ID to be OrthoGroup ID if necessary
                if not groupID.startswith("OG"):
                    groupID = "OG{0}{1}".format("0"*(7-len(groupID)), groupID) # OrthoGroup IDs are formatted like "OG0000000"
                # Store OrthoGroup ID in list
                families.append([groupID, direction])
    return families

def get_representative_from_orthogroups(sequencesDir, families, change_id):
    RIVAL_PERCENT = 0.6 # If the rival seq is at least 60% of the length of best, it's good
    repID_and_direction = []
    for family_and_direction in families:
        family = family_and_direction[0]
        direction = family_and_direction[1]
        # If we only want the family ID, immediately handle this now to save time
        if change_id:
            repID_and_direction.append([family, direction])
            continue
        # Handle other scenario where the specific gene ID matters
        # Derive the file name
        fileName = os.path.join(sequencesDir, "{0}.fa".format(family))
        assert os.path.isfile(fileName)
        # Read file as record dict
        records = custom_to_dict(fileName) # This indexes by description which is more stable
        # Find best representative with two metrics: M start, longest length
        bestID = None
        for key, value in records.items():
            # Initialise best finder
            if bestID == None:
                bestID = key
                continue
            # Check this against the best sequence
            bestSeq = str(records[bestID].seq)
            bestLen = len(bestSeq)
            rivalSeq = str(records[key].seq)
            rivalLen = len(rivalSeq)
            ## > Scenario 1: Both start with M
            if bestSeq.startswith("M") and rivalSeq.startswith("M"):
                ## > Scenario 1.1: ... and rival is longer
                if rivalLen > bestLen:
                    bestID = key
                    continue
            ## > Scenario 2: Neither start with M
            elif not bestSeq.startswith("M") and not rivalSeq.startswith("M"):
                ## > Scenario 2.1: ... and rival is longer
                if rivalLen > bestLen:
                    bestID = key
                    continue
            ## > Scenario 3: Only rival starts with M
            elif not bestSeq.startswith("M") and rivalSeq.startswith("M"):
                ## > Scenario 3.1: ... and current best is only a little longer
                if rivalLen >= bestLen * RIVAL_PERCENT: # Handicap the best seq
                    bestID = key
                    continue
            ## > Scenario 4: Only best starts with M
            elif bestSeq.startswith("M") and not rivalSeq.startswith("M"):
                ## > Scenario 4.1: ... and rival is only a little longer
                if bestLen <= rivalLen * RIVAL_PERCENT: # Handicap the rival seq
                    bestID = key
                    continue
        # Store the ID in our list
        repID_and_direction.append([bestID, direction])
    return repID_and_direction

def custom_to_dict(fastaFileName):
    outDict = {}
    with open(fastaFileName, "r") as fileIn:
        for record in SeqIO.parse(fileIn, "fasta"):
            outDict[record.description] = record
    return outDict

def main():
    # User input
    usage = """%(prog)s reads in the Base_change.tab file produced by CAFE
    and, provided the directory where OrthoFinder produced the individual orthogroup
    .fa FASTA files, will generate a FASTA output with a single sequence representing
    each orthogroup which experienced an expansion or contraction within the specified
    branch or node (-i). This branch or node should be represented by its ID as seen in
    the column header for the Base_change.tab file.
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-t", dest="tabFileName",
        help="Input Base_change.tab file name")
    p.add_argument("-d", dest="sequencesDir",
        help="Specify Orthogroup_Sequences directory produced by OrthoFinder")
    p.add_argument("-i", dest="id",
        help="Specify the branch/node ID to check for expansions/contractions")
    p.add_argument("-o", dest="outputPrefix",
        help="Output file prefix for the IDs")
    p.add_argument("--change_id", dest="change_id", action="store_true", default=False,
        help="Optionally produce a FASTA file where the representative ID is changed to be the Orthogroup ID")
    args = p.parse_args()
    validate_args(args)

    # Parse tab file
    families = parse_basechange_tab_for_change(args.tabFileName, args.id)

    # Extract representatives from family FASTA files
    repID_and_direction = get_representative_from_orthogroups(args.sequencesDir, families, args.change_id)

    # Produce output IDs for expanded families
    ## > Expanded
    with open("{0}.expanded.txt".format(args.outputPrefix), "w") as fileOut:
        for pair in repID_and_direction:
            family = pair[0]
            direction = pair[1]
            if direction == "expanded":
                fileOut.write("{0}\n".format(family))
    
    ## > Contracted
    with open("{0}.contracted.txt".format(args.outputPrefix), "w") as fileOut:
        for pair in repID_and_direction:
            family = pair[0]
            direction = pair[1]
            if direction == "contracted":
                fileOut.write("{0}\n".format(family))

if __name__ == "__main__":
    main()

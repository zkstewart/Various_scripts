#! python3
# get_cafe_family_representative.py
# Script to extract a single representative for each family
# identified as expanded/contracted within one or more lineages
# / nodes in the phylogenetic tree. It will produce a single
# FASTA file with sequences included.

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
    if os.path.isfile(args.outputFileName):
        print('File already exists at output location (' + args.outputFileName + ')')
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
            if num != 0:
                groupID = sl[0]
                # Update ID to be OrthoGroup ID if necessary
                if not groupID.startswith("OG"):
                    groupID = "OG{0}{1}".format("0"*(7-len(groupID)), groupID) # OrthoGroup IDs are formatted like "OG0000000"
                # Store OrthoGroup ID in list
                families.append(groupID)
    return families

def get_representative_from_orthogroups(sequencesDir, families):
    RIVAL_PERCENT = 0.6 # If the rival seq is at least 60% of the length of best, it's good
    repRecords = []
    for family in families:
        # Derive the file name
        fileName = os.path.join(sequencesDir, "{0}.fa".format(family))
        assert os.path.isfile(fileName)
        # Read file as record dict
        #try:
            #records = SeqIO.to_dict(SeqIO.parse(open(fileName, "r"), "fasta"))
        #except:
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
        # Change our sequence ID to be the family name
        records[bestID].description = ""
        records[bestID].id = family
        # Store the record in our list
        repRecords.append(records[bestID])
    return repRecords

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
    p.add_argument("-o", dest="outputFileName",
        help="Output file name for the FASTA")
    args = p.parse_args()
    validate_args(args)

    # Parse tab file
    families = parse_basechange_tab_for_change(args.tabFileName, args.id)

    # Extract representatives from family FASTA files
    repRecords = get_representative_from_orthogroups(args.sequencesDir, families)

    # Produce output FASTA
    with open(args.outputFileName, "w") as fileOut:
        for record in repRecords:
            fileOut.write(record.format("fasta"))

if __name__ == "__main__":
    main()

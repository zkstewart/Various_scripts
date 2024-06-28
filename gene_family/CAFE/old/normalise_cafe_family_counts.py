#! python3
# normalise_cafe_family_counts.py
# Script to assess the Base_count.tab file and
# "normalise" the values with respect to the predicted
# number of gene family members in the root node.

import os, argparse
from copy import deepcopy

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.tabFileName):
        print('I am unable to locate the Base_count.tab file (' + args.tabFileName + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print('File already exists at output location "{}"'.format(args.outputFileName))
        print('Make sure you specify a unique file name and try again.')
        quit()

def parse_basecount_tab(tabFileName):
    '''
    Simple function to parse a Base_count.tab file from CAFE into four
    outputs. Most of the outputs probably aren't relevant to what you need,
    but I don't abide by YAGNI so here, have all the possible ways you might
    want to parse this file.
    
    The nodeNames list comes with format:
        [node_1_name, node_2_name, ...]
    
    The familyNames list comes with format:
        [family_1_name, family_2_name, ...]
    
    The baseCountList list comes with format:
        [[node_1_count, node_2_count, ...], [...], ...]
    
    The baseCountDict dictionary comes with format:
        {node_name: node_int_count, ...}
    
    Returns:
        nodeNames -- a list containing each node name as in the header of the file.
        familyNames -- a list containing each family name as encountered
                       in file order.
        baseCountList -- a list containing each row's count values, with equivalent
                         order to familyNames list.
        baseCountDict -- a dict containing each node name as key, with the value
                         being a list of its counts ordered as in familyNames list.
    '''
    # nodeNames = [] # defined below
    familyNames = []
    baseCountList = []
    baseCountDict = {}
    with open(tabFileName, "r") as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n").split("\t")
            # Handle header line
            if sl[0] == "FamilyID":
                nodeNames = sl[1:]
                for name in nodeNames:
                    baseCountDict[name] = []
            # Handle other lines
            else:
                familyName, counts = sl[0], list(map(int, sl[1:]))
                familyNames.append(familyName)
                baseCountList.append(counts)
                for i in range(len(nodeNames)):
                    baseCountDict[nodeNames[i]].append(counts[i])
    return nodeNames, familyNames, baseCountList, baseCountDict

def normalise_baseCountList(baseCountList):
    '''
    Performs the operation (node_count - ancestral_count) for all count values
    in the provided baseCountList.
    
    Returns a new value without changing the original input.
    
    Parameters:
        baseCountList -- a list obtained through the parse_basecount_tab() function
                         containing rows of integers corresponding to gene counts.
                         It is assumed that the last value in each row is the ancestral
                         gene count.
    Returns:
        normalisedBaseCountList -- a similar list to the input parameter, but with
                                   the normalisation operation performed.
    '''
    normalisedBaseCountList = deepcopy(baseCountList) # don't change the input list
    
    for i in range(len(normalisedBaseCountList)):
        row = normalisedBaseCountList[i]
        for x in range(0, len(row)): # we'll also "normalised" the ancestral count itself to be 0
            row[x] = (row[x] - row[-1])
    return normalisedBaseCountList

def main():
    # User input
    usage = """%(prog)s reads in the Base_count.tab file produced by CAFE and
    normalises the number of members in each node with respect to the parent node.
    This helps to reveal how each node has changed over time with respect to the
    ancestral copy number of the gene family.
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-t", dest="tabFileName",
        help="Input Base_count.tab file name")
    p.add_argument("-o", dest="outputFileName",
        help="Output file name for the IDs")
    args = p.parse_args()
    validate_args(args)

    # Parse tab file
    nodeNames, familyNames, baseCountList, baseCountDict = parse_basecount_tab(args.tabFileName)

    # Normalise baseCountList
    normalisedBaseCountList = normalise_baseCountList(baseCountList)
    
    # Produce output file
    with open(args.outputFileName, "w") as fileOut:
        fileOut.write("FamilyID\t{0}\n".format("\t".join(nodeNames)))
        for i in range(len(familyNames)):
            name = familyNames[i]
            counts = normalisedBaseCountList[i]
            fileOut.write("{0}\t{1}\n".format(familyNames[i], "\t".join(map(str, counts))))

if __name__ == "__main__":
    main()

#! python3
# get_cafe_family_ids_expand_contract.py
# Script to pull out the IDs of families that
# expanded/contracted, using the single best representative
# for the family as an output ID

import os, argparse

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.tabFileName):
        print('I am unable to locate the Base_count.tab file (' + args.tabFileName + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.famResultsFileName):
        print('I am unable to locate the Base_family_results.tab file (' + args.famResultsFileName + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print('File already exists at output location ("{0}")'.format(args.outputFileName))
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
        {node_name: [node_int_count, ...]}
    
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

def parse_significant_familys(famResultsFileName):
    '''
    Function to parse the Base_family_results.txt file to generate a list of family
    names that were detected as being expanded in one or more lineages.
        
    Returns:
        significantFamilyNames -- a list containing strings of the family names that met
                                  CAFE's 0.05 significance threshold.
    '''    
    significantFamilyNames = []
    with open(famResultsFileName, "r") as fileIn:
        for line in fileIn:
            if line.startswith("#"):
                continue
            
            name, p, isSignificant = line.rstrip("\r\n").split("\t")
            if isSignificant == "y":
                significantFamilyNames.append(name)
    return significantFamilyNames


def get_expanded_family_ids(speciesNames, familyNames, significantFamilyNames, baseCountList):
    '''
    Performs some mathematical calculations to determine whether a family
    is expanded for each species.
    
    Parameters:
        familyNames -- a list containing ordered string values indicating the family
                       names.
        significantFamilyNames -- a list containing string values indicating the family
                                  names that were detected as being expanded in one or
                                  more lineages.
        baseCountList -- a list containing the integers from each row of the Base_count.tab
                         file as derived from parse_basecount_tab().
    Returns:
        familyIDsDict -- a dictionary with structure:
                            {
                                species_name: [
                                    [species_1_familyID_a, species_1_familyID_b],
                                    [species_2_familyID_c, ...],
                                    ...]
                                ]
                            }
    '''
    familyIDsDict = {sp: [] for sp in speciesNames if not sp.startswith("<")} # starting with < means it's an internal node, not tip
    for x in range(len(baseCountList)):
        # Check if this is significant or not
        familyName = familyNames[x]
        if familyName not in significantFamilyNames:
            continue
        # Get row values
        row = baseCountList[x][:-1] # drop the ancestral value
        # Index expanded IDs in their relevant list
        for i in range(len(speciesNames)):
            sp = speciesNames[i]
            if sp.startswith("<"): # skip internal nodes
                continue
            spValue = row[i]
            if spValue == max(row): # this allows tied values to be accepted as expanded in both lineages
                familyIDsDict[sp].append(familyName)
    return familyIDsDict

def main():
    # User input
    usage = """%(prog)s reads in the Base_change.tab file and Base_family_results.txt
    files produced by CAFE and produces an output file indicating which families have expanded
    and the lineage/species in which they expanded.
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i1", dest="tabFileName", required=True,
        help="Input Base_count.tab file name")
    p.add_argument("-i2", dest="famResultsFileName", required=True,
        help="Input Base_family_results.txt file name")
    p.add_argument("-o", dest="outputFileName",
        help="Output file name where tabular output will be written")
    args = p.parse_args()
    validate_args(args)
    
    # Parse Base_count.tab file
    nodeNames, familyNames, baseCountList, baseCountDict = parse_basecount_tab(args.tabFileName)
    
    # Parse Base_family_results.txt file
    significantFamilyNames = parse_significant_familys(args.famResultsFileName)
    
    # Get list of families that expanded for each species/tip node
    familyIDsDict = get_expanded_family_ids(nodeNames, familyNames, significantFamilyNames, baseCountList)
    
    # Produce output file for expanded families
    with open(args.outputFileName, "w") as fileOut:
        for speciesID, familyNames in familyIDsDict.items():
            speciesID = speciesID.split("<")[0] # remove the node identifier
            for family in familyNames:
                fileOut.write("{0}\t{1}\n".format(speciesID, family))
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

#! python3
# cafe_enrichment_analysis.py
# Script to assess the Base_count.tab file
# and performs an enrichment analysis for
# each tip node i.e., each species.

import os, argparse, statistics

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.tabFileName):
        print('I am unable to locate the Base_count.tab file (' + args.tabFileName + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.orthofinderFileName):
        print('I am unable to locate the OrthoFinder.tsv file (' + args.orthofinderFileName + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.annotationFileName):
        print('I am unable to locate the annotation file (' + args.annotationFileName + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate output file location
    if os.path.isfile("{0}.self.txt".format(args.outputPrefix)):
        print('File already exists at output location ("{0}.self.txt")'.format(args.outputPrefix))
        print('Make sure you specify a unique file name and try again.')
        quit()
    if os.path.isfile("{0}.others.txt".format(args.outputPrefix)):
        print('File already exists at output location ("{0}.others.txt")'.format(args.outputPrefix))
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

def parse_orthogroups_tsv(orthogroupsFile):
    '''
    Simple function to parse a Orthogroups.tsv file from OrthoFinder into
    a dictionary structure from which all relevant data can be obtained.
    
    The dictionary comes with format:
        {
            species_name: {
                orthogroup_name: [
                    gene_1,
                    gene_2,
                    ...
                ]
            }
        }
    
    Returns:
        orthoDict -- a dictionary with nested indexing, where the parent key corresponds
                     to the species name, and nested key is the orthogroup name which points
                     to the list of gene names.
    '''
    # species = [] # defined below
    orthoDict = {}
    firstLine = True
    with open(orthogroupsFile , "r") as orthoFile:
        for line in orthoFile:
            sl = line.rstrip("\r\n").replace('"', '').split("\t") # this won't support sequence IDs containing quotation marks
            # Handle the header line
            if firstLine == True:
                species = sl[1:] # drop the "Orthogroup" column text
                for sp in species:
                    orthoDict[sp] = {}
                firstLine = False
            # Handle all other lines
            else:
                    orthogroupName, genes = sl[0], sl[1:]
                    for i in range(len(species)):
                        sp = species[i]
                        spGenes = genes[i].split(", ") # orthofinder uses ", " as its gene name delimiter
                        spGenes = [] if spGenes == [''] else spGenes
                        orthoDict[sp][orthogroupName] = spGenes
    return orthoDict

def parse_annotation_file(annotationFileName):
    '''
    Simply parses a goseq style annotation table file
    to return a dictionary associating each gene
    to its annotation terms (GOs, PFAMs, whatever).
    
    The goseq style is defined by "; " separators between
    annotation terms, and "0" for values that are empty.
    
    Returns:
        annotDict -- a dictionary with structure like:
                        gene_name: [term1, term2, ...]
    '''
    annotDict = {}
    with open(annotationFileName, "r") as fileIn:
        for line in fileIn:
            model, terms = line.rstrip("\r\n").split("\t")
            if terms == "0":
                continue
            else:
                annotDict[model] = terms.split("; ")
    return annotDict

def merge_orthodict_and_annotdict(orthoDict, annotDict):
    '''
    Parameters:
        orthoDict -- a dictionary with nested indexing, where the parent key corresponds
                     to the species name, and nested key is the orthogroup name which points
                     to the list of gene names.
        annotDict -- a dictionary with structure like:
                        gene_name: [term1, term2, ...]
    Returns:
        orthoAnnotDict -- a dictionary with structure:
                            {
                                species_name: {
                                    orthogroup_name: [
                                        term_1,
                                        term_2,
                                        ...
                                    ]
                                }
                            }
    '''
    orthoAnnotDict = {}
    for _, groupDict in orthoDict.items(): # don't care about species name here
        for groupName, genes in groupDict.items():
            orthoAnnotDict.setdefault(groupName, set())
            for g in genes:
                if g in annotDict:
                    gAnnot = annotDict[g]
                    orthoAnnotDict[groupName] = orthoAnnotDict[groupName].union(set(gAnnot))
            orthoAnnotDict[groupName] = list(orthoAnnotDict[groupName]) # turn itself into list
    return orthoAnnotDict

def get_expanded_family_ids(speciesNames, familyNames, baseCountList):
    '''
    Performs some mathematical calculations involving median and st.dev
    to determine whether a family is expanded for each species.
    
    Parameters:
        familyNames -- a list containing ordered string values indicating the family
                       names.
        baseCountList -- a list containing the integers from each row of the Base_count.tab
                         file as derived from parse_basecount_tab().
    Returns:
        familyIDsList -- a list with structure:
                            [[species_1_familyID_a, species_1_familyID_b], [species_2_familyID_c, ...], ...]
    '''
    familyIDsList = [[] for sp in speciesNames if not sp.startswith("<")] # starting with < means it's an internal node, not tip
    for x in range(len(baseCountList)):
        # Get row values
        row = baseCountList[x][:-1] # drop the ancestral value
        median = statistics.median(row)
        stdev = statistics.stdev(row)
        # Drop outliers so we can re-calculate a better st.dev
        "We're going to define an outlier as a number that is not within a st.dev of any other number"
        outlierDropRow = []
        for j in range(len(row)):
            isWithinStdev = False
            for n in range(len(row)):
                if j == n:
                    continue
                elif (row[j]-stdev) <= row[n] and (row[j]+stdev) >= row[n]: # if row[n] "overlaps" this range
                    isWithinStdev = True
                    break
            if isWithinStdev:
                outlierDropRow.append(row[j])
        stdev = statistics.stdev(outlierDropRow) # re-calculate narrower st.dev value
        # Index expanded IDs in their relevant list
        ongoingCount = 0
        for i in range(len(speciesNames)):
            sp = speciesNames[i]
            if sp.startswith("<"): # skip internal nodes
                continue
            spValue = row[i]
            if spValue >= int(median + stdev) or spValue == max(row): # int() will "rescue" borderline values
                familyIDsList[ongoingCount].append(familyNames[x])
            ongoingCount += 1
    return familyIDsList

def main():
    # User input
    usage = """%(prog)s reads in 1) the Base_count.tab file produced by CAFE, 2) the 
    Orthogroups.tsv file produced by OrthoFinder, 3) the Base_branch_probabilities.tab
    (IMPORTANT: Make sure the column names
    in these 3 files match up). It will also require a final input i.e., the annotations
    file which should be formatted with lines "sequenceID\tterm1; term2; ...".
    
    A chi-square results table will be produced indicating the annotation terms which
    are overrepresented and underrepresented in the expanded families for each species.
    This will occur in two ways. First, the families expanded within a species will be compared
    to the families that did NOT expand within that same species. Second, the families expanded
    within a species will be compared to the families that expanded in OTHER species.
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i1", dest="tabFileName",
        help="Input Base_count.tab file name")
    p.add_argument("-i2", dest="orthofinderFileName",
        help="Input OrthoFinder.tsv file name")
    p.add_argument("-i3", dest="annotationFileName",
        help="Input annotation file name")
    p.add_argument("-o", dest="outputFilePrefix",
        help="Output file prefix for the chi-square results")
    args = p.parse_args()
    validate_args(args)
    
    # Parse Base_count.tab file
    nodeNames, familyNames, baseCountList, baseCountDict = parse_basecount_tab(args.tabFileName)
    
    # Parse OrthoFinder.tsv file
    orthoDict = parse_orthogroups_tsv(args.orthofinderFileName)
    
    # Parse annotation file
    annotDict = parse_annotation_file(args.annotationFileName)
    
    # Get annotations associated with each orthogroup
    orthoAnnotDict = merge_orthodict_and_annotdict(orthoDict, annotDict)
    
    # Get list of families that expanded for each species/tip node
    familyIDsList = get_expanded_family_ids(nodeNames, familyNames, baseCountList)
    
    print("Program completed successfully!")
#! python3
# cafe_enrichment_analysis.py
# Script to assess the Base_count.tab file
# and performs an enrichment analysis for
# each tip node i.e., each species.

import os, argparse, statistics
import numpy as np
from scipy.stats import chi2_contingency
from multipy.fdr import lsu

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
    if not os.path.isfile(args.orthofinderFileName):
        print('I am unable to locate the OrthoFinder.tsv file (' + args.orthofinderFileName + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.annotationFileName):
        print('I am unable to locate the annotation file (' + args.annotationFileName + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate output file location
    print("As part of argument validation, I can't be certain output files can be written")
    print("This program WILL NOT overwrite files, so make sure you're putting output somewhere fresh")

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
    Function to take the orthoDict from parse_orthogroups_tsv() and annotDict
    from parse_annotation_file() and smash them together to get a new dictionary
    which associates orthogroups (keys; string) to annotation terms (values; list).
    
    Parameters:
        orthoDict -- a dictionary with nested indexing, where the parent key corresponds
                     to the species name, and nested key is the orthogroup name which points
                     to the list of gene names.
        annotDict -- a dictionary with structure like:
                        gene_name: [term1, term2, ...]
    Returns:
        orthoAnnotDict -- a dictionary with structure:
                            {
                                orthogroup_name: [
                                    term_1,
                                    term_2,
                                    ...
                                ]
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
    for key in orthoAnnotDict.keys():
        orthoAnnotDict[key] = list(orthoAnnotDict[key])
    return orthoAnnotDict

def get_expanded_family_ids(speciesNames, familyNames, significantFamilyNames, baseCountList):
    '''
    Performs some mathematical calculations involving median and st.dev
    to determine whether a family is expanded for each species.
    
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
        # median = statistics.median(row)
        # stdev = statistics.stdev(row)
        # # Drop outliers so we can re-calculate a better st.dev
        # "We're going to define an outlier as a number that is not within a st.dev of any other number"
        # outlierDropRow = []
        # for j in range(len(row)):
        #     isWithinStdev = False
        #     for n in range(len(row)):
        #         if j == n:
        #             continue
        #         elif (row[j]-stdev) <= row[n] and (row[j]+stdev) >= row[n]: # if row[n] "overlaps" this range
        #             isWithinStdev = True
        #             break
        #     if isWithinStdev:
        #         outlierDropRow.append(row[j])
        # stdev = statistics.stdev(outlierDropRow) # re-calculate narrower st.dev value
        # Index expanded IDs in their relevant list
        for i in range(len(speciesNames)):
            sp = speciesNames[i]
            if sp.startswith("<"): # skip internal nodes
                continue
            spValue = row[i]
            # if spValue >= int(median + stdev) or spValue == max(row): # int() will "rescue" borderline values
            if spValue == max(row): # this allows tied values to be accepted as expanded in both lineages
                familyIDsDict[sp].append(familyName)
    return familyIDsDict

def enrichment_analysis(orthoAnnotDict, specialIDs, remainingIDs):
    '''
    Performs an enrichment analysis using chi-squared testing. Intended
    to be used as part of a project involving OrthoFinder and CAFE.
    The specialIDs and remainingIDs list can be redundant. If multiple
    IDs show up up, we'll count them each and every time. As such,
    make sure you want this behaviour to occur! If not, use list(set())
    on the data first!
    
    Downstream of this method, you MUST perform some sort of post-testing
    P-value correction like FDR correction.
    
    Note that this statistical method requires relatively large numbers
    of data points (e.g., greater than 5 in each group). Otherwise,
    we'd need to employ Fisher's exact test which isn't implemented here.
    
    Parameters:
        orthoAnnotDict -- a dictionary with structure:
                            {
                                orthogroup_name: [
                                    term_1,
                                    term_2,
                                    ...
                                ]
                            }
        specialIDs -- a list containing strings which map to keys in
                      orthoAnnotDict which are expanded/contracted/expressed/
                      "special" in some way.
        remainingIDs -- a list containing strings which map to keys in
                        orthoAnnotDict which are not "special" and hence are
                        the remaining IDs from an initially full list which
                        were separated into "special"/"non-special"
    Returns:
        chi2Results -- a list with format:
                            [[term, p_value, over_or_under_represented], ...]
                       ... where the over_or_under_represented value is a string
                       equal to "over" when the term is enriched, or "under" when
                       it is depleted relative to what we'd expect in a balanced
                       scenario. ALL results are provided regardless of their
                       significance.
    '''
    # Tally term presence in special versus remaining IDs sets
    specialCounts = {}
    remainingCounts = {}
    uniqueTerms = set()
    for groupName, annotTerms in orthoAnnotDict.items():
        if groupName in specialIDs:
            numIters = specialIDs.count(groupName)
            for _ in range(numIters):
                for term in annotTerms:
                    if term == '': continue
                    specialCounts.setdefault(term, 0)
                    specialCounts[term] += 1
                    uniqueTerms.add(term)
        if groupName in remainingIDs:
            numIters = remainingIDs.count(groupName)
            for _ in range(numIters):
                for term in annotTerms:
                    if term == '': continue
                    remainingCounts.setdefault(term, 0)
                    remainingCounts[term] += 1
                    uniqueTerms.add(term)
    
    # For each term, produce a contingency table and assess its enrichment
    chi2Results = []
    for term in uniqueTerms:
        # Get counts for the table
        specialOccurring = specialCounts[term] if term in specialCounts else 0
        specialNotOccurring = len(specialIDs) - specialOccurring
        remainingOccurring = remainingCounts[term] if term in remainingCounts else 0
        remainingNotOccurring = len(remainingIDs) - remainingOccurring
        
        # Format the 2x2 contingency table and run chi-squared analysis
        contingencyTable = [
            [specialOccurring, specialNotOccurring],
            [remainingOccurring, remainingNotOccurring]
        ]
        chi2, p, dof, expected = chi2_contingency(contingencyTable)
        
        # Figure out if it's over or underrepresented
        representation = "over" if expected[0][0] < specialOccurring else "under" if expected[0][0] > specialOccurring else "neither"
        
        # Store result
        chi2Results.append([term, p, representation])
    return chi2Results

def FDR_correction(testResults, ALPHA=0.05):
    '''
    Function to receive statistical testing results e.g., from enrichment_analysis()
    and determine which results are significant at the provided cut-off value.
    
    Parameters:
        testResults -- a list with format:
                            [[term, p_value, over_or_under_represented], ...]
                       ... where the over_or_under_represented value is a string
                       equal to "over" when the term is enriched, or "under" when
                       it is depleted relative to what we'd expect in a balanced
                       scenario.
        ALPHA -- a float value indicating what P-value to use as a significance threshold.
    Returns:
        significantResults -- the same list as the input parameter, but excluding
                              any values that are not significant after FDR correction.
    '''
    # Perform FDR correction
    significantPvals = lsu(np.array([p for _, p, _ in testResults]), q=ALPHA)
    
    # Retain only significant results
    significantResults = []
    for x in range(len(testResults)):
        if significantPvals[x] == True:
            significantResults.append(testResults[x])
    
    return significantResults

def write_chi2_results(chi2Results, outputFileName):
    if os.path.isfile(outputFileName):
        print("I refuse to overwrite '{0}'".format(outputFileName))
    else:
        chi2Results.sort(key = lambda x: (x[2], x[1])) # sort results for better interpretation
        with open(outputFileName, "w") as fileOut:
            fileOut.write("term\tP_value\trepresented\n")
            for result in chi2Results:
                fileOut.write("{0}\n".format("\t".join(map(str, result))))

def main():
    # User input
    usage = """%(prog)s reads in 1) the Base_count.tab file produced by CAFE, 2) the
    Base_family_results.txt file produced by CAFE, 3) the Orthogroups.tsv file produced
    by OrthoFinder, (IMPORTANT: Make sure the column names in these 3 files match up), and
    4) the annotations file which should be formatted with lines "sequenceID\tterm1;
    term2; ...".
    
    A chi-square results table will be produced indicating the annotation terms which
    are overrepresented and underrepresented in the expanded families for each species.
    This will occur in two ways. First, the families expanded within a species will be compared
    to the families that did NOT expand within that same species. Second, the families expanded
    within a species will be compared to the families that expanded in OTHER species.
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i1", dest="tabFileName", required=True,
        help="Input Base_count.tab file name")
    p.add_argument("-i2", dest="famResultsFileName", required=True,
        help="Input Base_family_results.txt file name")
    p.add_argument("-i3", dest="orthofinderFileName", required=True,
        help="Input OrthoFinder.tsv file name")
    p.add_argument("-i4", dest="annotationFileName", required=True,
        help="Input annotation file name")
    p.add_argument("-o", dest="outputFilePrefix", required=True,
        help="Output file prefix for the chi-square results")
    args = p.parse_args()
    validate_args(args)
    
    # Parse Base_count.tab file
    nodeNames, familyNames, baseCountList, baseCountDict = parse_basecount_tab(args.tabFileName)
    
    # Parse Base_family_results.txt file
    significantFamilyNames = parse_significant_familys(args.famResultsFileName)
    
    # Parse OrthoFinder.tsv file
    orthoDict = parse_orthogroups_tsv(args.orthofinderFileName)
    
    # Parse annotation file
    annotDict = parse_annotation_file(args.annotationFileName)
    
    # Get annotations associated with each orthogroup
    _orthoAnnotDict = merge_orthodict_and_annotdict(orthoDict, annotDict)
    
    # Cull orthogroups that weren't assessed by CAFE
    """
    This is important since CAFE will only assess gene families when all but one (?) species
    has at least one member present. Not culling the other orthogroups will bias our
    enrichment analysis in ways that are hard for me to predict.
    """
    orthoAnnotDict = {}
    for name in familyNames:
        orthoAnnotDict[name] = _orthoAnnotDict[name]
    
    # Get list of families that expanded for each species/tip node
    familyIDsDict = get_expanded_family_ids(nodeNames, familyNames, significantFamilyNames, baseCountList)
    
    # First test: self-expansion
    selfResults = {}
    for i in range(len(nodeNames)):
        sp = nodeNames[i]
        if sp.startswith("<"): # skip internal nodes which aren't special/tips
            continue
        
        # Get the IDs for this comparison
        specialIDs = set(familyIDsDict[sp])
        remainingIDs = set(orthoAnnotDict.keys()).difference(specialIDs)
        
        # Prevent biases by dropping orthogroups without terms
        """
        I'm not entirely sure of the implications, but I believe there's reason to suspect that
        our results will become biased if we keep orthogroups that lack annotations in the
        analysis. My reasoning is that smaller orthogroups are more likely to have not been
        studied, and hence they won't have received the same attention when it comes to
        annotation. These smaller groups, if they do expand in a species, will not be detected
        in our analysis as a positive result. However, they will bias our results to detect
        the more understood groups better due to how the chi-squared test works (i.e., the less
        understood orthogroups will all be counted into the notOccurring rows/columns). 
        """
        specialIDs = [id for id in specialIDs if orthoAnnotDict[id] != []]
        remainingIDs = [id for id in remainingIDs if orthoAnnotDict[id] != []]
        
        # Run the enrichment analysis
        chi2Results = enrichment_analysis(orthoAnnotDict, list(specialIDs), list(remainingIDs))
        selfResults[sp] = chi2Results
    
    # Second test: exclusive expansion
    exclusiveResults = {}
    for i in range(len(nodeNames)):
        sp = nodeNames[i]
        if sp.startswith("<"): # skip internal nodes which aren't special/tips
            continue
        
        # Get the IDs for this comparison
        specialIDs = set(familyIDsDict[sp])
        remainingIDs = []
        _ongoingCount = 0
        for x in range(len(nodeNames)):
            _sp = nodeNames[x]
            if _sp.startswith("<"): # skip internal nodes which aren't special/tips
                continue
            elif _sp == sp: # skip self comparison
                _ongoingCount += 1
                continue
            """
            We want to count orthogroups multiple times if they show up as expanded in more
            than one lineage. This will "counter-weight" the likelihood of any terms in that
            orthogroup being detected as enriched if the orthogroup's expansion is not itself
            unique to the lineage that our specialIDs are being drawn from.
            """
            remainingIDs += familyIDsDict[_sp]
            _ongoingCount += 1
        
        # Prevent biases by dropping orthogroups without terms
        specialIDs = [id for id in specialIDs if orthoAnnotDict[id] != []]
        remainingIDs = [id for id in remainingIDs if orthoAnnotDict[id] != []]
        
        # Run the enrichment analysis
        chi2Results = enrichment_analysis(orthoAnnotDict, list(specialIDs), list(remainingIDs))
        exclusiveResults[sp] = chi2Results
    
    # Perform FDR correction to obtain significant results from all comparisons
    correctedSelfResults = {}
    for i in range(len(nodeNames)):
        sp = nodeNames[i]
        if sp.startswith("<"):
            continue
        
        spSelfResult = selfResults[sp]
        correctedSpSelfResult = FDR_correction(spSelfResult)
        correctedSelfResults[sp] = correctedSpSelfResult
    
    correctedExclusiveResults = {}
    for i in range(len(nodeNames)):
        sp = nodeNames[i]
        if sp.startswith("<"):
            continue
        
        spExclusiveResult = exclusiveResults[sp]
        correctedSpExclusiveResult = FDR_correction(spExclusiveResult)
        correctedExclusiveResults[sp] = correctedSpExclusiveResult
    
    # Write outputs to file
    for i in range(len(nodeNames)):
        sp = nodeNames[i]
        if sp.startswith("<"):
            continue
        
        # Find the results again
        _selfResult = correctedSelfResults[sp]
        _exclusiveResult = correctedExclusiveResults[sp]
        
        # Get a nice file name and write it!
        sp = sp.split("<")[0] # drop the node identifier now
        write_chi2_results(_selfResult, "{0}.{1}.self.tsv".format(args.outputFilePrefix, sp))
        write_chi2_results(_exclusiveResult, "{0}.{1}.exclusive.tsv".format(args.outputFilePrefix, sp))
        
    print("Program completed successfully!")
    
if __name__ == "__main__":
    main()

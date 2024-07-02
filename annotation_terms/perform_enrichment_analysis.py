#! python3
# perform_enrichment_analysis.py
# Script to perform an enrichment analysis on an annotation
# TSV file using a set of orthogroups/gene IDs/families that
# are of interest.

import os, argparse
import numpy as np
from scipy.stats import fisher_exact
from multipy.fdr import lsu

from merge_annotation_tsvs import parse_annotation_file

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.annotationTSV):
        print('I am unable to locate the annotation TSV file (' + args.annotationTSV + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.selectedIDsFile):
        print('I am unable to locate the selected IDs file (' + args.selectedIDsFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if args.subsetFile != None:
        if not os.path.isfile(args.subsetFile):
            print('I am unable to locate the subset file (' + args.subsetFile + ')')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    # Validate numeric arguments
    if args.minimumOccurrence < 1:
        print("The minimum number of annotation term occurrences must be 1 or greater.")
        quit()
    # Handle file overwrites
    if os.path.exists(args.outputFileName):
        print(f"{args.outputDirectory} already exists, and I will not overwrite it.")
        print("Delete/move whatever exists here, or specify a different output name when you try again.")
        quit()

def parse_subset_ids(subsetFile):
    '''
    Function to parse a text file or TSV file containing IDs to retain.
    
    Parameters:
        subsetFile -- a string representing the path to the text or
                      TSV file containing the IDs to retain
    Returns:
        subsetIDs -- a set containing the IDs to retain
    '''
    subsetIDs = set()
    with open(subsetFile, "r") as fileIn:
        firstLine = True
        for line in fileIn:
            sl = line.strip("\r\n ").split("\t")
            if firstLine:
                firstLine = False
                if len(sl) > 1:
                    continue
            if sl[0] != "":
                subsetIDs.add(sl[0])
    return subsetIDs

def parse_selected_ids(selectedIDsFile):
    '''
    Function to parse a text file containing IDs of interest.
    
    Parameters:
        selectedIDsFile -- a string representing the path to the text file
                           containing the IDs of interest
    Returns:
        specialIDs -- a set containing the IDs of interest
    '''
    specialIDs = set()
    with open(selectedIDsFile, "r") as fileIn:
        for line in fileIn:
            l = line.strip("\r\n ")
            if l != "":
                specialIDs.add(l)
    return specialIDs

def enrichment_analysis(annotDict, selectedIDs):
    '''
    Performs an enrichment analysis using Fisher's exact testing. Intended
    to be used as part of a project involving OrthoFinder and CAFE, but
    it can be used for any annotation data.
    
    Downstream of this method, you MUST perform some sort of post-testing
    P-value correction like FDR correction.
    
    Parameters:
        annotDict -- a dictionary with structure:
                            {
                                "key1": set([
                                    "term_1",
                                    "term_2",
                                    ...
                                ]),
                                "key2": set([ ... ]),
                                ...
                            }
        selectedIDs -- a list containing strings which map to keys in
                       annotDict which are expanded/contracted/expressed/
                      "special" in some way.
    Returns:
        enrichmentResults -- a list with format:
                       [[term, p_value, selected_occurring, over_or_under_represented], ...]
                       where the over_or_under_represented value is a string
                       equal to "over" when the term is enriched, or "under" when
                       it is depleted relative to what we'd expect in a balanced
                       scenario. The p_value can be used for downstream FDR correction.
                       The selected_occurring value can be used to filter results to e.g.,
                       only consider over-representation/enrichment in terms that occur some
                       minimum number of times.
    '''
    # Drop all keys that have no annotations
    annotDict = { key: value for key, value in annotDict.items() if value != set() }
    assert len(annotDict) > 0, "No annotation keys were retained after dropping those without annotations."
    
    # Make sure that we still have some keys left to test
    selectedIDs = set(annotDict.keys()).intersection(selectedIDs)
    assert len(selectedIDs) > 0, "No selected keys were retained after dropping those without annotations."
    
    # Identify the remaining IDs
    remainingIDs = set(annotDict.keys()).difference(selectedIDs)
    
    # Get all unique annotation terms
    annotTerms = set().union(*annotDict.values())
    
    # Establish dictionary of terms and their presence in selected versus remaining IDs
    presenceDict = { x: [0, 0] for x in annotTerms } # [selected, remaining]
    
    # Tally term presence in selected versus remaining IDs sets
    for key, terms in annotDict.items():
        for term in terms:
            presenceDict[term][0] += 1 if key in selectedIDs else 0
            presenceDict[term][1] += 1
    
    # For each term, produce a contingency table and assess its enrichment
    testResults = []
    for term in annotTerms:
        # Get counts for the table
        selectedOccurring = presenceDict[term][0]
        selectedNotOccurring = len(selectedIDs) - selectedOccurring
        
        remainingOccurring = presenceDict[term][1]
        remainingNotOccurring = len(remainingIDs) - remainingOccurring
        
        # Format the 2x2 contingency table
        contingencyTable = [
            [selectedOccurring, selectedNotOccurring],
            [remainingOccurring, remainingNotOccurring]
        ]
        
        # Run Fisher's exact test
        statistic, p = fisher_exact(contingencyTable)
        
        # Figure out if it's over or underrepresented
        selectedRatio = selectedOccurring / len(selectedIDs)
        remainingRatio = remainingOccurring / len(remainingIDs)
        
        if selectedRatio > remainingRatio:
            representation = "over"
        else:
            representation = "under"
        
        # Store result
        testResults.append([term, p, selectedOccurring, representation])
    
    return testResults

def FDR_correction(testResults, pvalueIndex, ALPHA=0.05):
    '''
    Function to receive statistical testing results e.g., from enrichment_analysis()
    and determine which results are significant at the provided cut-off value.
    
    Parameters:
        testResults -- a list of sublists containing results of a statistical test
                       where at least one value is the P-value outcome.
        pvalueIndex -- an integer indicating which 0-based index in the sublists of
                       testResults contains the P-value.
        ALPHA -- a float value indicating what P-value to use as a significance threshold.
    Returns:
        significantResults -- the same list as the input parameter, but excluding
                              any values that are not significant after FDR correction.
    '''
    # Perform FDR correction
    significantPvals = lsu(np.array([p[pvalueIndex] for p in testResults]), q=ALPHA)
    
    # Retain only significant results
    significantResults = []
    for x in range(len(testResults)):
        if significantPvals[x] == True:
            significantResults.append(testResults[x])
    
    return significantResults

def main():
    # User input
    usage = """%(prog)s receives an annotation TSV pairing keys (e.g., orthogroups, gene families, DEGs, etc.)
    to annotation terms and a text file listing the keys of interest. It then performs an enrichment analysis
    on the annotation terms associated with the keys of interest versus the remaining keys. The output will be
    a TSV report on which annnotations were enriched or depleted in the keys of interest.
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="annotationTSV",
                   required=True,
                   help="Specify annotation TSV file")
    p.add_argument("-s", dest="selectedIDsFile",
                   required=True,
                   help="Specify the text file listing the IDs of interest")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file name for combined annotations")
    # Opts (behavioural)
    p.add_argument("--subset", dest="subsetFile",
                   required=False,
                   help="""Optionally, specify a text file listing the keys in the annotation TSV
                   to retain (all others will be dropped); this will also accept a TSV file with
                   the keys in the first column with all other columns being ignored (TSV file is
                   assumed to contain a header row but text file is not - make note of this!)""",
                   default=None)
    p.add_argument("--minimum", dest="minimumOccurrence",
                   required=False,
                   type=int,
                   help="""Optionally, specify how many times an annotation term must
                   occur in the selected keys to be considered for enrichment analysis;
                   default == 1""",
                   default=1)
    # Opts (file format)
    p.add_argument("--delimiter", dest="annotationDelimiter",
                   required=False,
                   help="""Optionally, specify the delimiter used in the 
                   annotation TSV file separating terms; default == '; '""",
                   default="; ")
    p.add_argument("--blank", dest="blankCharacter",
                   required=False,
                   help="""Optionally, specify the character used in the 
                   annotation TSV file to indicate no annotations; default == '0'""",
                   default="0")
    p.add_argument("--hasHeader", dest="hasHeader",
                   required=False,
                   action="store_true",
                   help="Optionally indicate if the annotation TSV file has a header",
                   default=False)
    args = p.parse_args()
    validate_args(args)
    
    # Parse the annotation TSV file
    annotDict = parse_annotation_file(args.annotationTSV, args.annotationDelimiter,
                                      args.blankCharacter, args.hasHeader)
    
    # Parse and use the subset file (if applicable)
    if args.subsetFile != None:
        subsetIDs = parse_subset_ids(args.subsetFile)
        annotDict = { key: value for key, value in annotDict.items() if key in subsetIDs }
        assert len(annotDict) > 0, \
            "No keys were retained from the subset file. Check that the file matches your annotation TSV."
    
    # Parse the selected IDs file
    selectedIDs = parse_selected_ids(args.selectedIDsFile)
    
    # Run the enrichment analysis
    enrichmentResults = enrichment_analysis(annotDict, selectedIDs)
    
    # Perform FDR correction to obtain significant results from all comparisons
    correctedResults = FDR_correction(enrichmentResults, 1)
    
    # Drop any results that don't meet the minimum occurrence threshold
    correctedResults = [
        x
        for x in correctedResults
        if x[3] == "under" or x[2] >= args.minimumOccurrence
    ]
    
    # Write outputs to file
    with open(args.outputFileName, "w") as fileOut:
        fileOut.write("term\tpvalue\tnum_occurrences\trepresentation\n")
        for result in correctedResults:
            fileOut.write("\t".join(map(str, result)) + "\n")
    
    print("Program completed successfully!")
    
if __name__ == "__main__":
    main()

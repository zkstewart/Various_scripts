#! python3
# plot_go_proportions.py
# Allows for the plotting of the proportion of genes associated with each
# GO term in a given annotation file

import os, argparse, sys, pickle, re
import numpy as np
import matplotlib.pyplot as plt
from goatools import obo_parser

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Function_packages import ZS_GO

colourPalette = ['#377eb8', '#ff7f00', '#4daf4a',
                 '#f781bf', '#a65628', '#984ea3',
                 '#999999', '#e41a1c', '#dede00']
categoryOrder = ["cellular_component", "molecular_function", "biological_process"]
goRegex = re.compile(r"^GO:\d{7}$")

# Define functions
def validate_args(args):
    # Validate input file locations
    for inputFile in args.inputFiles:
        if not os.path.isfile(inputFile):
            raise FileNotFoundError(f"I am unable to locate the -i file '{inputFile}'")
    if not os.path.isfile(args.goOboFile):
        raise FileNotFoundError(f"I am unable to locate the -g file '{args.goOboFile}'")
    if args.goIDsFile is not None and not os.path.isfile(args.goIDsFile):
        raise FileNotFoundError(f"I am unable to locate the --goIDsFile '{args.goIDsFile}'")
    
    # Validate that number of input files is supported
    if len(args.inputFiles) > len(colourPalette):
        raise ValueError(f"Number of input files ({len(args.inputFiles)}) exceeds the number " + 
                         f"of supported colours ({len(colourPalette)})")
    
    # Validate delimiter
    if args.columnDelimiter == "":
        raise ValueError("--columnDelimiter must not be empty")
    if args.annotationDelimiter == "":
        raise ValueError("--annotationDelimiter must not be empty")
    
    # Validate numeric arguments
    if args.minimumProportion < 0.0 or args.minimumProportion > 1.0:
        raise ValueError("--minimumProportion must be between 0.0 (==0%) and 1.0 (==100%)")
    if args.minimumNumber < 0:
        raise ValueError("--minimumNumber must be an integer >= 0")
    if args.width <= 0.0:
        raise ValueError("--width must be a positive float")
    if args.height <= 0.0:
        raise ValueError("--height must be a positive float")
    
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        raise FileExistsError(f"File already exists at -o location '{args.outputFileName}'")
    if not args.outputFileName.endswith(".png") and not args.outputFileName.endswith(".pdf"):
        raise ValueError("-o file file must end with '.png' or '.pdf'")

def parse_annotation_file(inputFile, columnDelimiter, annotationDelimiter, hasHeader):
    '''
    Parameters:
        inputFile -- a string pointing to the location of the input file
        columnDelimiter -- a string specifying the column columnDelimiter
        annotationDelimiter -- a string specifying the annotation delimiter
        hasHeader -- a boolean specifying whether the input file has a header row
    Returns:
        annotDict -- a dictionary with structure like:
                     {
                         "GO:0000001": { "geneID1", "geneID2", ... },
                         "GO:0000002": { "geneID2", "geneID4", ... },
                         ...
                     }
    '''
    annotDict = {}
    with open(inputFile, "r") as fileIn:
        skipHeader = hasHeader
        for line in fileIn:
            if skipHeader:
                skipHeader = False
                print(f"# Skipping header row: {line}")
                continue
            
            sl = line.rstrip("\r\n ").split(columnDelimiter)
            if len(sl) != 2:
                raise ValueError(f"Input file '{inputFile}' does not have 2 columns")
            
            geneID, annotations = sl
            for annot in annotations.split(annotationDelimiter):
                annotDict.setdefault(annot, set())
                annotDict[annot].add(geneID)
    return annotDict

def plot_annotation_proportions(annotDict, numGenesDict, outputFile, goDetails, width=12, height=8):
    '''
    Parameters:
        annotDict -- a dictionary with structure like:
                     {
                         "GO:0000001": { "geneID1", "geneID2", ... },
                         "GO:0000002": { "geneID2", "geneID4", ... },
                         ...
                     }
       numGenesDict -- a dictionary with structure like:
                       {
                           "filePrefix1": 100,
                           "filePrefix2": 200,
                           ...
                       }
         outputFile -- a string specifying the location of the output file
         goDetails -- a dictionary with structure like:
                      {
                          "GO:0000001": ["Name of GO 1", "Category of GO 1"],
                          "GO:0000002": ["Name of GO 2", "Category of GO 2"],
                          ...
                      }
    '''
    # Get consistently ordered annot IDs and sample names
    annotIDs = list(annotDict.keys())
    filePrefixes = list(numGenesDict.keys())
    
    # Order annotIDs according to category
    annotIDs.sort(key = lambda x: (categoryOrder.index(goDetails[x][1]), goDetails[x][0]) )
    categories = []
    for goID in annotIDs:
        category = goDetails[goID][1]
        if category not in categories:
            categories.append(category)
    
    # Calculate the bar width
    barWidth = 0.9 / len(numGenesDict) # leave padding between bars
    
    # Plot the data
    fig, (ax, axCategories) = plt.subplots(nrows=2, ncols=1,
                                           tight_layout=True,
                                           figsize=(width, height),
                                           height_ratios=[50, 1])
    ax2 = ax.twinx()
    
    categoryBoundaries = [ [ c, np.inf, 0] for c in categories ]
    for i, annot in enumerate(annotIDs):
        numAnnotatedDict = annotDict[annot]
        
        # Store the boundaries of the secondary X-axis
        thisCategory = goDetails[annot][1]
        thisIndex = [ z for z, boundary in enumerate(categoryBoundaries) if boundary[0] == thisCategory ][0]
        if i < categoryBoundaries[thisIndex][1]:
            categoryBoundaries[thisIndex][1] = i
        
        # Plot bar for each filePrefix
        for j, filePrefix in enumerate(filePrefixes):
            numAnnotated = numAnnotatedDict[filePrefix]
            
            # Get our x and y values
            x = i + (barWidth*j)
            y1 = (numAnnotated / numGenesDict[filePrefix]) * 100 # percentage
            y2 = numAnnotated
                        
            # Plot bar
            ax.bar(x, y1, color = colourPalette[j], width = barWidth,
                   label = filePrefix if i == 0 else None)
            ax2.bar(x, y2, color = colourPalette[j], width = barWidth)
            
            # Store the maximum value for the secondary X-axis
            categoryBoundaries[thisIndex][2] = x
    
    # Add labels
    textBoxes = ax.set_xticks(np.arange(len(annotDict)) + (barWidth / len(numGenesDict)),
                  labels=[ goDetails[x][0] for x in annotIDs ], rotation = 65, ha = 'right')
    ax.set_ylabel("Percentage of genes", fontweight = 'bold', fontsize = 12)
    ax2.set_ylabel("Number of genes", fontweight = 'bold', fontsize = 12)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=len(numGenesDict))
    
    # Add category boundary lines to bottom subplot
    categoryCentres = []
    for category, start, end in categoryBoundaries:
        axCategories.plot([start, end], [0.05, 0.05], color = 'black', lw=2,
                          marker = '|')
        categoryCentres.append((start + end) / 2)
    
    # Add category labels to bottom subplot
    for centre, category in zip(categoryCentres, categories):
        axCategories.text(centre, 0.1, category, ha = 'center', va = 'bottom', fontweight = 'bold', fontsize = 12)
    
    # Set axis limits of bottom subplot
    axCategories.set_xlim(ax.get_xlim())
    axCategories.set_ylim(0, 0.1)
    axCategories.axis('off')
    
    # Write plot to file
    fig.savefig(outputFile)#, bbox_inches="tight")

def main():
    usage = """%(prog)s receives a table pairing gene identifiers (left column)
    to annotation terms such as GOs (right column). It then plots the proportion
    of genes associated with each term.
    """
    
    p = argparse.ArgumentParser(description=usage)
    # Required
    p.add_argument("-i", dest="inputFiles",
                   required=True,
                   nargs="+",
                   help="Specify the location of one or more input annotation file(s)")
    p.add_argument("-g", dest="goOboFile",
                   required=True,
                   help="Specify the location of the go.obo file")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="""Specify the name for the output file; make sure to use
                   a suffix of .png or .pdf to determine the output format""")
    # Optional (input)
    p.add_argument("--goIDsFile", dest="goIDsFile",
                   required=False,
                   help="""Optionally, specify a text file containing GO IDs (one
                   per line) to include in the output plot; this will override
                   any other filtering options""")
    # Optional (plotting)
    p.add_argument("--height", dest="height",
                   type=float,
                   required=False,
                   help="""Optionally, specify the plot height in inches (default: 8.0)""",
                   default=8.0)
    p.add_argument("--width", dest="width",
                   type=float,
                   required=False,
                   help="""Optionally, specify the plot width in inches (default: 12.0)""",
                   default=12.0)
    p.add_argument("--minimumProportion", dest="minimumProportion",
                   type=float,
                   required=False,
                   help="""Optionally, specify the proportion of sequences that an
                   annotation must occur within on any given input file to be
                   included in the output plot (default: 0.0; no filtering done)""",
                   default=0.0)
    p.add_argument("--minimumNumber", dest="minimumNumber",
                   type=int,
                   required=False,
                   help="""Optionally, specify the number of of sequences that an
                   annotation must occur within on any given input file to be
                   included in the output plot (default: 0; no filtering done)""",
                   default=0)
    # Optional (file parsing)
    p.add_argument("--columnDelimiter", dest="columnDelimiter",
                   required=False,
                   help="""Optionally, specify the column delimiter; default is
                   tab-delimited ('\\t') but you can specify any value such
                   as a comma""",
                   default="\t")
    p.add_argument("--annotationDelimiter", dest="annotationDelimiter",
                   required=False,
                   help="""Optionally, specify the annotation delimiter; default is
                   comma-delimited (',') but you can specify any value such
                   as '; '""",
                   default=",")
    p.add_argument("--header", dest="hasHeader",
                   required=False,
                   action="store_true",
                   help="""Specify this flag if your input file has a header
                   row which should be skipped""",
                   default=False)
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse .obo file
    "Parse early to error out if format is incorrect"
    goObo = obo_parser.GODag(args.goOboFile)
    
    # Parse goIDsFile (if it exists)
    if args.goIDsFile is not None:
        goIDs = set()
        with open(args.goIDsFile, "r") as fileIn:
            for line in fileIn:
                goID = line.rstrip("\r\n ")
                if not goRegex.match(goID):
                    raise ValueError(f"GO ID '{goID}' does not match expected format")
                goIDs.add(goID)
    
    # Parse annotation files
    annotDicts = [
        [
            inputFile.rsplit(".", maxsplit=1)[0],
            parse_annotation_file(inputFile, args.columnDelimiter, args.annotationDelimiter, args.hasHeader)
        ]
        for inputFile in args.inputFiles
    ]
    
    # Combine and reformat annotDicts
    combinedAnnotDict = {}
    numGenesDict = {}
    for filePrefix, annotDict in annotDicts:
        numGenes = len( set().union(*annotDict.values()) )
        numGenesDict[filePrefix] = numGenes
        
        for annot, geneSet in annotDict.items():
            combinedAnnotDict.setdefault(annot, { v[0]: 0 for v in annotDicts })
            combinedAnnotDict[annot][filePrefix] = len(geneSet)
    
    # Filter annotations
    if args.goIDsFile is not None:
        filteredAnnotDict = {}
        for goID in goIDs:
            if goID in combinedAnnotDict:
                filteredAnnotDict[goID] = combinedAnnotDict[goID]
            else:
                raise KeyError(f"GO ID '{goID}' in '{args.goIDsFile}' not found in any -i file")
        
        filteredAnnotDict = {
            annot: annotDict
            for annot, annotDict in combinedAnnotDict.items()
            if annot in goIDs
        }
    else:
        filteredAnnotDict = {
            annot: annotDict
            for annot, annotDict in combinedAnnotDict.items()
            if
                any( numAnnotated >= args.minimumNumber for numAnnotated in annotDict.values() ) and
                any( numAnnotated / numGenesDict[filePrefix] >= args.minimumProportion for filePrefix, numAnnotated in annotDict.items() )
        }
    
    # Parse existing pickled data (if it exists)
    outputDirectory = os.path.dirname(os.path.abspath(args.outputFileName))
    pickleFileName = os.path.join(outputDirectory, ".plot_go_proportions.pickle")
    if os.path.isfile(pickleFileName):
        with open(pickleFileName, "rb") as fileIn:
            apiQueryDict = pickle.load(fileIn)
    else:
        apiQueryDict = {}
    
    # Resolve GO names and categories
    resolvedAnnotDict = {}
    goDetails = {}
    for goID, annotValue in filteredAnnotDict.items():
        # If we can find this GO ID already, use it
        if goID in goObo:
            newGOID = goID
        elif goID in apiQueryDict:
            newGOID = apiQueryDict[goID]
            
            if newGOID is None:
                print(f"WARNING: I am resuming a run and found that GO ID '{goID}' does not have a suitable replacement")
                continue
            elif not newGOID in goObo:
                print(f"WARNING: I am resuming a run and found that GO ID '{goID}' has a replacement '{newGOID}' that is not in your go.obo")
                continue
        # Try to locate a suitable replacement
        else:
            result = ZS_GO.query_go_api(goID)
            if result is None:
                apiQueryDict[goID] = None
                print(f"WARNING: GO ID '{goID}' not found in your go.obo with a suitable replacement not found")
                continue
            else:
                newGOID = result[0]
                apiQueryDict[goID] = newGOID
                if newGOID in goObo:
                    print(f"NOTE: GO ID '{goID}' not found in your go.obo; using '{newGOID}' as a replacement")
                else:
                    print(f"WARNING: GO ID '{goID}' not found in your go.obo with replacement '{newGOID}' also not found")
                    continue
        
        # Store the GO details
        goName = goObo[newGOID].name
        goCategory = goObo[newGOID].namespace
        goDetails[newGOID] = [goName, goCategory]
        resolvedAnnotDict[newGOID] = annotValue
    
    # Store the API query results
    with open(pickleFileName, "wb") as fileOut:
        pickle.dump(apiQueryDict, fileOut)
    
    # Combine any GO terms that, through replacement, end up being the same
    toCombine = []
    for goID, (goName, goCategory) in goDetails.items():
        occurrences = [
            otherGoID
            for otherGoID, (otherGoName, otherGoCategory) in goDetails.items()
            if goName == otherGoName and goCategory == otherGoCategory
        ]
        if len(occurrences) > 1:
            toCombine.append(occurrences)
    dontInclude = set([y for x in toCombine for y in x])
    
    finalAnnotDict = {
        goID: annotValue
        for goID, annotValue in resolvedAnnotDict.items()
        if goID not in dontInclude
    }
    for idsToCombine in toCombine:
        newGOID = idsToCombine[0]
        newDictValue = { filePrefix: 0 for filePrefix in numGenesDict.keys() }
        for goID in idsToCombine:
            for filePrefix, numAnnotated in resolvedAnnotDict[goID].items():
                newDictValue[filePrefix] += numAnnotated
        
        finalAnnotDict[newGOID] = newDictValue
    
    # Plot the data
    plot_annotation_proportions(finalAnnotDict, numGenesDict, args.outputFileName,
                                goDetails, args.width, args.height)
    
    # Done!
    print('Program completed successfully!')

if __name__ == '__main__':
    main()

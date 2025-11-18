#! python3
# tabulate_kasp_markers.py
# Reads in the outputs of 'kasp_check.py' to tabulate a report
# for easier interpretation of results

import os, argparse
import numpy as np
import pandas as pd
import scipy.spatial.distance as ssd
from scipy.cluster import hierarchy

class DirectoryNotFoundError(Exception):
    pass

def validate_args(args):
    # Validate input files
    args.markersDir = os.path.abspath(args.markersDir)
    if not os.path.isdir(args.markersDir):
        raise DirectoryNotFoundError(f"-i markers directory '{args.markersDir}' is not a directory!")
    
    # Locate kasp_check.py results
    args.kaspFiles = {}
    for _file in os.listdir(args.markersDir):
        fullPath = os.path.join(args.markersDir, _file)
        if os.path.isfile(fullPath) and _file.endswith(args.fileSuffix):
            # Get the sample name from the file prefix
            filePrefix = _file.rsplit(args.fileSuffix, maxsplit=1)[0]
            if args.isPaired:
                if filePrefix.endswith("1") or filePrefix.endswith("2"):
                    filePrefix = filePrefix[:-1].rstrip("_.-")
                else:
                    raise ValueError(f"The prefix for '{_file}' i.e., '{filePrefix}' should end with '1' or '2' since --paired was set")
            
            # Store the files
            args.kaspFiles.setdefault(filePrefix, [])
            args.kaspFiles[filePrefix].append(fullPath)
            args.kaspFiles[filePrefix].sort(key = lambda x: x[-1])
    
    # Validate output file location
    args.outputFileName = os.path.abspath(args.outputFileName)
    if os.path.isfile(args.outputFileName):
        raise FileExistsError(f"-o file '{args.outputFileName}' already exists!")

def tabulate_markers(kaspFiles):
    '''
    Parameters:
        kaspFiles -- a dictionary where keys are sample IDs and values are lists
                     with one or two files to check for marker occurrence
    Returns:
        markersDf -- a pandas DataFrame where rows are sample IDs and columns
                     are markers. Cell values are "." for absence, and
                     "Y" for marker presence.
    '''
    # Parse each file
    resultsDict = {}
    for sampleID, files in kaspFiles.items():
        resultsDict[sampleID] = {}
        for _file in files:
            with open(_file, "r") as fileIn:
                for line in fileIn:
                    l = line.strip()
                    if l != "":
                        resultsDict[sampleID][l] = 1
    
    # Find all unique markers that occurred
    markers = sorted(set([ key for values in resultsDict.values() for key in values.keys() ]))
    
    # Add in missing values
    for sampleID, markersDict in resultsDict.items():
        for marker in markers:
            if not marker in markersDict:
                markersDict[marker] = 0
    
    # Cluster samples by marker occurrence similarity
    ## See https://stackoverflow.com/questions/52612841/how-to-cluster-features-based-on-their-correlations-to-each-other-with-sklearn
    clusterDf = pd.DataFrame.from_dict(resultsDict)
    corrDf = clusterDf.corr(lambda x, y: sum([ 1 for i in range(0, x.size) if x[i] == y[i] ]) / x.size )
    distances = 1 - corrDf.abs().values  # pairwise distnces
    distArray = ssd.squareform(distances)  # scipy converts matrix to 1d array
    hier = hierarchy.linkage(distArray, method="ward")  # you can use other methods
    cluster_labels = hierarchy.fcluster(hier, 0, criterion="distance")
    
    # Get the ordered sample labels
    sampleSorting = sorted(zip(corrDf.columns.tolist(), cluster_labels), key = lambda x: x[1])
    samplesSorted = [ sampleID for sampleID, clustNum in sampleSorting ]
    
    # Convert to sorted dataframe for presentation
    markersDf = pd.DataFrame.from_dict(resultsDict, orient="index")
    markersDf = markersDf.reindex(samplesSorted)
    
    return markersDf

def main():
    usage = """%(prog)s provides a follow-up step to kasp_check.py after it has been run on
    multiple samples. The result files can be tabulated for an at-a-glance interpretation of
    which markers were found in each sample.
    """
    # Establish main parser
    p = argparse.ArgumentParser(description=usage)
    
    # Set arguments shared by subparsers
    p.add_argument("-i", dest="markersDir",
                   required=True,
                   help="Specify the directory containing kasp_check.py outputs")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="""Specify the location to write results to""")
    p.add_argument("--suffix", dest="fileSuffix",
                   required=False,
                   help="""Optionally specify the suffix of the kasp result files;
                   default == '.kasp.tsv'""",
                   default=".kasp.tsv")
    p.add_argument("--paired", dest="isPaired",
                   required=False,
                   action="store_true",
                   help="""Provide this flag if each read in a pair was checked; this means
                   we expect to find files with '1$SUFFIX' and '2$SUFFIX' and a common prefix""",
                   default=False)
    
    args = p.parse_args()
    validate_args(args) # sets args.kaspFiles
    
    # Combine marker results
    markersDf = tabulate_markers(args.kaspFiles)
    
    # Output results
    with open(args.outputFileName, "w") as fileOut: 
        markersDf.to_csv(fileOut, sep="\t", index=True)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

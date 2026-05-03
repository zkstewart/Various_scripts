#!/usr/bin/env python3
# tabulate_deseq2_results.py

import os, argparse
import pandas as pd
from statistics import median

class FileExistsError(Exception):
    pass

def validate_args(args):
    # Validate DESeq2 file locations depending on type of input
    args.files = []
    for inputLocation in args.inputLocations:
        inputLocation = os.path.abspath(inputLocation)
        
        # Handle directories
        if os.path.isdir(inputLocation):
            foundFiles = False
            for f in os.listdir(inputLocation):
                file = os.path.join(inputLocation, f)
                if os.path.isfile(file) and file.endswith(args.inputSuffix):
                    args.files.append(file)
                    foundFiles = True
            if not foundFiles:
                raise FileNotFoundError(f"'{inputLocation}' does not contain any files.")
        
        # Handle files
        elif os.path.isfile(inputLocation): # don't need to check suffix here
            args.files.append(inputLocation)
        else:
            raise FileNotFoundError(f"'{inputLocation}' is not a valid file or directory.")
    
    # Validate coldata and normalised count file locations
    args.coldataFileName = os.path.abspath(args.coldataFileName)
    if not os.path.isfile(args.coldataFileName):
        raise FileNotFoundError(f"-c '{args.outputFileName}' is not a file or does not exist!")
    
    args.normalisedCountsFileName = os.path.abspath(args.normalisedCountsFileName)
    if not os.path.isfile(args.normalisedCountsFileName):
        raise FileNotFoundError(f"-n '{args.normalisedCountsFileName}' is not a file or does not exist!")
    
    # Validate output file location
    args.outputFileName = os.path.abspath(args.outputFileName)
    if os.path.exists(args.outputFileName):
        raise FileExistsError(f"The -o output location '{args.outputFileName}' already exists and will not be overwritten!")
    
    parentDir = os.path.dirname(args.outputFileName)
    if not os.path.isdir(parentDir):
        raise ValueError(f"The -o output location '{args.outputFileName}' points to a non-existent parent directory '{parentDir}'")

def is_floatable(x):
    try:
        float(x)
        return True
    except:
        return False

def parse_coldata_file(fileName, groupingVariables):
    sampleGroups = {}
    with open(fileName, "r") as fileIn:
        firstLine = True
        for line in fileIn:
            sl = line.rstrip().split("\t")
            if firstLine:
                header = sl
                groupIndices = []
                for variable in groupingVariables:
                    try:
                        groupIndices.append(header.index(variable))
                    except:
                        raise ValueError(f"-g value '{variable}' does not occur in the coldata header ({header})")
                firstLine = False
            else:
                groupParts = tuple(( sl[i] for i in groupIndices ))
                
                sampleGroups.setdefault(groupParts, set())
                sampleGroups[groupParts].add(sl[0]) # first column is expected to be the sample
    
    # Figure out how we can sort the group parts in a logical way
    sortingTypes = []
    for parts in zip(*sampleGroups.keys()):
        if all([ p.isdigit() for p in parts ]):
            sortingTypes.append(int)
        elif all([ is_floatable(p) for p in parts ]):
            sortingTypes.append(float)
        else:
            sortingTypes.append(str)
    
    # Produce sorted groups that are not broken up into parts
    keys = list(sampleGroups.keys())
    keys.sort(key = lambda x: [sortingTypes[i](y) for i, y in enumerate(x)])
    
    # Produce a modified dict with ordered keys
    sortedSampleGroups = {}
    for groupParts in keys:
        group = "-".join(groupParts)
        sortedSampleGroups[group] = sampleGroups[groupParts]
    
    return sortedSampleGroups

def parse_deseq2_file(fileName, fileSuffix):
    IGNORE = ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue"]
    
    # Derive the comparison name
    comparison = os.path.basename(fileName).rsplit(fileSuffix, maxsplit=1)[0]
    
    # Parse the file
    d = {}
    with open(fileName, "r") as fileIn:
        firstLine = True
        for line in fileIn:
            sl = line.rstrip().split("\t")
            if firstLine:
                header = sl
                keys = [ (x, i) for i, x in enumerate(header) if ( i != 0 ) and ( not x in IGNORE ) ]
                firstLine = False
            else:
                geneID = sl[0]
                d[geneID] = { x:sl[i] for x, i in keys }
    return d, comparison

def parse_normalised_counts(fileName, sampleGroups, toParse=None):
    # Identify all samples that may be found within the counts file
    samples = set()
    for sgroup in sampleGroups.values():
        samples.update(sgroup)
    
    # Parse counts file into a dictionary for later DataFrame generation
    counts = {} # for later conversion into df
    with open(fileName, "r") as fileIn:
        firstLine = True
        for line in fileIn:
            sl = line.rstrip().split("\t")
            if firstLine:
                header = sl
                sampleIndices = {
                    s:header.index(s)
                    for s in samples
                    if s in header
                }
                firstLine = False
            else:
                geneID = sl[0]
                counts[geneID] = {}
                
                # Skip this gene if we don't need to parse it
                if toParse != None and (not geneID in toParse):
                    continue
                
                # Obtain and store the median of this group's expression
                for group, sgroup in sampleGroups.items():
                    groupCounts = []
                    for s in sgroup:
                        if s in sampleIndices:
                            groupCounts.append(float(sl[sampleIndices[s]]))
                    counts[geneID][group] = median(groupCounts)
    
    # Convert dict into pandas DataFrame
    countsDf = pd.DataFrame.from_dict(counts, orient="index")
    
    return countsDf

def main():
    ##### USER INPUT SECTION
    usage = """%(prog)s will modify DESeq2 output files to include their best BLAST
    name as found within a standard ZKS format annotation table (e.g., as produced by
    BINge). Existing files will NOT be overwritten, and warnings will be emitted where
    that is enforced.
    """
    # Required arguments
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="inputLocations",
                   required=True,
                   nargs="+",
                   help="Location of files or directories for DESeq2 results files to collate")
    p.add_argument("-c", dest="coldataFileName",
                   required=True,
                   help="Location of the coldata file used by DESeq2")
    p.add_argument("-g", dest="groupingVariables",
                   required=True,
                   nargs="+",
                   help="One or more coldata columns to identify biological replicates")
    p.add_argument("-n", dest="normalisedCountsFileName",
                   required=True,
                   help="Location of the normalised counts produced by DESeq2")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Location to write the output TSV file")
    # Optional arguments
    p.add_argument("--is", dest="inputSuffix",
                   required=False,
                   nargs="+",
                   help="""Optionally, specify the file ending of the DESeq2 output files
                   if a directory is given to -i (default='.edit.tsv')""",
                   default=".edit.tsv")
    p.add_argument("--is", dest="inputSuffix",
                   required=False,
                   help="""Optionally, specify the file ending of the DESeq2 output files
                   if a directory is given to -i (default='.edit.tsv')""",
                   default=".edit.tsv")
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse coldata to derive the way to group samples for expression averaging
    sampleGroups = parse_coldata_file(args.coldataFileName, args.groupingVariables)
    
    # Parse each DGE result file
    genes = {}
    comparisons = {}
    for fileName in args.files:
        genesDict, comparison = parse_deseq2_file(fileName, args.inputSuffix)
        genes.update({ g: {} for g in genesDict.keys() })
        comparisons[comparison] = genesDict
    
    # Format a dictionary for df conversion
    for comparison, genesDict in comparisons.items():
        for gene in genes.keys():
            if gene in genesDict:
                geneDetails = genesDict[gene]
                
                # Handle P-value
                genes[gene][comparison] = ["padj"]
                
                # Handle all non-Pvalue keys (expected to be redundant)
                for key, value in geneDetails.items():
                    if key == "padj":
                        genes[gene][comparison] = geneDetails["padj"]
                    else:
                        genes[gene][key] = value
            else:
                genes[gene][comparison] = "."
    df = pd.DataFrame.from_dict(genes, orient="index")
    
    # Parse the normalised counts
    countsDf = parse_normalised_counts(args.normalisedCountsFileName, sampleGroups, set(genes.keys()))
    
    # Join the DGE and counts dataframes together
    df = df.join(countsDf)
    
    # Write to file
    df.to_csv(args.outputFileName, sep="\t")
    
    # Notify user of successful completion
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

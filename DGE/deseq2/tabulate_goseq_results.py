#!/usr/bin/env python3
# tabulate_goseq_results.py

import os, argparse
import pandas as pd
from statistics import median
from goatools import obo_parser

class FileExistsError(Exception):
    pass

def validate_args(args):
    # Validate goseq file locations depending on type of input
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
    
    # Validate go.obo file location
    args.oboFileName = os.path.abspath(args.oboFileName)
    if not os.path.isfile(args.oboFileName):
        raise FileNotFoundError(f"-g '{args.oboFileName}' is not a file or does not exist!")
    
    # Validate output file location
    args.outputFileName = os.path.abspath(args.outputFileName)
    if os.path.exists(args.outputFileName):
        raise FileExistsError(f"The -o output location '{args.outputFileName}' already exists and will not be overwritten!")
    
    parentDir = os.path.dirname(args.outputFileName)
    if not os.path.isdir(parentDir):
        raise ValueError(f"The -o output location '{args.outputFileName}' points to a non-existent parent directory '{parentDir}'")

def parse_goseq_file(fileName, fileSuffix):
    # Derive the comparison name
    comparison = os.path.basename(fileName).rsplit(fileSuffix, maxsplit=1)[0]
    
    # Parse the file
    resultsDict = {}
    termDict = {}
    with open(fileName, "r") as fileIn:
        firstLine = True
        for line in fileIn:
            sl = line.rstrip().split("\t")
            if firstLine:
                firstLine = False
            else:
                category, over, under, numDE, num, term, ontology = sl
                direction = "over" if float(over) < float(under) else "under"
                resultsDict[category] = direction
                termDict[category] = term
    return resultsDict, termDict, comparison

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

def parse_table_to_dict(fileName, indexColumn):
    df = pd.read_csv(fileName, sep="\t")
    header = list(df.keys()[1:])
    tableDict = {
        k: list(v.values()) 
        for k, v 
        in df.set_index(indexColumn).to_dict('index').items()
    }
    return tableDict, header

def parse_candidates(fileName, valueColumn):
    tableDict = {}
    df = pd.read_csv(fileName, sep="\t")
    for index, row in df.iterrows():
        tableDict.setdefault(row["gene_id"], set())
        tableDict[row["gene_id"]].add(row[valueColumn])
    return tableDict

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
    p.add_argument("-g", dest="oboFileName",
                   required=True,
                   help="Location of the go.obo file")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Location to write the output TSV file")
    # Optional arguments
    p.add_argument("--is", dest="inputSuffix",
                   required=False,
                   help="""Optionally, specify the file ending of the DESeq2 output files
                   if a directory is given to -i (default='.edit.tsv')""",
                   default=".edit.tsv")
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse .obo file
    goObo = obo_parser.GODag(args.oboFileName)
    
    # Parse each goseq result file
    genes = {}
    comparisons = {}
    terms = {}
    for fileName in args.files:
        resultDict, termDict, comparison = parse_goseq_file(fileName, args.inputSuffix)
        genes.update({ g: {} for g in resultDict.keys() })
        terms.update(termDict)
        comparisons[comparison] = resultDict
    
    # Add any missing (NA) terms from the .obo
    for key in terms.keys():
        value = terms[key]
        if value == "NA":
            try:
                newValue = goObo[key].name
                terms[key] = newValue
            except:
                print(f"# Unable to resolve name for '{key}' which is currently specified as 'NA'")
    
    # Format a dictionary for df conversion
    for comparison, resultDict in comparisons.items():
        for gene in genes.keys():
            if gene in resultDict:
                direction = resultDict[gene]
                genes[gene][comparison] = direction
            else:
                genes[gene][comparison] = "."
    df = pd.DataFrame.from_dict(genes, orient="index")
    df.insert(0, "description", [ terms[geneID] for geneID in df.index ])
    
    # Write to file
    df.to_csv(args.outputFileName, sep="\t")
    
    # Notify user of successful completion
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

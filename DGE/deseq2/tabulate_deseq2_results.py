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
    
    # Validate optional WGCNA directory location
    if args.wgcnaDir != None:
        args.wgcnaDir = os.path.abspath(args.wgcnaDir)
        if not os.path.isdir(args.wgcnaDir):
            raise FileNotFoundError(f"--wgcna '{args.wgcnaDir}' is not a directory or does not exist!")
    
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

def make_sample_rcompatible(sample):
    sample = sample.replace("-", ".") # R will change hyphens to dots internally
    if sample[0].isdigit():
        sample = "X" + sample # R will also add an 'X' prior to a sample identifier that starts with a digit
    return sample

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
                sampleGroups[groupParts].add(make_sample_rcompatible(sl[0])) # first column is expected to be the sample
    
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
                   help="""Optionally, specify the file ending of the DESeq2 output files
                   if a directory is given to -i (default='.edit.tsv')""",
                   default=".edit.tsv")
    p.add_argument("--wgcna", dest="wgcnaDir",
                   required=False,
                   help="""Optionally, specify the directory where WGCNA files exist; if
                   unspecified, these results will not be included""",
                   default=None)
    
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
    
    # Parse the WGCNA files (if applicable)
    if args.wgcnaDir != None:
        wgcna = {}
        
        # Parse each file
        for file in os.listdir(args.wgcnaDir):
            # Parse gene_colours files
            if file.startswith("gene_colours") and file.endswith(".tsv"):
                suffix = file.split("gene_colours", maxsplit=1)[-1].split(".tsv")[0].strip(".")
                gcDict, _ = parse_table_to_dict(os.path.join(args.wgcnaDir, file), "gene_id")
                
                for geneID, colour in gcDict.items():
                    wgcna.setdefault(geneID, {})
                    wgcna[geneID].setdefault("colour", {})
                    wgcna[geneID]["colour"][suffix] = colour[0]
            
            # Parse gene_significance files
            if file.startswith("gene_significance") and file.endswith(".tsv"):
                suffix = file.split("gene_significance", maxsplit=1)[-1].split(".tsv")[0].strip(".")
                gsDict, _ = parse_table_to_dict(os.path.join(args.wgcnaDir, file), "Unnamed: 0")
                gsDict = { k:max(v) for k,v in gsDict.items() }
                
                for geneID, gs in gsDict.items():
                    wgcna.setdefault(geneID, {})
                    wgcna[geneID].setdefault("gene_significance", {})
                    wgcna[geneID]["gene_significance"][suffix] = gs
            
            # Parse module_kme files
            if file.startswith("module_kme") and file.endswith(".tsv"):
                suffix = file.split("module_kme", maxsplit=1)[-1].split(".tsv")[0].strip(".")
                kmeDict, kmeHeader = parse_table_to_dict(os.path.join(args.wgcnaDir, file), "Unnamed: 0")
                
                for geneID, kme in kmeDict.items():
                    orderedKme = sorted([ (kme[i], kmeHeader[i]) for i in range(len(kmeHeader)) ], key = lambda x: -x[0])
                    bestKme, bestModule = orderedKme[0]
                    
                    wgcna.setdefault(geneID, {})
                    wgcna[geneID].setdefault("colour_by_kme", {})
                    wgcna[geneID]["colour_by_kme"][suffix] = bestModule.split("MM.")[-1]
            
            # Parse central_candidates files
            if file.startswith("central_candidates") and file.endswith(".tsv"):
                suffix = file.split("central_candidates", maxsplit=1)[-1].split(".tsv")[0].strip(".")
                hubDict = parse_candidates(os.path.join(args.wgcnaDir, file), "trait")
                
                for geneID, traits in hubDict.items():
                    wgcna.setdefault(geneID, {})
                    wgcna[geneID].setdefault("central_candidates", {})
                    wgcna[geneID]["central_candidates"][suffix] = "; ".join(traits)
            
            # Parse network_screening_candidates files
            if file.startswith("network_screening_candidates") and file.endswith(".tsv"):
                suffix = file.split("network_screening_candidates", maxsplit=1)[-1].split(".tsv")[0].strip(".")
                networkDict = parse_candidates(os.path.join(args.wgcnaDir, file), "trait")
                
                for geneID, traits in networkDict.items():
                    wgcna.setdefault(geneID, {})
                    wgcna[geneID].setdefault("network_screening_candidates", {})
                    wgcna[geneID]["network_screening_candidates"][suffix] = "; ".join(traits)

        # Format into a dataframe-amenable dictionary
        dummyKey = list(wgcna.keys())[0]
        dummyDict = wgcna[dummyKey]
        wgcnaDfDict = {}
        for geneID in df.index:
            # Obtain the data for this gene (or use a dummy value for iteration)
            try:
                geneDict = wgcna[geneID]
                isDummy = False
            except:
                geneDict = dummyDict
                isDummy = True
            wgcnaDfDict[geneID] = {}
            
            # Store gene_colours value
            for suffix, value in geneDict["colour"].items():
                key = "colour" + f".{suffix}" if suffix != "" else ""
                
                if not isDummy:
                    wgcnaDfDict[geneID][key] = value
                else:
                    wgcnaDfDict[geneID][key] = "."

            # Store gene_significance value
            for suffix, value in geneDict["gene_significance"].items():
                key = "gene_significance" + f".{suffix}" if suffix != "" else ""
                
                if not isDummy:
                    wgcnaDfDict[geneID][key] = value
                else:
                    wgcnaDfDict[geneID][key] = "."
            
            # Store colour_by_kme value
            for suffix, value in geneDict["colour_by_kme"].items():
                key = "colour_by_kme" + f".{suffix}" if suffix != "" else ""
                
                if not isDummy:
                    wgcnaDfDict[geneID][key] = value
                else:
                    wgcnaDfDict[geneID][key] = "."
            
            # Store central_candidates value
            if "central_candidates" in geneDict:
                for suffix, value in geneDict["central_candidates"].items():
                    key = "central_candidates" + f".{suffix}" if suffix != "" else ""
                    
                    if not isDummy:
                        wgcnaDfDict[geneID][key] = value
                    else:
                        wgcnaDfDict[geneID][key] = "."
            else:
                for suffix in geneDict["colour"].keys():
                    key = "central_candidates" + f".{suffix}" if suffix != "" else ""
                    wgcnaDfDict[geneID][key] = "."
            
            # Store network_screening_candidates value
            if "network_screening_candidates" in geneDict:
                for suffix, value in geneDict["network_screening_candidates"].items():
                    key = "network_screening_candidates" + f".{suffix}" if suffix != "" else ""
                    
                    if not isDummy:
                        wgcnaDfDict[geneID][key] = value
                    else:
                        wgcnaDfDict[geneID][key] = "."
            else:
                for suffix in geneDict["colour"].keys():
                    key = "network_screening_candidates" + f".{suffix}" if suffix != "" else ""
                    wgcnaDfDict[geneID][key] = "."
        
        # Create a WGCNA dataframe
        wgcnaDf = pd.DataFrame.from_dict(wgcnaDfDict, orient="index")
        wgcnaDf.fillna(".", inplace=True)
        
        # Join the DGE and WGCNA dataframes together
        df = df.join(wgcnaDf)
    
    # Parse the normalised counts
    countsDf = parse_normalised_counts(args.normalisedCountsFileName, sampleGroups, set(genes.keys()))
    
    # Join the DGE(+WGCNA) and counts dataframes together
    df = df.join(countsDf)
    
    # Write to file
    df.to_csv(args.outputFileName, sep="\t")
    
    # Notify user of successful completion
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

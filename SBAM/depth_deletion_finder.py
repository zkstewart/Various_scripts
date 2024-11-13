#! python3
# depth_deletion_finder.py
# Script to detect homozygous deletion/non-deletion and heterozygous
# regions from a samtools depth histogram, as produced by plot_samtools_depth.py
# with the --tsv tag.

import os, argparse, re
import numpy as np
import pandas as pd

def validate_args(args):
    # Format the file pattern
    assert args.filePattern.count("()") == 2, \
        f"File pattern '{args.filePattern}' does not contain exactly two instances of '()'!"
    regexString = "^" + args.filePattern.replace(".", "\.").replace("()", "(.+)") + "$"
    args.fileRegex = re.compile(regexString)
    
    # Validate input data locations
    foundFiles = []
    for histogramFile in args.histogramFiles:
        if os.path.isfile(histogramFile):
            fileMatch = args.fileRegex.match(os.path.basename(histogramFile)) 
            if fileMatch:
                foundFiles.append(histogramFile)
            else:
                raise ValueError(f"Input histogram file '{histogramFile}' does not have the expected format '{regexString}'")
        elif os.path.isdir(histogramFile):
            foundAny = False
            for f in os.listdir(histogramFile):
                fileMatch = args.fileRegex.match(os.path.basename(f)) 
                if fileMatch:
                    foundFiles.append(os.path.join(histogramFile, f))
                    foundAny = True
            if not foundAny:
                raise FileNotFoundError(f"No histogram files found in directory '{histogramFile}' with format '{regexString}'")
        else:
            raise FileNotFoundError(f"Input histogram file '{histogramFile}' not found!")
    args.histogramFiles = foundFiles
    
    # Validate numeric arguments
    if args.binSize <= 0:
        raise ValueError(f"Bin size '{args.binSize}' must be a positive integer!")
    
    # Handle file output
    if os.path.exists(args.outputFileName):
        raise FileExistsError(f"Output file name '{args.outputFileName}' already exists!")

def parse_histogram_tsv(tsvFile):
    '''
    Parameters:
        tsvFile -- a string pointing to the TSV file with format:
                   contigID    position    depth
    Returns:
        histoDict -- a dictionary with structure like:
                     {
                         "contigID1": {
                             position1: depth1,
                             position2: depth2,
                             ...
                         },
                         "contigID2": { ... },
                         ...
                     }
    '''
    histoDict = {}
    with open(tsvFile, "r") as fileIn:
        firstLine = True
        for line in fileIn:
            sl = line.strip("\r\n ").split("\t")
            assert len(sl) == 3, f"Line '{line}' in histogram file '{tsvFile}' does not have 3 columns!"
            contigID, position, depth = sl
            
            # Skip header line (if present)
            if firstLine:
                firstLine = False
                if position.isdigit():
                    pass # if no header line, go ahead and process this line
                else:
                    continue # if there is a header line, skip this line
            
            # Validate format
            assert position.isdigit(), f"Position value '{position}' in histogram file '{tsvFile}' is not an integer!"
            assert depth.isdigit(), f"Depth value '{depth}' in histogram file '{tsvFile}' is not an integer!"
            
            # Store data
            histoDict.setdefault(contigID, {})
            histoDict[contigID][int(position)] = int(depth)
    return histoDict

def predict_deletions(binDict):
    '''
    Receives a histogram dictionary and predicts regions of homozygous deletion,
    homozygous presence, and heterozygous regions on the basis of depth coverage.
    Uses a simple heuristic approach.
    
    Parameters:
        binDict -- a dictionary with structure like:
                     {
                         position1: depth1,
                         position2: depth2,
                             ...
                    }
    Returns:
        alleles -- a numpy array with the same length as the input dictionary,
                   where 0 indicates homozygous deletion, 1 indicates heterozygous
                   deletion, and 2 indicates homozygous presence
    '''
    depths = np.array(list(binDict.values()))
    
    # Get breakpoints for deletions and presence
    medianDepth = np.median(depths)
    #duplicatedDepth = medianDepth * 2
    heteroDepth = medianDepth / 2.5
    homoDepth = medianDepth / 8
    
    # Predict deletions and presence
    hetero = np.where((depths > homoDepth) & (depths <= heteroDepth), 1, 0)
    homoPresent = np.where(depths > heteroDepth, 2, 0)
    alleles = hetero + homoPresent
    
    # Return the results
    return alleles

def convert_alleles_to_gt(alleles):
    '''
    Converts the allele predictions to a genotype array, where 0 indicates
    homozygous deletion, 1 indicates heterozygous deletion, and 2 indicates
    homozygous presence.
    
    Parameters:
        alleles -- a numpy array where 0 indicates homozygous deletion,
                   1 indicates heterozygous deletion, and 2 indicates
                   homozygous presence
    Returns:
        genotypes -- a list of strings with the same length as the input array,
                     where each string is a VCF-encoded genotype where '0/0'
                     indicates homozygous presence, '0/1' indicates heterozygous
                     deletion, and '1/1' indicates homozygous deletion
    '''
    genotypes = [
        "1/1" if allele == 0 # homozygous deletion
        else "0/1" if allele == 1 # heterozygous deletion
        else "0/0" # homozygous presence
        for allele in alleles
    ]
    return genotypes

def get_sorted_contig_ids(idsList):
    # Sort contig IDs by their numerical value (if possible)
    allHaveNumbers = all([ any([ c.isdigit() for c in contigID ]) for contigID in idsList ])
    if allHaveNumbers:
        numRegex = re.compile(r"\d+")
        return sorted(idsList, key=lambda x: int("".join(numRegex.findall(x))))
    else:
        return sorted(idsList)

def main():
    usage = """%(prog)s receives one or more histogram TSV files, where each file contains
    a single sample's depth coverage for a single contig. The script parses them in
    order to identify potential homozygous deletions, homozygous presence, and heterozygous
    regions. The output is a TSV file with VCF-like format, where each row represents a
    genomic bin and columns are VCF-encoded where 0 indicates presence and 1 indicates deletion.
    File names must conform to the pattern specified by the --filePattern flag, where two
    instances of '()' are used to denote the sample name and contig ID, respectively.
    """
    # Establish main parser
    p = argparse.ArgumentParser(description=usage)
    
    # Set arguments shared by subparsers
    ## Required arguments
    p.add_argument("-i", dest="histogramFiles",
                    required=True,
                    nargs="+",
                    help="""Specify one or more file names or directories containing
                    histogram TSV files""")
    p.add_argument("-b", dest="binSize",
                   type=int,
                   required=True,
                   help="Specify the bin size used when generating the histogram")
    p.add_argument("-o", dest="outputFileName",
                    required=True,
                    help="Output file name")
    ## Opts (plotting behaviour)
    p.add_argument("--filePattern", dest="filePattern",
                    required=False,
                    help="""Optionally, specify the pattern that uniquely identifies the
                    file names that should be considered as histogram files; must contain two
                    empty parentheses to denote the location of the same name and contig
                    ID (default="()_().histo.tsv")""",
                    default="()_().histo.tsv")
    
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse histogram files
    histoDict = {}
    samples = set()
    contigs = set()
    for histogramFile in args.histogramFiles:
        tsvDict = parse_histogram_tsv(histogramFile)
        assert len(tsvDict) == 1, \
            f"Expected a single contig in histogram file '{histogramFile}'!"
        
        sampleName, contigID = args.fileRegex.match(os.path.basename(histogramFile)).groups()
        
        histoDict.setdefault(sampleName, {})
        assert contigID not in histoDict[sampleName], \
            f"Contig ID '{contigID}' is not unique for sample '{sampleName}'!"
        histoDict[sampleName][contigID] = tsvDict[list(tsvDict.keys())[0]]
        
        samples.add(sampleName)
        contigs.add(contigID)
    
    # Notify user of the detected samples and contigs
    print(f"# Detected samples: {', '.join(samples)}")
    print(f"# Detected contigs: {', '.join(get_sorted_contig_ids(contigs))}")
    
    # Process each histogram file
    genotypesDict = {}
    for sampleName, contigDict in histoDict.items():
        genotypesDict[sampleName] = {}
        for contigID, binDict in contigDict.items():
            alleles = predict_deletions(binDict)
            genotypes = convert_alleles_to_gt(alleles)
            genotypesDict[sampleName][contigID] = genotypes
    
    # Convert to DataFrame with VCF-like format
    df = pd.DataFrame(genotypesDict)
    exploded_df = df.apply(lambda x: x.explode()).reset_index()
    exploded_df.rename(columns={"index": "#CHROM"}, inplace=True)    
    exploded_df["POS"] = (exploded_df.groupby("#CHROM").cumcount()+1) * args.binSize
    exploded_df["ID"] = "."
    exploded_df["REF"] = "N"
    exploded_df["ALT"] = "N"
    exploded_df["QUAL"] = "."
    exploded_df["FILTER"] = "."
    exploded_df["INFO"] = "."
    exploded_df["FORMAT"] = "GT"
    exploded_df = exploded_df[["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", *samples]]
    
    # Write to output file
    with open(args.outputFileName, "w") as fileOut:
        fileOut.write("##fileformat=VCF-like\n")
        fileOut.write("##depth_deletion_finder\n")
        exploded_df.to_csv(fileOut, sep="\t", index=False)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

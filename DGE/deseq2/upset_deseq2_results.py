#!/usr/bin/env python3
# upset_deseq2_results.py

import os, argparse, upsetplot
import pandas as pd
from matplotlib import pyplot

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
    
    # Validate output file suffix
    if not (args.outputFileName.endswith(".pdf") or args.outputFileName.endswith(".png") or args.outputFileName.endswith(".svg")):
        raise ValueError(f"-o output file '{args.outputFileName}' must end with '.pdf', '.png', or '.svg'!")
    
    # Validate output file location
    args.outputFileName = os.path.abspath(args.outputFileName)
    if os.path.exists(args.outputFileName):
        raise FileExistsError(f"The -o output location '{args.outputFileName}' already exists and will not be overwritten!")
    
    parentDir = os.path.dirname(args.outputFileName)
    if not os.path.isdir(parentDir):
        raise ValueError(f"The -o output location '{args.outputFileName}' points to a non-existent parent directory '{parentDir}'")

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

def parse_table_to_dict(fileName, indexColumn):
    df = pd.read_csv(fileName, sep="\t")
    header = list(df.keys()[1:])
    tableDict = {
        k: list(v.values()) 
        for k, v 
        in df.set_index(indexColumn).to_dict('index').items()
    }
    return tableDict, header

def main():
    ##### USER INPUT SECTION
    usage = """%(prog)s will receive DESeq2 output files and produce an UpSet plot
    showing the set membership of DEGs.
    """
    # Required arguments
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="inputLocations",
                   required=True,
                   nargs="+",
                   help="Location of files or directories for DESeq2 results files to collate")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Location to write the output plot file")
    # Optional arguments
    p.add_argument("--is", dest="inputSuffix",
                   required=False,
                   help="""Optionally, specify the file ending of the DESeq2 output files
                   if a directory is given to -i (default='.edit.tsv')""",
                   default=".edit.tsv")
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse each DGE result file
    contents = {}
    for fileName in args.files:
        genesDict, comparison = parse_deseq2_file(fileName, args.inputSuffix)
        contents[comparison] = list(genesDict.keys())
    
    # Load data into upsetplot object
    upset = upsetplot.from_contents(contents)
    
    # Plot to file
    upsetplot.plot(upset, sort_by="cardinality")
    pyplot.savefig(args.outputFileName)
    
    # Notify user of successful completion
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

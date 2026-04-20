#!/usr/bin/env python3
# annotate_deseq2_result.py
# Python script to take an annotation TSV and a raw DESeq2 output file
# and add data columns into the DESeq2 file (e.g., gene names)

import os, argparse
import pandas as pd

def validate_args(args):
    # Validate input file locations depending on type of input
    args.files = []
    for inputLocation in args.inputLocations:
        inputLocation = os.path.abspath(inputLocation)
        
        # Handle directories
        if os.path.isdir(inputLocation):
            foundFiles = False
            for f in os.listdir(inputLocation):
                file = os.path.join(inputLocation, f)
                if os.path.isfile(file) and file.endswith(args.inputSuffix) and (not file.endswith(args.outputSuffix)):
                    args.files.append(file)
                    foundFiles = True
            if not foundFiles:
                raise FileNotFoundError(f"'{inputLocation}' does not contain any files.")
        
        # Handle files
        elif os.path.isfile(inputLocation): # don't need to check suffix here
            args.files.append(inputLocation)
        else:
            raise FileNotFoundError(f"'{inputLocation}' is not a valid file or directory.")
    
    # Validate annotation file location
    args.annotFileName = os.path.abspath(args.annotFileName)
    if not os.path.isfile(args.annotFileName):
        raise FileNotFoundError(f"-a file '{args.annotFileName}' does not exist or is not a file!")
    
    # Validate numeric arguments
    if args.insertIndex < 0:
        raise ValueError("--insert value must be >= 0")

def parse_annot_table(annotTable, keyColumn="#Query", toIndex=["Gene_names"]):
    annot = {}
    with open(annotTable, "r") as fileIn:
        for line in fileIn:
            sl = line.rstrip().split("\t")
            if line.startswith("#"):
                header = sl
                keyIndex = header.index(keyColumn)
                dataIndices = [ [x, header.index(x)] for x in toIndex ]
            else:
                key = sl[keyIndex]
                data = { dkey:sl[dindex].split(" [")[0] for dkey, dindex in dataIndices }
                annot[key] = data
    return annot

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
                   help="Location of files or directories containing files to concatenate")
    p.add_argument("-a", dest="annotFileName",
                   required=True,
                   help="Location to the annotations TSV table file")
    # Optional arguments
    p.add_argument("--key", dest="keyColumn",
                   required=False,
                   help="""Optionally, specify the column header in the -a file where sequence
                   IDs are found (default='#Query')""",
                   default="#Query")
    p.add_argument("--data", dest="dataColumns",
                   required=False,
                   nargs="+",
                   help="""Optionally, specify one or more column headers in the -a file where data
                   column(s) to embed within the DESeq2 file(s) are found (default='Gene_names')""",
                   default=["Gene_names"])
    p.add_argument("--is", dest="inputSuffix",
                   required=False,
                   help="""Optionally, specify the file ending of the DESeq2 output files
                   if a directory is given to -i (default='.tsv')""",
                   default=".tsv")
    p.add_argument("--os", dest="outputSuffix",
                   required=False,
                   help="""Optionally, specify the file ending to give to the modified
                   DESeq2 files (default='.edit.tsv')""",
                   default=".edit.tsv")
    p.add_argument("--insert", dest="insertIndex",
                   required=False,
                   type=int,
                   help="""Optionally, specify the column index to insert data values
                   at (default=1)""",
                   default=1)
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse the annotation TSV
    annot = parse_annot_table(args.annotFileName, keyColumn=args.keyColumn, toIndex=args.dataColumns)
    
    # Modify each file
    for file in args.files:
        filePrefix = file.rsplit(args.inputSuffix, maxsplit=1)[0]
        outputFileName = f"{filePrefix}{args.outputSuffix}"
        if not os.path.exists(outputFileName):
            df = pd.read_csv(file, sep="\t")
            df.rename(columns={"Unnamed: 0": "gene_id"}, inplace=True)
            for col in args.dataColumns[::-1]:
                df.insert(args.insertIndex, col, [ annot[x][col] for x in df["gene_id"] ])
            df.to_csv(outputFileName, sep="\t", index=False)
        else:
            print(f"# WARNING: '{file}' appears to have an existing modified file at '{outputFileName}'; skipping ...")
    
    # Notify user of successful completion
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

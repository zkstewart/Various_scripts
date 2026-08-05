#!/usr/bin/env python3

import argparse
import os

import pandas as pd

ROW_LABELS = ["A", "B", "C", "D", "E", "F", "G", "H"]
COL_LABELS = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"]
EXPECTED_SHAPE = (8, 12)

def validate_args(args):
    # Validate input file
    args.inputExcelFile = os.path.abspath(args.inputExcelFile)
    if not os.path.isfile(args.inputExcelFile):
        raise FileNotFoundError(f"-i '{inputFile}' is not a file or does not exist.")
    
    # Validate output file
    args.outputFileName = os.path.abspath(args.outputFileName)
    if os.path.exists(args.outputFileName):
        raise FileExistsError(f"-o '{args.outputFileName}' already exists and will not be overwritten.")
    
    baseFile = os.path.basename(args.outputFileName)
    parentDir = os.path.dirname(args.outputFileName)
    if not os.path.isdir(parentDir):
        raise NotADirectoryError(f"-o cannot create '{baseFile}' in '{parentDir}' as " + 
                                 "that location does not exist or is not a directory.")

def main():
    ##### USER INPUT SECTION
    usage = """%(prog)s parses an Excel file with one or more sheets
    corresponding to microwell sequencing plates with 8 rows and 12 columns
    each. The output is a TSV file with three columns: plate well sample,
    where 'well' is row letter (from A to H) and column number (from 1 to 12).
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="inputExcelFile",
                   required=True,
                   help="Specify the input Excel file name")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Specify the output TSV file name")
    
    args = p.parse_args()
    validate_args(args)
    
    # Read the sheet names from the Excel file
    excel = pd.ExcelFile(args.inputExcelFile)
    sheets = excel.sheet_names
    
    formatted = [ f"'{name}' = Plate-{i+1}" for i, name in enumerate(sheets) ]
    print(f"# Sheets in '{args.inputExcelFile}' are being interpreted as:\n" + 
          "\n".join(formatted))
    
    # Iterate over each sheet for the well layouts
    metaDict = {}
    for plateNum, sheetName in enumerate(sheets):
        df = pd.read_excel(args.inputExcelFile, sheet_name=sheetName)
        
        # Identify and drop a label column (if applicable)
        firstColumn = df.iloc[:,0]
        if firstColumn.to_list() == ROW_LABELS:
            df = df.drop(firstColumn.name, axis=1)
        
        # Identify and drop a label row (if applicable)
        nrow, ncol = df.shape
        if nrow == 9:
            firstRow = df.iloc[0,:]
            if firstRow.to_list() == COL_LABELS:
                df = df.drop(firstRow.name, axis=0)
        
        # Make sure our table matches the expected shape
        nrow, ncol = df.shape
        expectedRows, expectedCols = EXPECTED_SHAPE
        if nrow != expectedRows:
            raise ValueError(f"Sheet '{sheetName}' has {nrow} rows when there should be {expectedRows}")
        if ncol != expectedCols:
            raise ValueError(f"Sheet '{sheetName}' has {ncol} cols when there should be {expectedCols}")
        
        # Parse the well layout into the metadata dict
        metaDict[plateNum+1] = { # make plateNum 1-based
            f"{ROW_LABELS[x]}{y+1}" : cell
            for x, row in enumerate(df.itertuples())
            for y, cell in enumerate(list(row)[1:]) # skip the index at [0]
            if not pd.isna(cell)
        }
    
    # Write the output file
    with open(args.outputFileName, "w") as fileOut:
        fileOut.write("plate\twell\tsample\n")
        for plateNum, layout in metaDict.items():
            for well, sample in layout.items():
                fileOut.write(f"{plateNum}\t{well}\t{sample}\n")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

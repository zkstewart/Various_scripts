#! python3
# table_merger.py
# Script to take two tables (file A and file B) and merge their
# contents using Pandas magic. Capable of handling files in three
# main output formats (csv, tsv, and Excel).

import os, argparse
import pandas as pd

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.fileA):
        print(f'I am unable to locate file A ({args.fileA})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.fileB):
        print(f'I am unable to locate file B ({args.fileB})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print(f'File already exists at output location ({args.outputFileName})')
        print('Make sure you specify a unique file name and try again.')
        quit()

def load_pd_df_from_file(tableFile, format):
    '''
    Utility function to load a Pandas dataframe
    depending on the file input format.
    
    Parameters:
        tableFile -- a string indicating a table file with format in
                     ["csv", "tsv", "xlsx"]
        format -- a string with value in ["csv", "tsv", "xlsx", "guess"]
                  which tells us how to read the tableFile
    Returns
        df -- a Pandas dataframe containing the loaded data
        format -- the file format which the file was loaded in from;
                  possible values are ["csv", "tsv", "xlsx"]
    '''
    assert format in ["csv", "tsv", "xlsx", "guess"], \
        f"Format '{format}' is not recognised as a valid argument!"
    
    # Guess the file type if possible
    if format == "guess":
        fileSuffix = tableFile.split(".")[-1].lower()
        if fileSuffix == "csv":
            format = "csv"
        elif fileSuffix == "tsv":
            format = "tsv"
        elif fileSuffix == "xlsx" or fileSuffix == "xls":
            format = "xlsx"
        elif fileSuffix == "txt" or fileSuffix == "text":
            # Check if csv or tsv
            with open(tableFile, "r") as fileIn:
                for line in fileIn:
                    if "\t" in line:
                        format = "tsv"
                    elif "," in line:
                        format = "csv"
                    else:
                        print(f"Unable to determine whether '{tableFile}' is TSV or CSV!")
                        print("(This means it doesn't contain a tab or comma in the first line)")
                        print("Make sure your file is formatted appropriately and try again.")
                        quit()
        else:
            print(f"Files with '.{fileSuffix}' ending aren't recognised.")
            print("Try to make sure your CSV, TSV, or Excel file has an appropriate file suffix.")
            print("(Or, specify the -aFormat and/or -bFormat arguments)")
            print("Program will exit now.")
            quit()
    
    # Load based on file type
    if format == "csv":
        df = pd.read_csv(tableFile)
    elif format == "tsv":
        df = pd.read_csv(tableFile, sep="\t")
    elif format == "xlsx":
        df = pd.read_excel(open(tableFile, "rb"), sheet_name=0)
    
    return df, format

def validate_primary_keys(dfA, aKey, dfB, bKey):
    '''
    Parameters:
        dfA -- a Pandas dataframe representing file A
        aKey -- a column header in file A
        (likewise for B values)
    '''
    # Check that the columns exist
    assert aKey in dfA, \
        f"'{aKey}' not found as a column in file A!"
    assert bKey in dfB, \
        f"'{bKey}' not found as a column in file B!"
    
    # Check that columns have overlapping keys
    keyIntersection = set(dfA[aKey]).intersection(dfB[bKey])
    if len(keyIntersection) == 0:
        print(f"WARNING: No values in file A's '{aKey}' column were found in file B's '{bKey}' column!")
        print("This indicates a possible incompatibility in your two files.")
        print("i.e., since there's no overlap, there's no possibility to carry over any values")
        print("If you expected this possibility, then don't worry; if this shouldn't happen, your results will be flawed!")

def carry_columns_from_A_to_B(dfA, aKey, dfB, bKey, carryColumns, howToMerge):
    '''
    Parameters:
        dfA -- a Pandas dataframe representing file A
        aKey -- a column header in file A
        (likewise for B values)
        carryColumns -- a list containing strings relating to column headings
                        in dfA
        howToMerge -- a string in the list of ["right", "inner"]
    '''
    # Validate carryColumns values
    for column in carryColumns:
        assert column in dfA, \
            f"'{column}' not found as a column in file A!"
    
    # Validate howToMerge value
    assert howToMerge in ["right", "inner"], \
        f"'{howToMerge}' not supported as table merging method"
    
    # Subset dfA to just the primary key and columns to carry over
    dfA = dfA.loc[:, dfA.columns.isin([aKey, *carryColumns])]
    
    # Merge dataframes together
    dfMerged = pd.merge(dfA, dfB, left_on=aKey, right_on=bKey, how=howToMerge)
    
    return dfMerged

## Main
def main():
    # User input
    usage = """%(prog)s reads in two table files and merges them together in a SQL-like
    fashion. To do so, you must specify "primary key" columns from
    each file which contain shared values. 
    
    Any values in file B which do not have
    a corresponding value in file A will have blank values indicated.
    
    Note 1: this script will write the output in the same format as fileB.
    Note 2: if loading in an Excel file, it's assumed to be the FIRST worksheet.
    """
    p = argparse.ArgumentParser(description=usage)
    ## Required file inputs
    p.add_argument("-a", dest="fileA",
                   required=True,
                   help="""Input file containing values that you want to carry over
                   into file B""")
    p.add_argument("-b", dest="fileB",
                   required=True,
                   help="Target file to receive values from file A")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="""Output file name representing file B with values
                   carried over from A""")
    ## Required behaviour parameters
    p.add_argument("--carryColumns", dest="carryColumns", nargs="+",
                   required=True,
                   help="""Indicate the column headings you want to carry over from A
                   to B. Heading values should ideally NOT have blank spaces in them;
                   if they do, make sure to enclose those values in quotation marks
                   e.g., '--carryColumns Gene_name "Toxin Family" details' is a valid way
                   to specify the three indicated columns in the presence of blank spaces""",)
    p.add_argument("--aKey", dest="aKey",
                   required=True,
                   help="""Indicate the column header containing IDs from file A""")
    p.add_argument("--bKey", dest="bKey",
                   required=True,
                   help="""Indicate the column header containing IDs from file B""")
    ## Optional
    p.add_argument("--blanks", dest="blanks",
                   required=False,
                   help="""Optionally, control what value you want when there are blank
                   values (default == '.')""",
                   default=".")
    p.add_argument("--aFormat", dest="aFormat",
                   required=False,
                   choices=["guess", "tsv", "csv", "excel"],
                   help="""Optionally provide explicit instruction for the file format
                   of file A (default == 'guess'); specify this argument if default behaviour
                   provides errors""",
                   default="guess")
    p.add_argument("--bFormat", dest="bFormat",
                   required=False,
                   choices=["guess", "tsv", "csv", "excel"],
                   help="""Optionally provide explicit instruction for the file format
                   of file B (default == 'guess'); specify this argument if default behaviour
                   provides errors""",
                   default="guess")
    p.add_argument("--howToMerge", dest="howToMerge",
                   required=False,
                   choices=["right", "inner"],
                   help="""Optionally indicate how the files should be merged together;
                   default == 'right' which means values from file A will be merged into
                   file B, with non-hits in B having blank values. Specifying 'inner' means
                   only the intersection of the two files (based on aKey and bKey) will be output""",
                   default="right")
    
    args = p.parse_args()
    validate_args(args)
    
    # Load in table files
    dfA, formatA = load_pd_df_from_file(args.fileA, args.aFormat)
    dfB, formatB = load_pd_df_from_file(args.fileB, args.bFormat)
    
    # Validate primary keys
    validate_primary_keys(dfA, args.aKey, dfB, args.bKey)
    
    # Merge tables
    dfMerged = carry_columns_from_A_to_B(dfA, args.aKey, dfB, args.bKey, args.carryColumns, args.howToMerge)
    
    ## TBD:
    # 1) With right join, make sure table is ordered the same as the original fileB
    # 2) Make sure the columns from fileA are the right-most in the output file
    
    # Write output file
    if formatB == "csv":
        dfMerged.to_csv(args.outputFileName, sep=",", index=False, na_rep=args.blanks)
    elif formatB == "tsv":
        dfMerged.to_csv(args.outputFileName, sep="\t", index=False, na_rep=args.blanks)
    elif formatB == "xlsx":
        dfMerged.to_excel(args.outputFileName, sep=",", index=False, na_rep=args.blanks)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

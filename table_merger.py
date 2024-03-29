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
    VALID_OUTPUT_SUFFIXES = ["csv", "tsv", "xlsx"]
    fileSuffix = args.outputFileName.split(".")[-1].lower()
    if not fileSuffix in VALID_OUTPUT_SUFFIXES:
        print(f"'{fileSuffix}' suffix is not recognised as a valid output file format")
        print("Your -o file name should end in .csv, .tsv, or .xlsx")
        print("Fix this up and try again.")
        quit()
    
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
    
    return df

def validate_primary_keys(dfA, aKey, dfB, bKey):
    '''
    Parameters:
        dfA -- a Pandas dataframe representing file A.
        aKey -- a column header in file A.
        (likewise for B values)
    '''
    # Check that the primary key columns exist
    assert aKey in dfA, \
        f"aKey value '{aKey}' not found as a column in file A!"
    assert bKey in dfB, \
        f"bKey value '{bKey}' not found as a column in file B!"
    
    # Check that columns are unique
    assert len(set(dfA[aKey])) == len(dfA[aKey]), \
        f"File A's '{aKey}' column does not contain unique values!"
    assert len(set(dfB[bKey])) == len(dfB[bKey]), \
        f"File B's '{bKey}' column does not contain unique values!"
    
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
        dfA -- a Pandas dataframe representing file A.
        aKey -- a column header in file A.
        (likewise for B values)
        carryColumns -- a list containing strings relating to column headings
                        in dfA.
        howToMerge -- a string in the list of ["right", "inner"] controlling how
                      dfA should be merged into dfB; these values correspond to
                      SQL-like merging methods, such that "right" indicates a
                      concatenation of values into dfB (blanks allowed),
                      whereas "inner" means only shared primary key rows will
                      be output (no blanks allowed).
    Returns:
        dfMerged -- a Pandas dataframe with carryColumns from dfA merged into
                    dfB.
    '''
    # Validate carryColumns values
    for column in carryColumns:
        assert column in dfA, \
            f"carryColumn value '{column}' not found as a column in file A!"
        assert column not in dfB, \
            f"carryColumn value '{column}' already exists as a column in file B! Carrying it over would replace it!"
    
    # Validate howToMerge value
    assert howToMerge in ["right", "inner"], \
        f"'{howToMerge}' not supported as table merging method"
    
    # Validate that bKey value is not redundantly shared by file A
    """If we're not merging by a shared common primary key, we have to assume
    the primary key from file B is not also in file A. If that's the case, we cannot
    rename the key below to prevent redundant column introduction. Also, we'd risk
    overwriting columns. At the end of the day this shouldn't happen if someone has a
    logically structured file..."""
    if bKey != aKey:
        assert bKey not in dfA, \
            f"'{bKey}' is already a column in file A, but we're not using it as the primary key? I can't handle this, sorry."
    
    # Set dfA key column to be identical to bKey value
    "This means we don't need to worry about a redundant aKey column being carried over, too"
    dfA.rename(columns = {aKey:bKey}, inplace=True)
    
    # Subset dfA to just primary key + columns to carry over
    dfA = dfA[[bKey] + carryColumns] # use bKey here since aKey doesn't exist anymore
    
    # Merge dataframes together
    """Have to logically flip left to right here since dfB (our right file)
    behaves as the left file when we directly merge into it like this with pandas"""
    dfMerged = dfB.merge(
        dfA,
        on=bKey, # both tables have the same key header values now, no need for left_on or right_on
        how="left" if howToMerge == "right" else howToMerge
    )
    
    return dfMerged

## Main
def main():
    # User input
    usage = """%(prog)s reads in two table files and merges them together in a SQL-like
    fashion. To do so, you must specify "primary key" columns from
    each file which contain shared values. Importantly, these columns are expected to
    have unique values.
    
    Any values in file B which do not have a corresponding key in file A will have
    blank values input; the --blanks argument controls what these blanks look like.
    
    Note 1: this script will write the output in the format indicated by the -o argument.
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
                   help="""Output file name representing file B with columns
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
    p.add_argument("--howToMerge", dest="howToMerge",
                   required=False,
                   choices=["right", "inner"],
                   help="""Optionally indicate how the files should be merged together;
                   default == 'right' which means values from file A will be merged into
                   file B, with non-hits in B having blank values. Specifying 'inner' means
                   only the intersection of the two files (i.e., shared aKey and bKey values)
                   will be output""",
                   default="right")
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
    
    args = p.parse_args()
    validate_args(args)
    
    # Load in table files
    dfA = load_pd_df_from_file(args.fileA, args.aFormat)
    dfB = load_pd_df_from_file(args.fileB, args.bFormat)
    
    # Validate primary keys
    validate_primary_keys(dfA, args.aKey, dfB, args.bKey)
    
    # Merge tables
    dfMerged = carry_columns_from_A_to_B(dfA, args.aKey, dfB, args.bKey, args.carryColumns, args.howToMerge)
    
    # Write output file
    fileSuffix = args.outputFileName.split(".")[-1].lower() # already validated in validate_args()
    
    if fileSuffix == "csv":
        dfMerged.to_csv(args.outputFileName, sep=",", index=False, na_rep=args.blanks)
    elif fileSuffix == "tsv":
        dfMerged.to_csv(args.outputFileName, sep="\t", index=False, na_rep=args.blanks)
    elif fileSuffix == "xlsx":
        dfMerged.to_excel(args.outputFileName, index=False, na_rep=args.blanks)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

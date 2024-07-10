#! python3
# filter_families_on_difference.py
# Script to handle the filtering of CAFE5 gene count files based on the difference
# between the maximum and minimum gene count per row. This can be necessary if CAFE 
# gives errors pertaining to not finding a sensible lambda value to initialise.

import os, argparse
import pandas as pd
import numpy as np

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.geneCountFile):
        print('I am unable to locate the GeneCount file (' + args.geneCountFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate numeric inputs
    if 0 >= args.percentile:
        print("--percentile should be greater than 0")
        quit()
    if 100 <= args.percentile:
        print("--percentile should be less than 100")
        quit()
    # Validate output file location
    if not os.path.isdir(args.outputDirectory):
        os.makedirs(args.outputDirectory)
        print(f"Output directory '{args.outputDirectory}' has been created as part of argument validation.")

def main():
    # User input
    usage = """%(prog)s performs filtration of OrthoFinder's Orthogroups.GeneCount.tsv file which
    is used by CAFE for analysis. The script filters out rows that have a high difference between
    the maximum and minimum gene count per row. This can be necessary if CAFE gives errors.
    For --percentile, the remaining value out of 100 is the portion that will be removed. Hence,
    if you provide a value close to 0, you'll filter almost everything, and a value close to 100
    will filter very little out. The default value is 99.95, which removes the top 0.5 percent of
    families and is likely to be sufficient to remove the families that are causing CAFE to error out.
    """
    p = argparse.ArgumentParser(description=usage)
    # Reqs
    p.add_argument("-i", dest="geneCountFile",
                   required=True,
                   help="Input Orthogroups.GeneCount.tsv file")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Specify the output file name")
    # Opt (behavioural)
    p.add_argument("--percentile", dest="percentile",
                   required=False,
                   type=float,
                   help="""Optionally, specify the percentile to use as the cutoff for filtering rows
                   based on the difference between the maximum and minimum gene count per row;
                   default == 99.95 to remove the top 0.5 percent of families""",
                   default=99.95)
    p.add_argument("--keepOneSpeciesRows", dest="keepOneSpeciesRows",
                   required=False,
                   action="store_true",
                   help="""Optionally provide this flag to NOT filter rows that only have one species;
                   you probably should do this if the output file is to be used by CAFE since these rows
                   are non-informative for CAFE analysis""",
                   default=False)
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse gene count file into pandas df
    df = pd.read_csv(args.geneCountFile, sep="\t", header=0)
    
    # Check that file meets expected format
    assert df.columns[0] == "Orthogroup", f"First column of '{args.geneCountFile}' should be 'Orthogroup'"
    assert df.columns[-1] == "Total", f"Last column of '{args.geneCountFile}' should be 'Total'"
    
    # Get the species prefixes from the df header
    "Everything inbetween 'Orthogroup' and 'Total' should be a species prefix"
    speciesPrefixes = [ col for col in df.columns if not col in ["Orthogroup", "Total"] ]
    
    # Filter rows with only one species represented (if applicable)
    preKeepRowNum = df.shape[0]
    if not args.keepOneSpeciesRows:
        df = df[df[speciesPrefixes].ne(0).sum(axis=1).ne(1)]
    
    # Get the difference between the maximum and minimum gene count per row
    df["minmax"] = df.loc[:, speciesPrefixes].apply(lambda g: g.max()-g.min(), axis=1)
    
    # Identify minmax value to use as cutoff when filtering rows 
    minmaxCutoff = np.percentile(df["minmax"], args.percentile)
    
    # Filter rows with high minmax
    preFilterRowNum = df.shape[0]
    df = df[df["minmax"] < minmaxCutoff]
    df = df.drop(columns=["minmax"])
    postFilterRowNum = df.shape[0]
    
    # Write output
    df.to_csv(args.outputFileName,
              sep="\t", index=False, lineterminator="\n")
    
    # Present statistics on what was removed
    print("# filter_families_on_difference.py output statistics:")
    print(f"# > '{args.geneCountFile}' originally had {preKeepRowNum} rows.")
    print(f"# > {preKeepRowNum-preFilterRowNum} rows were removed due to having only one species represented.")
    print(f"# > {preFilterRowNum-postFilterRowNum} rows were removed due to having a high difference " + 
          "between the maximum and minimum gene count per row.")
    print(f"# > The output '{args.outputFileName}' now has {postFilterRowNum} rows.")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

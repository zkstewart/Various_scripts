#! python3
# kreport_reformat.py
# Takes the result of 'combine_kreports.py' from KrakenTools
# and reformats it in a way that works better for me.

import os
import argparse
import pandas as pd

# Define functions
def validate_args(args):
    # Validate input file location
    args.kreportFile = os.path.abspath(args.kreportFile)
    if not os.path.isfile(args.kreportFile):
        raise FileNotFoundError(f"The -i location ('{args.kreportFile}') does not exist or is not a file")
    
    # Validate output file location
    if os.path.exists(args.outputFileName):
        raise FileExistsError(f"The -o location ('{args.outputFileName}') already exists; this script will not overwrite any files.")

## Main
def main():
    # User input
    usage = """%(prog)s accepts a kraken report file tabulated by KrakenTools' combine_kreports.py
    file and modifies the format to show a simpler breakdown of the percentage reads associated with
    each taxon across the samples.
    """
    p = argparse.ArgumentParser(description=usage)
    ## Required
    p.add_argument("-i", dest="kreportFile",
                   required=True,
                   help="Input kreport file")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file name for the modified repore")
    
    args = p.parse_args()
    validate_args(args)
    
    # Obtain the total sum of reads per sample
    df = pd.read_csv(args.kreportFile, sep="\t", comment="#", header=None)
    nrow, ncol = df.shape
    
    rootNodes = df[~ df[ncol-1].str.contains(" ")]
    sampleRootCols = rootNodes[[ n for n in range(3, ncol-4, 2)]] # ignore first 3 cols and last 4 cols; every 2nd row is the total count
    sampleTotalReads = sampleRootCols.sum().to_list()
    
    # Parse through the table to calculate each taxon level read percentage per sample    
    columns = []
    with open(args.kreportFile, "r") as fileIn, open(args.outputFileName, "w") as fileOut:
        for line in fileIn:
            # Handle header lines
            if line.startswith("#"):
                # Handle file identifying lines
                if line.startswith("#S"):
                    fileNum, indivFileName = line.strip("# \r\n").split("\t")
                    prefix = os.path.basename(indivFileName).split("_report.txt")[0]
                    columns.append(prefix)
                
                # Handle column identifying line
                elif line.startswith("#perc"):
                    fileOut.write("\t".join(columns) + "\ttaxon\n")
            
            # Handle body lines
            else:
                sl = line.rstrip().split("\t")
                row = map(str, [ round(int(x) / y, 4) for x, y in zip(sl[3:ncol-4:2], sampleTotalReads) ])
                fileOut.write("\t".join(row) + f"\t{sl[-1]}\n")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

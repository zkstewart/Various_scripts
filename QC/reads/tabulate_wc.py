#! python3
# tabulate_wc.py
# Simple script to read in .wc files
# and tabulate the results for easy understanding

import os, argparse

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isdir(args.wcDir):
        raise FileNotFoundError('I am unable to locate the parent directory where .wc files are (' + args.wcDir + '). ' + 
                                'Make sure you\'ve typed the file name or location correctly and try again.')
    # Validate numeric arguments
    if args.divideBy <= 0:
        raise ValueError('The divideBy argument must be a positive integer. Please provide a valid number and try again.')
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        raise FileExistsError('File already exists at output location (' + args.outputFileName + '). ' +
                              'Make sure you specify a unique file name and try again.')

def get_wc_from_files(wcFiles, suffix, divideBy=1):
    '''
    For FASTQ files, the read count is the number of lines divided by 4.
    This function reads in .wc files and returns a table of the read count
    values for each sample.
    
    Parameters:
        wcFiles -- a list containing strings that point to .wc
                   files containing read number values
        suffix -- a string that identifies the file suffix
        divideBy -- OPTIONAL; an integer to divide the read count
                    values by; default is 1 (no division)
    Returns:
        statsTable -- a list of strings with the following format:
                      [
                          ["sample\twc"], # with another header if divideBy != 1
                          [f"{sample1}\t{wc1}"], # with another column if divideBy != 1
                          ...
                      ]
    '''
    # Setup the output table format
    if divideBy != 1:
        statsTable=["sample\twc\tdivided_by_"+str(divideBy)]
    else:
        statsTable=["sample\twc"]
    
    # Parse all logs files
    uniqueNames = set()
    for file in wcFiles:
        # Get the sample name from the directory path
        "We assume the file name (sans suffixes) is the sample name"
        sampleName = os.path.basename(file).rsplit(suffix, maxsplit=1)[0]
        if sampleName in uniqueNames:
            raise ValueError("Duplicate sample name found: " + sampleName)
        uniqueNames.add(sampleName)
        
        # Parse the file
        with open(file, "r") as fileIn:            
            firstLine = True
            for line in fileIn:
                if not firstLine:
                    raise ValueError("File is not formatted correctly; has more than one line")
                
                sl = line.strip("\r\n ").split()
                readcount = sl[0]
                
                if divideBy != 1:
                    dividedReadCount = int(readcount) / divideBy
                    statsTable.append(f"{sampleName}\t{readcount}\t{dividedReadCount}")
                else:
                    statsTable.append(f"{sampleName}\t{readcount}")
                
                firstLine = False
    
    return statsTable

## Main
def main():
    # User input
    usage = """%(prog)s accepts a directory containing .wc files 
    and writes a condensed table of wc values to the output file. This is useful for
    FASTQ files to summarise the number of read counts in each sample. The output file
    will contain two columns: the sample name and the read count. Optionally, you can
    provide a number to divide the read count by; this will produce an additional column
    in the output file with the divided values (a value of 4 is useful here for FASTQ files).
    """
    p = argparse.ArgumentParser(description=usage)
    ## Required
    p.add_argument("-i", dest="wcDir",
                   required=True,
                   help="Input directory containing .wc files")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file name for the reads count table")
    ## Optional
    p.add_argument("--divide", dest="divideBy",
                   required=False,
                   type=int,
                   help="""Optionally, provide a number to divide the read count by; this
                   will produce an additional column in the output file with the divided values
                   but is not required and a value of 1 will be used if not provided (produces
                   no additional column) (default==1)""",
                   default=1)
    p.add_argument("--suffix", dest="suffix",
                   required=False,
                   help="""Optionally, indicate the file suffix that identifies files containing
                   wc results (default=='.wc')""",
                   default=".wc")
    
    args = p.parse_args()
    validate_args(args)
    
    # Get .wc files list
    flagFiles = []
    for file in os.listdir(args.wcDir):
        if file.endswith(args.suffix):
            flagFiles.append(os.path.join(args.wcDir, file))
    if len(flagFiles) == 0:
        raise FileNotFoundError(f"No '{args.suffix}' files found in the -i directory ({args.wcDir}). " + 
                                "Make sure you\'ve typed the file name or location correctly and try again")
    
    # Combine .wc files
    statsTable = get_wc_from_files(flagFiles, args.suffix, args.divideBy)
    
    # Write output
    with open(args.outputFileName, "w") as fileOut:
        fileOut.write("\n".join(statsTable))
    
    # Alert user to program success
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

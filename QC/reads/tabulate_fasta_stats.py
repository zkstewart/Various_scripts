#! python3
# tabulate_fasta_stats.py
# Script to crawl through all dirs, subdirs, and files nested under
# the provided location to identify output files from 'genome_stats.py'
# or 'annotarium.py fasta stats' to tabulate their results

import os, argparse, math
import pandas as pd

# Define functions
def validate_args(args):
    # Validate input file locations
    args.parentDir = os.path.abspath(args.parentDir)
    if not os.path.isdir(args.parentDir):
        raise FileNotFoundError(f"The -i location ('{args.parentDir}') does not exist or is not a directory")
    
    # Validate output file location
    if os.path.exists(args.outputFileName):
        raise FileExistsError(f"The -o location ('{args.outputFileName}') already exists; this script will not overwrite any files.")

def suffix_crawl(thisDir, suffix, foundFiles):
    '''
    Recursively locates all files with specified suffix nested under a parent directory.
    
    Parameters:
        thisDir -- a string indicating the full path of the location to crawl through
        suffix -- a string indicating the suffix that uniquely identifies all relevant files
        foundFiles -- a list storing any files found in a previous recursive loop, and to
                      be populated with any files found within thisDir.
    '''
    for location in os.listdir(thisDir):
        fullLocation = os.path.join(thisDir, location)
        if os.path.isfile(fullLocation) and fullLocation.endswith(suffix):
            foundFiles.append(fullLocation)
        elif os.path.isdir(fullLocation):
            suffix_crawl(fullLocation, suffix, foundFiles)

def parse_stats_file(fileName):
    '''
    Identifies whether the stats file should be handled by parse_genomestats_file()
    or parse_annotarium_file() and passes the result back. Each subfunction handles
    a differently formatted file
    '''
    with open(fileName, "r") as fileIn:
        for line in fileIn:
            break
    
    if line.startswith("Genome size (bp):"): # genome_stats.py
        return parse_genomestats_file(fileName)
    elif line.startswith("\tlength\tnum_lowercase"): # annotarium fasta stats
        return parse_annotarium_file(fileName)
    else: # unknown
        raise NotImplementedError(f"File '{fileName}' has a format which is unrecognised")

def parse_genomestats_file(fileName):
    '''
    Parses a statistic file with format like:
    Genome size (bp): 63,684,485,220
    Number of contigs: 4,894,702
    Shortest contig: 5,000
    Longest contig: 104,808

    N50: 15,463
    Median: 9,344
    Mean: 13,010
    
    Parameters:
        fileName -- a string indicating the stats file to parse
    Returns:
        statsDict -- a dictionary where keys are standardised metrics (strings)
                     and values are their corresponding value
    '''
    statsDict = {}
    with open(fileName, "r") as fileIn:
        for line in fileIn:
            l = line.rstrip()
            if l != "":
                if l.startswith("Genome size"):
                    statsDict["total_bp"] = int(l.split(": ")[-1].replace(",", ""))
                elif l.startswith("Number of contigs"):
                    statsDict["num_contigs"] = int(l.split(": ")[-1].replace(",", ""))
                elif l.startswith("Shortest contig"):
                    statsDict["shortest_contig"] = int(l.split(": ")[-1].replace(",", ""))
                elif l.startswith("Longest contig"):
                    statsDict["longest_contig"] = int(l.split(": ")[-1].replace(",", ""))
                elif l.startswith("N50"):
                    statsDict["n50_contig"] = int(l.split(": ")[-1].replace(",", ""))
                elif l.startswith("Median"):
                    statsDict["median_contig"] = int(l.split(": ")[-1].replace(",", ""))
                elif l.startswith("Mean"):
                    statsDict["mean_contig"] = int(l.split(": ")[-1].replace(",", ""))
                else:
                    raise ValueError(f"Line '{l}' is not recognised as a valid genome stats line")
    return statsDict

def parse_annotarium_file(fileName):
    '''
    Parses a statistic file with format like:
            length  num_lowercase   num_n
    chr1    1000    5               0
    chr2    ...     ...             ...
    #Total  9999
    #N50    5000
    #Median 3333
    #Mean   3300
    
    Parameters:
        fileName -- a string indicating the stats file to parse
    Returns:
        statsDict -- a dictionary where keys are standardised metrics (strings)
                     and values are their corresponding value
    '''
    statsDict = {}
    num = 0
    shortest = math.inf
    longest = -math.inf
    with open(fileName, "r") as fileIn:
        for line in fileIn:
            l = line.rstrip()
            
            if l.startswith("\tlength\tnum_lowercase"):
                pass # ignore the header line
            elif l.startswith("#"):
                if l.startswith("#Total"):
                    statsDict["total_bp"] = int(float(l.split("\t")[1]))
                elif l.startswith("#N50"):
                    statsDict["n50_contig"] = int(float(l.split("\t")[1]))
                elif l.startswith("#Median"):
                    statsDict["median_contig"] = int(float(l.split("\t")[1]))
                elif l.startswith("#Mean"):
                    statsDict["mean_contig"] = int(float(l.split("\t")[1]))
                else:
                    raise ValueError(f"Line '{l}' is not recognised as a valid annotarium stats line")
            else:
                sl = l.split("\t")
                try:
                    chrom, length, lowercase, n = sl
                    length, lowercase, n = int(float(length)), int(float(lowercase)), int(float(n)) # convert lowercase and n as format validation
                except:
                    raise ValueError(f"Line '{l}' is not recognised as a valid annotarium stats line")
                
                num += 1
                if length < shortest:
                    shortest = length
                if length > longest:
                    longest = length
    statsDict["num_contigs"] = num
    statsDict["shortest_contig"] = shortest
    statsDict["longest_contig"] = longest
    
    return statsDict

## Main
def main():
    # User input
    usage = """%(prog)s accepts a parent directory, within which all dirs and subdirs will be
    crawled through to locate files with an indicated suffix. These files are expected to be
    produced by 'genome_stats.py' or 'annotarium.py fasta stats'. The result is a table listing
    each statistic for all located samples/files.
    """
    p = argparse.ArgumentParser(description=usage)
    ## Required
    p.add_argument("-i", dest="parentDir",
                   required=True,
                   help="Input directory containing FASTA statistics files")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file name for the statistics")
    ## Optional
    p.add_argument("--suffix", dest="suffix",
                   required=False,
                   help="""Optionally, indicate the file suffix that identifies files containing
                   statistical results (default=='.stats')""",
                   default=".stats")
    
    args = p.parse_args()
    validate_args(args)
    
    # Locate all stats files
    statsFiles = []
    suffix_crawl(args.parentDir, args.suffix, statsFiles)
    if len(statsFiles) == 0:
        raise FileNotFoundError(f"No files ending in '{args.suffix}' found in the -i directory ({args.parentDir}).")
    
    # Parse each stats file
    statsDict = {}
    for fileName in statsFiles:
        baseName = os.path.basename(fileName)
        prefix = baseName.rsplit(args.suffix, maxsplit=1)[0]
        statsDict[prefix] = parse_stats_file(fileName)
    
    # Tabulate
    table = pd.DataFrame.from_dict(statsDict)
    table = table.T # samples as rows, statistics as columns
    table = table[["total_bp", "num_contigs", "shortest_contig", "longest_contig",
                   "n50_contig", "median_contig", "mean_contig"]] # reorder columns
    
    # Write output
    table.to_csv(args.outputFileName, sep="\t")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

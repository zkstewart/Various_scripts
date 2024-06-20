#! python3
# dorado_summary_stats.py
# A simple python program which reads a summary.tsv produced by dorado
# and calculates a handful of statistics including number of sequences,
# total basepairs sequenced, N50, median, and mean contig lengths.

import os, argparse, locale
from statistics import median, mean
locale.setlocale(locale.LC_ALL, '')

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.summaryTSV):
        print(f'I am unable to locate the summary TSV file ({args.representativesFasta})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Handle file overwrites
    if args.outputFileName != False:
        if os.path.exists(args.outputFileName):
            print(f'File already exists at output location ({args.outputFileName})')
            print('Make sure you specify a unique file name and try again.')
            quit()

def N50(numlist): 
    """ 
    Abstract: Returns the N50 value of the passed list of numbers. 
    Usage: N50(numlist) 

    Based on the definition from this SEQanswers post 
    http://seqanswers.com/forums/showpost.php?p=7496&postcount=4 
    (modified Broad Institute's definition 
    https://www.broad.harvard.edu/crd/wiki/index.php/N50) 

    See SEQanswers threads for details: 
    http://seqanswers.com/forums/showthread.php?t=2857 
    http://seqanswers.com/forums/showthread.php?t=2332 
    """ 
    numlist.sort(reverse = True) 
    s = sum(numlist) 
    limit = s * 0.5 
    for l in numlist: 
        s -= l 
        if s <= limit: 
            return l

def main():
    # Argparse handling
    usage = """%(prog)s reads in a dorado summary TSV file and calculates a handful
    of statistics, including the number of reads, the size distribution of reads
    including shortest, longest, median, mean, and N50 values, as well as the total
    amount of sequencing data in bp. These values are printed to terminal and,
    optionally, a text file may be produced.
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="summaryTSV",
                   help="Input dorado summary TSV file")
    p.add_argument("-o", dest="outputFileName",
                   default=False,
                   help="Optionally, produce an output statistics text with given file name")
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse the summary TSV file and obtain the read lengths
    statsList = []
    numSeqs = 0
    totalSize = 0
    firstLine = True
    with open(args.summaryTSV, "r") as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n ").split("\t")
            # Validate the header line
            if firstLine:
                assert sl[0] == "filename" and "sequence_length_template" in sl, \
                    "Invalid header line in summary TSV file"
                lengthIndex = sl.index("sequence_length_template")
                firstLine = False
            else:
                length = int(sl[lengthIndex])
                statsList.append(length)
                numSeqs += 1
                totalSize += length
    # Calculate additional statistics
    totalSize = locale.format_string("%d", totalSize, grouping=True)
    numSeqs = locale.format_string("%d", numSeqs, grouping=True)
    shortest = locale.format_string("%d", min(statsList), grouping=True)
    longest = locale.format_string("%d", max(statsList), grouping=True)
    n50 = locale.format_string("%d", N50(statsList), grouping=True)
    medianStat = locale.format_string("%d", median(statsList), grouping=True)
    meanStat = locale.format_string("%d", mean(statsList), grouping=True)
    # Print statistics
    print('Total amount of sequence (bp): ' + totalSize)
    print('Number of reads: ' + numSeqs)
    print('Shortest reads: ' + shortest)
    print('Longest reads: ' + longest)
    print('')
    print('N50: ' + n50)
    print('Median: ' + medianStat)
    print('Mean: ' + meanStat)
    # File output
    if args.outputFileName != False:
        with open(args.outputFileName, 'w') as output:
            output.write('Total amount of sequence (bp): ' + totalSize + '\n')
            output.write('Number of reads: ' + numSeqs + '\n')
            output.write('Shortest reads: ' + shortest + '\n')
            output.write('Longest reads: ' + longest + '\n')
            output.write('\n')
            output.write('N50: ' + n50 + '\n')
            output.write('Median: ' + medianStat + '\n')
            output.write('Mean: ' + meanStat + '\n')

if __name__ == '__main__':
    main()

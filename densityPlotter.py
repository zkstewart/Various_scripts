#! python3
# densityPlotter.py
# Script to create visualisations of various density metrics
# for assessing hypotheses usually pertaining to the clustering
# of biological markers along specific chromosomal segments

import os, argparse, math
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import zscore
from Bio import SeqIO
from scipy.ndimage.filters import gaussian_filter1d

def validate_args(args):
    # Validate input data locations
    if len(args.inputFiles) < 2:
        print("This function only makes sense if you have multiple file inputs")
        print("Make sure to specify 2 or more files, and try again.")
        quit()
    for file in args.inputFiles:
        if not os.path.isfile(file):
            print(f'I am unable to locate the input file ({file})')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    if not os.path.isfile(args.genomeFasta):
        print(f'I am unable to locate the input genome FASTA file ({args.genomeFasta})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Handle numeric parameters
    if args.windowSize < 1:
        print("windowSize must be a positive integer")
        quit()
    if args.minimumContigSize < (args.windowSize * 2):
        print("minimumContigSize must be at least 2x the window size")
        quit()
    if args.smoothingSigma < 0:
        print("smoothingSigma must be a float greater than or equal to zero")
        quit()
    # Handle headers format
    if args.headers == None:
        args.headers = [ "default" for _ in range(len(args.inputFiles)) ]
    updatedHeaders = []
    for header in args.headers:
        if header == "default":
            updatedHeaders.append(["contig", "position", "index"])
        else:
            h = header.split(",")
            if len(h) != 2 and len(h) != 3:
                print(f"'{header}' value when split by commas doesn't give 2 or 3 values.")
                print("Your formatting must be wrong; try to fix this parameter and try again.")
                quit()
            updatedHeaders.append(h)
    args.headers = updatedHeaders
    # Handle file output
    if os.path.isdir(args.outputDirectory):
        print('The specified output directory already exists. This program will attempt to resume an existing run where possible.')
        print('This program will not allowing overwriting. If you want to re-do an analysis, stop now and fix this issue.')

def get_densities(inputFileList, headersList, lengthsDict, windowSize=100000):
    '''
    Parameters:
        vcfFile -- a string pointing to the VCF file containing SNP annotation
        lengthsDict -- a dictionary with structure like:
                       {
                           'contig1': intLength1,
                           'contig2': intLength2,
                           ...
                       }
        windowSize -- an integer value indicating what size to bin genes in
    '''
    densityDicts = {}
    for i in range(len(inputFileList)):
        inputFile = inputFileList[i]
        firstLine = True
        densityDicts[i] = {}
        
        with open(inputFile, "r") as fileIn:
            for line in fileIn:
                sl = line.rstrip("\r\n ").split("\t")
                
                # Handle header line
                if line.startswith("#"):
                    continue
                if firstLine == True:
                    headers = headersList[i]
                    assert all([ header in sl for header in headers ]), \
                        f"Not all headers ({headers}) were found in '{inputFile}'!"
                    headerIndices = [ sl.index(header) for header in headers ]
                    firstLine = False
                    continue
                
                # Skip blank lines
                if sl == []:
                    continue
                
                # Extract information from content lines
                contigID = sl[headerIndices[0]]
                position = sl[headerIndices[1]]
                statistic = None if len(headerIndices) == 2 else float(sl[headerIndices[2]])
                windowChunkIndex = math.floor(int(position) / windowSize)
                
                # Handle occurrence content lines
                if len(headerIndices) == 2:
                    densityDicts[i].setdefault(contigID,
                        [ 0
                        for windowChunk in range(math.ceil(lengthsDict[contigID] / windowSize))
                        ]
                    )
                    densityDicts[i][contigID][windowChunkIndex] += 1
                
                # Handle statistic content lines
                elif len(headerIndices) == 3:
                    densityDicts[i].setdefault(contigID, {
                        "statistics": [ 0
                            for windowChunk in range(math.ceil(lengthsDict[contigID] / windowSize))
                        ],
                        "counts": [ 0
                            for windowChunk in range(math.ceil(lengthsDict[contigID] / windowSize))
                        ],
                        }
                    )
                    densityDicts[i][contigID]["statistics"][windowChunkIndex] += statistic
                    densityDicts[i][contigID]["counts"][windowChunkIndex] += 1
    
    # Average the SNP-index per window (if relevant)
    for i, densityDict in densityDicts.items():
        # Check if this is a statistic dict
        isStatistic = False
        for value in densityDict.values():
            if "statistics" in value and "counts" in value:
                isStatistic = True
            break
        
        # If it is a statistic dict, get the average over the window region
        if isStatistic:
            newDensityDict = {"dictType": "statistic"}
            for chrom in densityDict.keys():
                newDensityDict.setdefault(chrom, [])
                for windowChunkIndex in range(len(densityDict[chrom]["statistics"])):
                    try:
                        average = densityDict[chrom]["statistics"][windowChunkIndex] / densityDict[chrom]["counts"][windowChunkIndex]
                    except:
                        average = 0.0
                    newDensityDict[chrom].append(average)
            densityDicts[i] = newDensityDict
        else:
            densityDicts[i]["dictType"] = "occurrence"
    
    return densityDicts

def normalise_density(densityList):
    '''
    Using min-max normalisation
    
    Parameters:
        densityList -- a list containing integers of any length
    '''
    minValue, maxValue = min(densityList), max(densityList)
    normalisedDensityList = [ (x - minValue) / (maxValue - minValue) for x in densityList ]
    return normalisedDensityList

def plot_density(densityDicts, inputFiles, outputFileName, contigID, windowSize, yAxisText, smoothingSigma):
    '''
    Decides how to handle the plotting depending on the type of input data.
    Basically a Pythonic method overload.
    '''
    # Check if we're trying to plot multiple lines
    for key, value in densityDicts.items():
        if isinstance(value, dict):
            plot_density_multiple(densityDicts, inputFiles, outputFileName, contigID, windowSize, yAxisText, smoothingSigma)
        else:
            plot_density_single(densityDicts, outputFileName, contigID, windowSize, yAxisText, smoothingSigma)
        
        break # we only iterate to grab the value, we need to exit out after the above runs

def plot_density_multiple(densityDicts, inputFiles, outputFileName, contigID, windowSize, yAxisText, smoothingSigma):
    # Skip if we found nothing on this contig
    if not any([ contigID in densityDicts[i] for i in densityDicts.keys() ]):
        print(f"WARNING: '{contigID}' is in your genome file exceeding the minimum size but " +
            "has nothing associated with it; skipping...")
        return
    
    # Configure plot
    kbpWindowSize = round(windowSize / 1000, 2)
    fig = plt.figure(figsize=(10,6))
    
    ax = plt.axes()
    ax.set_xlabel(f"Chromosomal position ({kbpWindowSize} kbp windows)", fontweight="bold")
    ax.set_ylabel(yAxisText, fontweight="bold")
    ax.set_title(f"{contigID} density plot", fontweight="bold")
    
    # Plot all density values
    for i in densityDicts:
        try:
            densityList = densityDicts[i][contigID] # this may error out
        except:
            continue
        
        # Normalise and smooth density values
        if smoothingSigma > 0:
            densityList = gaussian_filter1d(densityList, sigma=smoothingSigma)
        
        if densityDicts[i]["dictType"] == "statistic":
            normalisedDensityList = densityList
        else:
            normalisedDensityList = normalise_density(densityList)
        
        # Plot the line
        filePrefix = inputFiles[i].rsplit(".", maxsplit=1)[0]
        ax.plot(normalisedDensityList, label=filePrefix, alpha=0.5)
    
    # Save output file
    plt.savefig(outputFileName)
    plt.close()

def plot_density_single(densityDict, outputFileName, contigID, windowSize, yAxisText, smoothingSigma):
    '''
    Parameters:
        densityDict -- a dictionary with structure like:
                       {
                           'contigID1': [ float1, float2, ... ],
                           'contigID2': [ ... ],
                           ...
                       }
    '''
    # Exit function if plotting is not possible
    if contigID not in densityDict:
        return
    
    # Configure plot
    kbpWindowSize = round(windowSize / 1000, 2)
    fig = plt.figure(figsize=(10,6))
    
    ax = plt.axes()
    ax.set_xlabel(f"Chromosomal position ({kbpWindowSize} kbp windows)", fontweight="bold")
    ax.set_ylabel(yAxisText, fontweight="bold")
    ax.set_title(f"{contigID} density plot", fontweight="bold")
    
    # Normalise and smooth density values
    densityList = densityDict[contigID]
    if smoothingSigma > 0:
        densityList = gaussian_filter1d(densityList, sigma=smoothingSigma)
    
    if densityDict["dictType"] == "statistic":
        normalisedDensityList = densityList
    else:
        normalisedDensityList = normalise_density(densityList)
    
    # Plot the line
    ax.plot(normalisedDensityList)
    
    # Save output file
    plt.savefig(outputFileName)
    plt.close()

def normalise_dict(dicts, contigs, normFunction=zscore):
    '''
    Parameters:
        dicts -- a dictionary with structure like:
                 {
                     index1: {
                         'dictType': 'statistic' OR 'occurrence',
                         'contig1': [ float1, float2, ... ],
                         'contig2': [ ... ],
                         ...
                     },
                     index2: {
                         ...
                     },
                     ...
                 }
    '''
    outDict = {}
    
    # Make a pandas df per contig
    for contigID, numWindows in contigs.items():
        dfList = []
        indices = []
        for index, indexDict in dicts.items():
            if contigID in indexDict:
                dfList.append(indexDict[contigID])
                indices.append(index)
        df = pd.DataFrame(dfList, index=indices, columns=[f"bin{i}" for i in range(numWindows)])
        
        # Normalise
        df = df.apply(normFunction, axis=0)
        
        # Convert back into a dict with the original format
        for index, indexDict in dicts.items():
            outDict.setdefault(index, {"dictType": indexDict["dictType"]})
            if contigID in dicts[index]:
                outDict[index][contigID] = df.loc[[index]].values.flatten().tolist()
    return outDict

def average_dict(dicts, contigs, dictType):
    '''
    Parameters:
        dicts -- a dictionary with structure like:
                 {
                     index1: {
                         'dictType': 'statistic' OR 'occurrence',
                         'contig1': [ float1, float2, ... ],
                         'contig2': [ ... ],
                         ...
                     },
                     index2: {
                         ...
                     },
                     ...
                 }
    '''
    outDict = {"dictType": dictType}
    
    # Make a pandas df per contig
    for contigID, numWindows in contigs.items():
        dfList = []
        indices = []
        for index, indexDict in dicts.items():
            if contigID in indexDict:
                dfList.append(indexDict[contigID])
                indices.append(index)
        df = pd.DataFrame(dfList, index=indices, columns=[f"bin{i}" for i in range(numWindows)])
        
        # Average
        averages = list(df.median())
        
        # Convert back into a dict with flattened structure
        outDict[contigID] = averages
    return outDict

def main():
    usage = """%(prog)s receives multiple TSV files containing at least two
    columns indicating the contig ID and position in nucleotides along said
    contig, with an optional third containing a test statistic that's likely
    to be on a 0->1 ratio.
    
    For each file, you can specify whether your input should have two or
    three columns read. In the first instance, you'll have an occurrence
    density line; in the second instance, you'll have the average of the
    test statistic plotted.
    
    This script will plot them all together to allow for visual assessment
    of whether these values cluster somehow along chromosomes.
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="inputFiles",
                   nargs="+",
                   required=True,
                   help="Specify the location of TSV files for combined plotting")
    p.add_argument("-f", dest="genomeFasta",
                   required=True,
                   help="Specify the location of the genome FASTA file")
    p.add_argument("-o", dest="outputDirectory",
                   required=True,
                   help="Output directory where plot files will be written")
    # Behavioural opts
    p.add_argument("--window_size", dest="windowSize",
                   type=int,
                   required=False,
                   help="""Optionally, specify the size of the window to sum
                   SNPs within (default=100000)""",
                   default=100000)
    p.add_argument("--minimum_contig", dest="minimumContigSize",
                   type=int,
                   required=False,
                   help="""Optionally, specify the minimum size of contigs which
                   should have plots created for (default=1000000)"; this value
                   must exceed window_size by at least 2x""",
                   default=1000000)
    p.add_argument("--smoothing", dest="smoothingSigma",
                   type=float,
                   required=False,
                   help="""Optionally, specify how much smoothing should be applied
                   to the plot; should be float value >= 0.0 (default=0.0)""",
                   default=0.0)
    # File parsing opts
    p.add_argument("--headers", dest="headers",
                   nargs="+",
                   required=False,
                   help="""Optionally, if all of your files do not fit within the default
                   header format of 'contig position index', then specify each files
                   header here with space separation of files, and comma separation of
                   header values; can use 'default' if a file conforms to default expectation""",
                   default=None)
    
    args = p.parse_args()
    validate_args(args)
    os.makedirs(args.outputDirectory, exist_ok=True)
    
    # Get contig lengths from genome FASTA
    genomeRecords = SeqIO.parse(open(args.genomeFasta, 'r'), "fasta")
    lengthsDict = { record.id:len(record) for record in genomeRecords }   
    
    # Parse each input file
    densityDicts = get_densities(args.inputFiles, args.headers, lengthsDict, args.windowSize)
    
    # Separate out density values by type
    statisticsDicts = { i : densityDicts[i] for i in densityDicts if densityDicts[i]["dictType"] == "statistic" }
    occurrenceDicts = { i : densityDicts[i] for i in densityDicts if densityDicts[i]["dictType"] == "occurrence" }
    
    # Get the average per dict in each window
    statisticsContigs = { contig:len(statisticsDicts[i][contig]) for i in statisticsDicts for contig in statisticsDicts[i] if contig != "dictType" }
    occurrenceContigs = { contig:len(occurrenceDicts[i][contig]) for i in occurrenceDicts for contig in occurrenceDicts[i] if contig != "dictType" }
    
    normStatisticsDicts = average_dict(statisticsDicts, statisticsContigs, "statistic")
    normOccurrenceDicts = average_dict(occurrenceDicts, occurrenceContigs, "occurrence")
    
    # Create plot per contig
    numContigsProcessed = 0
    numContigsPlotted = 0
    for contigID, length in lengthsDict.items():
        if length >= args.minimumContigSize:
            numContigsProcessed += 1
            
            # Plot statistics
            statFileOut = os.path.join(args.outputDirectory, f"{contigID}_statistic.png")
            if os.path.isfile(statFileOut):
                print(f"WARNING: Statistics plot for '{contigID}' already found in output directory; skipping...")
            else:
                plot_density(
                    normStatisticsDicts, args.inputFiles,
                    statFileOut, contigID, args.windowSize,
                    "Averaged statistic per window", args.smoothingSigma)
            
            # Plot occurrence
            occurrenceFileOut = os.path.join(args.outputDirectory, f"{contigID}_occurrence.png")
            if os.path.isfile(occurrenceFileOut):
                print(f"WARNING: Occurrence plot for '{contigID}' already found in output directory; skipping...")
            else:
                plot_density(
                    normOccurrenceDicts, args.inputFiles,
                    occurrenceFileOut, contigID, args.windowSize,
                    "Min-max normalised number per window", args.smoothingSigma)
            
            numContigsPlotted += 1
    
    # Raise relevant warnings
    if numContigsProcessed == 0:
        print(f"WARNING: We didn't find any contigs which exceeded {args.minimumContigSize}bp in size")
        print("Hence, no output files have been generated! Maybe you should fix your --minimum_contig value?")
    elif numContigsPlotted == 0:
        print("WARNING: We ended up skipping every contig! This means the program has already run to completion previously.")
        print("Hence, no new output files have been generated! Maybe you should delete the existing folder to restart?")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

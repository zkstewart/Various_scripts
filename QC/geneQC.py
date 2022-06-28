#! python 3
# geneQC.py
# A script/program designed to compute QC metrics
# using only a BAM file and a GFF3. It will be flexible
# and easy to use unlike every other approach available.

import os, argparse, sys
import matplotlib.pyplot as plt
import numpy as np

sys.path.append(os.path.dirname(os.path.dirname(__file__))) # 2 dirs up is where we find dependencies
from Function_packages import ZS_GFF3IO, ZS_BAMIO


def validate_args(args):
    # Validate input data locations
    if not os.path.isfile(args.bamFile):
        print('I am unable to locate the input protein BAM file (' + args.bamFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.gff3File):
        print('I am unable to locate the input GFF3 file (' + args.gff3File + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if args.mappingFile != None:
        if not os.path.isfile(args.mappingFile):
            print('I am unable to locate the input ID mapping file (' + args.mappingFile + ')')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    # Handle file output
    if os.path.isdir(args.outputDirectory):
        print('The specified output directory already exists. This program will attempt to resume an existing run where possible.')
        print('This program will not allowing overwriting. If you want to re-do an analysis, stop now and fix this issue.')

def parse_mapping_file(mappingFile):
    mappingDict = {}
    with open(mappingFile, "r") as fileIn:
        for line in fileIn:
            line = line.strip("\r\n ")
            if line.startswith("#") or line == "":
                continue
            
            sl = line.split("\t")
            mappingDict[sl[0]] = sl[1]
    return mappingDict

def plot_genebody_coverage(bamObj, outputFileName):
    '''
    Parameters:
        bamObj -- a ZS_BAMIO.BAM object which has the .gbc field set
                  with a Pandas DataFrame that has had all that good
                  min-max normalisation done on it.
        outputFileName -- where to write the plot figure to.
    '''
    # Summarise genebody coverage statistics into per-bin values
    summarisedBinValues = []
    for columnNumber in bamObj.gbc.columns:
        summarisedBinValues.append(np.average(bamObj.gbc[columnNumber]))
    
    # Plot as line chart
    fig = plt.figure()
    ax = plt.axes()
    
    ax.set_xlabel("Gene coverage bins (%)", fontweight="bold")
    ax.set_ylabel("Min-max normalised coverage", fontweight="bold")
    ax.set_title("{0} genebody coverage".format(os.path.basename(bamObj.fileLocation)),
                 fontweight="bold"
    )
    
    binSizes = int(100 / len(summarisedBinValues))
    ax.set_xticklabels(np.arange(-binSizes, 100 + binSizes, 2*binSizes))
    
    ax.plot(summarisedBinValues)
    
    plt.savefig(outputFileName)

def main():
    usage = """%(prog)s receives a BAM file in addition to possible two more files,
    and performs some post-alignment QC checks of the BAM file. These include 1)
    a genebody coverage calculation, and 2) ... forthcoming?
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-b", dest="bamFile", required=True,
                help="Specify the location of the input BAM file")
    p.add_argument("-g", dest="gff3File", required=True,
                help="Specify the location of the GFF3 with gene annotations")
    p.add_argument("-o", dest="outputDirectory", required=True,
                help="Output directory where QC files will be written")
    # Opts
    p.add_argument("-m", dest="mappingFile", required=False,
                help="Optionally, specify a TSV file mapping BAM IDs (left) to GFF3 IDs (right)")
    args = p.parse_args()
    validate_args(args)
    
    class Eg:
        def __init__(self):
            self.gff3File = r"F:\flies\chloe_2022\1_qc\GCF_016617805.1_CSIRO_BtryS06_freeze2_genomic.fixed.gff"
            self.bamFile = r"F:\flies\chloe_2022\1_qc\development\bowtie2\2A_output.sorted.bam"
            self.mappingFile = r"F:\flies\chloe_2022\1_qc\bowtie2_mapping_dict.txt"
    
    args = Eg()
    
    # Parse GFF3
    gff3Obj = ZS_GFF3IO.GFF3(args.gff3File)
    
    # Parse BAM
    bamObj = ZS_BAMIO.BAM(args.bamFile)
    
    # Parse mapping dict (if needed)
    if args.mappingFile != None:
        mappingDict = parse_mapping_file(args.mappingFile)
    
    # Obtain genebody coverage statistics
    bamObj.compute_coverage()
    bamObj.summarise_coverage_into_histogram()
    bamObj.qc_genebody_coverage() # bamObj.gbc is now set with a Pandas DataFrame
    
    # Plot it!
    os.makedirs(args.outputDirectory, exist_ok=True)
    plotFileName = os.path.join(
        args.outputDirectory,
        "{0}.png".format(os.path.basename(bamObj.fileLocation))
    )
    plot_genebody_coverage(bamObj, plotFileName)
    

if __name__ == "__main__":
    main()
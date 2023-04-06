#! python 3
# geneQC.py
# A script/program designed to compute QC metrics
# using only a BAM file and a GFF3. It will be flexible
# and easy to use unlike every other approach available.

import os, argparse, sys
import matplotlib.pyplot as plt
import numpy as np

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) # 2 dirs up is where we find dependencies
from Function_packages import ZS_GFF3IO, ZS_BAMIO

def validate_args(args):
    # Validate input data locations
    if not os.path.isfile(args.bamFile):
        print('I am unable to locate the input protein BAM file (' + args.bamFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if args.gff3File != None:
        if not os.path.isfile(args.gff3File):
            print('I am unable to locate the input GFF3 file (' + args.gff3File + ')')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    if args.mappingFile != None:
        if not os.path.isfile(args.mappingFile):
            print('I am unable to locate the input ID mapping file (' + args.mappingFile + ')')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    # Handle program behaviour parameters
    if args.bamType != None:
        if args.isChromosomes != False:
            print("When specifying --bamType, you should not specify --isChromosomes")
            quit()
        if args.needsFlip != False:
            print("When specifying --bamType, you should not specify --needsFlip")
            quit()
        
        # Set isChromosomes and needsFlip depending on the input type
        if args.bamType == "genome":
            args.isChromosomes = True
            args.needsFlip = True
        elif args.bamType == "rna":
            args.needsFlip = True
        elif args.bamType == "mrna":
            pass
    
    # Make sure the GFF3 argument is set (if relevant)
    if (args.isChromosomes != False or args.needsFlip != False) and args.gff3File == None:
        print("You've specified either --isChromosomes or --needsFlip, but haven't given a GFF3 file")
        print("These behaviours require a GFF3 to work; fix this and try again.")
        quit()
    
    # Make sure the GFF3 argument is NOT set (if relevant)
    if args.bamType == "mrna" and args.gff3File != None:
        print("You've specified --bamType as mrna, and have given a GFF3 file as input")
        print(
            "This behaviour does not require a GFF3 input. Out of an abundance of caution," +
            " I'm going to stop here so you can double check what you want this program to do."
        )
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
        summarisedBinValues.append(np.median(bamObj.gbc[columnNumber]))
    
    # Plot as line chart
    fig = plt.figure()
    ax = plt.axes()
    
    ax.set_xlabel("Gene coverage bins (%)", fontweight="bold")
    ax.set_ylabel("Median min-max normalised coverage", fontweight="bold")
    ax.set_title("{0} genebody coverage".format(os.path.basename(bamObj.fileLocation)),
                 fontweight="bold"
    )
    
    binSizes = int(100 / len(summarisedBinValues))
    ax.set_xticklabels(np.arange(-binSizes, 100 + binSizes, 2*binSizes))
    
    ax.plot(summarisedBinValues)
    plt.savefig(outputFileName)

def write_genebody_coverage(bamObj, perGeneFileName, summarisedFileName):
    '''
    Parameters:
        bamObj -- a ZS_BAMIO.BAM object which has the .gbc field set
                  with a Pandas DataFrame that has had all that good
                  min-max normalisation done on it.
        perGeneFileName -- where to write the statistics per-gene table.
        summarisedFileName -- where to write the summarised statistics table.
    '''
    # Summarise genebody coverage statistics into per-bin values
    summarisedBinValues = []
    for columnNumber in bamObj.gbc.columns:
        summarisedBinValues.append(np.median(bamObj.gbc[columnNumber]))
    
    # Write as tables
    bamObj.gbc.to_csv(perGeneFileName, sep="\t")
    with open(summarisedFileName, "w") as fileOut:
        fileOut.write("\t{0}\n".format("\t".join([str(num) for num in bamObj.gbc.columns]))) # header
        fileOut.write("median_cov_pct\t{0}\n".format("\t".join([str(num) for num in summarisedBinValues]))) # result

def main():
    usage = """%(prog)s receives a BAM file in addition to possibly two more files,
    and performs some post-alignment QC checks of the BAM file. At this time, only
    genebody coverage calculation is performed.
    
    For the --bamType parameter, "genome" should be used if your reads are mapped
    directly against the genome; this sets --isChromosomes and --needsFlip, and requires
    a GFF3 input.
    
    "mrna" should be used if your reads are aligned to mRNA models which means that
    each sequence is in the 5' to 3' orientation; this omits both --isChromosomes and
    --needsFlip, and does NOT require a GFF3 input.
    
    "rna" should be used instead when the sequences which were mapped against are
    in the genomic strand's orientation, rather than in the reading orientation of
    the ORF (i.e., 5' -> 3'); this sets --needsFlip, and requires a GFF3 input.
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-b", dest="bamFile",
                   required=True,
                   help="Specify the location of the input BAM file")
    p.add_argument("-o", dest="outputDirectory",
                   required=True,
                   help="Output directory where QC files will be written")
    # Opts
    p.add_argument("-g", "--gff3", dest="gff3File",
                   required=False,
                   help="Optionally, specify the location of the GFF3 with gene annotations")
    p.add_argument("-m", "--mapping", dest="mappingFile",
                   required=False,
                   help="""Optionally, specify a TSV file mapping BAM contig IDs (left) to
                   GFF3 contig IDs (right) if these are not identical for some reason""")
    p.add_argument("--bamType", dest="bamType",
                   required=False,
                   choices=["genome", "mrna", "rna"],
                   help="""Rather than setting the below two parameters yourself,
                   indicate what type of BAM file you are working with to have those
                   parameters picked for you""",
                   default=None)
    p.add_argument("--isChromosomes", dest="isChromosomes",
                   required=False,
                   action="store_true",
                   help="""Optionally, if your reads were aligned to the genome,
                   specify this flag to extract reads that overlap predicted gene
                   regions in the GFF3""",
                   default=False)
    p.add_argument("--needsFlip", dest="needsFlip",
                   required=False,
                   action="store_true",
                   help="""Optionally if you aligned to the genome OR your reads
                   were aligned to RNA models which are in an orientation dictated
                   by the genomic strand (and not all in the 5'->3' orientation),
                   specify this flag so we can check the GFF3 and perform flipping
                   of the counts""",
                   default=False)
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse BAM
    bamObj = ZS_BAMIO.BAM(args.bamFile)
    
    # Parse GFF3 (if needed)
    if args.gff3File != None:
        gff3Obj = ZS_GFF3IO.GFF3(args.gff3File, strict_parse=False)
    
    # Parse mapping dict (if needed)
    if args.mappingFile != None:
        mappingDict = parse_mapping_file(args.mappingFile)
    
    # Obtain genebody coverage statistics
    bamObj.compute_coverage(
        gff3Obj = None if args.gff3File == None else gff3Obj if args.isChromosomes else None
    )
    bamObj.summarise_coverage_into_histogram(
        mappingDict = None if args.mappingFile == None else mappingDict,
        gff3Obj = None if args.gff3File == None else gff3Obj if args.needsFlip else None
    )
    bamObj.qc_genebody_coverage() # bamObj.gbc is now set with a Pandas DataFrame
    
    os.makedirs(args.outputDirectory, exist_ok=True)
    # Plot the genebody coverage
    plotFileName = os.path.join(
        args.outputDirectory,
        "{0}.png".format(os.path.basename(bamObj.fileLocation))
    )
    if not os.path.isfile(plotFileName):
        plot_genebody_coverage(bamObj, plotFileName)
    
    # Write an output table file
    perGeneFileName = os.path.join(
        args.outputDirectory,
        "{0}.genebody.perGene.tsv".format(os.path.basename(bamObj.fileLocation))
    )
    summaryFileName = os.path.join(
        args.outputDirectory,
        "{0}.genebody.summary.tsv".format(os.path.basename(bamObj.fileLocation))
    )
    if not os.path.isfile(perGeneFileName) and not os.path.isfile(summaryFileName):
        write_genebody_coverage(bamObj, perGeneFileName, summaryFileName)
        print("Program completed successfully!")

if __name__ == "__main__":
    main()

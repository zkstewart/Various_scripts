#! python3
# plot_diffratio.py
# Script to create visualisations of the difference ratio statistics
# for assessing hypotheses of variant segregation along chromosomes.

import os, argparse, sys, re, pickle, math
import numpy as np
from Bio import SeqIO

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from euclidean_dist.plot_ed import get_statistics_for_dotting, get_sorted_contig_ids, \
    lineplot_per_contig, lineplot_horizontal, lineplot_regions, \
    histo_per_contig, histo_horizontal, histo_regions

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def validate_args(args):
    # Validate input data locations
    if not os.path.isfile(args.diffratioFile):
        eprint(f'I am unable to locate the difference ratio file ({args.diffratioFile})')
        eprint('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.genomeFasta):
        eprint(f'I am unable to locate the input genome FASTA file ({args.genomeFasta})')
        eprint('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Handle numeric parameters
    if args.width < 0:
        eprint("width must be a positive integer")
        quit()
    if args.height < 0:
        eprint("height must be a positive integer")
        quit()
    if args.bulkAlleles != []:
        if len(args.bulkAlleles) != 2:
            eprint("bulkAlleles must be a list of two integers")
            quit()
        if any([ not val >= 2 for val in args.bulkAlleles ]):
            eprint("bulkAlleles values must be integers >= 2")
            quit()
    if args.bulkOccurrence != None:
        if args.bulkOccurrence < 0 or args.bulkOccurrence > 1:
            eprint("bulkOccurrence must be a float value >0 and <=1")
    if args.reportAboveCutoff != None:
        if 0 > args.reportAboveCutoff:
            eprint("reportAboveCutoff must be a float or int value >0")
            quit()
    # Check for conflicting arguments
    if args.onePlot and args.regions != []:
        eprint("You can't provide both --onePlot and --regions; please choose one and try again.")
        quit()
    if args.bulkAlleles != [] and args.bulkOccurrence == None:
        eprint("You must provide a --bulkOccurrence value if you provide --bulkAlleles; please try again.")
        quit()
    if args.bulkAlleles == [] and args.bulkOccurrence != None:
        eprint("You must provide --bulkAlleles if you provide --bulkOccurrence; please try again.")
        quit()
    # Handle regions
    for region in args.regions:
        if not re.match(r"^.+:\d+:\d+$", region):
            eprint(f"Region '{region}' is not in the expected format (contig:start:end)")
            eprint("Please provide regions in the format 'contig:start:end' and try again.")
            quit()
    # Handle file output
    if os.path.isdir(args.outputDirectory) and os.listdir(args.outputDirectory) != []:
        eprint(f"Output directory '{args.outputDirectory}' already exists; I'll write output files here.")
        eprint("But, I won't overwrite any existing files, so beware that if a previous run had issues, " +
              "you may need to delete/move files first.")
    if not os.path.isdir(args.outputDirectory):
        os.makedirs(args.outputDirectory)
        eprint(f"Output directory '{args.outputDirectory}' has been created as part of argument validation.")
    # Handle mode-specific arguments
    if args.mode == "line":
        if args.wmaSize < 1:
            eprint("wmaSize must be an integer >= 1")
            quit()
        if args.linewidth < 1:
            eprint("linewidth must be an integer >= 1")
            quit()
    elif args.mode == "histogram":
        if args.binSize < 1:
            eprint("binSize must be an integer >= 1")
            quit()
        if args.binThreshold <= 0:
            eprint("binThreshold must be a float value >0")
            quit()

def get_sorted_contig_ids(idsList):
    # Sort contig IDs by their numerical value (if possible)
    allHaveNumbers = all([ any([ c.isdigit() for c in contigID ]) for contigID in idsList ])
    if allHaveNumbers:
        numRegex = re.compile(r"\d+")
        return sorted(idsList, key=lambda x: int("".join(numRegex.findall(x))))
    else:
        return sorted(idsList)

def main():
    usage = """%(prog)s receives a difference ratio TSV file and creates
     difference ratio plots in various formats (line plots and histograms)
    in different manners (per contig, horizontally, or for specific regions).
    Through calculating a weighted moving average (for line plots), it can
    help to visualise where in the genome regions that segregate between
    bulks occur.
    """
    # Establish main parser
    p = argparse.ArgumentParser(description=usage)
    
    # Set arguments shared by subparsers
    ## Required arguments
    p.add_argument("-d", dest="diffratioFile",
                    required=True,
                    help="Specify the location of the input difference ratio file")
    p.add_argument("-f", dest="genomeFasta",
                    required=True,
                    help="Specify the location of the genome FASTA file")
    p.add_argument("-o", dest="outputDirectory",
                    required=True,
                    help="Output directory where plot files will be written")
    ## Opts (metadata behaviour)
    p.add_argument("--bulkAlleles", dest="bulkAlleles",
                    required=False,
                    nargs="+",
                    type=int,
                    help="""Optionally, indicate the number of maximum possible alleles
                    in each bulk in order to calculate the occurrence fraction for
                    filtering""",
                    default=[])
    p.add_argument("--bulkOccurrence", dest="bulkOccurrence",
                    type=float,
                    required=False,
                    help="""Optionally, specify the minimum fraction of occurrence
                    for one of the two bulks to be considered for plotting""",
                    default=None)
    ## Opts (plotting behaviour)
    p.add_argument("--width", dest="width",
                    type=int,
                    required=False,
                    help="""Optionally, specify the output plot width (default=10)""",
                    default=10)
    p.add_argument("--height", dest="height",
                    type=int,
                    required=False,
                    help="""Optionally, specify the output plot height (default=6)""",
                    default=6)
    ## Opts (statistical behaviour)
    p.add_argument("--minimum_contig", dest="minimumContigSize",
                    type=int,
                    required=False,
                    help="""Optionally, specify the minimum size of contig to
                    create plots for (default=200000 i.e., 2Mb)""",
                    default=200000)
    ## Opts (output)
    p.add_argument("--regions", dest="regions",
                    required=False,
                    nargs="+",
                    help="""Optionally, indicate one or more regions to plot in greater detail
                    by providing the contig ID and start and end positions in bp (e.g.
                    contig1:10000:20000)""",
                    default=[])
    p.add_argument("--onePlot", dest="onePlot",
                    required=False,
                    action="store_true",
                    help="""Optionally, provide this flag if you want a single plot to be
                    produced with all chromosomes positioned horizontally""",
                    default=False)
    p.add_argument("--pdf", dest="plotPDF",
                    required=False,
                    action="store_true",
                    help="""Optionally, provide this flag if you want outputs to be
                    in PDF format instead of PNG format""",
                    default=False)
    p.add_argument("--reportAboveCutoff", dest="reportAboveCutoff",
                    required=False,
                    type=float,
                    help="""Optionally, indicate a cutoff value for which any variant statistic
                    ('line' mode) or bin quantity ('histogram' mode) above this value will
                    be reported to the console""",
                    default=None)
    
    # Establish subparsers
    subParentParser = argparse.ArgumentParser()
    subparsers = subParentParser.add_subparsers(dest="mode",
                                                required=True)
    
    lineparser = subparsers.add_parser("line",
                                       parents=[p],
                                       add_help=False,
                                       help="Create line plots of difference ratios")
    lineparser.set_defaults(func=linemain)
    
    histoparser = subparsers.add_parser("histogram",
                                        aliases=["histo"],
                                        parents=[p],
                                        add_help=False,
                                        help="Create histograms of difference ratios")
    histoparser.set_defaults(func=histomain)
    
    # Line-subparser arguments
    lineparser.add_argument("--wmaSize", dest="wmaSize",
                            type=int,
                            required=False,
                            help="""Optionally, specify the number of previous values to consider
                            during weighted moving average calculation (default=5)""",
                            default=5)
    lineparser.add_argument("--linewidth", dest="linewidth",
                            type=int,
                            required=False,
                            help="""Optionally, specify the line width (default=1)""",
                            default=1)
    lineparser.add_argument("--showDots", dest="showDots",
                            required=False,
                            action="store_true",
                            help="""Optionally, provide this flag if you want plots to show dots
                            of each data point in addition to the line plot""",
                            default=False)
    
    # Histogram-subparser arguments
    histoparser.add_argument("--binSize", dest="binSize",
                            type=int,
                            required=False,
                            help="""Optionally, specify the bin size to count variants
                            within (default=10000)""",
                            default=10000)
    histoparser.add_argument("--binThreshold", dest="binThreshold",
                            type=float,
                            required=False,
                            help="""Optionally, specify the difference ratio threshold
                            to set for counting a SNP (default=0.4)""",
                            default=0.4)
    
    args = subParentParser.parse_args()
    validate_args(args)
    
    # Figure out what our pickle file should be called
    pickleFile = os.path.join(
        args.outputDirectory,
        f"{os.path.basename(args.diffratioFile)}.al{'_'.join(map(str, args.bulkAlleles))}.oc{args.bulkOccurrence}.pkl"
    )
    
    # Load pickle if it exists to skip computation
    if os.path.isfile(pickleFile):
        with open(pickleFile, "rb") as fileIn:
            dotsX, dotsY = pickle.load(fileIn)
    
    # Otherwise ...
    else:
        # Parse difference ratio data
        dotsX, dotsY = get_statistics_for_dotting(args.diffratioFile, args.bulkAlleles, args.bulkOccurrence,
                           HEADER_VALUES = ["CHROM", "POSI", "differenceRatio", "bulk1_alleles", "bulk2_alleles"])
        
        # Save data
        with open(pickleFile, "wb") as fileOut:
            pickle.dump([dotsX, dotsY], fileOut)
    
    # Get contig lengths from genome FASTA
    genomeRecords = SeqIO.parse(open(args.genomeFasta, 'r'), "fasta")
    lengthsDict = { record.id:len(record) for record in genomeRecords }
    
    # Drop any contigs which don't meet our length cutoff
    for contigID, length in lengthsDict.items():
        if length < args.minimumContigSize:
            print(f"NOTE: '{contigID}' is below the minimum contig size and will be skipped")
            try:
                del dotsX[contigID]
                del dotsY[contigID]
            except:
                raise ValueError(f"ERROR: '{contigID}' was not found in the difference ratio file but " +
                                 "was found in the genome FASTA file; this is unexpected and suggests " +
                                 "a mismatch between the two files")
    
    # Check that we still have contigs to plot
    if dotsX == {}:
        raise ValueError("ERROR: We didn't find any contigs which matched or exceeded the minimum size. " +
                         "Hence, no output files have been generated! Maybe you should fix your " +
                         "--minimum_contig value?")
    
    # Split into mode-specific functions
    if args.mode == "line":
        linemain(args, dotsX, dotsY)
    elif args.mode in ["histogram", "histo"]:
        histomain(args, dotsX, dotsY, lengthsDict)

def linemain(args, dotsX, dotsY):
    # Report any variants above the cutoff
    if args.reportAboveCutoff != None:
        print(f"# Difference ratio >= {args.reportAboveCutoff} report:")
        for contigID in dotsX.keys():
            for x, y in zip(dotsX[contigID], dotsY[contigID]):
                if y >= args.reportAboveCutoff:
                    print(f"# {contigID}:{x} = {y}")
    
    # Create plots
    if args.onePlot:
        numContigsPlotted = lineplot_horizontal(dotsX, dotsY, args.wmaSize,
                                                args.width, args.height, 1, # power=1
                                                args.outputDirectory, args.plotPDF,
                                                args.showDots, args.linewidth)
    elif args.regions != []:
        numContigsPlotted = lineplot_regions(dotsX, dotsY, args.regions, args.wmaSize,
                                             args.width, args.height, 1, # power=1
                                             args.outputDirectory, args.plotPDF,
                                             args.showDots, args.linewidth)
    else:
        numContigsPlotted = lineplot_per_contig(dotsX, dotsY, args.wmaSize,
                                                args.width, args.height, 1, # power=1
                                                args.outputDirectory, args.plotPDF,
                                                args.showDots, args.linewidth)
    
    # Raise errors if necessary
    if numContigsPlotted == 0:
        raise ValueError("ERROR: We ended up skipping every contig! This means the program has " + 
                         "already run to completion previously. Hence, no new output files have been " +
                         "generated! Maybe you should delete the existing files to restart?")
    
    print("Program completed successfully!")

def histomain(args, dotsX, dotsY, lengthsDict):
    # Bin data into histograms
    histoDict = {}
    for contigID in get_sorted_contig_ids(dotsX.keys()):
        histoDict[contigID] = np.array([ 0
                for windowChunk in range(math.ceil(lengthsDict[contigID] / args.binSize))
            ])
        for x, y in zip(dotsX[contigID], dotsY[contigID]):
            binIndex = x // args.binSize
            if y >= args.binThreshold:
                histoDict[contigID][binIndex] += 1
    
    # Report any bins above the cutoff
    if args.reportAboveCutoff != None:
        print(f"# Bin containing difference ratio >= {args.binThreshold} report:")
        for contigID, binCounts in histoDict.items():
            for i, count in enumerate(binCounts):
                if count >= args.reportAboveCutoff:
                    print(f"# {contigID}:{i*args.binSize}-{(i+1)*args.binSize} = {count}")
    
    # Create plots
    if args.onePlot:
        numContigsPlotted = histo_horizontal(histoDict,
                                             args.width, args.height, 1, # power=1
                                             args.outputDirectory, args.binSize,
                                             args.binThreshold, args.plotPDF)
    elif args.regions != []:
        numContigsPlotted = histo_regions(histoDict, args.regions,
                                          args.width, args.height, 1, # power=1
                                          args.outputDirectory, args.binSize,
                                          args.binThreshold, args.plotPDF)
    else:
        numContigsPlotted = histo_per_contig(histoDict,
                                             args.width, args.height, 1, # power=1
                                             args.outputDirectory, args.binSize,
                                             args.binThreshold, args.plotPDF)
    
    # Raise errors if necessary
    if numContigsPlotted == 0:
        raise ValueError("ERROR: We ended up skipping every contig! This means the program has " + 
                         "already run to completion previously. Hence, no new output files have been " +
                         "generated! Maybe you should delete the existing files to restart?")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

#! python3
# markers_from_ed.py
# Script to locate informative marker SNPs based on
# Euclidean distance statistics.

import os, argparse, sys, pickle, math
import numpy as np
from Bio import SeqIO
from plot_ed import get_statistics_for_dotting, get_sorted_contig_ids

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))
from Function_packages import ZS_VCFIO

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def validate_args(args):
    # Validate input data locations
    if not os.path.isfile(args.edistFile):
        eprint(f'I am unable to locate the Euclidean distance file ({args.edistFile})')
        eprint('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.vcfFile):
        eprint(f'I am unable to locate the input VCF file ({args.vcfFile})')
        eprint('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.genomeFasta):
        eprint(f'I am unable to locate the input genome FASTA file ({args.genomeFasta})')
        eprint('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if args.pickleFile != None:
        if not os.path.isfile(args.vcfFile):
            eprint(f'I am unable to locate the input VCF file ({args.vcfFile})')
            eprint('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    # Handle numeric parameters
    if args.binSize < 1:
        eprint("binSize must be an integer >= 1")
        quit()
    if args.binThreshold <= 0:
        eprint("binThreshold must be a float value >0")
        quit()
    if args.reportThreshold <= 0:
        eprint("reportThreshold must be a float value >0")
        quit()
    if args.power < 1:
        eprint("power must be an integer >= 1")
        quit()
    if args.minimumContigSize < 1:
        eprint("minimumContigSize must be an integer >= 1")
        quit()
    # Handle file output
    if os.path.isfile(args.outputFileName):
        eprint(f'File already exists at output location ({args.outputFileName})')
        eprint('Make sure you specify a unique file name and try again.')
        quit()
    elif not os.path.isdir(os.path.dirname(os.path.abspath(args.outputFileName))):
        eprint(f"Output directory '{os.path.dirname(os.path.abspath(args.outputFileName))}' doesn't exist")
        eprint("Please create the directory and try again.")
        quit()    

def main():
    usage = """%(prog)s receives a Euclidean distance TSV file and a VCF
    file and produces a table of markers and their associated information for
    manual inspection.
    """
    # Establish main parser
    p = argparse.ArgumentParser(description=usage)
    
    # Set arguments shared by subparsers
    ## Required arguments
    p.add_argument("-d", dest="edistFile",
                    required=True,
                    help="Specify the location of the input Euclidean distance file")
    p.add_argument("-f", dest="genomeFasta",
                    required=True,
                    help="Specify the location of the genome FASTA file")
    p.add_argument("-v", dest="vcfFile",
                    required=True,
                    help="Specify the location of the VCF file")
    p.add_argument("-o", dest="outputFileName",
                    required=True,
                    help="Output name for the marker selection file")
    p.add_argument("--binThreshold", dest="binThreshold",
                    type=True,
                    required=False,
                    help="""Specify the Euclidean distance threshold
                    to set for binning a SNP (default=0.4)""",
                    default=0.4)
    p.add_argument("--reportThreshold", dest="reportThreshold",
                    type=True,
                    required=False,
                    help="""Specify the Euclidean distance threshold
                    to set for outputting a SNP (default=0.1)""",
                    default=0.4)
    ## Opts
    p.add_argument("--minimum_contig", dest="minimumContigSize",
                    type=int,
                    required=False,
                    help="""Optionally, specify the minimum size of contig to
                    create plots for (default=200000 i.e., 2Mb)""",
                    default=200000)
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
    p.add_argument("--binSize", dest="binSize",
                    type=int,
                    required=False,
                    help="""Optionally, specify the bin size to count variants
                    within (default=10000)""",
                    default=10000)
    p.add_argument("--power", dest="power",
                    type=int,
                    required=False,
                    help="""Optionally, specify the power to raise Euclidean distances to
                    reduce noise (default=4)""",
                    default=4)
    p.add_argument("--pickle", dest="pickleFile",
                    required=False,
                    help="""Optionally, specify the location of a pickle file""",
                    default=None)
    
    args = p.parse_args()
    validate_args(args)
    
    # Load pickle if it exists to skip computation
    if args.pickleFile != None:
        with open(args.pickleFile, "rb") as fileIn:
            dotsX, dotsY = pickle.load(fileIn)
    
    # Otherwise, parse Euclidean distance data
    else:
        dotsX, dotsY = get_statistics_for_dotting(args.edistFile, args.bulkAlleles, args.bulkOccurrence)
    
    # Power-transform values
    powerY = {}
    for contigID in dotsX.keys():
        y = np.array(dotsY[contigID])**args.power
        powerY[contigID] = y
    
    # Get contig lengths from genome FASTA
    genomeRecords = SeqIO.parse(open(args.genomeFasta, 'r'), "fasta")
    lengthsDict = { record.id:len(record) for record in genomeRecords }
    
    # Drop any contigs which don't meet our length cutoff
    for contigID, length in lengthsDict.items():
        if length < args.minimumContigSize:
            print(f"NOTE: '{contigID}' is below the minimum contig size and will be skipped")
            try:
                del dotsX[contigID]
                del powerY[contigID]
            except:
                raise ValueError(f"ERROR: '{contigID}' was not found in the Euclidean distance file but " +
                                 "was found in the genome FASTA file; this is unexpected and suggests " +
                                 "a mismatch between the two files")
    
    # Check that we still have contigs to plot
    if dotsX == {}:
        raise ValueError("ERROR: We didn't find any contigs which matched or exceeded the minimum size. " +
                         "Hence, no output files has been generated! Maybe you should fix your " +
                         "--minimum_contig value?")
    
    # Bin data into histograms
    histoDict = {}
    for contigID in get_sorted_contig_ids(dotsX.keys()):
        histoDict[contigID] = np.array([ 0
                for windowChunk in range(math.ceil(lengthsDict[contigID] / args.binSize))
            ])
        for x, y in zip(dotsX[contigID], powerY[contigID]):
            binIndex = x // args.binSize
            if y >= args.binThreshold:
                histoDict[contigID][binIndex] += 1
    
    # Subset data based on report threshold
    subsetDict = {}
    contigsReported = False
    for contigID in histoDict.keys():
            for position, Y in zip(dotsX[contigID], powerY[contigID]):
                if Y < args.reportThreshold:
                    continue
                
                binIndex = position // args.binSize
                subsetDict.setdefault(contigID, {})
                subsetDict[contigID][position] = [Y, histoDict[contigID][binIndex] ]
                contigsReported = True
    
    # Raise errors if necessary
    if contigsReported == False:
        raise ValueError("ERROR: We ended up skipping every contig! This probably means " +
                         "that your --reportThreshold is too high. Try lowering it and " +
                         "running the program again.")
    
    # Parse through VCF for relevant genotypes and write output report
    with open(args.outputFileName, "w") as fileOut:
        vcfIterator = ZS_VCFIO.SimpleGenotypeIterator(args.vcfFile)
        samples = None
        for value in vcfIterator:
            # Write header line
            if samples == None:
                samples = value
                fileOut.write("{0}\n".format("\t".join([
                    "contig", "position", f"ed**{args.power}", "binSNPs",
                    *samples
                ])))
                continue
            
            contig, position, ref, alt, snpDict = value
            ref_alt = [ref, *alt]
            
            # Skip if we're not interested in this contig+pos
            if contig not in subsetDict or position not in subsetDict[contig]:
                continue
            
            # Format genotypes for each sample
            genotypes = [
                f"{ref_alt[snpDict[sample][0]]}/{ref_alt[snpDict[sample][1]]}" if sample in snpDict
                else "./."
                for sample in samples ]
            
            # Get ED data
            Y, binSNPs = subsetDict[contig][position]
            
            # Format and write output line
            fileOut.write("{0}\n".format("\t".join([
                contig, str(position), str(Y), str(binSNPs),
                *genotypes
            ])))
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

#! python3
# markers_matching_ed.py
# Script to locate genotyped positions from a VCF that
# match chromosome/position pairs from a marker set

import os, argparse, sys, pickle, math
import numpy as np
from Bio import SeqIO
from plot_ed import get_sorted_contig_ids

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))
from Function_packages import ZS_VCFIO

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def validate_args(args):
    # Validate input data locations
    if not os.path.isfile(args.markerFile):
        eprint(f'I am unable to locate the marker TSV file ({args.markerFile})')
        eprint('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.vcfFile):
        eprint(f'I am unable to locate the input VCF file ({args.vcfFile})')
        eprint('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Handle column headers
    if not len(args.columnHeaders) == 2:
        eprint('Please specify exactly two column headers for the chromosome and position columns')
        eprint('Make sure to separate the two headers with a space and try again.')
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
    usage = """%(prog)s receives a TSV file as produced by 'markers_from_ed.py'
    and a VCF file of a new experiment. It then locates the genotyped positions
    from the VCF that match chromosome/position pairs from the marker set. The
    output is a TSV file with the matching genotypes for each sample at any
    positions that occur as variants in the new experiment VCF.
    """
    p = argparse.ArgumentParser(description=usage)
    ## Reqs
    p.add_argument("-i", dest="markerFile",
                   required=True,
                   help="Specify the location of the input marker TSV file")
    p.add_argument("-v", dest="vcfFile",
                   required=True,
                   help="Specify the location of the VCF file")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output name for the matching markers file")
    ## Opts
    p.add_argument("--columns", dest="columnHeaders",
                   required=False,
                   help="""Optionally, specify the column headers for the
                   chromosome and position columns in the marker file - in
                   that order! Default is 'contig position'""",
                   default=["contig", "position"])
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse marker file
    markers = {}
    with open(args.markerFile, "r") as fileIn:
        firstLine = True
        for line in fileIn:
            if firstLine:
                sl = line.strip().split("\t")
                assert [x in sl for x in args.columnHeaders], \
                    "Column headers not found in marker file!"
                columnIndices = [sl.index(x) for x in args.columnHeaders]
                firstLine = False
            else:
                sl = line.strip().split("\t")
                contig, position = [sl[x] for x in columnIndices]
                markers.setdefault(contig, set()).add(int(position))
    
    # Parse through VCF for relevant genotypes and write output report
    with open(args.outputFileName, "w") as fileOut:
        vcfIterator = ZS_VCFIO.SimpleGenotypeIterator(args.vcfFile)
        samples = None
        for value in vcfIterator:
            # Write header line
            if samples == None:
                samples = value
                fileOut.write("{0}\n".format("\t".join([
                    "contig", "position", *samples
                ])))
                continue
            
            contig, position, ref, alt, snpDict = value
            ref_alt = [ref, *alt]
            
            # Skip if we're not interested in this contig+pos
            if contig not in markers or position not in markers[contig]:
                continue
            
            # Format genotypes for each sample
            genotypes = [
                f"{ref_alt[snpDict[sample][0]]}/{ref_alt[snpDict[sample][1]]}" if sample in snpDict
                else "./."
                for sample in samples ]
            
            # Format and write output line
            fileOut.write("{0}\n".format("\t".join([
                contig, str(position), *genotypes
            ])))
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

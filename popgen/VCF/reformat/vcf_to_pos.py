#! python3
# vcf_to_pos.py
# A script to produce pos format from a VCF file.

import os, argparse, gzip
from contextlib import contextmanager

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.vcfFile):
        print('I am unable to locate the VCF file (' + args.vcfFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print('File already exists at output location (' + args.outputFileName + ')')
        print('Make sure you specify a unique file name and try again.')
        quit()

@contextmanager
def open_vcf_file(filename):
    if filename.endswith(".gz"):
        with gzip.open(filename, "rt") as f:
            yield f
    else:
        with open(filename) as f:
            yield f

## Main
def main():
    # User input
    usage = """%(prog)s reads in a VCF to produce a pos format file
    for use with ngsLD. It will skip multiallelic and indel calls by
    default when generating the output file, but you can choose to retain
    them with the appropriate flag arguments. The output is a TSV where
    the first column indicates the chromosome and the second column indicates
    the position of the SNP. A header can be added to the output file with the
    --header flag. Make sure you specify the same --keepM and --keepI flags as
    you did when running vcf_to_geno.py
    """
    p = argparse.ArgumentParser(description=usage)
    # Reqs
    p.add_argument("-i", dest="vcfFile",
                   required=True,
                   help="Input VCF file for filtering")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="""Output name for the pos file""")
    # Opts
    p.add_argument("--header", dest="header",
                   required=False,
                   action="store_true",
                   help="Optionally, add a header to the output file",
                   default=False)
    p.add_argument("--keepM", dest="keepMultiallelic",
                   required=False,
                   action="store_true",
                   help="Optionally, keep multiallelic sites in the output file",
                   default=False)
    p.add_argument("--keepI", dest="keepIndels",
                   required=False,
                   action="store_true",
                   help="Optionally, keep indels in the output file",
                   default=False)
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse VCF while writing output to file
    with open_vcf_file(args.vcfFile) as fileIn, open(args.outputFileName, "w") as fileOut:
        for line in fileIn:
            if line.startswith("#"):
                continue
            sl = line.rstrip("\r\n ").split("\t")
            
            # Extract relevant details of the SNP
            chrom = sl[0]
            pos = sl[1]
            ref = sl[3]
            alt = sl[4].split(",")
            
            # Skip depending on --keepM and --keepI flags
            if (not args.keepMultiallelic) and len(alt) > 1:
                continue
            if (not args.keepIndels) and ( ref == "." or any([ x == "." for x in alt ]) or any([len(x) != len(ref) for x in alt]) ):
                continue
            
            # Write the pos line to the output file
            fileOut.write(f"{chrom}\t{pos}\n")
    
    print("Program finished successfully!")

if __name__ == "__main__":
    main()

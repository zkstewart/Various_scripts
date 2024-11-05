#! python3
# vcf_to_geno.py
# A script to produce geno format from a VCF file.

import os, argparse, gzip
from contextlib import contextmanager

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.vcfFile):
        print('I am unable to locate the VCF file (' + args.vcfFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate uncalled character argument
    if len(args.uncalledCharacter) < 1:
        print('Uncalled character must be one or more characters.')
        print('Please correct this and try again.')
        quit()
    # Validate output file location
    if not args.outputFileName.endswith(".gz"):
        args.outputFileName += ".gz"
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
    usage = """%(prog)s reads in a VCF to produce a geno format file
    for use with ngsLD. It will skip multiallelic and indel calls by
    default when generating the output file, but you can choose to retain
    them with the appropriate flag arguments. The output is a gzip'd TSV where each
    row is a SNP and each column is a sample whose value is 0, 1, or 2 indicating
    how many reference alleles are present in the sample at that SNP. Uncalled
    genotypes are represented by the character you specify with the --uncalled
    argument.
    """
    p = argparse.ArgumentParser(description=usage)
    # Reqs
    p.add_argument("-i", dest="vcfFile",
                   required=True,
                   help="Input VCF file for filtering")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="""Output file name for the filtered SNPs; files that 
                   don't end with .gz will have it appended to the end""")
    # Opts
    p.add_argument("--uncalled", dest="uncalledCharacter",
                   required=False,
                   help="""Optionally, set the character to use for uncalled
                   genotypes. Default is -1.""",
                   default="-1")
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
    with open_vcf_file(args.vcfFile) as fileIn, gzip.open(args.outputFileName, "wt") as fileOut:
        for line in fileIn:
            if line.startswith("#"):
                continue
            sl = line.rstrip("\r\n ").split("\t")
            
            # Extract relevant details of the SNP
            ref = sl[3]
            alt = sl[4].split(",")
            
            # Determine which field position we're extracting to get our GT value
            fieldsDescription = sl[8]
            if ":" not in fieldsDescription:
                gtIndex = -1
            else:
                gtIndex = fieldsDescription.split(":").index("GT")
            
            # Skip depending on --keepM and --keepI flags
            if (not args.keepMultiallelic) and len(alt) > 1:
                continue
            if (not args.keepIndels) and ( ref == "." or any([ x == "." for x in alt ]) or any([len(x) != len(ref) for x in alt]) ):
                continue
            
            # Format geno line for this variant
            geno = []
            for sampleResult in sl[9:]:
                # Grab our genotype
                if gtIndex != -1:
                    genotype = sampleResult.split(":")[gtIndex]
                else:
                    genotype = sampleResult
                
                # Edit genotype to have a consistently predictable separator
                "We don't care if the VCF is phased or not for this function"
                if genotype == ".":
                    genotype = "./."
                genotype = genotype.replace("/", "|")
                genotype = genotype.split("|")
                
                # Get the genotype character for this sample
                genoNum = sum([ 1 if x == "0" else 0 for x in genotype ])
                geno.append(str(genoNum))
            
            # Write the geno line to the output file
            fileOut.write("\t".join(geno) + "\n")
    
    print("Program finished successfully!")

if __name__ == "__main__":
    main()

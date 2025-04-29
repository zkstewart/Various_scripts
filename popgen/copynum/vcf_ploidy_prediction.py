#! python3
# vcf_ploidy_prediction.py
# Script to parse VCF files and identify the proportion of sites where
# the depth of coverage (DP) is greater than the number of alleles (AD).
# This is a proxy for identifying non-diploidy in the sample.

import os, argparse, gzip, codecs
from contextlib import contextmanager

def validate_args(args):
    # Validate input file locations
    if not os.path.isdir(args.parentDirectory):
        raise FileNotFoundError(f"I am unable to locate the parent directory '{args.parentDirectory}'")
    args.parentDirectory = os.path.abspath(args.parentDirectory)
    
    # Locate VCF files
    args.vcfFiles = []
    for root, dirs, files in os.walk(".", topdown=False):
        for name in files:
            if name.endswith(args.vcfSuffix):
                args.vcfFiles.append(os.path.join(root, name))
    if not args.vcfFiles:
        raise FileNotFoundError(f"No VCF files with suffix '{args.vcfSuffix}' found after " + 
                                f"recursive traversal of '{args.parentDirectory}'")
    
    # Validate numeric arguments
    if args.cutoff <= 0.0 or args.cutoff > 1.0:
        raise ValueError(f"Cutoff value '{args.cutoff}' is not > 0.0 and <= 1.0")
    
    # Validate output file names
    if os.path.exists(args.outputFileName):
        raise FileExistsError(f"Cannot write -o '{CALLING_SCRIPT}' since file already exists.")

def get_codec(fileName):
    try:
        f = codecs.open(fileName, encoding='utf-8', errors='strict')
        for line in f:
            break
        return "utf-8"
    except:
        try:
            f = codecs.open(fileName, encoding='utf-16', errors='strict')
            for line in f:
                break
            return "utf-16"
        except UnicodeDecodeError:
            print(f"Can't tell what codec '{fileName}' is!!")

@contextmanager
def read_gz_file(filename):
    if filename.endswith(".gz"):
        with gzip.open(filename, "rt") as f:
            yield f
    else:
        with open(filename, "r", encoding=get_codec(filename)) as f:
            yield f

def parse_vcf_alleles(vcfFile):
    '''
    Parses a VCF file to tally the number of SNPs and multiallelic sites.
    
    Parameters:
        vcfFile -- a string indicating the path to a VCF file
    Returns:
        dpDiff -- a dictionary with sample names as keys and lists of differences between DP and AD as values
    '''
    dpDiff = {}
    with read_gz_file(vcfFile) as fileIn:
        for line in fileIn:
            l = line.strip('\r\n\t "') # remove quotations to help with files opened by Excel
            
            # Skip blank lines
            if l == "":
                continue
            
            # Handle header lines
            if l.startswith("#CHROM"):
                try:
                    samples = l.split("\t")[9:]
                    sampleAlleles = {s: [0, 0, 0] for s in samples} # initialize dictionary with sample names as keys
                except:
                    raise ValueError(f"#CHROM header line is malformed; offending line is '{l}'")
            elif l.startswith("#"):
                continue
            
            # Handle variant lines
            else:
                # Split line based on expected delimiter
                try:
                    sl = l.split("\t")
                except:
                    raise ValueError(f"VCF file lacks tab-delimited formatting; offending line is '{l}'")
                
                # Validate line length
                if len(sl) < 10: # 9 fixed columns + >=1 genotype columns
                    raise ValueError(f"VCF file has too few columns; offending line is '{l}'")
                
                # Determine which field position we're extracting to get our GT value
                fieldsDescription = sl[8]
                dpIndex = fieldsDescription.split(":").index("DP")
                adIndex = fieldsDescription.split(":").index("AD")
                
                # Format a dictionary to store the dpDiff for each sample
                ongoingCount = 0 # This gives us the index for our samples header list 
                for sampleResult in sl[9:]: # This gives us the results for each sample as per fieldsDescription
                    dpDiff.setdefault(samples[ongoingCount], [])
                    
                    # Grab our DP and AD
                    dp = sampleResult.split(":")[dpIndex]
                    ad = sampleResult.split(":")[adIndex]
                    
                    # Skip if either DP or AD is missing
                    if dp == "." or ad == ".":
                        continue
                    
                    # Calculate the ratio difference between DP and AD
                    diff = sum(map(int, ad.split(","))) / int(dp)
                    
                    # Store the difference in the dictionary
                    dpDiff[samples[ongoingCount]].append(diff)
                    
                    ongoingCount += 1
    
    # Return results
    return dpDiff

####

## Main
def main():
    # User input
    usage = """%(prog)s will receive a directory to crawl for VCF files and parse them
    to identify the proportion of sites where the depth of coverage (DP) is greater than the number of alleles (AD).
    This is a proxy for identifying non-diploidy in the sample.
    """
    
    p = argparse.ArgumentParser(description=usage)
    # Reqs
    p.add_argument("-d", dest="parentDirectory",
                   required=True,
                   help="Input directory to crawl for VCF file(s)")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file name")
    # Opts
    p.add_argument("--vcfSuffix", dest="vcfSuffix",
                   required=False,
                   help="""Indicate the suffix of the VCF files to look for;
                   default == '.vcf.gz'""",
                   default=".vcf.gz")
    p.add_argument("--cutoff", dest="cutoff",
                   required=False,
                   type=float,
                   help="""Optional cutoff value for the ratio of DP to AD;
                   default == 1.0""",
                   default=1.0)
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse the VCF files for each sample
    combinedDpDiff = {}
    dupes = {}
    for vcfFile in args.vcfFiles:
        # Parse the VCF file
        dpDiff = parse_vcf_alleles(vcfFile)
        
        # Combine the results into a single dictionary
        for sample, diffs in dpDiff.items():
            # Handle duplicate sample names
            if sample in combinedDpDiff:
                dupes.setdefault(sample, 0)
                dupes[sample] += 1
                print(f"Duplicate sample name '{sample}' found in VCF files; " +
                      f"sample will be called f'{sample}_{dupes[sample]}'")
                sample = f"{sample}_{dupes[sample]}"
            # Add the sample to the combined dictionary
            combinedDpDiff[sample] = diffs
    
    # Write result stats to the output file
    with open(args.outputFileName, "w") as outFile:
        outFile.write("sample\tpct_with_dp_ad_diff\n")
        for sample, diffs in combinedDpDiff.items():
            d = [ x < args.cutoff for x in diffs ]
            p = sum(d) / len(d) * 100
            outFile.write(f"{sample}\t{p}\n")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

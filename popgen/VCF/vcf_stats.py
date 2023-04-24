#! python3
# vcf_stats.py
# A script to tabulate statistics for a VCF

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

def parse_vcf_to_statistics(vcfFile):
    statsDict = {
        "biallelic": 0,
        "multiallelic": 0,
        "snp": 0,
        "indel": 0,
        "calls": {}
    }
    
    with open_vcf_file(vcfFile) as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n").split("\t")
            
            # Handle header lines
            if line.startswith("#CHROM"):
                samples = sl[9:] # This gives us the ordered sample IDs
            if line.startswith("#"):
                continue
            
            # Extract details from line
            chrom, pos, id, ref, alt, qual, filter, \
                info, format = sl[0:9]
            format = format.split(":")
            
            # Tally variant types
            if "," in alt:
                statsDict["multiallelic"] += 1
            else:
                statsDict["biallelic"] += 1
            
            # Tally indel / SNP
            if any([ len(allele) != len(ref) for allele in alt.split(",") ]):
                statsDict["indel"] += 1
            else:
                statsDict["snp"] += 1
            
            # Store sample details in compact format
            varDict = { samples[i]: {} for i in range(len(samples)) }
            for i in range(len(samples)):
                sampleDetails = sl[9+i].split(":")
                varDict[samples[i]] = {
                    format[x]: sampleDetails[x] for x in range(len(format))
                }
            
            # Tally sample calls
            for sampleID, detailsDict in varDict.items():
                statsDict["calls"].setdefault(sampleID, {
                    "called": 0,
                    "missed": 0
                })
                
                sampleGT = list(set(detailsDict["GT"].replace("|", "/").split("/")))
                if sampleGT == ["."]:
                    statsDict["calls"][sampleID]["missed"] += 1
                else:
                    statsDict["calls"][sampleID]["called"] += 1
    
    return statsDict

## Main
def main():
    # User input
    usage = """%(prog)s reads in a VCF and produces a tabular TSV
    output summarising some relevant statistics.
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-v", dest="vcfFile", required=True,
        help="Input VCF file for filtering")
    p.add_argument("-o", dest="outputFileName", required=True,
        help="Output file name for the statistics TSV file")
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse VCF as genotypes dictionary
    statsDict = parse_vcf_to_statistics(args.vcfFile)
    
    # Produce TSV file
    with open(args.outputFileName, "w") as fileOut:
        # Write header
        fileOut.write("{0}\n".format("\t".join([
            f"# vcf_stats.py output for {args.vcfFile}"
        ])))
        
        # Write overall statistics
        fileOut.write("{0}\n".format("\t".join([
            "# SNPs that are:", "biallelic",
            "multiallelic", "snp", "indel"
        ])))
        fileOut.write("{0}\n\n".format("\t".join([
            "", 
            str(statsDict["biallelic"]), str(statsDict["multiallelic"]), 
            str(statsDict["snp"]), str(statsDict["indel"])
        ])))
        
        # Write per-sample statistics
        totalSNPs = statsDict["biallelic"] + statsDict["multiallelic"]
        
        fileOut.write("{0}\n".format("\t".join([
            "# sample ID", 
            "called num", "called percent",
            "missed num", "missed percent"
        ])))
        for sampleID, sampleStatsDict in statsDict["calls"].items():
            pctCalled = sampleStatsDict["called"] / totalSNPs
            pctMissed = sampleStatsDict["missed"] / totalSNPs
            fileOut.write("{0}\n".format("\t".join([
                sampleID, 
                str(sampleStatsDict["called"]), str(pctCalled),
                str(sampleStatsDict["missed"]), str(pctMissed)
            ])))
    
    print("Program finished successfully!")

if __name__ == "__main__":
    main()

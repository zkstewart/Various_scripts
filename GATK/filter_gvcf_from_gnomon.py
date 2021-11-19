#! python3
# filter_gvcf_from_gnomon.py
# Script to parse a combined gvcf and remove
# SNPs within introns

import os, argparse

# Define functions
def validate_args(args):
    # Validate that all arguments have been provided
    for key, value in vars(args).items():
        if value == None:
            print(key + ' argument was not specified. Fix this and try again.')
            quit()
    # Validate input file locations
    if not os.path.isfile(args.vcf):
        print('I am unable to locate the VCF file (' + args.vcf + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.gff3):
        print('I am unable to locate the GFF3 file (' + args.gff3 + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print('File already exists at output location (' + args.outputFileName + ')')
        print('Make sure you specify a unique file name and try again.')
        quit()

def cds_regions(gff3File):
    cdsRegions = {}
    with open(gff3File, "r") as fileIn:
        for line in fileIn:
            # Skip irrelevant lines
            if line.startswith("#"):
                continue
            sl = line.rstrip("\r\n ").split("\t")
            if sl[2] != "CDS":
                continue
            if "protein_id=" not in sl[8]:
                continue
            # Extract relevant details
            contig = sl[0]
            start = sl[3]
            end = sl[4] + 1 # 1-based range for checking within
            # Add relevant CDS coordinates to dict
            if contig not in cdsRegions:
                cdsRegions[contig] = []
            cdsRegions[contig].append(range(int(start), int(end)))
    return cdsRegions

def filter_gvcf_qual_introns(vcfFile, outputFileName, cdsRegions):
    with open(vcfFile, "r") as fileIn, open(outputFileName, "w") as fileOut:
        for line in fileIn:
            if line.startswith("#"):
                fileOut.write(line)
                continue
            # Get relevant details
            sl = line.split("\t")
            contig = sl[0]
            pos = int(sl[1])
            # Make filtration checks
            isIntron = True
            if contig not in cdsRegions: # Skip contigs which have no CDS regions
                continue
            for region in cdsRegions[contig]:
                if pos in region:
                    isIntron = False
                    break
            if isIntron:
                continue
            # Write good line to file
            fileOut.write(line)

def main():
    # User input
    usage = """%(prog)s reads in a VCF with multiple samples formatted by combining
    GVCF's using GATK's recommended methodology. It also reads in a GFF3 containing the gene
    that corresponds to the VCF region. It can filter variants by their QUAL field, and it will
    also remove variants that are not coding i.e., within a CDS region annotated in the GFF3.
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-v", dest="vcf",
        help="Input VCF file name")
    p.add_argument("-g", dest="gff3",
        help="Input GFF3 file name")
    p.add_argument("-o", dest="outputFileName",
        help="Output file name for the filtered VCF")
    args = p.parse_args()
    validate_args(args)

    # Simple GFF3 load for CDS regions
    cdsRegions = cds_regions(args.gff3)

    # Perform filtration
    filter_gvcf_qual_introns(args.vcf, args.outputFileName, cdsRegions)

if __name__ == "__main__":
    main()

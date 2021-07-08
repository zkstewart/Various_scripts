#! python3
# filter_gvcf.py
# Script to parse a combined gvcf and perform various
# filtrations including setting a minimum qual value
# as well was removing regions corresponding to introns

import os, argparse
import gzip
from Bio import SeqIO

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
    # Validate float inputs
    if args.qualMin < 0:
        print('Quality minimum value must be greater than 0.')
        quit()

def cds_regions_for_gene_id(gff3File, targetID):
    cdsRegions = []
    with open(gff3File, "r") as fileIn:
        for line in fileIn:
            # Skip irrelevant lines
            if line.startswith("#"):
                continue
            sl = line.rstrip("\r\n").split("\t")
            if sl[2] != "CDS":
                continue
            # Check attributes for gene ID
            attributes = sl[8].split(";")
            geneID = None
            parentID = None
            for a in attributes:
                if a.startswith("ID="):
                    geneID = a[3:]
                elif a.startswith("Parent="):
                    parentID = a[7:]
            if targetID not in geneID and targetID not in parentID:
                # Rescue 1: remove the .suffix and try again
                if ".path" or ".mrna" in targetID:
                    targetAlt = targetID.rsplit(".", maxsplit=1)[0]
                    if targetAlt not in geneID and targetAlt not in parentID:
                        continue
                else:
                    continue
            # Add relevant CDS coordinates to dict
            cdsRegions.append(range(int(sl[3]), int(sl[4])+1)) # 1-based range for checking within
    return cdsRegions

def filter_gvcf_qual_introns(vcfFile, outputFileName, cdsRegions, qualMin):
    with open(vcfFile, "r") as fileIn, open(outputFileName, "w") as fileOut:
        for line in fileIn:
            if line.startswith("#"):
                fileOut.write(line)
                continue
            # Get relevant details
            sl = line.split("\t")
            pos = int(sl[1])
            qual = float(sl[5])
            # Make filtration checks
            if qual < qualMin:
                continue
            isIntron = True
            for region in cdsRegions:
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
    p.add_argument("-id", dest="geneID",
        help="Input geneID to check for within the GFF3")
    p.add_argument("-q", dest="qualMin", type=float, default=0,
        help="Optional: Input geneID to check for within the GFF3")
    p.add_argument("-o", dest="outputFileName",
        help="Output file name for the filtered VCF")
    args = p.parse_args()
    validate_args(args)

    # Simple GFF3 load for CDS regions
    cdsRegions = cds_regions_for_gene_id(args.gff3, args.geneID)
    print(cdsRegions)

    # Perform filtration
    filter_gvcf_qual_introns(args.vcf, args.outputFileName, cdsRegions, args.qualMin)

if __name__ == "__main__":
    main()

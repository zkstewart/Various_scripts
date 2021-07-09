#! python3
# convert_gvcf_to_phase.py
# Script to parse a gvcf that has been appropriately filtered
# and convert it to a file format that PHASE will accept for
# reconstructing haplotypes

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
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print('File already exists at output location (' + args.outputFileName + ')')
        print('Make sure you specify a unique file name and try again.')
        quit()

def parse_gvcf_to_phase(vcfFile, outputFileName):
    COLS_BEFORE_GENOTYPES = 9 # As of 8/07/21, GATK produces a VCF with "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" as leadup to species
    with open(vcfFile, "r") as fileIn:
        for line in fileIn:
            # Skip unnecessary comment lines
            if line.startswith("##"):
                continue
            # Handle the first comment line
            sl = line.split("\t")
            if line.startswith("#"):
                numSpecies = len(sl) - COLS_BEFORE_GENOTYPES 
                speciesGenotypes = [ [] for i in range(numSpecies) ] # Produce container for genotype
                positions = [] # Produce container for chromosome positions
                continue
            # Handle variant lines
            positions.append(sl[1]) # Add POS to positions list
            for i in range(numSpecies):
                colIndex = i + COLS_BEFORE_GENOTYPES
                gt = sl[colIndex].split(":")[0] # gt = genotype, which is GT field in VCF; it's always the first field in GATK VCF
                if "|" in gt:
                    gt = gt.split("|")
                else:
                    gt = gt.split("/")
                if gt[0] == "." or gt[1] == ".":
                    gt = ["?", "?"]
                speciesGenotypes[i].append(gt) # Add GT to genotypes list
    # Produce output file
    with open(outputFileName, "w") as fileOut:
        # Write individuals / loci number header
        fileOut.write(str(numSpecies) + "\n")
        fileOut.write(str(len(positions)) + "\n")
        # Write positions header
        fileOut.write(" ".join(positions) + "\n")
        # Write locus type header
        ## Note: Assume all alleles are bi-allelic; if not, this is bad
        fileOut.write("S"*len(positions) + "\n")
        # Loop and write each species' data to file
        for i in range(numSpecies):
            fileOut.write("#{0}\n".format(i+1))
            gt1 = []
            gt2 = []
            for x in range(len(speciesGenotypes[i])):
                gt1.append(speciesGenotypes[i][x][0])
                gt2.append(speciesGenotypes[i][x][1])
            fileOut.write(" ".join(gt1) + "\n")
            fileOut.write(" ".join(gt2) + "\n")

def main():
    # User input
    usage = """%(prog)s reads in a VCF with multiple samples formatted by combining
    GVCF's using GATK's recommended methodology. It should have been filtered already
    by filter_gvcf.py so it's ready to be used. Using this file, an output will be generated
    that is compatible with PHASE for haplotype reconstruction.
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-v", dest="vcf",
        help="Input VCF file name")
    p.add_argument("-o", dest="outputFileName",
        help="Output file name for the filtered VCF")
    args = p.parse_args()
    validate_args(args)

    # Perform conversion
    parse_gvcf_to_phase(args.vcf, args.outputFileName)

if __name__ == "__main__":
    main()

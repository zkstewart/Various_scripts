#! python3
# convert_vcf_to_abh.py
# Script to enable conversion from VCF to ABH genotype format.
# That format is highly simplified, being A for 0/0 genotype,
# B for 1/1 genotype, and H for 0/1 genotype.

import os, argparse

# Define functions
def validate_args(args):
    # Validate input file location
    if not os.path.isfile(args.vcf):
        print('I am unable to locate the QTLseq VCF file (' + args.vcf + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Handle file overwrites
    if os.path.isfile(args.outputFileName):
        print(args.outputFileName + ' already exists. Delete/move/rename this file and run the program again.')
        quit()

def main():
    # User input
    usage = """%(prog)s receives a VCF and converts it to ABH format. Output
    will be written in CSV format.
    """
    # Required
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-v", dest="vcf",
                   required=True,
                   help="Input VCF file produced by Freebayes")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Specify the output file name")
    args = p.parse_args()
    validate_args(args)
    
    # Parse VCF and write to file
    with open(args.vcf, "r") as fileIn, open(args.outputFileName, "w") as fileOut:
        for line in fileIn:
            sl = line.rstrip("\r\n ").split("\t")
            # Skip comment lines
            if line.startswith("#") and not line.startswith("#CHROM"):
                continue
            # Handle the sample header line
            elif line.startswith("#CHROM"):
                fileOut.write(",{0}\n".format(
                    ",".join(sl[9:])
                ))
            # Handle content lines
            else:
                chrom, pos, id, ref, alt, \
                    qual, filter, info, format = sl[0:9]
                sampleFormat = format.split(":")
                
                # Skip non-biallelic sites
                if "," in alt:
                    continue
                
                # Get genotype per sample (as A/B/H encoded)
                genotypes = []
                for sampleData in sl[9:]:
                    sampleDetails = sampleData.split(":")
                    gt = sampleDetails[sampleFormat.index("GT")].replace(".", "0")
                    if gt == "0/0" or gt == "0|0":
                        genotypes.append("A")
                    elif gt == "1/1" or gt == "1|1":
                        genotypes.append("B")
                    else:
                        genotypes.append("H")
                
                # Write to file
                newLine = "{0}_{1},{2}\n".format(
                    chrom, pos,
                    ",".join(genotypes)
                )
                
                fileOut.write(newLine)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

#! python3
# find_vcf_snps_in_gff3.py
# Script to take in a VCF and GFF3 file and locate any genes
# that contain a SNP identified in the VCF.

import os, argparse

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.vcfFile):
        print('I am unable to locate the VCF file (' + args.vcfFile + ')')
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

def get_cds_regions(gff3File):
    cdsRegions = {}
    geneIDs = [] # we use this to find CDS' that map directly to genes rather than mRNAs
    with open(gff3File, "r") as fileIn:
        for line in fileIn:
            # Handle header lines
            if line.startswith("#"):
                continue
            
            # Extract relevant details
            l = line.rstrip("\r\n").split("\t")
            chrom = l[0]
            annotType = l[2]
            start = int(l[3])
            end = int(l[4])
            details = l[8].strip("\"").split(';')
            detail_dict = {}
            for i in range(len(details)):
                if details[i] == '':
                    continue
                split_details = details[i].split('=', maxsplit=1)
                detail_dict[split_details[0]] = split_details[1]
            
            # Set up storage structure at the start of each chromosome
            if chrom not in cdsRegions:
                cdsRegions[chrom] = {}
            
            # Work out the chrom -> gene -> mRNA -> CDS pairings 
            if annotType == "gene":
                cdsRegions[chrom][detail_dict["ID"]] = {}
                geneIDs.append(detail_dict["ID"])
            elif annotType == "mRNA":
                cdsRegions[chrom][detail_dict["Parent"]][detail_dict["ID"]] = [] # top-down indexing
                cdsRegions[chrom][detail_dict["ID"]] = detail_dict["Parent"] # index back up a level for CDS storing
            elif annotType == "CDS":
                # Get coordinates range
                coords = range(start, end+1)
                # Index into dict structure
                mrnaID = detail_dict["Parent"]
                if mrnaID in geneIDs: # this means its a CDS which maps direct to a gene i.e., it's not protein coding
                    continue
                geneID = cdsRegions[chrom][mrnaID]
                cdsRegions[chrom][geneID][mrnaID].append(coords) # store the CDS coordinates range in the 
            
    return cdsRegions

def get_snp_positions_from_vcf(vcfFile):
    snpPositions = {}
    with open(vcfFile, "r") as fileIn:
        for line in fileIn:
            # Handle header lines
            if line.startswith("#"):
                continue
            
            # Extract relevant details
            l = line.rstrip("\r\n").split("\t")
            chrom = l[0]
            pos = int(l[1])
            
            # Store in dictionary
            if chrom not in snpPositions:
                snpPositions[chrom] = []
            snpPositions[chrom].append(pos)
            
    return snpPositions

def get_genes_with_snp_in_cdsRegions(cdsRegions, snpPositions):
    genes = []
    for chrom, positions in snpPositions.items():
        for p in positions:
            found = False
            for geneID, mrnaDict in cdsRegions[chrom].items():
                if type(mrnaDict).__name__ == "str": # This is a back indexing value
                    continue
                for mrnaID, rangesList in mrnaDict.items():
                    for range in rangesList:
                        if p in range:
                            found = True
                            break
                if found: # This will keep geneID at the value where found was set == True
                    break
            if found:
                genes.append(geneID)
    return genes

## Main
def main():
    # User input
    usage = """%(prog)s reads in a VCF and GFF3 to locate SNPs
    that reside within protein-coding gene sequences i.e., within
    regions marked as CDS within the GFF3. It will output
    a text file containing any relevant gene IDs.
    """
    p = argparse.ArgumentParser(description=usage)
    ## Required
    p.add_argument("-v", dest="vcfFile", required=True,
        help="Input VCF file containing SNPs")
    p.add_argument("-g", dest="gff3", required=True,
        help="Input GFF3 file")
    p.add_argument("-o", dest="outputFileName", required=True,
        help="Output text file name")
    
    args = p.parse_args()
    validate_args(args)
    
    # Simple parse of GFF3 for CDS regions and their corresponding gene ID
    cdsRegions = get_cds_regions(args.gff3)
    
    # Parse VCF as genotypes dictionary
    snpPositions = get_snp_positions_from_vcf(args.vcfFile)
    
    # Compare CDS and VCF SNPs to find relevant genes
    genes = get_genes_with_snp_in_cdsRegions(cdsRegions, snpPositions)
    
    # Produce output file0
    with open(args.outputFileName, "w") as fileOut:
        fileOut.write("\n".join(genes))
    
    # Let user know everything went swimmingly
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

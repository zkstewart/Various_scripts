#! python3
# put_variants_into_context.py
# Script to take in a GFF3 file and a VCF,
# predicting the effect of variants upon
# gene CDS.

import os, argparse, sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))) # 3 dirs up is where we find GFF3IO
from Function_packages import ZS_GFF3IO, ZS_VCFIO

# Define functions
def validate_args(args):
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
        print(f'File already exists at output location ({args.outputFileName})')
        print('Make sure you specify a unique file name and try again.')
        quit()

def determine_snp_location_from_matches(pos, matches):
    '''
    Parameters:
        pos -- an integer indicating the position of a SNP in the genome
        matches -- list containing ZS_GFF3IO.Feature objects that are
                   the result of a NCLS search for a SNP position; the
                   NCLS search is assumed to have been performed over
                   the "gene" feature type.
    Returns:
        location -- a string indicating where the SNP is located
                    i.e., "CDS", "UTR", "intron", or "intergenic".
        geneID -- a string indicating the gene ID of the gene that
                  the SNP is located in; if the SNP is intergenic,
                  this value is None.
    '''
    if matches == []:
        return ["intergenic", None]
    else:
        bestMatch = ["intron", None]
        for geneFeature in matches:
            # Default condition for gene feature
            "If we don't find a CDS or exon overlap in this gene, that means it must be an intron"
            if bestMatch[0] == "intron" and bestMatch[1] == None:
                bestMatch = ["intron", geneFeature.ID]
            
            # Updating condition based on overlap position
            for childFeature in geneFeature.retrieve_all_children():
                if int(pos) >= childFeature.start and int(pos) <= childFeature.end:
                    # Update best match according to feature type
                    if childFeature.type == "CDS":
                        bestMatch = ["CDS", geneFeature.ID]
                    elif childFeature.type == "exon":
                        if bestMatch[0] != "CDS":
                            bestMatch = ["UTR", geneFeature.ID]
    return bestMatch

## Main
def main():
    # User input
    usage = """%(prog)s reads in GFF3 and a VCF file and generates a tabular file
    indicating the genomic context for each variant i.e., whether they are intergenic,
    intron, UTR, or CDS variants.
    """
    p = argparse.ArgumentParser(description=usage)
    ## Required
    p.add_argument("-v", dest="vcf",
                   required=True,
                   help="Input phased VCF")
    p.add_argument("-g", dest="gff3",
                   required=True,
                   help="Input GFF3 file")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Specify file name to write output to")
    args = p.parse_args()
    validate_args(args)
    
    # Parse GFF3 with NCLS indexing
    gff3Obj = ZS_GFF3IO.GFF3(args.gff3, strict_parse=False)
    gff3Obj.create_ncls_index(typeToIndex="gene")
    
    # Parse VCF and iteratively produce output
    with ZS_VCFIO.open_vcf_file(args.vcf) as fileIn, open(args.outputFileName, "w") as fileOut:
        # Write header line
        fileOut.write("#contig\tposition\tref\talt\tlocation\tgene_id\n")
        # Iterate through VCF
        for line in fileIn:
            # Skip header lines
            if line.startswith("#"):
                continue
            else:
                # Extract details from line
                chrom, pos, id, ref, alt, qual, filter, info, \
                    format = line.rstrip("\r\n ").split("\t")[0:9]
                pos = int(pos)
                
                # Identify SNP location in genomic context
                matches = gff3Obj.ncls_finder(pos, pos, "contig", chrom)
                snpLocation, geneID = determine_snp_location_from_matches(pos, matches)
                
                # Write output
                fileOut.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(
                    chrom, pos, ref, alt, snpLocation, geneID if geneID != None else "."
                ))
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

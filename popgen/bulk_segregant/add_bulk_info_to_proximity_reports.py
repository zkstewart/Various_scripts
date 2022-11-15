#! python3
# filter_vcf_to_match_snp_index.py
# Script to enable easy filtering of the qtlseq.vcf file so as to
# only contain the SNPs identified in the snp_index.p##.tsv file.

import os, argparse

# Define functions
def validate_args(args):
    # Validate input file location
    for suffix in [".gene_proximity.tsv", ".snp_proximity.tsv"]:
        if not os.path.isfile(args.proximityReportPrefix + suffix):
            print(f'I am unable to locate one of the proximity report files ({os.path.isfile(args.proximityReportPrefix + suffix)})')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    if not os.path.isfile(args.vcf):
        print('I am unable to locate the QTLseq VCF file (' + args.vcf + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Handle file overwrites
    for suffix in [".gene_proximity.tsv", ".snp_proximity.tsv"]:
        if os.path.isfile(args.outputFilePrefix + suffix):
            print(f'"{args.outputFilePrefix + suffix}" already exists. Delete/move/rename this file and run the program again.')
            quit()

def parse_snp_proximity_file(snpProxFile):
    '''
    Parameters:
        snpProxFile -- a string indicating the location of the .snp_proximity.tsv file
    Returns:
        snpProxDict -- a dictionary with structure like:
                       {
                           'contigID1':
                           {
                               'pos1':
                               {
                                   'at location': 'UTR',
                                   'left of': '.',
                                   ...
                               },
                               'pos2':
                               {
                                   ...
                               },
                               ...
                           },
                           ...
                       }
    '''
    snpProxDict = {}
    with open(snpProxFile, "r") as fileIn:
        for line in fileIn:
            if line.startswith("#"):
                header = line[1:].rstrip("\r\n ").split("\t")
                continue
            sl = line.rstrip("\r\n ").split("\t")
            
            # Index values
            geneID, pos = sl[0], sl[1]
            snpProxDict.setdefault(geneID, {})
            thisSnpDict = {}
            for i in range(2, len(sl)):
                thisSnpDict[header[i]] = sl[i]
            snpProxDict[geneID][pos] = thisSnpDict
    return snpProxDict

def parse_gene_proximity_file(geneProxFile):
    '''
    Parameters:
        geneProxFile -- a string indicating the location of the .gene_proximity.tsv file
    Returns:
        geneProxDict -- a dictionary with structure like:
                       {
                           ##TBD
                       }
    '''
    geneProxDict = {}
    with open(geneProxFile, "r") as fileIn:
        for line in fileIn:
            if line.startswith("#"):
                header = line[1:].rstrip("\r\n ").split("\t")
                continue
            sl = line.rstrip("\r\n ").split("\t")
            
            # Index values
            geneID = sl[0]
            thisSnpDict = {}
            for i in range(1, len(sl)):
                thisSnpDict[header[i]] = sl[i]
            geneProxDict[geneID] = thisSnpDict
    return geneProxDict

def get_genotypes_from_qtlseq_vcf(vcfFile, snpProxDict):
    '''
    Reads in a VCF file and obtains information on SNP genotypes
    for SNPs present in snpProxDict.
    
    It's assumed the VCF is formatted like QTL-seq would. So that means
    there's 3 samples; the first is the parent, followed by bulk1,
    then bulk2.
    '''
    genotypesDict = {}
    with open(vcfFile, "r") as fileIn:
        firstLine = True
        for line in fileIn:
            l = line.rstrip("\r\n").split("\t")
            
            # Handle header lines
            if line.startswith("#"):
                firstLine = False
                continue
            elif firstLine is True:
                firstLine = False
                if not l[1].isdigit():
                    continue
            
            # Extract details
            chrom = l[0]
            pos = l[1]
            
            # Skip if this SNP is not relevant to us
            if not chrom in snpProxDict and not pos in snpProxDict[chrom]:
                continue
            
            # Determine which field position we're extracting to get our GT value
            fieldsDescription = l[8]
            if ":" not in fieldsDescription:
                gtIndex = -1
            else:
                gtIndex = fieldsDescription.split(":").index("GT")
            
            # Get genotypes for parent and bulks
            parentGT = l[9].split(":")[gtIndex]
            bulk1GT = l[10].split(":")[gtIndex]
            bulk2GT = l[11].split(":")[gtIndex]
            
            # Store in dictiomnary
            genotypesDict.setdefault(chrom, {})
            genotypesDict[chrom][pos] = {
                "parent": parentGT,
                "bulk1": bulk1GT,
                "bulk2": bulk2GT
            }
    
    return genotypesDict

def get_tally_for_genes(snpProxDict, genotypesDict):
    '''
    Using the snpProxDict and genotypesDict, we can get the genotype
    for any SNPs contained within coding regions.
    '''
    geneProxTally = {}
    for contigID, posDict in snpProxDict.items():
        for pos, rowDict in posDict.items():
            # Skip SNPs that are irrelevant to us here
            geneID = rowDict["inside"]
            if geneID == "." or rowDict["at location"] != "CDS":
                continue
            
            # Get the genotype for this SNP and tally it
            thisGtDict = genotypesDict[contigID][pos]
            parent, bulk1, bulk2 = thisGtDict["parent"], thisGtDict["bulk1"], thisGtDict["bulk2"]
            geneProxTally.setdefault(geneID, {
                "bulk1_ratio": [0, 0], # [same, total]
                "bulk2_ratio": [0, 0]
            })
            
            geneProxTally[geneID]["bulk1_ratio"][1] += 1
            if bulk1 == parent:
                geneProxTally[geneID]["bulk1_ratio"][0] += 1
            
            geneProxTally[geneID]["bulk2_ratio"][1] += 1
            if bulk2 == parent:
                geneProxTally[geneID]["bulk2_ratio"][0] += 1
    return geneProxTally

def main():
    # User input
    usage = """%(prog)s receives the proximity report files generated by
    snp_proximity_report.py and adds bulk and parent information into the file.
    Specifically, it will note the genotype of the parent and two bulks for the
    snp_proximity.tsv file, and for the gene_proximity.tsv it will tally the
    proportion of SNPs in the CDS that are shared with the parent.
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-p", dest="proximityReportPrefix",
                   required=True,
                   help="""Specify the prefix to the two proximity reports e.g.,
                   '/path/to/qtlseq_comparison' which should then match to two files
                   with '.gene_proximity.tsv' and '.snp_proximity.tsv' suffixes""")
    p.add_argument("-v", dest="vcf",
                   required=True,
                   help="Input VCF file produced by QTL-seq")
    p.add_argument("-o", dest="outputFilePrefix",
                   required=True,
                   help="Specify the output file prefix (will create two new prox report files)")
    args = p.parse_args()
    validate_args(args)
    
    # Get our proximity report file names
    snpProxFile = os.path.join(args.proximityReportPrefix + ".snp_proximity.tsv")
    geneProxFile = os.path.join(args.proximityReportPrefix + ".gene_proximity.tsv")
    
    # Parse the proximity files
    snpProxDict = parse_snp_proximity_file(snpProxFile)
    # geneProxDict = parse_gene_proximity_file(geneProxFile) ## not actually needed
    
    # Parse the VCF to obtain genotypes for relevant SNPs
    genotypesDict = get_genotypes_from_qtlseq_vcf(args.vcf, snpProxDict)
    
    # Tally genotypes for each gene
    geneProxTally = get_tally_for_genes(snpProxDict, genotypesDict)
    
    # Write new SNP proximity file
    with open(snpProxFile, "r") as fileIn, open(args.outputFilePrefix + ".snp_proximity.tsv", "w") as fileOut:
        for line in fileIn:
            l = line.rstrip("\r\n ")
            
            # Update header
            if l.startswith("#"):
                l += "\tparent_GT\tbulk1_GT\tbulk2_GT"
                fileOut.write(f"{l}\n")
            # Update content lines
            else:
                contig, pos = l.split("\t")[0:2]
                thisGtDict = genotypesDict[contig][pos]
                l += "\t_{0}_\t_{1}_\t_{2}_".format(
                    thisGtDict["parent"], thisGtDict["bulk1"], thisGtDict["bulk2"]
                )
                fileOut.write(f"{l}\n")
    
    # Write new gene proximity file
    with open(geneProxFile, "r") as fileIn, open(args.outputFilePrefix + ".gene_proximity.tsv", "w") as fileOut:
        for line in fileIn:
            l = line.rstrip("\r\n ")
            
            # Update header
            if l.startswith("#"):
                l += "\tbulk1 different from parent (%)\tbulk2 different from parent (%)"
                fileOut.write(f"{l}\n")
            # Update content lines
            else:
                geneID = l.split("\t")[0]
                if geneID not in geneProxTally:
                    fileOut.write(f"{l}\t.\t.\n")
                else:
                    thisTally = geneProxTally[geneID]
                    bulk1Ratio = round(thisTally["bulk1_ratio"][0]/thisTally["bulk1_ratio"][1], 2)
                    bulk2Ratio = round(thisTally["bulk2_ratio"][0]/thisTally["bulk2_ratio"][1], 2)
                    fileOut.write(f"{l}\t{bulk1Ratio}\t{bulk2Ratio}\n")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

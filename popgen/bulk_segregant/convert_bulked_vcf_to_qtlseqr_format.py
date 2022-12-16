#! python3
# convert_bulked_vcf_to_qtlseqr_format.py
# Script to produce a re-formatted output amenable to loading in
# to the QTLseqr package.

import os, argparse

# Define functions
def validate_args(args):
    # Validate input file location
    if not os.path.isfile(args.vcf):
        print('I am unable to locate the VCF file (' + args.vcf + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Handle file overwrites
    if os.path.isfile(args.outputFileName):
        print(args.outputFileName + ' already exists. Delete/move/rename this file and run the program again.')
        quit()

def parse_vcf(vcfFile):
    '''
    This function will parse a VCF file and return a dictionary which can be used
    for various purposes
    
    Parameters:
        vcfFile -- a string indicating the file location of the VCF file
    Returns:
        vcfDict -- a dictionary with structure like:
                      {
                          'chrom1': {
                              pos1: {"REF": ___, "ALT": ___, "GT": ___, "AD: ___, ...},
                              pos2: { ... },
                              ...
                          },
                          'chrom2': ...,
                          ...
                      }
    '''
    vcfDict = {}
    with open(vcfFile, "r") as fileIn:
        firstLine = True
        for line in fileIn:
            # Skip comment lines
            if line.startswith("#") and not line.startswith("#CHROM"):
                continue
            
            l = line.rstrip("\r\n ")
            sl = l.split("\t")
            
            # Handle the header line
            if firstLine == True:
                assert line.startswith("#CHROM") or line.startswith("CHROM"), \
                    "First non-comment line isn't the expected VCF header line format"
                
                sampleIDs = sl[9:]
                firstLine = False
            # Handle content lines
            else:
                chrom, pos, id, ref, alt, \
                    qual, filter, info, format = sl[0:9]
                vcfDict.setdefault(chrom, {})
                vcfDict[chrom][pos] = {
                    "REF": ref,
                    "ALT": alt,
                }
                
                # Parse DP from info in case imputing is required below
                if "DP=" in info:
                    infoDP = [x.split("=")[1] for x in info.split(";") if x.startswith("DP=")][0]
                else:
                    infoDP = None
                
                # Parse sample details according to FORMAT 
                for i in range(len(sampleIDs)):
                    sampleFormat = format.split(":")
                    sampleDetails = sl[9+i].split(":")
                    
                    sampleDict = {}
                    for x in range(len(sampleFormat)):
                        sampleDict[sampleFormat[x]] = sampleDetails[x]
                    
                    # Store result
                    vcfDict[chrom][pos][sampleIDs[i]] = sampleDict # sampleIDs[i] == sampleName
                
                # Roughly impute DP per sample if not already existing (and if it's even possible)
                if "DP" not in sampleDict and infoDP != None and "AD" in sampleDict: # most recently used sampleDict is fine
                    sampleADs = [vcfDict[chrom][pos][id]["AD"] for id in sampleIDs]
                    adsRatio = [sum(map(int, ad.split(","))) / int(infoDP) for ad in sampleADs]
                    dpsImpute = [int(adR*int(infoDP)) for adR in adsRatio]
                    for i in range(len(sampleIDs)):
                        vcfDict[chrom][pos][sampleIDs[i]]["DP"] = str(dpsImpute[i])
    return vcfDict

def main():
    # User input
    usage = """%(prog)s receives a VCF produced by e.g., QTL-seq which contains at least
    2 sample columns for the bulked progeny. It will reformat this so that QTLseqr can load
    it in.
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="vcf",
                   required=True,
                   help="Input bulked VCF file")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Specify the output file name")
    p.add_argument("-b1", dest="bulk1Column",
                   required=True,
                   help="Specify the column name for bulk 1")
    p.add_argument("-b2", dest="bulk2Column",
                   required=True,
                   help="Specify the column name for bulk 2")
    args = p.parse_args()
    validate_args(args)
    
    # Parse VCF file
    vcfDict = parse_vcf(args.vcf)
    
    # Validate that bulk columns exist, and check which fields we are writing
    FIELDS_TO_WRITE = ["GT", "AD", "DP", "GQ", "PL"]
    for chromDict in vcfDict.values():
        for posDict in chromDict.values():
            assert args.bulk1Column in posDict and args.bulk2Column in posDict, \
                "Provided bulk column IDs not found in {0}".format(posDict.keys())
            fieldsThatExist = [x for x in FIELDS_TO_WRITE if x in posDict[args.bulk1Column]]
            assert "AD" in fieldsThatExist, \
                "AD value not found in VCF; it must exist for QTLseqr to work!"
            break
    
    # Write new output
    with open(args.outputFileName, "w") as fileOut:
        fileOut.write("CHROM\tPOS\tREF\tALT\t{0}\t{1}\n".format(
            "\t".join([f"{args.bulk1Column}.{field}" for field in fieldsThatExist]),
            "\t".join([f"{args.bulk2Column}.{field}" for field in fieldsThatExist]),
        ))
        for chrom, chromDict in vcfDict.items():
            for pos, posDict in chromDict.items():
                # Update GT
                ref_alts = [posDict["REF"], *posDict["ALT"].split(",")]
                posDict[args.bulk1Column]["GT"] = posDict[args.bulk1Column]["GT"].replace(".", "0")
                posDict[args.bulk2Column]["GT"] = posDict[args.bulk2Column]["GT"].replace(".", "0")
                
                posDict[args.bulk1Column]["GT"] = "/".join([
                    ref_alts[int(gtNum)] for gtNum in posDict[args.bulk1Column]["GT"].split("/")
                ])
                posDict[args.bulk2Column]["GT"] = "/".join([
                    ref_alts[int(gtNum)] for gtNum in posDict[args.bulk2Column]["GT"].split("/")
                ])
                
                # Write line
                line = "{chrom}\t{pos}\t{ref}\t{alt}\t{b1Values}\t{b2Values}\n".format(
                    chrom=chrom,
                    pos=pos,
                    ref=posDict["REF"],
                    alt=posDict["ALT"],
                    b1Values="\t".join([posDict[args.bulk1Column][field] for field in fieldsThatExist]),
                    b2Values="\t".join([posDict[args.bulk2Column][field] for field in fieldsThatExist])
                )
                fileOut.write(line)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

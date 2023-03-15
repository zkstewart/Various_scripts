#! python3
# bulked_VariantsToTable.py
# Script which emulates GATK's VariantsToTable function, but with
# imputation of missing fields so as to be more appropriate for the
# QTLseqr package, or for PyBSASeq if that's actually working.

import os, argparse

# Define functions
def validate_args(args):
    # Validate input file location
    if not os.path.isfile(args.vcf):
        print(f'I am unable to locate the VCF file ({args.vcf})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Handle file overwrites
    if os.path.isfile(args.outputFileName):
        print(f'"{args.outputFileName}" already exists. Delete/move/rename this file and run the program again.')
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
                # Extract relevant details from line
                chrom, pos, id, ref, alt, \
                    qual, filter, info, format = sl[0:9]
                sampleFormat = format.split(":")
                sampleData = sl[9:]
                
                # Set up data storage structure
                vcfDict.setdefault(chrom, {})
                vcfDict[chrom][pos] = {
                    "REF": ref,
                    "ALT": alt,
                    "QUAL": qual
                }
                
                # Parse sample details according to FORMAT 
                for i in range(len(sampleIDs)):
                    sampleDetails = sampleData[i].split(":")
                    
                    # Set up sample dict for this VCF row
                    sampleDict = {}
                    for x in range(len(sampleFormat)):
                        sampleDict[sampleFormat[x]] = sampleDetails[x]
                    
                    # Index it in the overall dictionary
                    vcfDict[chrom][pos][sampleIDs[i]] = sampleDict # sampleIDs[i] == sampleName
                
                # Impute DP if required
                impute_DP(vcfDict[chrom][pos], info, "DP" in sampleFormat, "AD" in sampleFormat)
                
                # Impute PL if required
                impute_PL(vcfDict[chrom][pos], "PL" in sampleFormat, "GL" in sampleFormat)
                
                # Impute GQ if required
                impute_GQ(vcfDict[chrom][pos], "GQ" in sampleFormat, "PL" in sampleDict) # most recent sampleDict may have been imputed
                
    return vcfDict

def _derive_sample_keys(dictKeys):
    '''
    Helper function for getting keys out of the VCF dictionary.
    
    Parameters:
        dictKeys -- a dict_keys() object which should have REF, ALT keys
                    removed
    Returns:
        sampleIDs -- a list containing just the keys representing sample IDs
    '''
    return [ key for key in dictKeys if not key in ["REF", "ALT", "QUAL"]]

def impute_DP(vcfDictValue, info, dpFieldExists, adFieldExists):
    '''
    Parameters:
        vcfDictValue -- a dictionary indexed under the vcfDict of parse_vcf(),
                        which should correspond to vcfDict[chrom][pos] where
                        chrom == a string of the chromosome/contig the variant
                        is located at, and pos == a string of a digit indicating
                        where the variant resides.
        info -- a string corresponding to the 8th column of a VCF file.
                It should have a structure akin to "DP=285;VDB=0;MQ=60;...".
        dpFieldExists -- a boolean indicating whether we can expect to find
                         the DP field in this VCF entry.
        adFieldExists -- as above, but for the AD field.
    '''
    # Exit out of function if there's no need to impute DP; it already exists!
    if dpFieldExists:
        return
    
    # Exit out of function if AD field doesn't exist; impossible to impute!
    if not adFieldExists:
        return
    
    # Roughly impute DP per sample if it's possible
    if "DP=" in info:
        infoDP = [x.split("=")[1] for x in info.split(";") if x.startswith("DP=")][0]
        sampleIDs = _derive_sample_keys(vcfDictValue.keys())
        
        sampleADs = [vcfDictValue[id]["AD"] for id in sampleIDs]
        adsRatio = [sum(map(int, ad.split(","))) / int(infoDP) for ad in sampleADs]
        dpsImpute = [int(adR*int(infoDP)) for adR in adsRatio]
        
        for i in range(len(sampleIDs)):
            vcfDictValue[sampleIDs[i]]["DP"] = str(dpsImpute[i])
    
    # If DP is not provided over the entry VCF row, there's nothing we can do
    else:
        return

def impute_PL(vcfDictValue, plFieldExists, glFieldExists):
    '''
    Parameters:
        vcfDictValue -- a dictionary indexed under the vcfDict of parse_vcf(),
                        which should correspond to vcfDict[chrom][pos] where
                        chrom == a string of the chromosome/contig the variant
                        is located at, and pos == a string of a digit indicating
                        where the variant resides.
        plFieldExists -- a boolean indicating whether we can expect to find
                         the PL field in this VCF entry.
        glFieldExists -- as above, but for the GL field.
    '''
    # Exit out of function if there's no need to impute PL; it already exists!
    if plFieldExists:
        return
    
    # Exit out of function if GL field doesn't exist; impossible to impute!
    if not glFieldExists:
        return
    
    # Impute PL for each sample
    for sampleID in _derive_sample_keys(vcfDictValue.keys()):
        # First, grab the GL values as floats
        glValues = list(map(float, vcfDictValue[sampleID]["GL"].split(",")))
        
        # Transform it to PL values
        plValues = [str(int(abs(-gl*10))) for gl in glValues] # must use int!
        
        # Store the PL values
        vcfDictValue[sampleID]["PL"] = plValues

def impute_GQ(vcfDictValue, gqFieldExists, plFieldExists):
    '''
    Parameters:
        vcfDictValue -- a dictionary indexed under the vcfDict of parse_vcf(),
                        which should correspond to vcfDict[chrom][pos] where
                        chrom == a string of the chromosome/contig the variant
                        is located at, and pos == a string of a digit indicating
                        where the variant resides.
        gqFieldExists -- a boolean indicating whether we can expect to find
                         the GQ field in this VCF entry.
        plFieldExists -- as above, but for the PL field.
    '''
    # Exit out of function if there's no need to impute GQ; it already exists!
    if gqFieldExists:
        return
    
    # Exit out of function if PL field doesn't exist; impossible to impute!
    if not plFieldExists:
        return
    
    # Impute GQ for each sample
    for sampleID in _derive_sample_keys(vcfDictValue.keys()):
        # First, grab the PL values as integers
        plValues = list(map(int, vcfDictValue[sampleID]["PL"].split(",")))
        
        # Sort it from lowest to highest
        plValues.sort()
        
        # Store the second lowest value; that's the GQ
        vcfDictValue[sampleID]["GQ"] = str(plValues[1])

def main():
    # User input
    usage = """%(prog)s receives a VCF produced by e.g., QTL-seq which contains at least
    2 sample columns for the bulked progeny. It will reformat this so that QTLseqr can load
    it in.
    """
    ## Required
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
    ## Optional
    p.add_argument("--defaultHeader", dest="defaultHeader",
                   required=False,
                   action="store_true",
                   help="""Optionally indicate if you'd prefer the header to be
                   written out as bulk1.GT, bulk2.GT, etc, regardless of whatever
                   your bulk columns are""",
                   default=False)
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
        # Write header
        if args.defaultHeader:
            fileOut.write("CHROM\tPOS\tREF\tALT\tQUAL\t{0}\t{1}\n".format(
                "\t".join([f"bulk1.{field}" for field in fieldsThatExist]),
                "\t".join([f"bulk2.{field}" for field in fieldsThatExist]),
            ))
        else:
            fileOut.write("CHROM\tPOS\tREF\tALT\tQUAL\t{0}\t{1}\n".format(
                "\t".join([f"{args.bulk1Column}.{field}" for field in fieldsThatExist]),
                "\t".join([f"{args.bulk2Column}.{field}" for field in fieldsThatExist]),
            ))
        
        # Write content lines
        for chrom, chromDict in vcfDict.items():
            for pos, posDict in chromDict.items():
                ref_alts = [posDict["REF"], *posDict["ALT"].split(",")]
                
                # Handle GT depending on whether the value exists or not
                posDict[args.bulk1Column]["GT"] = "/".join(
                    [
                        ref_alts[int(gtNum)] if gtNum != "." else "."
                            for gtNum in posDict[args.bulk1Column]["GT"].split("/")
                    ]
                )
                posDict[args.bulk2Column]["GT"] = "/".join(
                    [
                        ref_alts[int(gtNum)] if gtNum != "." else "."
                            for gtNum in posDict[args.bulk2Column]["GT"].split("/")
                    ]
                )
                
                # Write line
                line = "{chrom}\t{pos}\t{ref}\t{alt}\t{qual}\t{b1Values}\t{b2Values}\n".format(
                    chrom=chrom,
                    pos=pos,
                    ref=posDict["REF"],
                    alt=posDict["ALT"],
                    qual=posDict["QUAL"],
                    b1Values="\t".join([posDict[args.bulk1Column][field] for field in fieldsThatExist]),
                    b2Values="\t".join([posDict[args.bulk2Column][field] for field in fieldsThatExist])
                )
                fileOut.write(line)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

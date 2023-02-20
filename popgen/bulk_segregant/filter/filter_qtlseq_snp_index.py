#! python3
# filter_qtlseq_snp_index.py
# Script to enable easy filtering of the snp_index.p##.tsv
# file produced by QTL-seq

import os, argparse
from collections import OrderedDict

# Define functions
def validate_args(args):
    # Validate input file location
    if not os.path.isfile(args.snpIndexFile):
        print('I am unable to locate the SNP index p## file (' + args.snpIndexFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate numeric inputs
    if args.bulk1_depthFilter != None:
        if 0 > args.bulk1_depthFilter:
            print("bulk1_depthFilter argument only makes sense to be 0 or greater")
            quit()
    if args.bulk2_depthFilter != None:
        if 0 > args.bulk2_depthFilter:
            print("bulk2_depthFilter argument only makes sense to be 0 or greater")
            quit()
    if args.bulk1_SNPindexFilter != None:
        if 0 > args.bulk1_SNPindexFilter:
            print("bulk1_SNPindexFilter argument only makes sense to be 0 or greater")
            quit()
        if 1 < args.bulk1_SNPindexFilter:
            print("bulk1_SNPindexFilter argument only makes sense to be 1 or less")
            quit()
    if args.bulk2_SNPindexFilter != None:
        if 0 > args.bulk2_SNPindexFilter:
            print("bulk2_SNPindexFilter argument only makes sense to be 0 or greater")
            quit()
        if 1 < args.bulk2_SNPindexFilter:
            print("bulk2_SNPindexFilter argument only makes sense to be 1 or less")
            quit()
    if args.delta_SNPindexFilter != None:
        if 0 > args.delta_SNPindexFilter:
            print("delta_SNPindexFilter argument only makes sense to be 0 or greater")
            quit()
        if 1 < args.delta_SNPindexFilter:
            print("delta_SNPindexFilter argument only makes sense to be 1 or less")
            quit()
    # Handle file overwrites
    if os.path.isfile(args.outputFileName):
        print(args.outputFileName + ' already exists. Delete/move/rename this file and run the program again.')
        quit()

def parse_qtlseq_snp_index_file(snpIndexFile):
    '''
    This function will parse a a snp_index.p##.tsv file produced by QTL-seq
    into an OrderedDict containing all relevant values.
    
    Parameters:
        snpIndexFile -- a string indicating the file location of the QTL-seq
                        snp_index.p##.tsv file.
    Returns:
        snpIndexDict -- an OrderedDict with structure like:
                      {
                          'contig1': OrderedDict{
                              pos1: {"variant": ___, "bulk1_depth": ___, ...},
                              pos2: { ... },
                              ...
                          },
                          'contig2': ...,
                          ...
                      }
    '''
    snpIndexDict = OrderedDict()
    with open(snpIndexFile, "r") as fileIn:
        firstLine = True
        for line in fileIn:
            # Handle header line
            if firstLine is True:
                firstLine = False
                continue
            # Handle content lines
            else:
                chrom, pos, variantType, bulk1_depth, bulk2_depth, \
                    p99, p95, bulk1_SNPindex, bulk2_SNPindex, \
                    delta_SNPindex = line.rstrip("\r\n ").split("\t")
                
                # Store in dictionary
                snpIndexDict.setdefault(chrom, OrderedDict())
                snpIndexDict[chrom][pos] = {
                    "variant": variantType,
                    "bulk1_depth": int(bulk1_depth),
                    "bulk2_depth": int(bulk2_depth),
                    "p99": float(p99),
                    "p95": float(p95),
                    "bulk1_SNPindex": float(bulk1_SNPindex),
                    "bulk2_SNPindex": float(bulk2_SNPindex),
                    "delta_SNPindex": float(delta_SNPindex)
                }
    return snpIndexDict

def filter_snp_index(snpIndexDict, variantFilter,
                     bulk1_depthFilter, bulk2_depthFilter,
                     bulk1_SNPindexFilter, bulk2_SNPindexFilter,
                     delta_SNPindexFilter):
    '''
    This function will filter a snpIndexDict as created by parse_qtlseq_snp_index_file()
    and return a new snpIndexDict object (note: no changes made to original) with only
    variants that pass all provided filters. The sub-dictionaries are merely copied,
    so the objects will be linked to each other.
    
    Parameters:
        snpIndexDict -- an OrderedDict with structure like:
                      {
                          'contig1': OrderedDict{
                              pos1: {"variant": ___, "bulk1_depth": ___, ...},
                              pos2: { ... },
                              ...
                          },
                          'contig2': ...,
                          ...
                      }
        ... [I don't want to comment everything else. You can figure it out based
        on the main() argparse values]
    '''
    newSnpIndexDict = OrderedDict()
    for contigID, contigSnpsDict in snpIndexDict.items():
        for pos, variantDict in contigSnpsDict.items():
            filtersPassed = True
            
            if variantFilter != None:
                if variantDict["variant"] != variantFilter:
                    filtersPassed = False
            
            if bulk1_depthFilter != None:
                if variantDict["bulk1_depth"] < bulk1_depthFilter:
                    filtersPassed = False
            
            if bulk2_depthFilter != None:
                if variantDict["bulk2_depth"] < bulk2_depthFilter:
                    filtersPassed = False
            
            if bulk1_SNPindexFilter != None:
                lowerFilter = 0 + bulk1_SNPindexFilter
                upperFilter = 1 - bulk1_SNPindexFilter
                
                if not (variantDict["bulk1_SNPindex"] <= lowerFilter or variantDict["bulk1_SNPindex"] >= upperFilter):
                    filtersPassed = False
            
            if bulk2_SNPindexFilter != None:
                lowerFilter = 0 + bulk2_SNPindexFilter
                upperFilter = 1 - bulk2_SNPindexFilter
                
                if not (variantDict["bulk2_SNPindex"] <= lowerFilter or variantDict["bulk2_SNPindex"] >= upperFilter):
                    filtersPassed = False
            
            if delta_SNPindexFilter != None:
                if abs(variantDict["delta_SNPindex"]) < delta_SNPindexFilter:
                    filtersPassed = False
            
            if filtersPassed:
                newSnpIndexDict.setdefault(contigID, OrderedDict())
                newSnpIndexDict[contigID][pos] = variantDict
    return newSnpIndexDict

def main():
    # User input
    usage = """%(prog)s filters 
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="snpIndexFile",
                   required=True,
                   help="Input SNP index p## file produced by QTL-seq")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Specify the output file name")
    # Opts
    p.add_argument("--variant", dest="variantFilter",
                   required=False,
                   help="Optionally, filter to retain the indicated variant type",
                   choices=["snp", "indel"])
    p.add_argument("--bulk1_depth", dest="bulk1_depthFilter",
                   type=int,
                   required=False,
                   help="Optionally, filter to retain variants >= this bulk 1 depth")
    p.add_argument("--bulk2_depth", dest="bulk2_depthFilter",
                   type=int,
                   required=False,
                   help="Optionally, filter to retain variants >= this bulk 2 depth")
    p.add_argument("--bulk1_SNPindex", dest="bulk1_SNPindexFilter",
                   type=float,
                   required=False,
                   help="Optionally, filter to retain variants <= (0 + this) or >= (1 - this)")
    p.add_argument("--bulk2_SNPindex", dest="bulk2_SNPindexFilter",
                   type=float,
                   required=False,
                   help="Optionally, filter to retain variants <= (0 + this) or >= (1 - this)")
    p.add_argument("--delta_SNPindex", dest="delta_SNPindexFilter",
                   type=float,
                   required=False,
                   help="Optionally, filter to retain variants with >= absolute delta value")
    args = p.parse_args()
    validate_args(args)
    
    # Parse SNP index file
    snpIndexDict = parse_qtlseq_snp_index_file(args.snpIndexFile)
    
    # Filter SNP index dictionary based on indicated parameters
    snpIndexDict = filter_snp_index(snpIndexDict,
                                    args.variantFilter,
                                    args.bulk1_depthFilter, args.bulk2_depthFilter,
                                    args.bulk1_SNPindexFilter, args.bulk2_SNPindexFilter,
                                    args.delta_SNPindexFilter)
    
    # Write new output file
    with open(args.outputFileName, "w") as fileOut:
        fileOut.write("CHROM\tPOSI\tvariant\tbulk1_depth\tbulk2_depth\tp99\tp95\tbulk1_SNPindex\tbulk2_SNPindex\tdelta_SNPindex\n")
        for contigID, contigSnpsDict in snpIndexDict.items():
            for pos, variantDict in contigSnpsDict.items():
                fileOut.write(
                    "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n".format(
                        contigID, pos, variantDict["variant"], variantDict["bulk1_depth"],
                        variantDict["bulk2_depth"], variantDict["p99"], variantDict["p95"],
                        variantDict["bulk1_SNPindex"], variantDict["bulk2_SNPindex"],
                        variantDict["delta_SNPindex"]
                    )
                )
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

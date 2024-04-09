#! python3
# filter_persample_diffratio.py
# Script to enable filtering of the persample difference ratio
# table file produced by calculate_persample_diffratio.py

import os, argparse

# Define functions
def validate_args(args):
    # Validate input file location
    if not os.path.isfile(args.diffratioFile):
        print('I am unable to locate the difference ratio TSV file (' + args.diffratioFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate numeric inputs
    foundAFilter = False
    if args.bulk1_alleles != None:
        foundAFilter = True
        if 0 > args.bulk1_alleles:
            print("bulk1_alleles argument must be 0 (no filtering) or greater")
            quit()
    if args.bulk2_alleles != None:
        foundAFilter = True
        if 0 > args.bulk2_alleles:
            print("bulk2_alleles argument must be 0 (no filtering) or greater")
            quit()
    ##
    if args.bulk1_refIndex != None:
        foundAFilter = True
        if 0 > args.bulk1_refIndex or 1 < args.bulk1_refIndex:
            print("bulk1_refIndex must be in the range 0 (no filtering) to 1 " + 
                  "(retain only variants completely different to the reference)")
            quit()
    if args.bulk2_refIndex != None:
        foundAFilter = True
        if 0 > args.bulk2_refIndex or 1 < args.bulk2_refIndex:
            print("bulk2_refIndex must be in the range 0 (no filtering) to 1 " + 
                  "(retain only variants completely different to the reference)")
            quit()
    ##
    if args.delta_refIndex != None:
        foundAFilter = True
        if 0 > args.delta_refIndex or 1 < args.delta_refIndex:
            print("delta_refIndex must be in the range 0 (no filtering) to 1 " + 
                  "(retain only variants where one population is completely different to " + 
                  "the reference, and the other is identical to the reference)")
            quit()
    ##
    if args.differenceRatio != None:
        foundAFilter = True
        if 0 > args.differenceRatio or 1 < args.differenceRatio:
            print("differenceRatio must be in the range 0 (no filtering) to 1 " + 
                  "(retain only variants that completely segregate between the two populations)")
            quit()
    # If no filters were provided, print an error and quit
    if not foundAFilter:
        print("You must provide at least one filter to apply to the difference ratio file!")
        quit()
    # Handle file overwrites
    if os.path.isfile(args.outputFileName):
        print(args.outputFileName + ' already exists. Delete/move/rename this file and run the program again.')
        quit()

def parse_persample_diffratio_file(diffratioFile):
    '''
    This function will parse a TSV file produced by calculate_per_sample_diffratio.py
    into a dictionary containing all relevant values.
    
    Parameters:
        diffratioFile -- a string indicating the file location to be parsed.
    Returns:
        diffratioDict -- a dict with structure like:
                      {
                          'contig1': {
                              pos1: {"variant": ___, "bulk1_depth": ___, ...},
                              pos2: { ... },
                              ...
                          },
                          'contig2': ...,
                          ...
                      }
    '''
    snpIndexDict = {}
    with open(diffratioFile, "r") as fileIn:
        firstLine = True
        for line in fileIn:
            # Handle header line
            if firstLine is True:
                firstLine = False
                continue
            # Handle content lines
            else:
                chrom, pos, variantType, bulk1_alleles, \
                    bulk2_alleles, bulk1_refIndex, \
                    bulk2_refIndex, delta_refIndex, \
                    differenceRatio = line.rstrip("\r\n ").split("\t")
                
                # Store in dictionary
                snpIndexDict.setdefault(chrom, {})
                snpIndexDict[chrom][pos] = {
                    "variant": variantType,
                    "bulk1_alleles": int(bulk1_alleles),
                    "bulk2_alleles": int(bulk2_alleles),
                    "bulk1_refIndex": float(bulk1_refIndex) if bulk1_refIndex != "." else ".",
                    "bulk2_refIndex": float(bulk2_refIndex) if bulk2_refIndex != "." else ".",
                    "delta_refIndex": float(delta_refIndex) if delta_refIndex != "." else ".",
                    "differenceRatio": float(differenceRatio)
                }
    return snpIndexDict

def filter_diffratio(diffratioDict, variantFilter,
                     bulk1_alleles, bulk2_alleles,
                     bulk1_refIndex, bulk2_refIndex,
                     delta_refIndex, differenceRatio):
    '''
    This function will filter a snpIndexDict as created by parse_qtlseq_snp_index_file()
    and return a new snpIndexDict object (note: no changes made to original) with only
    variants that pass all provided filters. The sub-dictionaries are merely copied,
    so the objects will be linked to each other.
    
    Parameters:
        diffratioDict -- a dict with structure like:
                      {
                          'contig1': {
                              pos1: {"variant": ___, "bulk1_alleles": ___, ...},
                              pos2: { ... },
                              ...
                          },
                          'contig2': ...,
                          ...
                      }
        ... [I don't want to comment everything else. You can figure it out based
        on the main() argparse values]
    '''
    newDiffratioDict = {}
    for contigID, contigSnpsDict in diffratioDict.items():
        for pos, variantDict in contigSnpsDict.items():
            filtersPassed = True
            
            if variantFilter != None:
                if variantDict["variant"] != variantFilter:
                    filtersPassed = False
            
            if bulk1_alleles != None:
                if variantDict["bulk1_alleles"] < bulk1_alleles:
                    filtersPassed = False
            
            if bulk2_alleles != None:
                if variantDict["bulk2_alleles"] < bulk2_alleles:
                    filtersPassed = False
            
            if bulk1_refIndex != None:
                if variantDict["bulk1_refIndex"] != "." and (variantDict["bulk1_refIndex"] < bulk1_refIndex):
                    filtersPassed = False
            
            if bulk2_refIndex != None:
                if variantDict["bulk2_refIndex"] != "." and (variantDict["bulk2_refIndex"] < bulk2_refIndex):
                    filtersPassed = False
            
            if delta_refIndex != None:
                if variantDict["delta_refIndex"] != "." and (variantDict["delta_refIndex"] < delta_refIndex):
                    filtersPassed = False
            
            if differenceRatio != None:
                if variantDict["differenceRatio"] < differenceRatio:
                    filtersPassed = False
            
            if filtersPassed:
                newDiffratioDict.setdefault(contigID, {})
                newDiffratioDict[contigID][pos] = variantDict
    return newDiffratioDict

def main():
    # User input
    usage = """%(prog)s filters a TSV file produced by calculate_persample_diffratio.py
    based on any combination of its results columns.
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="diffratioFile",
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
    p.add_argument("--bulk1_alleles", dest="bulk1_alleles",
                   type=int,
                   required=False,
                   help="Optionally, filter to retain variants >= bulk1_alleles value")
    p.add_argument("--bulk2_alleles", dest="bulk2_alleles",
                   type=int,
                   required=False,
                   help="Optionally, filter to retain variants >= bulk2_alleles value")
    p.add_argument("--bulk1_refIndex", dest="bulk1_refIndex",
                   type=int,
                   required=False,
                   help="Optionally, filter to retain variants >= bulk1_refIndex value")
    p.add_argument("--bulk2_refIndex", dest="bulk2_refIndex",
                   type=float,
                   required=False,
                   help="Optionally, filter to retain variants >= bulk2_refIndex value")
    p.add_argument("--delta_refIndex", dest="delta_refIndex",
                   type=float,
                   required=False,
                   help="Optionally, filter to retain variants >= delta_refIndex value")
    p.add_argument("--differenceRatio", dest="differenceRatio",
                   type=float,
                   required=False,
                   help="Optionally, filter to retain variants with >= difference ratio value")
    args = p.parse_args()
    validate_args(args)
    
    # Parse difference ratio file
    diffratioDict = parse_persample_diffratio_file(args.diffratioFile)
    
    # Filter difference ratio dictionary based on indicated parameters
    diffratioDict = filter_diffratio(diffratioDict,
                                    args.variantFilter,
                                    args.bulk1_alleles, args.bulk2_alleles,
                                    args.bulk1_refIndex, args.bulk2_refIndex,
                                    args.delta_refIndex, args.differenceRatio)
    
    # Write new output file
    with open(args.outputFileName, "w") as fileOut:
        # Write header line
        fileOut.write("{0}\n".format("\t".join([
            "CHROM", "POSI", "variant", "bulk1_alleles",
            "bulk2_alleles", "bulk1_refIndex",
            "bulk2_refIndex", "delta_refIndex",
            "differenceRatio"
        ])))
        
        # Write content lines
        for contigID, contigSnpsDict in diffratioDict.items():
            for pos, variantDict in contigSnpsDict.items():
                fileOut.write(
                    "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(
                        contigID, pos, variantDict["variant"],
                        variantDict["bulk1_alleles"], variantDict["bulk2_alleles"],
                        variantDict["bulk1_refIndex"], variantDict["bulk2_refIndex"],
                        variantDict["delta_refIndex"], variantDict["differenceRatio"]
                    )
                )
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

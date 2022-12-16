#! python3
# convert_vcf_to_mstmap.py
# Script to enable easy filtering of the qtlseq.vcf file so as to
# only contain the SNPs identified in the snp_index.p##.tsv file.

import os, argparse

# Define functions
def validate_args(args):
    # Validate input file location
    if not os.path.isfile(args.vcf):
        print('I am unable to locate the QTLseq VCF file (' + args.vcf + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.pops):
        print('I am unable to locate the population map file (' + args.pops + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Handle file overwrites
    if os.path.isfile(args.outputFileName):
        print(args.outputFileName + ' already exists. Delete/move/rename this file and run the program again.')
        quit()

def parse_mstmap_file(mstPopsFile, columnsList):
    '''
    The provided file is expected to conform to a format that is:
        - TSV
        - First line (regardless of #) is the header line
        - Multiple columns can exist, but there must be two columns as indicated
          by the provided columnsList
    '''
    popDict = {}
    with open(mstPopsFile, "r") as fileIn:
        firstLine = True
        for line in fileIn:
            l = line.strip("#\r\n ")
            if l == "":
                continue
            sl = l.split("\t")
            
            # Handle header line
            if firstLine == True:
                for column in columnsList:
                    assert column in sl, \
                        f"{column} not found in pop file header!"
                columnIndices = [sl.index(c) for c in columnsList]
                firstLine = False
            # Handle content lines
            else:
                sampleID, pop = [sl[index] for index in columnIndices]
                popDict[sampleID] = pop
                popDict.setdefault(pop, [])
                popDict[pop].append(sampleID)
    
    # Validations
    # >> Check that all expected keys are found
    for expectedKey in ["p1", "p2", "c1", "c2"]:
        assert expectedKey in popDict, \
            f"{expectedKey} not found in pop file!"
    
    # >> Check that there's only one parent each
    for parentKey in ["p1", "p2"]:
        assert len(popDict[parentKey]) == 1, \
            f"{parentKey} has {len(popDict[parentKey])} values! Expected 1..."
    
    # >> Check that there's more than one child each
    for childKey in ["c1", "c2"]:
        assert len(popDict[childKey]) > 1, \
            f"{childKey} has {len(popDict[childKey])} values! Expected more than 1..."
        
        # >> And check that there's no duplicate IDs
        assert len(popDict[childKey]) == len(set(popDict[childKey])), \
            f"{childKey} has duplicate sample IDs!"
    
    return popDict

def main():
    # User input
    usage = """%(prog)s receives a Freebayes VCF (genotypes predicted per-sample) as well as a population
    map TSV file which indicates the bulk1 parent (p1), bulk2 parent (p2), and their children (c1 / c2).
    
    Note: when referring to bulk1 parent, I mean the parent that has the extreme phenotype that comprises bulk1 (c1),
    and the other parent should have the opposite extreme that is characteristic of bulk2 (c2).
    """
    # Required
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-v", dest="vcf",
                   required=True,
                   help="Input VCF file produced by Freebayes")
    p.add_argument("-p", dest="pops",
                   required=True,
                   help="Input population map file in MSTmap format")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Specify the output file name")
    # Optional
    p.add_argument("--columns", dest="columns", nargs="+",
                   required=False,
                   help="""Optionally, specify the two columns that indicate 1) the sample identifier,
                   and 2) the mapping codes""",
                   default=["sampleid", "map"])
    args = p.parse_args()
    validate_args(args)
    
    # Parse population map file
    popDict = parse_mstmap_file(args.pops, args.columns)
    
    # Parse VCF and write to file
    with open(args.vcf, "r") as fileIn, open(args.outputFileName, "w") as fileOut:
        # Write new header line
        fileOut.write("locus_name\t{0}\t{1}\t{2}\t{3}\n".format(
            "p1", # popDict["p1"][0]
            "p2",
            "\t".join([f"bone{i}" for i in range(len(popDict["c1"]))]),
            "\t".join([f"btwo{i}" for i in range(len(popDict["c2"]))])
        ))
        
        # Go through VCF lines
        for line in fileIn:
            # Skip comment lines
            if line.startswith("#"):
                continue
            # Handle content lines
            else:
                '''
                Let's lay out the rules for how this works. We need to distill the 
                '''
                fileOut.write(line)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

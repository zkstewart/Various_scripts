#! python3
# persample_experiment_comparison.py
# Script to ...

# Load normal/pip packages
import os, argparse, sys, re

# Load ZS_IO Class code
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))) # 4 dirs up is where we find Function_packages
from Function_packages import ZS_GFF3IO, ZS_SeqIO

# Load functions from other scripts
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))) # 3 dirs up is where we find windows
import haplotypes

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.vcfFile):
        print(f'I am unable to locate the VCF file ({args.vcfFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.metadataFile):
        print(f'I am unable to locate the metadata file ({args.metadataFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate that behavioural parameters make sense
    numBehaviours = sum([ param for param in [args.homoParents, args.provideInfo]])
    if numBehaviours == 0:
        print("Either --homoParents or --provideInfo must be specified!")
        quit()
    elif numBehaviours > 1:
        print("Only one of --homoParents or --provideInfo should be specified!")
        quit()
    
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print(f'File already exists at output location ({args.outputFileName})')
        print('Make sure you specify a unique file name and try again.')
        quit()

def _split_tsv_or_csv(line):
    if "\t" in line:
        return line.rstrip("\r\n ").split("\t")
    else:
        assert "," in line, \
            f"File line '{line}' doesn't conform to CSV or TSV expectations!"
        return line.rstrip("\r\n ").split(",")

def parse_pops_metadata(metadataFile, ancestors, parents):
    COLUMNS_TO_PARSE = ["sample", "experiment", "bulk"]
    assert len(ancestors) == 2
    assert len(parents) == 2
    
    # Parse file
    metadataDict = {}
    foundAncestors = 0
    foundParents = 0
    foundExp1 = 0
    foundExp2 = 0
    
    with open(metadataFile, "r") as fileIn:
        header = _split_tsv_or_csv(next(fileIn))
        headerIndices = [COLUMNS_TO_PARSE.index(h) for h in header if h in COLUMNS_TO_PARSE]
        assert len(headerIndices) == 3, \
            f"Wasn't able to find {COLUMNS_TO_PARSE} in header line '{header}'"
        
        for line in fileIn:
            sl = _split_tsv_or_csv(line)
            if sl == []:
                continue
            else:
                sample, experiment, bulk = [ sl[x] for x in headerIndices ]
                if sample in ancestors[0] or sample in ancestors[1]:
                    metadataDict[sample] = "ancestor"
                    foundAncestors += 1
                elif sample in parents[0] or sample in parents[1]:
                    metadataDict[sample] = "parent"
                    foundParents += 1
                elif int(experiment) == 2: # will raise error if experiment column malformatted
                    metadataDict[sample] = "comparison"
                    foundExp2 += 1
                else:
                    metadataDict[sample] = int(bulk[-1:]) # will raise error if bulk column malformatted
                    foundExp1 += 1
    
    # Peform some quick validations
    if foundAncestors < 2:
        print("Did not find all your ancestors in the metadata file!")
        quit()
    if foundParents < 2:
        print("Did not find all your parents in the metadata file!")
        quit()
    if foundExp1 == 0:
        print("Did not find any experiment 1 samples in the metadata file!")
        quit()
    if foundExp2 == 0:
        print("Did not find any experiment 2 samples in the metadata file!")
        quit()
    
    # Print some QC information that a human might recognise as being wrong
    print("# Metadata information")
    print(f"# Found {foundExp1} samples from experiment 1 (original bulks)")
    print(f"# Found {foundExp2} samples from experiment 2 (comparison group)")
    
    return metadataDict

def is_homo(gtList):
    if not isinstance(gtList[0], list):
        flatList = gtList
    else:
        flatList = [ gt for sublist in gtList for gt in sublist ]
    return True if len(set(flatList)) == 1 else False

def is_het(gtList):
    if not isinstance(gtList[0], list):
        flatList = gtList
    else:
        flatList = [ gt for sublist in gtList for gt in sublist ]
    return True if len(set(flatList)) == 2 else False

def is_same(gtList1, gtList2):
    if not isinstance(gtList1[0], list):
        flatList1 = gtList1
    else:
        flatList1 = [ gt for sublist in gtList1 for gt in sublist ]
    if not isinstance(gtList2[0], list):
        flatList2 = gtList2
    else:
        flatList2 = [ gt for sublist in gtList2 for gt in sublist ]
    return True if set(flatList1) == set(flatList2) else False

def find_relevant_heterozygotes(snpGenotypes, metadataDict, ancestorGood, ancestorPoor,
                                parentGood, parentPoor):
    '''
    Parameters:
        snpGenotypes -- a dictionary with structure like:
                        {
                            'contigID1':
                            {
                                pos1:
                                {
                                    'sampleID1': [0, 1], # genotypes
                                    'sampleID2': [0, 0]
                                },
                                pos2: { ... },
                                ...
                            },
                            'contigID2': { ... },
                            ...
                        }
        metadataDict -- a dictionary with structure like:
                        {
                            'sampleID1': 'comparison', # each example is a
                            'sampleID2': 1,            # different category of
                            'sampleID3': 2,            # sample type; all should
                            'sampleID4': 'parent',     # present in the
                            'sampleID5': 'ancestor',   # metadataDict object
                            ...
                        }
        ancestorGood, ancestorPoor -- strings indicating the sample ID of the
                                      ancestors which have a good and bad phenotype
        parentGood, parentPoor -- strings indicating the sample ID of the
                                  parents which have a good and bad phenotype
    '''
    GOOD_PCT = 0.75
    # BAD_PCT = 0.5
    
    numBulk1 = sum([ 1 for sample, group in metadataDict.items() if group == 1 ])
    # numBulk2 = sum([ 1 for sample, group in metadataDict.items() if group == 2 ])
    
    # Iterate through SNPs looking for ones which match our aims
    selectedSNPs = {}
    for contig, posDict in snpGenotypes.items():
        for pos, sampleGTDict in posDict.items():
            # Obtain ancestor and parent GTs
            try:
                ancestorGoodGT = sampleGTDict[ancestorGood]
                ancestorPoorGT = sampleGTDict[ancestorPoor]
                
                parentGoodGT = sampleGTDict[parentGood]
                parentPoorGT = sampleGTDict[parentPoor]
            except:
                continue
            
            # See if this meets a pattern that makes sense
            pattern1, pattern2, pattern3 = False, False, False
            '''
            The pattern we want to find is when our good parent and ancestor are both the
            same homozygote. The bad parent and ancestor are both different to this homozygote.
            The precocious bulk should also be this homozygote. And our selected plants should
            be heterozygotes.
            '''
            
            isPattern = False
            if (is_homo(ancestorGoodGT) and not is_same(ancestorGoodGT, ancestorPoorGT)
            and (is_homo(parentGoodGT) and not is_same(parentGoodGT, parentPoorGT))
            and is_same(ancestorGoodGT, parentGoodGT)):
                isPattern = True
            
            if not isPattern:
                continue
            
            # Since this has a good pattern, figure out whether our bulk1 supports it
            combinedGT = set(ancestorGoodGT).union(set(ancestorPoorGT))
            assert len(combinedGT) == 2, "Unhandled genotype situation"
            
            desiredGT = sorted(combinedGT)
            numGoodBulk1 = sum([
                1
                for sample, group in metadataDict.items()
                if sample in sampleGTDict and group == 1
                and sampleGTDict[sample] == ancestorGoodGT
            ])
            # numBadBulk2 = sum([
            #     1
            #     for sample, group in metadataDict.items()
            #     if sample in sampleGTDict and group == 2
            #     and sampleGTDict[sample] != ancestorGoodGT
            # ])
            
            supportedPctGood = numGoodBulk1 / numBulk1
            if supportedPctGood < GOOD_PCT:
                continue
            
            # supportedPctBad = numBadBulk2 / numBulk2
            # if supportedPctBad < BAD_PCT:
            #     continue
            
            # Find samples which meet our wildest desires
            selectedSNPs.setdefault(contig, {})
            selectedSNPs[contig].setdefault(pos, [])
            
            for sample, group in metadataDict.items():
                if group == "comparison":
                    if sample in sampleGTDict:
                        sampleGT = sampleGTDict[sample]
                        if sampleGT == desiredGT:
                            selectedSNPs[contig][pos].append(sample)
    
    return selectedSNPs

def find_relevant_snps(snpGenotypes, metadataDict, ancestorGood, ancestorPoor,
                                parentGood, parentPoor):
    '''
    See the method header for the other function
    '''
    GOOD_PCT = 0.75
    
    numBulk1 = sum([ 1 for sample, group in metadataDict.items() if group == 1 ])
    
    # Get a consistent order for the samples
    sampleOrder = list(metadataDict.keys())
    
    # Iterate through SNPs looking for ones which match our aims
    selectedSNPs = {}
    for contig, posDict in snpGenotypes.items():
        for pos, sampleGTDict in posDict.items():
            # Obtain ancestor and parent GTs
            try:
                ancestorGoodGTs = [ sampleGTDict[ag] for ag in ancestorGood ]
                ancestorPoorGTs = [ sampleGTDict[ap] for ap in ancestorPoor ]
                
                parentGoodGTs = [ sampleGTDict[pg] for pg in parentGood ]
                parentPoorGTs = [ sampleGTDict[pp] for pp in parentPoor ]
            except:
                continue
            
            # See if this meets a pattern that makes sense
            pattern1, pattern2, pattern3 = False, False, False
            '''
            The pattern we want to find is when our good parent and ancestor are both the
            same homozygote. The bad parent and ancestor are both different to this homozygote.
            The precocious bulk should also be this homozygote. And our selected plants should
            have ANY genotype. It doesn't matter.
            '''
            isPattern = False
            if (is_homo(ancestorGoodGTs) and not is_same(ancestorGoodGTs, ancestorPoorGTs)
            and (is_homo(parentGoodGTs) and not is_same(parentGoodGTs, parentPoorGTs))
            and is_same(ancestorGoodGTs, parentGoodGTs)):
                isPattern = True
            
            if not isPattern:
                continue
            
            # Store the genotypes of each sample
            selectedSNPs.setdefault(contig, {})
            selectedSNPs[contig].setdefault(pos, [])
            
            #for sample, group in metadataDict.items():
            for sample in sampleOrder:
                #if group == "comparison":
                specifiedSample = False
                if sample in sampleGTDict:
                    # Figure out what type of SNP this is
                    sampleGT = sampleGTDict[sample]
                    if is_homo(sampleGT) and is_same(parentGoodGTs, sampleGT):
                        snpType = "goodHom"
                    elif is_homo(sampleGT):
                        snpType = "badHom"
                    else:
                        snpType = "het"
                    # Store it
                    selectedSNPs[contig][pos].append(snpType)
                    specifiedSample = True
                
                # Store a blank for this sample if we didn't find it
                if specifiedSample == False:
                    selectedSNPs[contig][pos].append(".")
    
    return selectedSNPs, sampleOrder

def main():
    # User input
    usage = """%(prog)s receives several files ...
    """
    
    # Required
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-v", dest="vcfFile",
                   required=True,
                   help="Specify the VCF file name")
    p.add_argument("-m", dest="metadataFile",
                   required=True,
                   help="Specify the metadata file name")
    # p.add_argument("-w", dest="windowsFile",
    #                required=True,
    #                help="Specify the file name containing window regions")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Specify the output file name")
    p.add_argument("-ancestor_poor", dest="ancestorPoor",
                   required=True,
                   nargs="+",
                   help="""Indicate which sample corresponds to an original ancestor
                   with poor phenotype""")
    p.add_argument("-ancestor_good", dest="ancestorGood",
                   required=True,
                   nargs="+",
                   help="""Indicate which sample corresponds to an original ancestor
                   with GOOD phenotype""")
    p.add_argument("-parent_poor", dest="parentPoor",
                   required=True,
                   nargs="+",
                   help="""Indicate which sample corresponds to an original parent
                   with poor phenotype""")
    p.add_argument("-parent_good", dest="parentGood",
                   required=True,
                   nargs="+",
                   help="""Indicate which sample corresponds to an original parent
                   with GOOD phenotype""")
    # Optional
    p.add_argument("--homoParents", dest="homoParents",
                   required=False,
                   action="store_true",
                   help="""Specify this flag to detect SNPs that meet the pattern
                   where the good parent and ancestor are homozygotes, and the good
                   offspring are heterozygotes""")
    p.add_argument("--provideInfo", dest="provideInfo",
                   required=False,
                   action="store_true",
                   help="""Specify this flag to simply generate a table of SNPs in
                   which the good parent and ancestors are both homozygotes, and the
                   bad parent is different; this table should be manually inspected
                   to make decisions""")
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse metadata dict
    metadataDict = parse_pops_metadata(args.metadataFile,
                                       [args.ancestorGood, args.ancestorPoor],
                                       [args.parentGood, args.parentPoor])
    
    # Parse VCF data for outlier SNP genotypes
    snpGenotypes = haplotypes.get_genotypes_from_vcf(args.vcfFile,
                                                     snpPositions=None,
                                                     imputeMissing=False)
    
    # Find SNPs that match our expectations for heterozygote selection
    '''
    This script is being designed specifically for secret business with a goal
    that can be summarised as 'we want to find crosses that are heterozygous,
    where the original ancestors and parents were not both heterozygous'. This
    approach could be useful for sample selection.
    '''
    if args.homoParents:
        selectedSNPs = find_relevant_heterozygotes(snpGenotypes, metadataDict,
                                                  args.ancestorGood, args.ancestorPoor,
                                                  args.parentGood, args.parentPoor)
    elif args.provideInfo:
        selectedSNPs, sampleOrder = find_relevant_snps(snpGenotypes, metadataDict,
                                                       args.ancestorGood, args.ancestorPoor,
                                                       args.parentGood, args.parentPoor)
    
    # Format results for human interpretation
    comparisonSamples = [ sample for sample, group in metadataDict.items() if group == "comparison" ]
    
    with open(args.outputFileName, "w") as fileOut:
        # Write header
        if args.homoParents:
            fileOut.write("contig\tposition\t{0}\n".format("\t".join(comparisonSamples)))
        elif args.provideInfo:
            fileOut.write("contig\tposition\t{0}\n".format("\t".join(sampleOrder)))
        # Write contents
        for contig, posDict in selectedSNPs.items():
            for pos, sampleList in posDict.items():
                if sampleList != []:
                    if args.homoParents:
                        sampleFormat = "\t".join([ "." if c not in sampleList else "Y" for c in comparisonSamples ])
                        fileOut.write(f"{contig}\t{pos}\t{sampleFormat}\n")
                    elif args.provideInfo:
                        sampleFormat = "\t".join(sampleList)
                        fileOut.write(f"{contig}\t{pos}\t{sampleFormat}\n")
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

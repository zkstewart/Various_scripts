#! python3
# filter_vcf.py
# Simple filtration options that aren't easily accessible
# through things like vcftools

import os, sys, argparse

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from Function_packages import ZS_VCFIO

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.vcfFile):
        print('I am unable to locate the VCF file (' + args.vcfFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print('File already exists at output location (' + args.outputFileName + ')')
        print('Make sure you specify a unique file name and try again.')
        quit()
    # Validate numeric inputs
    if args.missingPerPopulation < 0:
        print("missingPerPopulation should be 0 or greater")
        quit()
    elif args.missingPerPopulation > 1:
        print("missingPerPopulation should be 1 or less")
        quit()
    elif args.minimumSnps < 0:
        print("minimumSnps should be 0 or greater")
        quit()
    elif args.sampleDP < 0:
        print("sampleDP should be 0 or greater")
        quit()

## Data parsing and filtering
def parse_pops_file(popsFile):
    pops = {}
    with open(popsFile, "r") as fileIn:
        for line in fileIn:
            if line.startswith("#"):
                continue
            
            sample, population = line.rstrip("\r\n").split("\t")
            # Bidirectional indexing
            pops[sample] = population
            if population not in pops:
                pops[population] = []
            pops[population].append(sample)
    return pops

def filter_vcf_by_population_missingness(vcf, pops, missingPerPopulation=0.5):
    '''
    Receives a ZS_VCFIO.VCF object and filters it (i.e., removes sites) that
    don't have at least >= missingPerPopulation presence in each population.
    
    Parameters:
        vcf -- a ZS_VCFIO.VCF object
        pops -- a dictionary which, at minimum, contains keys (population IDs)
                which index lists that contain the samples associated with
                the population
        missingPerPopulation -- a float from 0->1 (inclusive) indicating what
                                proportion of missingness will be tolerated before
                                excluding the variant
    '''
    numPopGroups = sum([1 for key, value in pops.items() if isinstance(value, list)])
    
    # Check VCF to find variants that fail filtration
    toDrop = {}
    for contigID, contigDict in vcf.items():
        for pos, posDict in contigDict.items():
            gtIndex = posDict["FORMAT"].index("GT")
            
            # Tally presence across populations
            popsCount = {}
            for sampleID in vcf.samples:
                # Skip non-called sites
                "This should handle formats produced by Freebayes and VCFtools"
                if set(posDict[sampleID]) == {'.', './.'} or set(posDict[sampleID]) == {'.'}:
                    continue
                
                sampleGT = posDict[sampleID][gtIndex]
                samplePop = pops[sampleID]
                popsCount.setdefault(samplePop, 0)
                
                if sampleGT != "." and sampleGT != "./.":
                    popsCount[samplePop] += 1
            
            # Compute the missingness per population
            failed = False
            if popsCount == {}: # happens if ALL samples are blank after we dropped samples
                failed = True
            elif len(popsCount) != numPopGroups: # happens if ALL OF A POPULATION'S samples are blank
                failed = True
            else:
                for popID, popCount in popsCount.items():
                    pctMissing = 1 - (popCount / len(pops[popID]))
                    if pctMissing > missingPerPopulation: # if this FAILS the filter
                        failed = True
                        break
            if failed:
                toDrop.setdefault(contigID, set())
                toDrop[contigID].add(pos)
    
    # Eliminate variants if they failed the filter
    for contigID, posSet in toDrop.items():
        for pos in posSet:
            vcf.del_variant(contigID, pos)

def filter_vcf_by_sample_dp(vcf, dpMinimum=5):
    '''
    Receives a ZS_VCFIO.VCF object and modifies sample variant calls to be ambiguous
    (all "."s) if they have a DP value smaller than the given cutoff.
    
    Note that this function assumes each sample has a DP value, which should be safe
    in most instances since the VCF class automatically imputes this value if possible.
    If this wasn't able to happen, this function simply cannot work and it'll error
    out ungracefully (but still informatively).
    
    Parameters:
        vcf -- a ZS_VCFIO.VCF object
        dpMinimum -- an integer specifying the minimum DP required for a sample to have its
                     genotype called
    '''
    # Check VCF to find variants that fail filtration
    for contigID, contigDict in vcf.items():
        for pos, posDict in contigDict.items():
            gtIndex = posDict["FORMAT"].index("GT")
            dpIndex = posDict["FORMAT"].index("DP")
            
            # Check DP per sample
            for sampleID in vcf.samples:
                if posDict[sampleID] == ["."]: # freebayes does non-called like this
                    continue
                
                sampleDP = posDict[sampleID][dpIndex]
                if sampleDP == ".": # vcftools does non-called like this
                    continue
                
                elif int(sampleDP) < dpMinimum:
                    for i in range(0, len(posDict[sampleID])):
                        if i == gtIndex:
                            posDict[sampleID][i] = "./."
                        else:
                            posDict[sampleID][i] = "."

def filter_vcf_by_contig_snps(vcf, minimumSnps):
    '''
    Receives a ZS_VCFIO.VCF object and drops contigs and any SNPs annotated on them
    if the contig has fewer than the indicated SNPs located on it.
    
    Parameters:
        vcf -- a ZS_VCFIO.VCF object
        minimumSnps -- an integer specifying the minimum number of SNPs required for us
                       to retain it in the VCF
    '''
    toDrop = []
    # Check VCF to find variants that fail filtration
    for contigID, contigDict in vcf.items():
        if len(contigDict) < minimumSnps:
            toDrop.append(contigID)
    
    # Remove any contigs flagged through this
    for contigID in toDrop:
        vcf.del_contig(contigID)

def filter_vcf_where_no_alt_allele(vcf):
    '''
    Simple filter to remove variants where no samples with the variant are
    called.
    '''
    toDrop = {}
    for contigID, contigDict in vcf.items():
        for pos, posDict in contigDict.items():
            gtIndex = posDict["FORMAT"].index("GT")
            
            # Check GT per sample
            foundAlt = False
            for sampleID in vcf.samples:
                if posDict[sampleID] == ["."]: # freebayes does non-called like this
                    continue
                
                sampleGT = posDict[sampleID][gtIndex]
                if sampleGT == "." or sampleGT == "./." or sampleGT == "0/0":
                    continue
                else:
                    foundAlt = True
                    break
            
            # Fail this site if we didn't find any alt alleles
            if foundAlt == False:
                toDrop.setdefault(contigID, set())
                toDrop[contigID].add(pos)
    
    # Eliminate variants if they failed the filter
    for contigID, posSet in toDrop.items():
        for pos in posSet:
            vcf.del_variant(contigID, pos)

def filter_vcf_where_sample_is_untyped(vcf, mustBeTyped):
    '''
    Simple filter to remove variants where the samples in mustBeTyped
    are untyped.
    '''
    toDrop = {}
    for contigID, contigDict in vcf.items():
        for pos, posDict in contigDict.items():
            gtIndex = posDict["FORMAT"].index("GT")
            
            # Check GT for mustBeTyped samples
            wasNotTyped = False
            for sampleID in mustBeTyped:
                if posDict[sampleID] == ["."]: # freebayes does non-called like this
                    wasNotTyped = True
                    break
                
                sampleGT = posDict[sampleID][gtIndex]
                if sampleGT == "." or sampleGT == "./.":
                    wasNotTyped = True
                    break
            
            # Fail this site if any samples were untyped
            if wasNotTyped == True:
                toDrop.setdefault(contigID, set())
                toDrop[contigID].add(pos)
    
    # Eliminate variants if they failed the filter
    for contigID, posSet in toDrop.items():
        for pos in posSet:
            vcf.del_variant(contigID, pos)

## Main
def main():
    # User input
    usage = """%(prog)s reads in a VCF and performs various filtrations
    that aren't easily doable with other VCF programs e.g., vcftools.
    
    The population file is assumed to be tab-separated into two columns
    i.e., sample_name : population_ID.
    """
    p = argparse.ArgumentParser(description=usage)
    ## Required
    p.add_argument("-v", dest="vcfFile",
                   required=True,
                   help="Input VCF file for filtering")
    p.add_argument("-p", dest="popsFile",
                   required=True,
                   help="Input population file")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file name for the filtered SNPs")
    ## Optional
    p.add_argument("--mbt", dest="mustBeTyped", nargs="+",
                   required=False,
                   help="""Optionally, provide one or more sample IDs that must be typed
                   for us to retain the site overall""",
                   default=[])
    p.add_argument("--mpp", dest="missingPerPopulation", type=float,
                   required=False,
                   help="""This number is the minimum proportion of samples per population
                   where ambiguity will be tolerated before dropping the site; default=0.5
                   (range 0 -> 1)""",
                   default=0.5)
    p.add_argument("--min", dest="minimumSnps", type=int,
                   required=False,
                   help="""Optionally, specify how many SNPs must be on a contig
                   for it to be retained (default == 1)""",
                   default=0)
    p.add_argument("--dp", dest="sampleDP", type=int,
                   required=False,
                   help="""Optionally, specify the per-sample DP required for the
                   genotype call in that sample to not be converted to ambiguous (default=5)""",
                   default=5)
    p.add_argument("--geno", dest="genoOutput", action="store_true",
                   required=False,
                   help="""Optionally, produce a .geno formatted output rather than .vcf""",
                   default=False)
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse VCF file
    vcf = ZS_VCFIO.VCF(args.vcfFile)
    
    # Parse populations into bidirectional dictionary structure
    pops = parse_pops_file(args.popsFile)
    
    # Remove any samples that aren't in our population file
    samplesNotFound = set(vcf.samples).difference(pops.keys())
    for sampleID in samplesNotFound:
        vcf.del_sample(sampleID)
    
    # Modify VCF by per-sample DP to set sample as ambiguous if its DP is too low
    filter_vcf_by_sample_dp(vcf, args.sampleDP)
    
    # Filter VCF by population missingness/ambiguity
    filter_vcf_by_population_missingness(vcf, pops, args.missingPerPopulation)
    
    # Filter VCF by contig SNP count
    filter_vcf_by_contig_snps(vcf, args.minimumSnps)
    
    # Filter VCF to remove sites that no longer have a variant predicted
    "Might happen if the only sample(s) that had a ALT allele predicted were removed/made ambiguous"
    filter_vcf_where_no_alt_allele(vcf)
    
    # Filter VCF to remove sites where specific samples have not been typed
    if args.mustBeTyped != []:
        filter_vcf_where_sample_is_untyped(vcf, args.mustBeTyped)
    
    # Write filtered output (if relevant)
    if len(vcf) == 0:
        print("In the process of filtering this VCF, we ended up with 0 SNPs remaining!")
        print("You should check to make sure your parameters make sense...")
        print("(Since there's no SNPs left, there's nothing to write to an output file)")
    else:
        vcf.comments["footer"].append("##filter_vcf {0} {1} {2} {3}".format(
            f"mpp={args.missingPerPopulation}",
            f"min={args.minimumSnps}",
            f"dp={args.sampleDP}",
            "mbt={0}".format(",".join(args.mustBeTyped)) if args.mustBeTyped != [] else "mb=NA"
        ))
        if args.genoOutput:
            vcf.write_geno(args.outputFileName)
        else:
            vcf.write_vcf(args.outputFileName)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

#! python3
# filter_vcf.py
# Simple filtration options that aren't easily accessible
# through things like vcftools

import os, argparse

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

## Data parsing and filtering
def parse_pops_file(popsFile):
    pops = {}
    with open(popsFile, "r") as fileIn:
        for line in fileIn:
            sample, population = line.rstrip("\r\n").split("\t")
            # Bidirectional indexing
            pops[sample] = population
            if population not in pops:
                pops[population] = []
            pops[population].append(sample)
    return pops

def filter_vcf(vcfFile, pops, missingPerPopulation=0.5):
    outputLines = []
    with open(vcfFile, "r") as fileIn:
        for line in fileIn:
            l = line.rstrip("\r\n").split("\t")
            
            # Handle header lines
            if line.startswith("#CHROM"):
                samples = l[9:] # This gives us the ordered sample IDs
            if line.startswith("#"):
                outputLines.append("\t".join(l)) # Store all header lines for output
                continue
            
            # Determine which field position we're extracting for each filtration
            fieldsDescription = l[8]
            if ":" not in fieldsDescription:
                gtIndex = -1
            else:
                gtIndex = fieldsDescription.split(":").index("GT")
            
            # Filter 1: Per-population missing samples
            popsCount = {}
            ongoingCount = 0 # This gives us the index for our samples header list 
            for sampleResult in l[9:]: # This gives us the results for each sample as per fieldsDescription
                # Grab our genotype
                if gtIndex != -1:
                    genotype = sampleResult.split(":")[gtIndex]
                else:
                    genotype = sampleResult
                # Assess data presence/absence for this sample
                present, absent = 0, 0
                for g in genotype.replace("|", "/").split("/"): # Temporarily ensure all genotype fields can be iterated through
                    if g == ".":
                        absent += 1
                    else:
                        present += 1
                # Tally presence/absence data across populations
                samplePopulation = pops[samples[ongoingCount]] # Samples is our column headers; ongoingCount the sample index for this loop
                if samplePopulation not in popsCount:
                    popsCount[samplePopulation] = [0, 0]
                popsCount[samplePopulation][0] += present
                popsCount[samplePopulation][1] += absent
                
                ongoingCount += 1
            # Check if filter 1 passes
            skip = False
            for key, value in popsCount.items():
                total = sum(value)
                if (value[1] / total) > missingPerPopulation:
                    skip = True
                    break
            if skip:
                continue
            
            # Store results if all filter(s) pass
            outputLines.append("\t".join(l))
            
    return outputLines

## File out
def write_vcf_file(vcfLines, outputFileName):
    with open(outputFileName, "w") as fileOut:
        fileOut.write("\n".join(vcfLines))

def write_geno_file(vcfLines, outputFileName):
    # Find header from vcfLines
    for line in vcfLines:
        if not line.startswith("#CHROM"):
            continue
        else:
            l = line.rstrip("\r\n").split("\t")
            samples = l[9:] # This gives us the ordered sample IDs
            break
    
    # Write output geno
    with open(outputFileName, "w") as fileOut:
        # Write header to file
        fileOut.write("#CHROM\tPOS\t{0}\n".format("\t".join(samples)))
        # Write everything else
        for line in vcfLines:
            if line.startswith("#"): continue
            
            # Extract relevant info
            l = line.split("\t")
            chrom = l[0]
            pos = l[1]
            ref = l[3]
            alt = l[4]
            gtOptions = [ref, *alt.split(",")] # This enables us to deal with multi-allelic calling
            fieldsDescription = l[8]
            
            # Determine which field position we're extracting to get genotype
            if ":" not in fieldsDescription:
                gtIndex = -1
            else:
                gtIndex = fieldsDescription.split(":").index("GT")
            
            # Parse genotype per sample
            genotypes = []
            ongoingCount = 0 # This gives us the index for the sample in order
            for sampleResult in l[9:]: # This gives us the results for each sample as per fieldsDescription
                if gtIndex != -1:
                    genotype = sampleResult.split(":")[gtIndex]
                else:
                    genotype = sampleResult
                # Multi-allelic compatible genotype finding
                for i in range(0, len(gtOptions)):
                    genotype = genotype.replace(str(i), gtOptions[i]).replace(".", "N") # need to switch to N to satisfy the .geno requirements
                # Store results
                genotypes.append(genotype)
                ongoingCount += 1
            
            # Write to file
            fileOut.write("{0}\t{1}\t{2}\n".format(chrom, pos, "\t".join(genotypes)))

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
    p.add_argument("-v", dest="vcfFile", required=True,
        help="Input VCF file for filtering")
    p.add_argument("-p", dest="popsFile", required=True,
        help="Input population file")
    p.add_argument("-o", dest="outputFileName", required=True,
        help="Output file name for the filtered SNPs")
    ## Optional
    p.add_argument("--mpp", dest="missingPerPopulation", type=float,
        help="""This number if the minimum proportion of samples per population
        that is tolerated; default=0.5 (range 0 -> 1)""",
        default=0.5)
    p.add_argument("--geno", dest="genoOutput", action="store_true",
        help="""Optionally, produce a .geno formatted output rather than .vcf""",
        default=False)
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse populations into bidirectional dictionary structure
    pops = parse_pops_file(args.popsFile)
    
    # Get filtered VCF lines
    vcfLines = filter_vcf(args.vcfFile, pops, args.missingPerPopulation)
    
    # Write output write_vcf_file file
    if not args.genoOutput:
        write_vcf_file(vcfLines, args.outputFileName)
    else:
        write_geno_file(vcfLines, args.outputFileName)

if __name__ == "__main__":
    main()

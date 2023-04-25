#! python3
# vcf_to_cf.py
# A script to produce countsfile (cf) format from a
# biallelic, SNPs only VCF. Only making this script to
# handle a bug in the cflib script which is impeding
# the progress of my own projects.

import os, argparse, gzip
from Bio import SeqIO
from contextlib import contextmanager

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.vcfFile):
        print('I am unable to locate the VCF file (' + args.vcfFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.genomeFile):
        print('I am unable to locate the genome FASTA file (' + args.genomeFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print('File already exists at output location (' + args.outputFileName + ')')
        print('Make sure you specify a unique file name and try again.')
        quit()

@contextmanager
def open_vcf_file(filename):
    if filename.endswith(".gz"):
        with gzip.open(filename, "rt") as f:
            yield f
    else:
        with open(filename) as f:
            yield f

def vcf_parse_as_genotypeDict(vcfFile):
    genotypeDict = {}
    with open_vcf_file(vcfFile) as fileIn:
        for line in fileIn:
            l = line.rstrip("\r\n").split("\t")
            
            # Handle header lines
            if line.startswith("#CHROM"):
                samples = l[9:] # This gives us the ordered sample IDs
            if line.startswith("#"):
                continue
            
            # Extract relevant details
            chrom = l[0]
            pos = int(l[1])
            ref = l[3]
            alt = l[4]
            
            # Skip multiallelic variants
            if "," in alt or "." in alt:
                continue
            
            # Skip indels
            if len(ref) != len(alt):
                continue
            
            # Determine which field position we're extracting to get our genotype
            fieldsDescription = l[8]
            if ":" not in fieldsDescription:
                gtIndex = -1
            else:
                gtIndex = fieldsDescription.split(":").index("GT")
            
            # Get genotype per sample
            genotypes = []
            for sampleResult in l[9:]: # This gives us the results for each sample as per fieldsDescription
                if gtIndex != -1:
                    genotype = sampleResult.split(":")[gtIndex]
                else:
                    genotype = sampleResult
                genotype = genotype.replace("|", "/").replace(".", "0") # Impute missing data as reference allele ## Ensure all genotypes can be split by "/"
                genotypes.append(genotype)
            
            # Store in dictionary
            genotypeDict.setdefault(chrom, {})
            genotypeDict[chrom][pos] = [ref, alt, genotypes]
            
    return genotypeDict, samples

def count_genome_length(records):
    nsites = 0
    for record in records:
        seq = str(record.seq).upper()
        nsites += sum([ 1 for nuc in seq if nuc != "N" ])
    return nsites

def genotypeDict_to_cf(genotypeDict, samplesList, genomeFile, outputFileName, onlySNPs):
    '''
    cf format is like:
        A, C, G, T
        0, 0, 0, 0
    
    Where if no SNP is present, and the ref allele is T, you get:
        0, 0, 0, 2
    '''
    # Calculate nsites for use in COUNTSFILE header format
    if onlySNPs:
        genomeRecords = SeqIO.parse(open(genomeFile, 'r'), "fasta")
        nsites = count_genome_length(genomeRecords)
    else:
        nsites = sum([ 1 for chrom in genotypeDict for pos in genotypeDict[chrom] ])
    
    # Load genome records in for parsing
    genomeRecords = SeqIO.parse(open(genomeFile, 'r'), "fasta")
        
    with open(outputFileName, "w") as fileOut:
        # Write header to file
        fileOut.write("COUNTSFILE NPOP {0} NSITES {1}\n".format(len(samplesList), nsites))
        fileOut.write("CHROM POS {0}\n".format(" ".join(samplesList)))
        
        # Write everything else
        for record in genomeRecords:
            ongoingCount = 1 # This is the genomic position
            chrom = record.id
            for nucleotide in record.seq:
                # Speed up onlySNPs operation
                if onlySNPs and ((chrom not in genotypeDict) or (ongoingCount not in genotypeDict[chrom])):
                    ongoingCount += 1
                    continue
                # Generate reference genotype
                if nucleotide.lower() == "a":
                    genotype = "2,0,0,0"
                elif nucleotide.lower() == "c":
                    genotype = "0,2,0,0"
                elif nucleotide.lower() == "g":
                    genotype = "0,0,2,0"
                elif nucleotide.lower() == "t":
                    genotype = "0,0,0,2"
                elif nucleotide.lower() == "n":
                    ongoingCount += 1
                    continue # skip gaps
                else:
                    print("Unknown nucleotide encountered '{0}'".format(nucleotide))
                    print("Program must exit to prevent erroneous behaviour")
                    quit()
                
                isSNP = False
                # Produce genotypes for non-variant position
                if (chrom not in genotypeDict) or (ongoingCount not in genotypeDict[chrom]):
                    cfGenotypes = [genotype for i in range(len(samplesList))]
                # Produce genotypes for variant position
                else:
                    print("Found SNP in chrom {0}, pos {1}".format(chrom, ongoingCount))
                    ref, alt, genotypes = genotypeDict[chrom][ongoingCount]
                    if len(ref) > 1: # counts file format does not support block substitution; we must ignore this SNP
                        print("Skipping block substitution at above position")
                        cfGenotypes = [genotype for i in range(len(samplesList))] # do the same as for non-variant positions
                    else:
                        isSNP = True
                        cfGenotypes = []
                        for g in genotypes:
                            assert "|" not in g # ensure formatting compliance
                            
                            newGenotype = [0,0,0,0]
                            g = g.replace("0", ref).replace("1", alt)
                            for allele in g.split("/"):
                                if allele.lower() == "a":
                                    newGenotype[0] += 1
                                elif allele.lower() == "c":
                                    newGenotype[1] += 1
                                elif allele.lower() == "g":
                                    newGenotype[2] += 1
                                elif allele.lower() == "t":
                                    newGenotype[3] += 1
                                elif allele.lower() == ".":
                                    newGenotype[0] += 1 # impute gaps as reference allele
                                else:
                                    print("Unknown allele encountered '{0}'".format(allele))
                                    print("Program must exit to prevent erroneous behaviour")
                                    quit()
                            cfGenotypes.append(str(newGenotype).replace("[", "").replace("]", "").replace(" ", ""))
                
                # Write line to file
                if onlySNPs == False or isSNP == True:
                    outLine = "{0} {1} {2}\n".format(chrom, ongoingCount, " ".join(cfGenotypes))
                    fileOut.write(outLine)
                
                # Iterate position counter
                ongoingCount += 1

## Main
def main():
    # User input
    usage = """%(prog)s reads in a VCF and genome FASTA file to produce
    a cflib-specification .cf file for use with IQ-Tree. It will
    skip multiallelic and indel calls by default when generating the output file.
    """
    p = argparse.ArgumentParser(description=usage)
    ## Required
    p.add_argument("-v", dest="vcfFile", required=True,
        help="Input VCF file for filtering")
    p.add_argument("-g", dest="genomeFile", required=True,
        help="Input population file")
    p.add_argument("-o", dest="outputFileName", required=True,
        help="Output file name for the filtered SNPs")
    ## Optional
    p.add_argument("--snps", dest="onlySNPs", required=False, action="store_true",
        help="Optionally, produce a .cf with only SNPs present",
        default=False)
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse VCF as genotypes dictionary
    genotypeDict, samplesList = vcf_parse_as_genotypeDict(args.vcfFile)
    
    # Produce CF file
    genotypeDict_to_cf(genotypeDict, samplesList, args.genomeFile, args.outputFileName, args.onlySNPs)

    print("Program finished successfully!")

if __name__ == "__main__":
    main()

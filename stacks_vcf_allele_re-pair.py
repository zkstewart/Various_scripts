#! python3
# stacks_vcf_allele_re-pair

import os, argparse

# Define functions for later use
## Validate arguments
def validate_args(args):
        # Ensure no None arguments exist
        for key, value in vars(args).items():
                if value == None:
                        print(key + ' argument was not specified. Fix this and try again.')
                        quit()
        # Validate input file locations
        if not os.path.isfile(args.vcfFile):
                print('I am unable to locate the VCF file (' + args.vcfFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if not os.path.isfile(args.tagsFile):
                print('I am unable to locate the catalog tags file (' + args.tagsFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if not os.path.isfile(args.idsFile):
                print('I am unable to locate the IDs text file (' + args.idsFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()   
        # Ensure that integer inputs are sensible
        if args.numTrim < 0:
                print('numTrim cannot be a negative integer; fix this and try again.')
                quit()
        # Handle file overwrites
        if os.path.isfile(args.outputFileName):
                print(args.outputFileName + ' already exists. Delete/move/rename this file and run the program again.')
                quit()

def text_file_to_list(textFile):
        outList = []
        with open(textFile, 'r') as fileIn:
                for line in fileIn:
                        outList.append(line.rstrip('\r\n'))
        return outList

def vcf_stacks_parse(vcfFile):
        # Setup
        alleleDict = {}
        ongoingCount = 1        # Most outlier loci programs are 1-based w/r/t number indexing
        # Main function
        with open(vcfFile, 'r') as fileIn:
                for line in fileIn:
                        # Ensure lines can be handled as Linux-formatted
                        line = line.replace('\r', '')
                        # Skip unnecessary lines
                        if line.startswith('#') or line == '\n':
                                continue
                        # Get details
                        sl = line.split('\t')
                        alleleID = sl[2]
                        # Add to dict then iterate ongoingCount
                        alleleDict[str(ongoingCount)] = alleleID
                        ongoingCount += 1
        return alleleDict

def tags_stacks_parse(tagsFile, numTrim):
        # Setup
        consensusDict = {}
        # Main function
        with open(tagsFile, 'r') as fileIn:
                for line in fileIn:
                        # Ensure lines can be handled as Linux-formatted
                        line = line.replace('\r', '')
                        # Skip unnecessary lines
                        if line.startswith('#') or line == '\n':
                                continue
                        # Get details
                        sl = line.split('\t')
                        alleleID = sl[2]
                        consensus = sl[9][numTrim:]
                        # Add to dict
                        consensusDict[alleleID] = consensus
        return consensusDict

def stacks_pair_out(snpIDs, alleleDict, consensusDict, outputFileName):
        # Setup
        alreadyWrote = []
        with open(outputFileName, 'w') as fileOut, open(outputFileName + '_idmap', 'w') as idmapOut:
                for snp in snpIDs:
                        alleleID = alleleDict[snp]
                        if alleleID not in alreadyWrote:
                                consensus = consensusDict[alleleID]
                                fileOut.write('>' + alleleID + '\n' + consensus + '\n')
                                alreadyWrote.append(alleleID)
                        idmapOut.write(snp + '\t' + alleleID + '\n')

##### USER INPUT SECTION
usage = """%(prog)s reads in a STACKS .vcf file and corresponding .catalog.alleles file
and, using a text file of SNP IDs, generates a FASTA file corresponding to those SNPs. This
is useful when you have a list of outlier loci and you want to retrieve their allelic sequence
for downstream investigation.
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-v", "-vcf", dest="vcfFile",
                  help="Specify the .vcf file.")
p.add_argument("-t", "-tags", dest="tagsFile",
                  help="Specify the .catalog.tags.tsv file.")
p.add_argument("-n", "-numTrim", dest="numTrim", type=int,
                  help="Specify the number of base pairs to trim off the front of each consensus stack.")
p.add_argument("-i", "-ids", dest="idsFile",
               help="Specify the text file listing SNP IDs.")
p.add_argument("-o", "-outputFile", dest="outputFileName",
                   help="Output file name.")

args = p.parse_args()
validate_args(args)

# Parse text IDs file
snpIDs = text_file_to_list(args.idsFile)

# Parse VCF file to get the allele IDs for these SNPs
alleleDict = vcf_stacks_parse(args.vcfFile)

# Parse tags file to retrieve consensus sequences
consensusDict = tags_stacks_parse(args.tagsFile, args.numTrim)

# Write nonredundant FASTA file output
stacks_pair_out(snpIDs, alleleDict, consensusDict, args.outputFileName)

# All done!
print('Program completed successfully!')

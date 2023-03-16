#! python3
# predict_variant_effects.py
# Script to take in a GFF3 file and a VCF,
# predicting the effect of variants upon
# gene CDS.

import os, argparse, sys

sys.path.append(os.path.dirname(os.path.dirname(__file__))) # 2 dirs up is where we find windows
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))) # 3 dirs up is where we find GFF3IO
from windows import snp_proximity_report
from Function_packages import ZS_GFF3IO, ZS_SeqIO

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.vcf):
        print('I am unable to locate the VCF file (' + args.vcf + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.gff3):
        print('I am unable to locate the GFF3 file (' + args.gff3 + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print(f'File already exists at output location ({args.outputFileName})')
        print('Make sure you specify a unique file name and try again.')
        quit()

def get_genotypes_from_vcf(vcfFile, snpPositions=None, imputeMissing=True):
    '''
    This function will simply read in a VCF and store data in a dictionary so as to
    allow easy query operations.
    
    Parameters:
        vcfFile -- a string indicating the location of a VCF file to be parsed.
        snpPositions -- optional; providing a dictionary with structure indicated
                        below will only index SNPs at the indicated positions:
                        {
                            'contig1': set([
                                pos1, pos2, pos3, ...
                            ]),
                            'contig2': set([ ... ]),
                            ...
                        }
        imputeMissing -- optional; boolean indicating whether missing data should
                         be imputed as 0/0 or skipped.
    Returns:
        snpGenotypes -- a dictionary with structure like:
                        {
                            'contig1': {
                                pos1: {
                                    "ref_alt": [refAllele, [altAllele1, ...]],
                                    'sample1': genotype,
                                    'sample2': genotype,
                                    ...
                                },
                                ...
                            },
                            'contig2': {
                                "ref_alt": [...],
                                pos1: { ... },
                                ...
                            }
                            ...,
                        }
    '''
    snpGenotypes = {}
    with open(vcfFile, "r") as fileIn:
        for line in fileIn:
            l = line.rstrip("\r\n").replace('"', '').split("\t") # replace quotations may help with Excel files
            
            # Handle header lines
            if line.startswith("#CHROM"):
                samples = l[9:] # This gives us the ordered sample IDs
            if line.startswith("#"):
                continue
            
            # Extract relevant details of the SNP
            chrom = l[0]
            pos = int(l[1])
            ref = l[3]
            alt = l[4].split(",")
            
            # Skip indexing this line if snpPositions is provided
            if snpPositions != None:
                if not (chrom in snpPositions and pos in snpPositions[chrom]):
                    continue
            
            # Determine which field position we're extracting to get our GT value
            fieldsDescription = l[8]
            if ":" not in fieldsDescription:
                gtIndex = -1
            else:
                gtIndex = fieldsDescription.split(":").index("GT")
            
            # Format a dictionary to store sample genotypes for this position
            posGenotypeDict = {}
            ongoingCount = 0 # This gives us the index for our samples header list 
            for sampleResult in l[9:]: # This gives us the results for each sample as per fieldsDescription
                # Grab our genotype
                if gtIndex != -1:
                    genotype = sampleResult.split(":")[gtIndex]
                else:
                    genotype = sampleResult
                
                # Edit genotype to have a consistently predictable separator
                "We don't care if the VCF is phased or not for this function"
                genotype = genotype.replace("/", "|")
                
                # Impute empty genotypes (if applicable) or skip otherwise
                if imputeMissing == True:
                    genotype = genotype.replace(".", "0")
                else:
                    if "." in genotype:
                        ongoingCount += 1
                        continue
                
                # Parse and store genotype
                samplePopulation = samples[ongoingCount]
                posGenotypeDict[samplePopulation] = list(map(int, genotype.split("|")))
                
                ongoingCount += 1
            
            # Check to see if this genotype, after imputation, still has a variant allele
            alleles = set([allele for key, allelePair in posGenotypeDict.items() if key != "ref_alt" for allele in allelePair])
            if alleles == {0}: # skip if it doesn't
                continue
            
            # Store in our overall dictionaries
            snpGenotypes.setdefault(chrom, {})
            snpGenotypes[chrom][pos] = posGenotypeDict
            snpGenotypes[chrom][pos]["ref_alt"] = [ref, *alt]
    
    return snpGenotypes

def get_genotyped_snps_for_genes(gff3Obj, snpGenotypes):
    '''
    Parameters:
        gff3Obj -- a ZS_GFF3IO.GFF3 object with NCLS indexing of the gene features
        snpGenotypes -- a dictionary with structure like:
                        {
                            'contig1': {
                                "ref_alt": [refAllele, [altAllele1, ...]],
                                pos1: {
                                    'sample1': genotype,
                                    'sample2': genotype,
                                    ...
                                },
                                ...
                            },
                            'contig2': {
                                "ref_alt": [...],
                                pos1: { ... },
                                ...
                            }
                            ...,
                        }
    Returns:
        geneSnpDict -- a dictionary with structure like:
                       {
                           'geneID1': {
                               "ref_alt": [refAllele, [altAllele1, ...]],
                               pos1: {
                                   'sample1': genotype,
                                   'sample2': genotype,
                                   ...
                               }
                           },
                           ...
                       }
    '''
    geneSnpDict = {}
    
    for contig, positionsDict in snpGenotypes.items():
        for pos, genotypesDict in positionsDict.items():
            # Skip indexed values that don't point to positions
            if not isinstance(pos, int):
                continue
            
            # If the SNP is located within a gene
            matches = gff3Obj.ncls_finder(pos, pos, "contig", contig)
            if matches != []:
                # Get the longest/representative mRNA for each gene
                mrnaFeatures = [ZS_GFF3IO.GFF3.longest_isoform(m) for m in matches]
                
                # For each mRNA, locate the SNP in the gene
                for feature in mrnaFeatures:
                    snpLocation, geneID = snp_proximity_report.locate_snp_within_gene([feature], pos) # hacky use of existing func
                    
                    # If we found the SNP in the CDS, index it
                    if snpLocation == "CDS":
                        geneSnpDict.setdefault(geneID, {})
                        geneSnpDict[geneID][pos] = genotypesDict
        
    return geneSnpDict

def get_variant_effects(gff3Obj, genomeFASTA_obj, geneSnpDict):
    '''
    Parameters:
        gff3Obj -- a ZS_GFF3IO.GFF3 object containing all the genes
                   referred to in geneSnpDict
        genomeFASTA_obj -- a ZS_SeqIO.FASTA object containing the contigs
                   that are referred to in our gff3Obj
        geneSnpDict -- a dictionary with structure like:
                       {
                           'geneID1': {
                               pos1: {
                                   "ref_alt": [refAllele, [altAllele1, ...]],
                                   'sample1': genotype,
                                   'sample2': genotype,
                                   ...
                               }
                           },
                           ...
                       }
        minSnps -- an integer indicating how many SNPs must be
                   present in a gene for us to make a haplotype of it.
                   Genes without this many SNPs will not have any
                   sequences reported.
    Returns:
        effectDict -- a dictionary with structure like:
                        {
                            "geneID1": {
                                pos1: {
                                    "ref": "nucleotide (amino_acid)",
                                    "alts": [
                                        "nucleotide (amino acid)",
                                        "nucleotide (amino acid)",
                                        ...
                                    ],
                                    "homozygotes": [
                                        "sample1 (nucleotide)",
                                        "sample2 (nucleotide)",
                                        ...
                                    ],
                                    "heterozygotes": [
                                        "sample3 (nucleotide / nucleotide)",
                                        "sample4 (nucleotide / nucleotide)",
                                        ...
                                    ]
                                },
                                pos2: {
                                    ...
                                }
                            },
                            "geneID2": {
                                ...
                            }
                        }
    '''
    effectDict = {}
    for geneID, snpDict in geneSnpDict.items():
        # Get the representative mRNA feature and sequence
        mrnaFeature = ZS_GFF3IO.GFF3.longest_isoform(gff3Obj[geneID])
        cds_FastASeq_obj, cds_featureType, cds_startingFrame = gff3Obj.retrieve_sequence_from_FASTA(genomeFASTA_obj, mrnaFeature.ID, "CDS")
        
        # Update snpDict coordinates to be localised to the extracted CDS (rather than genomic coordinates)
        newSnpDict = convert_vcf_snps_to_cds_snps(mrnaFeature, snpDict)
        
        # Skip if we've found no variants actually within the CDS
        if newSnpDict == {}:
            continue
        
        # Loop through this gene's variants and ...
        for pos, genotypeDict in newSnpDict.items():
            # Predict variant details
            formattedRef, formattedAlts, homozygotes, heterozygotes \
                = predict_variant(genotypeDict, pos, cds_FastASeq_obj, mrnaFeature.strand)
            
            # Store data in our dictionary
            effectDict.setdefault(geneID, {})
            effectDict[geneID][pos] = {
                "ref": formattedRef,
                "alts": formattedAlts,
                "homo": homozygotes,
                "hetero": heterozygotes
            }
    return effectDict

def predict_variant(genotypeDict, pos, cds_FastASeq_obj, strand):
    '''
    Parameters:
        genotypeDict -- a dictionary with structure like:
                        {
                            "ref_alt": [refAllele, [altAllele1, ...]],
                            'sample1': genotype,
                            'sample2': genotype,
                            ...
                        }
        pos -- an integer indicating the location of this SNP in CDS.
        cds_FastASeq_obj -- a ZS_SeqIO.FastASeq object containing the reference
                            coding sequence for the gene in question.
        strand -- a value in the list ["+", "-"] indicating whether the gene
                  is encoded on the +ve or -ve strand.
    Returns:
        formattedRef -- 
        formattedAlts -- 
        homozygotes -- 
        heterozygotes -- 
    '''
    # Figure out what amino acid our reference allele codes for
    ref = genotypeDict["ref_alt"][0]
    refRange = range(pos, pos + len(ref))
    
    refCodons = [
        cds_FastASeq_obj.seq[x:x+3]
            for x in range(0, len(cds_FastASeq_obj.seq), 3)
                if set(range(x,x+3)).intersection(refRange) != set()
    ]
    refAminos = [
        ZS_SeqIO.FastASeq.dna_to_protein(codon)
            for codon in refCodons
    ]
    
    formattedRef = "{0}({1})".format(ref, ",".join(refAminos))
    
    # Figure out what amino acid our alt alleles code for
    formattedAlts = []
    for alt in genotypeDict["ref_alt"][1:]:
        altRange = range(pos, pos + len(alt))
        altSequence = edit_reference_to_alt_sequence(cds_FastASeq_obj.seq[:], pos, ref, alt, strand)
        
        altCodons = [
            altSequence[x:x+3]
                for x in range(0, len(altSequence), 3)
                    if set(range(x,x+3)).intersection(altRange) != set()
        ]
        altAminos = [
            ZS_SeqIO.FastASeq.dna_to_protein(codon)
                for codon in altCodons
        ]
        
        formattedAlt = "{0}({1})".format(alt, ",".join(altAminos))
        formattedAlts.append(formattedAlt)
    
    # Partition samples into homozygotes and heterozygotes
    homozygotes = [
        sample + "({0})".format(genotypeDict["ref_alt"][genotypeDict[sample][0]]) # any index will do
        
            for sample in genotypeDict.keys()
                if sample != "ref_alt"
                and len(set(genotypeDict[sample])) == 1
    ]
    homozygotes.sort(key = lambda x: 
        (
            list(map(len, x.split("_"))), # try to sort by longest value in _ separated list; order by species name??
            x.split("(")[-1].rstrip(")"), # sort by the variant
            list(map(int, [
                thing for thing in x.split("_") if thing.isdigit()
            ])) # sort by any digits in the _ separated list
        )
    )
    
    heterozygotes = [
        sample + "({0})".format(genotypeDict["ref_alt"][genotypeDict[sample][0]] + "/" + genotypeDict["ref_alt"][genotypeDict[sample][1]])
        
        for sample in genotypeDict.keys()
            if sample != "ref_alt"
            and len(set(genotypeDict[sample])) != 1
    ]
    heterozygotes.sort(key = lambda x: 
        (
            list(map(len, x.split("_"))),
            x.split("(")[-1].rstrip(")"), # sort by the variant,
            list(map(int, [
                thing for thing in x.split("_") if thing.isdigit()
            ]))
        )
    )
    return formattedRef, formattedAlts, homozygotes, heterozygotes

def convert_vcf_snps_to_cds_snps(mrnaFeature, snpDict, embedOriginalPos=False):
    '''
    Receives a mRNA feature and a dictionary indicating SNP locations as interpreted
    from a VCF, and alters the positions to point to locations in the CDS where edits
    should be made. This function handles +ve and -ve stranded mRNA features differenly
    to give the appropriate coordinates in a 5' -> 3' reading direction.
    
    Parameters:
        mrnaFeature -- a ZS_GFF3IO.Feature object representing a mRNA
        snpDict -- a dictionary with structure like:
                   {
                       pos1: { ... } # contents of dictionary don't matter
                   }
        embedOriginalPos -- optional; boolean to indicate whether the output dict
                            should also index the original SNP index. Provided as
                            optional since I'm adding this functionality into
                            legacy code and I don't want to break stuff.
    '''
    newSnpDict = {}
    ongoingCount = 0
    for cdsFeature in mrnaFeature.CDS:
        # Check each position to see if we need to localise it
        for pos, genotypeDict in snpDict.items():
            # If the position overlaps this CDS section
            if pos >= cdsFeature.start and pos <= cdsFeature.end:
                # Get the strand-adjusted position
                if mrnaFeature.strand == "+":
                    newPos = (pos - cdsFeature.start) + ongoingCount
                else:
                    newPos = cdsFeature.end - pos + ongoingCount
                # Handle splice site variants
                '''
                Originally, I was trying to make this function do fancy things with splice
                sites like allowing variants to run up to the closest splice site motif.
                But, it gets complex because then we'd need to re-adjust the whole
                gene model just in case theres an intron variant that changes the acceptor
                site. I'd need to dedicate an entire program to just accounting for this
                kind of scenario. Ultimately, I'm deciding to make a trade off here against
                complexity in favour of simplicity because, if we find a haplotype that
                is inducing alterations to splice donor/acceptor sites that has a position
                in the CDS responsible for this, we should still detect it. Selecting for
                that SNP will still give us our favourable haplotype, theoretically, so long
                as the acceptor variant is in linkage disequilibrium which is likely.
                '''
                refAllele = genotypeDict["ref_alt"][0]
                skipThisPos = False
                if (pos + len(refAllele) - 1) > cdsFeature.end:
                    # Modify the ref allele to be contained within the exon
                    allowedAlleleLength = cdsFeature.end - pos + 1
                    newRefAllele = refAllele[0:allowedAlleleLength]
                    # Modify alt allele(s)
                    newAltAlleles = [allele[0:allowedAlleleLength] for allele in genotypeDict["ref_alt"][1:]]
                    # If an alt allele is identical to our reference, eliminate it now
                    deleteIndices = []
                    for x in range(len(newAltAlleles)):
                        if newAltAlleles[x] == newRefAllele:
                            deleteIndices.append(x)
                            for sampleID, genotype in genotypeDict.items():
                                if sampleID != "ref_alt":
                                    for z in range(len(genotype)):
                                        if genotype[z] == x+1: # x+1 gives the index of our alt allele in ["ref_alt"]
                                            genotype[z] = 0 # set it to the ref allele index
                                    genotypeDict[sampleID] = genotype
                    for index in deleteIndices[::-1]:
                        del newAltAlleles[index] # remove it from our alt alleles values
                    # If we no longer have any variants contained within the CDS region, eliminate this variant position
                    if newAltAlleles == []:
                        skipThisPos = True
                    # Otherwise, update the ref_alt allele in our dictionary
                    else:
                        genotypeDict["ref_alt"] = [newRefAllele, *newAltAlleles]
                # Handle normal scenarios / index the modified alleles if relevant
                if skipThisPos is False:
                    newSnpDict[newPos] = genotypeDict
                    if embedOriginalPos is True:
                        newSnpDict[newPos]["originalPos"] = pos
        ongoingCount += cdsFeature.end - cdsFeature.start + 1 # feature coords are 1-based inclusive, so 1->1 is a valid coord
    return newSnpDict

def edit_reference_to_alt_sequence(referenceSeq, pos, refAllele, varAllele, strand):
    '''
    This function will take in a reference nucleotide sequence, typically representing
    a CDS for a gene in either +ve or -ve strand, and generates the haplotype version
    of that sequence.
    
    Parameters:
        referenceSeq -- the sequence as a string prior to any editing
        pos -- an integer indicating what position in the referenceSeq should have its allele replaced
        refAllele -- a string including one or more nucleotides starting at the given position
        varAllele -- a string including one or more nucleotides for substitution at the given position
        strand -- a string equal to "+" if +ve stranded, or "-" if -ve stranded
    Returns:
        editedSeq -- an edited version of the input referenceSeq with all variations made
    '''
    # Get strand-appropriate alleles
    if strand == "-":
        refAllele = ZS_SeqIO.FastASeq.get_reverse_complement(None, refAllele)
        varAllele = ZS_SeqIO.FastASeq.get_reverse_complement(None, varAllele)
    
    # Validate that our position is correct and edit the sequence
    if strand == "+":
        assert referenceSeq[pos:pos+len(refAllele)].upper() == refAllele.upper(), \
            "Zac, you need to fix your +ve haplotype positioning code!"
        referenceSeq = referenceSeq[:pos] + varAllele + referenceSeq[pos+len(refAllele):]
    else:
        assert referenceSeq[pos-len(refAllele)+1:pos+1].upper() == refAllele.upper(), \
            "Zac, you need to fix your -ve haplotype positioning code!"
        referenceSeq = referenceSeq[:pos-len(refAllele)+1] + varAllele + referenceSeq[pos+1:]
    
    editedSeq = referenceSeq # just for clarity since this sequence is a new, modified object
    return editedSeq.upper() # we'd like all sequences to be upper cased

## Main
def main():
    # User input
    usage = """%(prog)s reads in GFF3 and a VCF file and generates a tabular file
    indicating what coding sequence changes exist relative to the reference coding
    sequence. Outputs will be written to the provided file as TSV.
    """
    p = argparse.ArgumentParser(description=usage)
    ## Required
    p.add_argument("-v", dest="vcf",
                   required=True,
                   help="Input phased VCF")
    p.add_argument("-g", dest="gff3",
                   required=True,
                   help="Input GFF3 file")
    p.add_argument("-f", dest="fasta",
                   required=True,
                   help="Input genome FASTA file")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Specify file name to write output to")
    args = p.parse_args()
    validate_args(args)
    
    # Parse GFF3 with NCLS indexing
    gff3Obj = ZS_GFF3IO.GFF3(args.gff3, strict_parse=False)
    gff3Obj.sort_CDS() # makes sure our CDS children lists work properly
    gff3Obj.create_ncls_index(typeToIndex="gene")
    
    # Parse VCF to get SNP genotypes
    snpGenotypes = get_genotypes_from_vcf(args.vcf)
    
    # Parse genome FASTA
    genomeFASTA_obj = ZS_SeqIO.FASTA(args.fasta)
    
    # Get the SNPs that are present within gene model CDS regions
    geneSnpDict = get_genotyped_snps_for_genes(gff3Obj, snpGenotypes)
    
    # Predict variant effects for each gene
    effectDict = get_variant_effects(gff3Obj, genomeFASTA_obj, geneSnpDict)
    
    # Create output files
    with open(args.outputFileName, "w") as fileOut:
        # Write header
        fileOut.write("{0}\n".format("\t".join([
            "#gene_ID", "ref_allele", "alt_alleles",
            "homozygotes", "heterozygotes"
        ])))
        
        # Write content lines
        for geneID, snpDict in effectDict.items():
            for pos, effect in snpDict.items():
                fileOut.write("{0}\n".format("\t".join([
                    geneID, effect["ref"],
                    ",".join(effect["alts"]),
                    ";".join(effect["homo"]),
                    ";".join(effect["hetero"])
                ])))
            
    # Let user know everything went swimmingly
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

#! python3
# predict_gene_haplotypes.py
# Script to take in a GFF3 file and a phased VCF,
# predicting haplotypes for gene sequences (CDS only)

import os, argparse, sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) # 2 dirs up is where we find windows
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))) # 3 dirs up is where we find GFF3IO
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
    # Validate numeric inputs
    if 0 > args.minSnps:
        print("minSnps should be 0 or greater")
        quit()
    # Validate output file location
    if os.path.isdir(args.outputDirectory) and os.listdir(args.outputDirectory) != []:
        print('Specified output directory already exists and contains files.')
        print('This script would rather not mess with the possibility of overwriting stuff.')
        print('Make sure you specify a non-existing directory or delete the files in the directory and try again.')
        quit()
    if not os.path.isdir(args.outputDirectory):
        os.makedirs(args.outputDirectory)
        print(f"Output directory '{args.outputDirectory}' has been created as part of argument validation.")

def get_phased_genotypes_from_vcf(vcfFile, assumePhased=False):
    '''
    This function expects the provided VCF file to contain phased genotypes
    within the GT field separated by |, and genotypes that were unphased as \.
    
    Missing genotypes will be imputed as 0|0, since this may be the case when
    merging VCFs as in the workflow of using Clair3.
    
    Ambiguous genotypes (e.g., 0/1) will simply be imputed as 0|0; if doing so
    makes the variant "disappear" such that all samples are now 0|0, we'll
    just ignore the variant call.
    
    Returns:
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
            
            # Determine which field position we're extracting to get our GT value
            fieldsDescription = l[8]
            if ":" not in fieldsDescription:
                gtIndex = -1
            else:
                gtIndex = fieldsDescription.split(":").index("GT")
            
            # Extract relevant details of the SNP
            chrom = l[0]
            pos = int(l[1])
            ref = l[3]
            alt = l[4].split(",")
            
            # Format a dictionary to store sample genotypes for this position
            posGenotypeDict = {}
            ongoingCount = 0 # This gives us the index for our samples header list 
            for sampleResult in l[9:]: # This gives us the results for each sample as per fieldsDescription
                # Grab our genotype
                if gtIndex != -1:
                    genotype = sampleResult.split(":")[gtIndex]
                else:
                    genotype = sampleResult
                
                # Edit genotype if "/" ambiguity is given for a homozygous allele
                if ("/" in genotype and len(set(genotype.split("/"))) == 1) or assumePhased == True:
                    genotype = genotype.replace("/", "|")
                
                # Impute empty genotypes
                genotype = genotype.replace(".", "0")
                
                # Impute ambiguous genotypes
                if "/" in genotype:
                    genotype = "0|0"
                assert "/" not in genotype
                
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

def get_haplotype_sequences(gff3Obj, genomeFASTA_obj, geneSnpDict, minSnps=0):
    '''
    Parameters:
        gff3Obj -- a ZS_GFF3IO.GFF3 object containing all the genes
                   referred to in geneSnpDict
        genomeFASTA_obj -- a ZS_SeqIO.FASTA object containing the contigs
                   that are referred to in our gff3Obj
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
        minSnps -- an integer indicating how many SNPs must be
                   present in a gene for us to make a haplotype of it.
                   Genes without this many SNPs will not have any
                   sequences reported.
    Returns:
        haploSeqDict -- a dictionary with structure like:
                        {
                            "geneID1": {
                                "reference_sequence": 'ATCG...',
                                "sample_haplotypes": {
                                    "sampleID1": [[hapCode, hapCode, ...], [hapCode, hapCode, ...]],
                                    "sampleID2": [[...], [...]],
                                    ...
                                },
                                "haplotype_sequences": {
                                    (hapCode, hapCode, ...): 'ATCG...',
                                    (hapCode, hapCode, ...): 'AGCG...',
                                    ...
                                }
                            },
                            "geneID2": { ... },
                            ...
                        }
    '''
    assert isinstance(minSnps, int), \
        "minSnps must be an integer value"
    
    haploSeqDict = {}
    for geneID, snpDict in geneSnpDict.items():
        # Skip genes which fail minSnps threshold
        if len(snpDict) < minSnps: # this len gives us the number of SNPs in this gene
            continue
        
        # Get the representative mRNA feature and sequence
        mrnaFeature = ZS_GFF3IO.GFF3.longest_isoform(gff3Obj[geneID])
        cds_FastASeq_obj, cds_featureType, cds_startingFrame = gff3Obj.retrieve_sequence_from_FASTA(genomeFASTA_obj, mrnaFeature.ID, "CDS")
        
        # Update snpDict coordinates to be localised to the extracted CDS (rather than genomic coordinates)
        newSnpDict = convert_vcf_snps_to_cds_snps(mrnaFeature, snpDict)
        
        # Skip if we've found no variants actually within the CDS
        if newSnpDict == {}:
            continue
        
        # Get haplotype codes for all samples
        haploCodeDict = {sampleID:[[], []] for genotypeDict in newSnpDict.values() for sampleID in genotypeDict.keys() if sampleID != "ref_alt"}
        
        orderedPositions = list(newSnpDict.keys())
        orderedPositions.sort()
        for pos in orderedPositions:
            genotypeDict = newSnpDict[pos]
            for sampleID, genotype in genotypeDict.items():
                if sampleID == "ref_alt":
                    continue
                
                haploCodeDict[sampleID][0].append(genotype[0])
                haploCodeDict[sampleID][1].append(genotype[1])
        
        # Get just the unique haplotypes so we can avoid redundant sequence extraction
        uniqueHaplotypes = set()
        for genotypes in haploCodeDict.values():
            uniqueHaplotypes.add(tuple(genotypes[0]))
            uniqueHaplotypes.add(tuple(genotypes[1]))
        uniqueHaplotypes = list(uniqueHaplotypes)
        
        # Extract haplotype sequences for all unique haplotypes
        haplotypeSequences = []
        for haplotype in uniqueHaplotypes:
            haplotypeSequence = edit_reference_to_haplotype_sequence(cds_FastASeq_obj.seq[:], haplotype, orderedPositions, newSnpDict, mrnaFeature.strand)
            haplotypeSequences.append(haplotypeSequence)
        
        # Store in our overall dictionary
        haploSeqDict[geneID] = {
            "reference_sequence": cds_FastASeq_obj.seq[:],
            "sample_haplotypes": haploCodeDict,
            "haplotype_sequences": {uniqueHaplotypes[i]:haplotypeSequences[i] for i in range(len(haplotypeSequences))}
        }
    return haploSeqDict

def convert_vcf_snps_to_cds_snps(mrnaFeature, snpDict):
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
        ongoingCount += cdsFeature.end - cdsFeature.start + 1 # feature coords are 1-based inclusive, so 1->1 is a valid coord
    return newSnpDict

def edit_reference_to_haplotype_sequence(referenceSeq, haplotype, orderedPositions, snpDict, strand):
    '''
    This function will take in a reference nucleotide sequence, typically representing
    a CDS for a gene in either +ve or -ve strand, and generates the haplotype version
    of that sequence.
    
    Parameters:
        referenceSeq -- the sequence as a string prior to any editing
        haplotype -- a list or tuple containing integers which are ordered with
                     respect to orderedPositions and refAltList
        orderedPositions -- a list or tuple containing integers which point to positions
                            in the referenceSeq that should be edited
        snpDict -- a dictionary with (at least) a structure like:
                   {
                       pos1: { "ref_alt": ['variantNucleotides1', 'variantNucleotides2', ...]},
                       pos2: { ... },
                       ...
                   }
    Returns:
        editedSeq -- an edited version of the input referenceSeq with all variations made
    '''
    for i in range(len(haplotype)-1, -1, -1): # iterate backwards through positions and variants
        pos = orderedPositions[i]
        variant = haplotype[i]
        
        # Get strand-appropriate alleles
        if strand == "+":
            refAllele = snpDict[pos]["ref_alt"][0] # ref allele might be more than 1 character long
            varAllele = snpDict[pos]["ref_alt"][variant]
        else:
            refAllele = ZS_SeqIO.FastASeq.get_reverse_complement(None, snpDict[pos]["ref_alt"][0])
            varAllele = ZS_SeqIO.FastASeq.get_reverse_complement(None, snpDict[pos]["ref_alt"][variant])
        
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

def sequence_function_alteration_inference(referenceSequence, modifiedSequence):
    '''
    The purpose of this function is to make a quick heuristic guess as to whether
    a sequence modification is likely to result in any significant change at the
    protein level.
    
    We look at a few things:
        1) Are there amino acid changes?
        2) Is there an internal stop?
        3) Did the final stop codon get removed?
    
    Parameters:
        referenceSequence -- a string indicating a nucleotide sequence that corresponds
                             to the reference genome sequence (for example)
        modifiedSequence -- a string indicating a modified sequence of the original reference
    Returns:
        modIsIdentical -- a boolean indicating whether changes have occurred or not in the
                          coding portion of the sequence
        internalStopAddition -- a boolean indicating whether an internal stop codon has been
                                added in
        finalStopRemoval -- a boolean indicating whether the final stop codon has been
                            removed
    '''
    ref_FastASeq_obj = ZS_SeqIO.FastASeq("ref", referenceSequence)
    mod_FastASeq_obj = ZS_SeqIO.FastASeq("mod", modifiedSequence)
    
    # Get translations
    refProteinSeq, _, _ = ref_FastASeq_obj.get_translation(findBestFrame=False, strand=1, frame=0)
    modProteinSeq, _, _ = mod_FastASeq_obj.get_translation(findBestFrame=False, strand=1, frame=0)
    
    # Get heuristic features from the sequence translations
    refHasNoInternalStops = True if "*" not in refProteinSeq[:-1] else False
    modHasNoInternalStops = True if "*" not in modProteinSeq[:-1] else False
    
    refHasStop = True if refProteinSeq[-1] == "*" else False
    modHasStop = True if modProteinSeq[-1] == "*" else False
    
    # Make conclusions from heuristics
    modIsIdentical = True if modProteinSeq[:-1] == refProteinSeq[:-1] else False
    internalStopAddition = True if (refHasNoInternalStops and not modHasNoInternalStops) else False
    finalStopRemoval = True if (refHasStop and not modHasStop) else False
    
    return modIsIdentical, internalStopAddition, finalStopRemoval

## Main
def main():
    # User input
    usage = """%(prog)s reads in GFF3 and a phased VCF file and generates haplotype
    predictions for genes containing a variant within its CDS. Outputs will be written
    to the provided directory.
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
    p.add_argument("-o", dest="outputDirectory",
                   required=True,
                   help="Specify location to write output files to")
    ## Optional
    p.add_argument("--minSnps", dest="minSnps",
                   type=int,
                   required=False,
                   help="""Optionally, specify how many SNPs are needed
                   to be considered a haplotype (default==2)""",
                   default=2)
    p.add_argument("--assumePhased", dest="assumePhased",
                   action="store_true",
                   required=False,
                   help="""Optionally, specify this flag if your VCF isn't
                   phased but you want to treat it like it is anyway""",
                   default=False)
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse GFF3 with NCLS indexing
    gff3Obj = ZS_GFF3IO.GFF3(args.gff3, strict_parse=False)
    gff3Obj.sort_CDS() # makes sure our CDS children lists work properly
    gff3Obj.create_ncls_index(typeToIndex="gene")
    
    # Parse VCF to get SNP genotypes
    snpGenotypes = get_phased_genotypes_from_vcf(args.vcf, args.assumePhased)
    
    # Parse genome FASTA
    genomeFASTA_obj = ZS_SeqIO.FASTA(args.fasta)
    
    # Get the SNPs that are present within gene model CDS regions
    geneSnpDict = get_genotyped_snps_for_genes(gff3Obj, snpGenotypes)
    
    # Get haplotypes for each gene and sample
    haploSeqDict = get_haplotype_sequences(gff3Obj, genomeFASTA_obj, geneSnpDict, args.minSnps)
    
    # Create output files
    with open(os.path.join(args.outputDirectory, "haplotype_predictions.tsv"), "w") as fileOut:
        fileOut.write("gene_ID\thaplotype_code\tfrequency_in_population\tsamples_with_haplotype\tlikely_change\n")
        
        for geneID, haplotypesDict in haploSeqDict.items():
            ## Output 1: Tabulation file
            # Calculate haplotype frequencies
            haplotypeCount = {"".join(map(str, genotype)):0 for genotype in list(map(list, haplotypesDict["haplotype_sequences"].keys()))}
            for genotypes in haplotypesDict["sample_haplotypes"].values():
                for gt in genotypes:
                    haplotypeCount["".join(map(str, gt))] += 1
            haplotypeFrequency = [[gt, count / sum(haplotypeCount.values())] for gt, count in haplotypeCount.items()]
            haplotypeFrequency.sort(key = lambda x: -x[1])
            
            # Get samples associated to each haplotype
            haplotypeSampleAssoc = {gtPair[0]: [] for gtPair in haplotypeFrequency}
            for sampleID, genotypes in haplotypesDict["sample_haplotypes"].items():
                if genotypes[0] == genotypes[1]:
                    haplotypeSampleAssoc["".join(map(str, genotypes[0]))].append(f"{sampleID}_(2/2)")
                else:
                    for gt in genotypes:
                        haplotypeSampleAssoc["".join(map(str, gt))].append(f"{sampleID}_(1/2)")
            
            # Write gene rows
            for i in range(len(haplotypeFrequency)):
                row = []
                
                # Handle first row for the gene
                if i == 0:
                    row.append(geneID)
                else:
                    row.append("")
                
                # Add details to row
                row.append(f"H{haplotypeFrequency[i][0]}")
                row.append(str(haplotypeFrequency[i][1]))
                row.append(", ".join(sorted(haplotypeSampleAssoc[haplotypeFrequency[i][0]])))
                
                # Figure out if any relevant changes have occurred
                modIsIdentical, internalStopAddition, finalStopRemoval = \
                    sequence_function_alteration_inference(haplotypesDict["reference_sequence"],
                        haplotypesDict["haplotype_sequences"][tuple(map(int,tuple(haplotypeFrequency[i][0])))]) # good god
                if modIsIdentical:
                    row.append("no_change")
                elif internalStopAddition:
                    row.append("internal_stop")
                elif finalStopRemoval:
                    row.append("final_stop_deletion")
                else:
                    row.append("amino_acid_change")
                
                # Write to file
                fileOut.write("\t".join(row) + "\n")
            
            ## Output 2: FASTA files
            fastaFilePrefix = os.path.join(args.outputDirectory, f"{geneID}_haplotypes")
            refSeqID = f"{geneID}_seq{i+1} haplotypeCode={'0'*len(haplotypeFrequency[0][0])} frequency=REFERENCE"
            with open(fastaFilePrefix + ".nucl.fasta", "w") as nuclFileOut, open(fastaFilePrefix + ".prot.fasta", "w") as protFileOut:
                # Write reference sequences
                refNuclSeq = haplotypesDict["reference_sequence"]
                refProtSeq, _, _ = ZS_SeqIO.FastASeq("id", refNuclSeq).get_translation(findBestFrame=False, strand=1, frame=0)
                
                nuclFileOut.write(f">{refSeqID}\n{refNuclSeq}\n")
                protFileOut.write(f">{refSeqID}\n{refProtSeq}\n")
                
                # Write haplotype sequences
                for i in range(len(haplotypeFrequency)):
                    code, frequency = haplotypeFrequency[i]
                    seqID = f"{geneID}_seq{i+1} haplotypeCode={code} frequency={frequency}"
                    
                    nuclSeq = haplotypesDict["haplotype_sequences"][tuple(map(int,tuple(code)))]
                    protSeq, _, _ = ZS_SeqIO.FastASeq("id", nuclSeq).get_translation(findBestFrame=False, strand=1, frame=0)
                    
                    nuclFileOut.write(f">{seqID}\n{nuclSeq}\n")
                    protFileOut.write(f">{seqID}\n{protSeq}\n")
    
    # Let user know everything went swimmingly
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

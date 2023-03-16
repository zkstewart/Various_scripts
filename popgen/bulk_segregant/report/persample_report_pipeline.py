#! python3
# persample_report_pipeline.py
# Script to combine several input files into a cohesive report file
# to understand the results of a SNP prediction project where the
# goals are to obtain information on segregating populations, but where
# SNP prediction has occurred per-sample.

# Load normal/pip packages
import os, argparse, sys, re
from copy import deepcopy
from goatools import obo_parser

# Load ZS_IO Class code
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))) # 4 dirs up is where we find Function_packages
from Function_packages import ZS_GFF3IO, ZS_SeqIO

# Load functions from other scripts
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))) # 3 dirs up is where we find windows
import windows # pulls relevant functions from snp_proximity_report.py
import haplotypes

# Hard-coded GO replacements
replacedGOs = {"GO:0000453": "GO:0006364", "GO:0031224": "GO:0016020",
               "GO:0031225": "GO:0016020", "GO:0031300": "GO:0031090",
               "GO:0030176": "GO:0005789", "GO:0031227": "GO:0005789",
               "GO:0031301": "GO:0031090"}

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.annotFile):
        print(f'I am unable to locate the GOextended annotation file ({args.annotFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.vcfFile):
        print(f'I am unable to locate the VCF file ({args.vcfFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.gffFile):
        print(f'I am unable to locate the GFF3 file ({args.gffFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.fastaFile):
        print(f'I am unable to locate the FASTA file ({args.fastaFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.goOboFile):
        print(f'I am unable to locate the GO .obo file ({args.goOboFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.metadataFile):
        print(f'I am unable to locate the metadata file ({args.metadataFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print(f'File already exists at output location ({args.outputFileName})')
        print('Make sure you specify a unique file name and try again.')
        quit()

def parse_pops_metadata(metadataFile):
    # Parse file
    metadataDict = {}
    with open(metadataFile, "r") as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n ").split("\t")
            if sl == []:
                continue
            else:
                sample, pop = sl
                metadataDict.setdefault(pop, set())
                metadataDict[pop].add(sample)
    
    # Make sure it's valid for this script
    if len(metadataDict) == 2:
        assert all([ pop in ['bulk1', 'bulk2'] for pop in metadataDict.keys() ]), \
            "ERROR: Metadata file should have 2 or 3 populations: bulk1, bulk2, and (optionally) parent!"
    elif len(metadataDict) == 3:
        assert all([ pop in ['bulk1', 'bulk2', 'parent'] for pop in metadataDict.keys() ]), \
            "ERROR: Metadata file should have 2 or 3 populations: bulk1, bulk2, and (optionally) parent!"
    else:
        print("ERROR: Metadata file should have 2 or 3 populations: bulk1, bulk2, and (optionally) parent!")
        quit()
    
    return metadataDict

def generate_bsa_proximity_dict(gff3Obj, genomeFASTA_obj, snpGenotypes):
    '''
    Reads in a GFF3 and locates SNPs in relation to gene features.
    Unlike the functions in snp_proximity_report.py, this function
    returns a different output which is somewhat of an amalgamation
    of the SNP and gene-centric reports. In reality, if this works
    I'll likely move it over there and just import it here!
    
    Parameters:
        gff3Obj -- a ZS_GFF3IO.GFF3 object with NCLS indexing of the gene features
        snpGenotypes -- a dictionary as produced by get_genotyped_snps_for_genes()
                        with structure like:
                        {
                            'contig1': {
                                pos1: {
                                    'ref_alt': ['nucleotide1', 'nucleotide2'],
                                    'sample1': [allele1, allele2],
                                    'sample2': [...],
                                    ...
                                },
                                pos2: { ... },
                                ...
                            }
                            'contig2': ...,
                            ...
                        }
    Returns:
        proximityDict -- a dictionary with structure like:
                         {
                             'gene1': {
                                 'contig': 'contigID',
                                 'coords': 'start-end',
                                 'within': [ snpData1, snpData2, ... ],
                                 'left_of': [  snpData1, snpData2, ... ],
                                 'right_of': [  snpData1, snpData2, ... ],
                                 'samples': [ sample1_GT, sample2_GT, ...]
                             },
                             'gene2': {
                                 ...
                             },
                             ...
                         }
    '''
    proximityDict = {}
    
    for contig, snpDict in snpGenotypes.items():
        for pos in snpDict.keys():
            # Scenario 1: SNP located within a gene
            matches = gff3Obj.ncls_finder(pos, pos, "contig", contig)
            if matches != []: # i.e., if the SNP is located within a gene
                snpLocation, geneID = windows.locate_snp_within_gene(matches, pos)
                
                # Set default value for gene
                _setdefault_proximity(proximityDict, geneID, gff3Obj)
                
                # For non-CDS locations, just store data in dict
                if snpLocation != "CDS":
                    proximityDict[geneID][pos] = {
                        "location": snpLocation
                    }
                # For CDS locations, predict its effect (if any)
                else:
                    # Get the mRNA feature for this gene
                    mrnaFeature = _get_mrnaFeature_for_cds_match(matches, pos)
                    newSnpDict = haplotypes.convert_vcf_snps_to_cds_snps(mrnaFeature, snpDict, embedOriginalPos=True)
                    cds_FastASeq_obj, cds_featureType, cds_startingFrame = gff3Obj.retrieve_sequence_from_FASTA(genomeFASTA_obj, mrnaFeature.ID, "CDS")
                    
                    # Just use the variant for this SNP position
                    for newPos, genotypeDict in newSnpDict.items():
                        if genotypeDict["originalPos"] == pos:
                            break # break so our variables are set
                    del genotypeDict["originalPos"] # need this to disappear to not break code! ... ugh, this is gross
                    
                    # Get data
                    formattedRef, formattedAlts, homozygotes, heterozygotes \
                        = haplotypes.predict_variant(genotypeDict, newPos, cds_FastASeq_obj, mrnaFeature.strand)
                    
                    # Store data in our dictionary
                    proximityDict[geneID][pos] = {
                        "location": "CDS",
                        "ref": formattedRef,
                        "alts": formattedAlts,
                        "homo": homozygotes,
                        "hetero": heterozygotes
                    }
            
            # Scenario 2: SNP located not inside, but near a gene
            else:
                nearestLeft, nearestRight = windows.locate_genes_near_snp(gff3Obj, contig, pos, mrnaIndexing=False) # NCLS has gene indexing
                if nearestLeft[1] != None:
                    distance, geneID = nearestLeft
                    
                    # Set default value for gene
                    _setdefault_proximity(proximityDict, geneID, gff3Obj)
                    
                    # Store data in dict
                    proximityDict[geneID][pos] = {
                        "location": "right", # gene is left of snp; hence snp is right of gene
                        "distance": distance
                    }
                
                if nearestRight[1] != None:
                    distance, geneID = nearestRight
                    
                    # Set default value for gene
                    _setdefault_proximity(proximityDict, geneID, gff3Obj)
                    
                    # Store data in dict
                    proximityDict[geneID][pos] = {
                        "location": "left", # gene is right of snp; hence snp is left of gene
                        "distance": distance
                    }
    
    return proximityDict

def _setdefault_proximity(proximityDict, geneID, gff3Obj):
    # Get gene feature from GFF3
    feature = gff3Obj[geneID]
    
    # Get representative mRNA feature
    mrnaFeature = ZS_GFF3IO.GFF3.longest_isoform(feature)
    
    # Setdefault in dict
    proximityDict.setdefault(geneID, {
        "contig": feature.contig,
        "strand": feature.strand,
        "mRNA": mrnaFeature.ID,
        "start": feature.start,
        "end": feature.end
    })

def _get_mrnaFeature_for_cds_match(matches, snpPos):
    for feature in matches:
        # Look for CDS overlap
        if hasattr(feature, "mRNA"):
            for mrnaFeature in feature.mRNA:
                for cdsFeature in mrnaFeature.CDS:
                    if cdsFeature.start <= snpPos and snpPos <= cdsFeature.end: # if it overlaps...
                        return mrnaFeature
        elif hasattr(feature, "CDS"):
            for cdsFeature in feature.CDS:
                if cdsFeature.start <= snpPos and snpPos <= cdsFeature.end: # if it overlaps...
                    return feature
    return None

def parse_annot_file(annotFile, idsDict, proximityDict):
    '''
    Simple function to parse the GO extended annotation file and grab some details
    from it.
    '''
    HEADER_VALUES = ["#Query", "Gene_names", "Best_mapped_GOs_+_parents"]
    
    annotDict = {}
    with open(annotFile, "r") as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n ").split("\t")
            
            # Handle header lines
            if line.startswith("#"):
                idIndex = sl.index(HEADER_VALUES[0])
                nameIndex = sl.index(HEADER_VALUES[1])
                goIndex = sl.index(HEADER_VALUES[2])
                continue
            
            # Handle content lines
            else:
                # Grab basic details
                id = sl[idIndex]
                names = sl[nameIndex]
                gos = sl[goIndex]
                
                # Skip if we don't care about this gene
                if id in idsDict:
                    geneID = idsDict[id]
                elif id in proximityDict:
                    geneID = id
                else:
                    continue
                
                # Get best hit name
                name = names.split(" [")[0].rstrip(" ")
                
                # Store it
                annotDict[geneID] = {
                    "name": name,
                    "gos": gos
                }
    return annotDict

def format_snp_types(snpDict):
    '''
    Helper function for formatting SNP data into within/left/right categories.
    
    Parameters:
        snpDict -- a dictionary with structure like:
                   {
                       'contig': 'contig_id',
                       'strand': '+', # or '-'
                       ... ,
                       snpPos1: { 'location': 'left/right/CDS', 'distance': intValue },
                       snpPos2: { ... },
                       ...
                   }
    Returns:
        snpDetails -- a dictionary with structure like:
                      {
                          'within_pos': '', # "; " separated list of values in [sCDS, UTR, T>G, etc.]
                          'within_location': '.', # numeric list of SNP positions
                          'left': '4095701', # numeric list of SNP positions
                          'left_dist': '23043', # numeric list of distances of SNP from nearest gene
                          'right': '.', # numeric list of SNP positions
                          'right_dist': '.' # numeric list of distances of SNP from nearest gene
                      }
    '''
    orderedSNPs = sorted([x for x in snpDict.keys() if isinstance(x, int)])
    within_pos, within_location, left, left_dist, right, right_dist = [], [], [], [], [], []
    
    for snp in orderedSNPs:
        posDict = snpDict[snp]
        # Handle CDS
        if posDict["location"] == "CDS":
            within_pos.append(str(snp))
            
            # Handle neutral variants
            refAmino = posDict["ref"].split("(")[1].rstrip(")")
            altAmino = ",".join([ alt.split("(")[1].rstrip(")") for alt in posDict["alts"] ])
            
            if refAmino == altAmino:
                within_location.append("sCDS")
            # Handle non-synonymous mutations
            else:
                within_location.append(f"{refAmino}>{altAmino}")
        
        # Handle intron/UTR
        elif posDict["location"] in ["UTR", "intron"]:
            within_pos.append(str(snp))
            within_location.append(posDict["location"])
        
        # Handle left/right
        elif posDict["location"] == "left":
            left.append(str(snp))
            left_dist.append(str(posDict["distance"]))
        elif posDict["location"] == "right":
            right.append(str(snp))
            right_dist.append(str(posDict["distance"]))
    
    return {
        "within_pos": "; ".join(within_pos) if within_pos != [] else ".",
        "within_location": "; ".join(within_location) if within_location != [] else ".",
        "left": "; ".join(left) if left != [] else ".",
        "left_dist": "; ".join(left_dist) if left_dist != [] else ".",
        "right": "; ".join(right) if right != [] else ".",
        "right_dist": "; ".join(right_dist) if right_dist != [] else "."
    }

def format_genotypes(snpGenotypes, snpDetails, contigID, metadataDict):
    '''
    Parameters:
        snpGenotypes -- a dictionary as produced by get_genotyped_snps_for_genes()
                        with structure like:
                        {
                            'contig1': {
                                pos1: {
                                    'ref_alt': ['nucleotide1', 'nucleotide2'],
                                    'sample1': [allele1, allele2],
                                    'sample2': [...],
                                    ...
                                },
                                pos2: { ... },
                                ...
                            }
                            'contig2': ...,
                            ...
                        }
        snpDetails -- a dictionary as produced by format_snp_types() with structure like:
                      {
                          'within_pos': '', # "; " separated list of values in [sCDS, UTR, T>G, etc.]
                          'within_location': '.', # numeric list of SNP positions
                          'left': '4095701', # numeric list of SNP positions
                          'left_dist': '23043', # numeric list of distances of SNP from nearest gene
                          'right': '.', # numeric list of SNP positions
                          'right_dist': '.' # numeric list of distances of SNP from nearest gene
                      }
        contigID -- a string indicating which contig the SNP is located on.
        metadataDict -- a dictionary with structure like:
                        {
                            'parent': set(['parent1']),
                            'bulk1': set(['b1sample1', 'b1sample2', ...]),
                            'bulk2': set(['b2sample1', b2sample2', ...])
                        }
    Returns:
        snpGtFormat -- a list with three values, with structure like:
                       [
                           'parental GT per SNP',
                           'bulk1 GT per SNP',
                           'bulk2 GT per SNP'
                       ]
                       Note that the format for each string value is akin to:
                            "snpPosition: 0/0=11,1/1=5,0/1=2"
                       The snpPosition is a numeric value. All present genotypes
                       and their tallies are presented as shown above with comma
                       separation.
    '''
    KEYS_TO_SKIP = ["ref_alt", "originalPos"]
    
    # Get SNP positions from snpDetails dict
    snpNumbers = [
        int(position)
            for key in ["within_pos", "left", "right"]
                for position in snpDetails[key].split("; ")
                    if position != "."
    ]
    
    # Format data pertaining to this SNP
    snpGtDict = {}
    for pos in snpNumbers:
        # Figure out what genotypes exist for each SNP position
        gtDict = snpGenotypes[contigID][pos]
        gts = set([ "/".join(map(str, v)) for k,v in gtDict.items() if not k in KEYS_TO_SKIP ])
        
        # Set up for data storage of this SNP
        snpGtDict[pos] = {}
        for pop in metadataDict.keys():
            snpGtDict[pos][pop] = { gt:0 for gt in gts}
        
        # Tally genotypes for each population
        anyFound = False
        for key, value in gtDict.items():
            if key in KEYS_TO_SKIP:
                continue
            else:
                # Figure out what population this sample belongs to
                pop = [ k for k,v in metadataDict.items() if key in v ]
                
                if len(pop) != 1:
                    continue
                else:
                    anyFound = True
                
                pop = pop[0]
                
                # Figure out and store this genotype
                gt = "/".join(map(str, value))
                snpGtDict[pos][pop][gt] += 1
        if anyFound == False:
            print("ERROR: Issue in format_genotypes; your metadata file doesn't match your VCF")
            print("In other words, we failed to find ANY samples from your metadata file in the VCF")
            print("This is an irreconcilable error and we must exit out now...")
            quit()
    
    # Reformat in sorted order of allele frequency
    """
    This format will directly correspond to something that can be tabulated. Hence,
    it will have 3 columns for each population group (parents, bulk1, bulk2).
    """
    snpGtFormat = [[], [], []]
    for pos, gtDict in snpGtDict.items():
        # Handle parent
        if "parent" in gtDict:
            parentGTs = sorted(gtDict["parent"].items(), key = lambda item: -item[1])
            parentGTs = ",".join([ f"{gt}={count}" for gt,count in parentGTs if count > 0])
        else:
            parentGTs = "."
        
        # Handle bulk1
        bulk1GTs = sorted(gtDict["bulk1"].items(), key = lambda item: -item[1])
        bulk1GTs = ",".join([ f"{gt}={count}" for gt,count in bulk1GTs if count > 0])
        bulk1GTs = "." if bulk1GTs == "" else bulk1GTs
        
        # Handle bulk2
        bulk2GTs = sorted(gtDict["bulk2"].items(), key = lambda item: -item[1])
        bulk2GTs = ",".join([ f"{gt}={count}" for gt,count in bulk2GTs if count > 0])
        bulk2GTs = "." if bulk2GTs == "" else bulk2GTs
        
        # Store into formatting list
        snpGtFormat[0].append(f"{pos}: {parentGTs}")
        snpGtFormat[1].append(f"{pos}: {bulk1GTs}")
        snpGtFormat[2].append(f"{pos}: {bulk2GTs}")
    
    snpGtFormat[0] = "; ".join(snpGtFormat[0])
    snpGtFormat[1] = "; ".join(snpGtFormat[1])
    snpGtFormat[2] = "; ".join(snpGtFormat[2])
    
    return snpGtFormat

def format_difference_ratio(snpGtFormat):
    '''
    Parameters:
        snpGtFormat -- a list with three values, with structure like:
                       [
                           'parental GT per SNP',
                           'bulk1 GT per SNP',
                           'bulk2 GT per SNP'
                       ]
                       Note that the format for each string value is akin to:
                            "snpPosition: 0/0=11,1/1=5,0/1=2"
                       The snpPosition is a numeric value. All present genotypes
                       and their tallies are presented as shown above with comma
                       separation.
    Returns:
        differenceRatio -- a string representing a ratio from 0->1, where 1 indicates
                           entirely different allele composition, and 0 indicates
                           identical alleles in the two populations. In short, close to 1
                           gives us interesting regions, and close to 0 should be
                           ignored when doing a BSA-like inspection.
    '''
    GT_REGEX = re.compile(r"\d/\d")
    
    # Extract relevant bits of data from snpGtFormat
    b1 = snpGtFormat[1] # we don't care about the parent for this index
    b2 = snpGtFormat[2]
    
    b1Alleles = b1.split("; ") # can be multiple SNP positions separated by "; "
    b2Alleles = b2.split("; ")
    
    # See if we have any multiallelic variant
    genotypes = set(GT_REGEX.findall(b1)).union(set(GT_REGEX.findall(b2)))
    alleles = set([ a for gt in genotypes for a in gt.split("/") ])
    thisAlleles = {
        "b1": [0 for _ in range(len(alleles))],
        "b2": [0 for _ in range(len(alleles))]
    }
    
    # Handle identical allele compositions
    "If everything is 0/0 or 1/1, we know the ratio == 0.0; stop here to prevent bugs"
    if len(alleles) == 1:
        return "0.0"
    
    # Tally alleles at this locus
    for x in range(len(b1Alleles)):
        # Get alleles at this position
        b1x = b1Alleles[x]
        b2x = b2Alleles[x]
        
        # Skip tallying if it's not possible at this locus
        if "." in b1x or "." in b2x:
            continue
        
        # Tally for bulk 1
        for b1xAllele in b1x.split(","):
            b1xAllele = b1xAllele.split(": ")[-1]
            gt, count = b1xAllele.split("=")
            for gtValue in gt.split("/"):
                thisAlleles["b1"][int(gtValue)] += int(count)
        
        # Tally for bulk 2
        for b2xAllele in b2x.split(","):
            b2xAllele = b2xAllele.split(": ")[-1]
            gt, count = b2xAllele.split("=")
            for gtValue in gt.split("/"):
                thisAlleles["b2"][int(gtValue)] += int(count)
    
    # Calculate the proportions of each allele for this locus
    b1Sum = sum(thisAlleles["b1"])
    b2Sum = sum(thisAlleles["b2"])
    
    if b1Sum == 0 or b2Sum == 0: # if this happens, we can't meaningfully calculate anything
        return "."
    
    thisProportion = {
        "b1": [ alleleCount / b1Sum for alleleCount in thisAlleles["b1"] ],
        "b2": [ alleleCount / b2Sum for alleleCount in thisAlleles["b2"] ]
    }
    
    # Derive our difference ratio value
    proportionCommon = sum([
        min(thisProportion["b1"][x], thisProportion["b2"][x])
        for x in range(len(thisProportion["b1"]))
    ])
    differenceRatio = 1 - proportionCommon
    
    return str(differenceRatio)

def main():
    # User input
    usage = """%(prog)s receives several files associated with a per-sample SNP
    prediction project dealing with two segregating populations. It combines relevant
    information into a cohesive and comprehensive report TSV file.
    
    Required input files are:
    1) GOextended annotation table from annotation_table_extend_GOs.py.
    2) VCF containing variant calls for multiple samples belonging to two segregating phenotype
    groups (optionally, can include parents).
    3) Gene annotation GFF3 file.
    4) Genome FASTA file.
    5) GO .obo file.
    6) Metadata file indicating the groups each sample belongs to.
    
    Alongside these, you should specify several parameters to control the behaviour
    of this script.
    """
    
    # Required
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-annot", dest="annotFile",
                   required=True,
                   help="Specify the GOextended annotation file name")
    p.add_argument("-vcf", dest="vcfFile",
                   required=True,
                   help="Specify the VCF file name")
    p.add_argument("-gff", dest="gffFile",
                   required=True,
                   help="Specify the genome annotation GFF3 file name")
    p.add_argument("-fasta", dest="fastaFile",
                   required=True,
                   help="Specify the genome FASTA file name")
    p.add_argument("-obo", dest="goOboFile",
                   required=True,
                   help="Specify the GO .obo file name")
    p.add_argument("-meta", dest="metadataFile",
                   required=True,
                   help="Specify the metadata file name")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Specify the output file name")
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse GO.obo
    go = obo_parser.GODag(args.goOboFile)
    
    # Parse metadata file
    metadataDict = parse_pops_metadata(args.metadataFile)
    
    # Parse GFF3 with NCLS indexing
    gff3Obj = ZS_GFF3IO.GFF3(args.gffFile, strict_parse=False)
    gff3Obj.create_ncls_index(typeToIndex="gene")
    
    # Parse FASTA file
    genomeFASTA_obj = ZS_SeqIO.FASTA(args.fastaFile)
    
    # Parse VCF data for outlier SNP genotypes
    snpGenotypes = haplotypes.get_genotypes_from_vcf(args.vcfFile, snpPositions=None, imputeMissing=False)
        
    # Get the proximity reporting dict (including effects prediction)
    proximityDict = generate_bsa_proximity_dict(gff3Obj, genomeFASTA_obj, snpGenotypes)
    
    # Parse annotation file for relevant entries
    mrnaDict = { proximityDict[key]["mRNA"]: key for key in proximityDict.keys() }
    annotDict = parse_annot_file(args.annotFile, mrnaDict, proximityDict)
    
    # Combine everything into a table output
    needsReplace = set()
    with open(args.outputFileName, "w") as fileOut:
        # Write header line
        fileOut.write("{0}\n".format("\t".join([
            "#contig", "gene_ID", "mRNA_ID", "strand",
            "coords", "gene_name", "GO_IDs", "GO_names",
            "within_pos", "within_location", "left_pos",
            "left_distance", "right_pos", "right_distance",
            "parentGT", "bulk1GT", "bulk2GT", "difference_ratio"
            #"desirable_delta"
        ])))
        
        # Write content lines
        for geneID, snpDict in proximityDict.items():
            annot = annotDict[geneID]
            
            # Get GO names from IDs
            if annot["gos"] == ".":
                goNames = "."
            else:
                goNames = []
                for term in annot["gos"].split("; "):
                    if term in replacedGOs:
                        goNames.append(go.get(replacedGOs[term]).name)
                    else:
                        try:
                            goNames.append(go.get(term).name)
                        except:
                            needsReplace.add(term)
                goNames = "; ".join(goNames)
            
            # Format SNP details
            snpDetails = format_snp_types(snpDict)
            
            # Format genotype info
            snpGtFormat = format_genotypes(snpGenotypes, snpDetails, snpDict["contig"], metadataDict)
            
            # Calculate difference ratio
            differenceRatio = format_difference_ratio(snpGtFormat)
            
            # Format output line
            outputLine = "{contig}\t{geneID}\t{mrnaID}\t{strand}\t{coords}\t{geneName}\
\t{gos}\t{goNames}\t{within_pos}\t{within_location}\t{left}\
\t{left_dist}\t{right}\t{right_dist}\
\t{gtFormat}\t{differenceRatio}\n".format(
                contig = snpDict["contig"],
                geneID = geneID,
                mrnaID = snpDict["mRNA"],
                strand = snpDict["strand"],
                coords = "{0}-{1}".format(snpDict["start"], snpDict["end"]),
                geneName = annot["name"],
                gos = annot["gos"],
                goNames = goNames,
                within_pos = snpDetails["within_pos"],
                within_location = snpDetails["within_location"],
                left = snpDetails["left"],
                left_dist = snpDetails["left_dist"],
                right = snpDetails["right"],
                right_dist = snpDetails["right_dist"],
                gtFormat = "\t".join(snpGtFormat),
                differenceRatio = differenceRatio
            )
            
            # Write output row as line
            fileOut.write(outputLine)
    
    if len(needsReplace) != 0:
        print("Some GOs need replacing for this to work correctly!")
        print(f"These are: {needsReplace}")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

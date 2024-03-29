#! python3
# bsa_report_pipeline.py
# Script to combine several input files into a cohesive report file
# to understand the results of a bulked segregant analysis.

# Load normal/pip packages
import os, argparse, sys
from goatools import obo_parser

# Load ZS_IO Class code
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))) # 4 dirs up is where we find Function_packages
from Function_packages import ZS_GFF3IO, ZS_SeqIO

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))) # 3 dirs up is where we find windows
import windows # pulls relevant functions from snp_proximity_report.py
import haplotypes

# Hard-coded GO replacements
replacedGOs = {"GO:0000453": "GO:0006364", "GO:0031224": "GO:0016020",
               "GO:0031225": "GO:0016020", "GO:0031300": "GO:0031090",
               "GO:0030176": "GO:0005789", "GO:0031227": "GO:0005789",
               "GO:0031301": "GO:0031090", 'GO:0055076': 'GO:0030003',
               'GO:0006875': 'GO:0030003', 'GO:0046916': 'GO:0030003',
               'GO:0055067': 'GO:0055080', 'GO:0030004': 'GO:0030003',
               'GO:0005779': 'GO:0005778', 'GO:0032592': 'GO:0031966',
               'GO:0098573': 'GO:0031966', 'GO:0031231': 'GO:0005778',
               'GO:0055065': 'GO:0030003', 'GO:0018196': 'GO:0018193',
               'GO:0006471': 'GO:1990404', 'GO:0008022': 'GO:0005515',
               'GO:0071458': 'GO:0098554', 'GO:0031304': 'GO:0005743',
               'GO:0071556': 'GO:0098553', 'GO:0140323': 'GO:0022853',
               'GO:0015299': 'GO:0015078', 'GO:0015298': 'GO:0022853',
               'GO:0044728': 'GO:0006304', 'GO:0031226': 'GO:0005886',
               'GO:0046658': 'GO:0005886', 'GO:0031305': 'GO:0005743',
               'GO:0005451': 'GO:0022890'}



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
    if not os.path.isfile(args.snpFile):
        print(f'I am unable to locate the SNP index-like file ({args.snpFile})')
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
    # Validate optional input file locations
    if args.scanoneFile != None:
        if not os.path.isfile(args.scanoneFile):
            print(f'I am unable to locate the R/qtl scanone file ({args.scanoneFile})')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    if args.scantwoFile != None:
        if not os.path.isfile(args.scantwoFile):
            print(f'I am unable to locate the R/qtl scantwo file ({args.scantwoFile})')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    # Validate column header parameter
    if len(args.snpColumns) != 3:
        print("snpColumns should provide three values corresponding to header labels")
        print(f"You provided '{args.snpColumns}', which contains {len(args.snpColumns)} values")
        quit()
    # Validate numeric arguments
    if args.kbRadius < 0:
        print(f"kbRadius should be an integer greater than 0")
        quit()
    if args.pvalue < 0:
        print(f"pvalue should be an integer greater than 0")
        quit()
    if args.lod < 0:
        print(f"lod should be an integer greater than 0")
        quit()
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print(f'File already exists at output location ({args.outputFileName})')
        print('Make sure you specify a unique file name and try again.')
        quit()

def validate_snpGenotypes(snpGenotypes):
    '''
    This function exists to detect problems in the SNP parsing early into the process.
    I'm making it because I found problems slipped through the cracks in my own development.
    Luckily I noticed it, but it was otherwise creating a silent problem for the code.
    
    ACTUALLY - this is dumb and doesn't work except for QTL-seq files. Freebayes per-sample
    prediction may have many more than 4 samples present. We're just gonna have to allow
    the silent error to slip through; it is visible in the output file.
    '''
    for contigID, posDict in snpGenotypes.items():
        for pos, gtDict in posDict.items():
            assert len(gtDict) == 4, \
                "validate_snpGenotypes: detected problem with snpGenotypes; your SNP index is bad??"

def parse_pops_metadata(metadataFile):
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
    return metadataDict


def parse_scanone_file(scanoneFile, pvalueFilter=0.05, lodFilter=10.0):
    '''
    Expects a TSV format produced by me (ZKS) from R/qtl. The file has a single
    header row with split by tab looking like:
    
    ["", "chr", "pos", "lod", "pval"]
    
    The first column itself is formatted like:
    
        S{chr}_{bp}
    
    Where {chr} == the value in the next column over, and {bp} is the numeric position
    in the chromosome measured in basepairs.
    
    If the file hasn't had a P-value filter enforced, this function will do that here.
    
    Parameters:
        scanoneFile -- a string indicating the file location of the scanone output
                       from R/qtl.
        pvalueFilter -- a float indicating the P-value significance to enforce.
    Returns:
        scanoneDict -- a dictionary with structure like:
            {
                "chr1": [bp1, bp2] # bp == nucleotide bp as int value
            }
    '''
    scanoneDict = {}
    
    firstLine = True
    with open(scanoneFile, "r") as fileIn:
        for line in fileIn:
            if firstLine == True:
                firstLine = False
                continue
            else:
                # Get relevant information from this line
                chr_bp, chr, pos, lod, pval = line.rstrip("\r\n ").split("\t")
                bp = int(chr_bp.split("_")[-1])
                lod = float(lod)
                pval = float(pval)
                
                # Skip if pval or lod don't pass cutoffs
                if pval > pvalueFilter or lod < lodFilter:
                    continue
                
                # Store results otherwise
                scanoneDict.setdefault(chr, [])
                scanoneDict[chr].append(bp)
    
    return scanoneDict

def parse_scantwo_file(scantwoFile, pvalueFilter=0.05, lodFilter=10.0):
    '''
    Expects a TSV format produced by me (ZKS) from R/qtl. The file has a single
    header row with split by tab looking like:
    
    ["", "chr1", "chr2", "pos1f", "pos2f", "lod.full", "pval", ... "chr1_pos", "chr2_pos"]
    (intervening columns don't matter to us)
    
    If the file hasn't had a P-value filter enforced, this function will do that here.
    
    Parameters:
        scantwoFile -- a string indicating the file location of the scantwo output
                       from R/qtl.
        pvalueFilter -- a float indicating the P-value significance to enforce.
    Returns:
        scantwoDict -- a dictionary with structure like:
            {
                n: [ # == integer indexing this pair of QTL values
                    "chr1", "pos1", "chr2", "pos2"
                ],
                ...
            }
    '''
    scantwoDict = {}
    
    firstLine = True
    ongoingCount = 0
    with open(scantwoFile, "r") as fileIn:
        for line in fileIn:
            if firstLine == True:
                firstLine = False
                continue
            else:
                # Get relevant information from this line
                sl = line.rstrip("\r\n ").split("\t")
                chr1, chr2, lod, pval, pos1, pos2 = sl[1], sl[2], sl[5], sl[6], sl[-2], sl[-1]
                pos1, pos2 = int(pos1), int(pos2)
                lod = float(lod)
                pval = float(pval)
                
                # Skip if pval or lod don't pass cutoffs
                if pval > pvalueFilter or lod < lodFilter:
                    continue
                
                # Store results otherwise
                scantwoDict[ongoingCount] = [
                    chr1, pos1, chr2, pos2
                ]
                ongoingCount += 1
    
    return scantwoDict

def get_snp_positions_from_indexlike(snpFile, snpColumns):
    '''
    Parameters:
        snpFile -- a string indicating the location of the SNP index-like
                   file from which we want to parse contig>position pairs.
        snpColumns -- a list containing two values indicating the header
                      values for the relevant columns we can extract
                      the pairs from. Should be contig first and position
                      second in the list.
    Returns:
        snpPositions -- a dictionary with structure like:
                        {
                            "contig1": set([
                                pos1, pos2, pos3, ...
                            ]),
                            "contig2": set([ ... ]),
                            ...
                        }
    '''
    assert len(snpColumns) == 2, \
        "get_snp_positions_from_indexlike needs a snpColumns list with only 2 values!"
    
    snpPositions = {}
    with open(snpFile, "r") as fileIn:
        firstLine = True
        for line in fileIn:
            sl = line.rstrip("\r\n ").split("\t")
            # Handle header line
            if firstLine == True:
                contigIndex = sl.index(snpColumns[0])
                posIndex = sl.index(snpColumns[1])
                firstLine = False
            # Handle content lines
            else:
                contig, pos = sl[contigIndex], sl[posIndex]
                snpPositions.setdefault(contig, set())
                snpPositions[contig].add(int(pos))
    return snpPositions

def get_delta_snpindex_from_indexlike(snpFile, snpColumns):
    '''
    Parameters:
        snpFile -- a string indicating the location of the SNP index-like
                   file from which we want to parse contig>position pairs.
        snpColumns -- a list containing three values indicating the header
                      values for the relevant columns we can extract
                      the pairs from. Should be contig first, position
                      second, and delta SNP index third in the list.
    Returns:
        deltaSnpIndex -- a dictionary with structure like:
                        {
                            "contig1": {
                                pos1: delta_SNPindex1,
                                pos2: delta_SNPindex2,
                                pos3: ...,
                                ...
                            }
                            "contig2": {
                                ...
                            },
                            ...
                        }
    '''
    assert len(snpColumns) == 3, \
        "get_delta_snpindex_from_indexlike needs a snpColumns list with only 2 values!"
    
    deltaSnpIndex = {}
    with open(snpFile, "r") as fileIn:
        firstLine = True
        for line in fileIn:
            sl = line.rstrip("\r\n ").split("\t")
            # Handle header line
            if firstLine == True:
                contigIndex = sl.index(snpColumns[0])
                posIndex = sl.index(snpColumns[1])
                deltaIndex = sl.index(snpColumns[2])
                firstLine = False
            # Handle content lines
            else:
                contig, pos, delta = sl[contigIndex], sl[posIndex], sl[deltaIndex]
                deltaSnpIndex.setdefault(contig, {})
                deltaSnpIndex[contig][int(pos)] = float(delta)
    return deltaSnpIndex

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

def parse_annot_file(annotFile, idsDict):
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
                if not id in idsDict:
                    continue
                geneID = idsDict[id]
                
                # Get best hit name
                name = names.split(" [")[0].rstrip(" ")
                
                # Store it
                annotDict[geneID] = {
                    "name": name,
                    "gos": gos
                }
    return annotDict

def convert_scanone_to_windows(scanoneDict, windowSize):
    '''
    Parameters:
        scanoneDict -- a dictionary with structure like:
                       {
                           'contig1': [ qtl1, qtl2, ... ],
                           'contig2': [ qtl1, qtl2, ... ],
                           ...
                       }
        windowSize -- an integer representing the value in kb
                      to look to the left and right of a qtl
                      when forming a window
    '''
    scanoneWindows = {}
    for contig, qtls in scanoneDict.items():
        contig = contig.lower()
        
        scanoneWindows.setdefault(contig, [])
        for qtl in qtls:
            start = qtl - (windowSize*1000)
            end = qtl + (windowSize*1000)
            scanoneWindows[contig].append(range(start, end+1))
    return scanoneWindows

def convert_scantwo_to_windows(scantwoDict, windowSize):
    '''
    Parameters:
        scantwoDict -- a dictionary with structure like:
                       {
                           qtlIndex1: [ contig1_id, contig1_pos, contig2_id, contig2_pos ],
                           qtlIndex2: [ ... ],
                           ...
                       }
        windowSize -- an integer representing the value in kb
                      to look to the left and right of a qtl
                      when forming a window
    '''
    scantwoWindows = {}
    for contig1_id, contig1_pos, contig2_id, contig2_pos in scantwoDict.values():
        contig1_id = contig1_id.lower()
        contig2_id = contig2_id.lower()
        
        scantwoWindows.setdefault(contig1_id, [])
        scantwoWindows.setdefault(contig2_id, [])
        
        qtl1Start = contig1_pos - (windowSize*1000)
        qtl2Start = contig2_pos - (windowSize*1000)
        
        qtl1End = contig1_pos + (windowSize*1000)
        qtl2End = contig2_pos + (windowSize*1000)
        
        if range(qtl1Start, qtl1End+1) not in scantwoWindows[contig1_id]:
            scantwoWindows[contig1_id].append(range(qtl1Start, qtl1End+1))
        if range(qtl2Start, qtl2End+1) not in scantwoWindows[contig2_id]:
            scantwoWindows[contig2_id].append(range(qtl2Start, qtl2End+1))
    return scantwoWindows

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
        for key, value in gtDict.items():
            if key in KEYS_TO_SKIP:
                continue
            else:
                # Figure out what population this sample belongs to
                pop = [ k for k,v in metadataDict.items() if key in v ]
                assert len(pop) == 1, \
                    "Issue in format_genotypes; your metadata file doesn't match your VCF"
                pop = pop[0]
                
                # Figure out and store this genotype
                gt = "/".join(map(str, value))
                snpGtDict[pos][pop][gt] += 1
    
    # Reformat in sorted order of allele frequency
    """
    This format will directly correspond to something that can be tabulated. Hence,
    it will have 3 columns for each population group (parents, bulk1, bulk2).
    """
    snpGtFormat = [[], [], []]
    for pos, gtDict in snpGtDict.items():
        # Handle parent
        parentGTs = sorted(gtDict["parent"].items(), key = lambda item: -item[1])
        parentGTs = ",".join([ f"{gt}={count}" for gt,count in parentGTs if count > 0])
        
        # Handle bulk1
        bulk1GTs = sorted(gtDict["bulk1"].items(), key = lambda item: -item[1])
        bulk1GTs = ",".join([ f"{gt}={count}" for gt,count in bulk1GTs if count > 0])
        
        # Handle bulk2
        bulk2GTs = sorted(gtDict["bulk2"].items(), key = lambda item: -item[1])
        bulk2GTs = ",".join([ f"{gt}={count}" for gt,count in bulk2GTs if count > 0])
        
        # Store into formatting list
        snpGtFormat[0].append(f"{pos}: {parentGTs}")
        snpGtFormat[1].append(f"{pos}: {bulk1GTs}")
        snpGtFormat[2].append(f"{pos}: {bulk2GTs}")
    
    snpGtFormat[0] = "; ".join(snpGtFormat[0])
    snpGtFormat[1] = "; ".join(snpGtFormat[1])
    snpGtFormat[2] = "; ".join(snpGtFormat[2])
    
    return snpGtFormat

def format_delta(deltaSnpIndex, snpDetails, contigID):
    '''
    Parameters:
        deltaSnpIndex -- a dictionary as produced by get_delta_snpindex_from_indexlike()
                        with structure like:
                        {
                            'contig1': {
                                pos1: delta1,
                                pos2: delta2,
                                ...
                            },
                            'contig2': {
                                pos1: delta1,
                                ...
                            },
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
    snpDeltaFormat = []
    for pos in snpNumbers:
        delta = deltaSnpIndex[contigID][pos]
        snpDeltaFormat.append(f"{pos}: {delta}")
    
    snpDeltaFormat = "; ".join(snpDeltaFormat)
    return snpDeltaFormat

def determine_desirable_delta(snpGtFormat, parentValue):
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
        parentValue -- a string equal to "high" if the parent's phenotype is high
                       value, or "low" for the opposite
    Returns:
        desirableDelta -- a string pairing the SNP positions to a possible value
                          in the range of:
                          [ "1", "-1", "." ]
                          "." represents uncertainty which occurs when the
                          parental genotype is a heterozygote or has ambiguity
                          itself in terms of its genotype.
                          
                          For multiple SNPs, it will provide these values with "; "
                          separation akin to:
                          "snpPosition1: -1; snpPositions2: 1; ..."
    '''
    assert parentValue.lower() in ["high", "low"], \
        f"determine_desirable_delta: '{parentValue}' must be equal to 'high' or 'low'!"
    
    parentalSnpGT = snpGtFormat[0]
    desirableDelta = []
    for snpGT in parentalSnpGT.split("; "):
        position, gtCount = snpGT.split(": ")
        gt = gtCount.split("=")[0]
        gtSet = set(gt.replace("|", "/").split("/")) # replace just for safety purposes
        # Handle ambiguity
        if "." in gt:
            desirableDelta.append(".")
        # Handle homozygotes [note that we're assuming biallelic variants]
        elif len(gtSet) == 1:
            if parentValue.lower() == "high":
                "We want the genotype to be EQUAL to the parent"
                if "0" in gtSet: # if the desirable GT is the reference allele...
                    desirableDelta.append(f"{position}: -1") # the good bulk SNP index == 0, and bad == 1, hence delta == -1
                else:
                    desirableDelta.append(f"{position}: 1") # the good bulk SNP index == 1, and bad == 0, hence delta == 1
            else:
                "We want the genotype to be OPPOSITE to the parent"
                if "0" in gtSet: # if the desirable GT is the reference allele...
                    desirableDelta.append(f"{position}: 1") # the good bulk SNP index == 1, and bad == 0, hence delta == 1
                else:
                    desirableDelta.append(f"{position}: -1") # the good bulk SNP index == 0, and bad == 1, hence delta == -1
        # Handle heterozygotes
        else:
            desirableDelta.append(f"{position}: .")
    
    return "; ".join(desirableDelta)

def has_desirable_delta(deltaFormat, desirableDelta):
    '''
    Returns a boolean indicating whether there's any correspondance between
    our actual deltas, and our desired deltas.
    
    Note: This method isn't planned for use since QTL-seq doesn't appear to
    calculate delta SNP index according to the formula they use in their
    own paper. I'd have to dig into things to figure out what's going on...
    
    Parameters:
        deltaFormat -- simialr to the below, but instead of having "1", "-1", or ".",
                       just presents the actual delta value.
        desirableDelta -- a string pairing the SNP positions to a possible value
                          in the range of:
                          [ "1", "-1", "." ]
                          "." represents uncertainty which occurs when the
                          parental genotype is a heterozygote or has ambiguity
                          itself in terms of its genotype.
                          
                          For multiple SNPs, it will provide these values with "; "
                          separation akin to:
                          "snpPosition1: -1; snpPositions2: 1; ..."
    '''
    deltas = deltaFormat.split("; ")
    desireds = desirableDelta.split("; ")
    
    return NotImplementedError()

def main():
    # User input
    usage = """%(prog)s receives several files associated with a BSA
    project and combines that information into a cohesive and comprehensive
    report TSV file.
    
    Required input files are:
    1) GOextended annotation table from annotation_table_extend_GOs.py.
    2) VCF containing variant calls that BSA was based upon.
    3) SNP index-like file containing at least two columns linking chromosome
    and position for variants detected as being outliers.
    4) Gene annotation GFF3 file.
    5) Genome FASTA file.
    6) GO .obo file.
    
    Optional input files are:
    1) R/qtl scanone file.
    2) R/qtl scantwo file.
    
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
    p.add_argument("-snp", dest="snpFile",
                   required=True,
                   help="Specify the SNP index-like file name")
    p.add_argument("--snpColumns", dest="snpColumns", nargs="+",
                   required=False,
                   help="""Optionally, specify the header values of the three columns
                   to parse from the SNP index-like file
                   (default == 'CHROM POSI delta_SNPindex' """,
                   default=["CHROM", "POSI", "delta_SNPindex"])
    p.add_argument("-obo", dest="goOboFile",
                   required=True,
                   help="Specify the GO .obo file name")
    p.add_argument("-meta", dest="metadataFile",
                   required=True,
                   help="Specify the metadata file name")
    # p.add_argument("-pheno", dest="parentalPhenotype",
    #                required=True,
    #                choices=["high", "low"],
    #                help="""For the PARENT, specify whether it has a high or low
    #                value phenotype""")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Specify the output file name")
    # Optional
    p.add_argument("--s1", dest="scanoneFile",
                   required=False,
                   help="Input scanone file produced by Rqtl")
    p.add_argument("--s2", dest="scantwoFile",
                   required=False,
                   help="Input scantwo file produced by Rqtl")
    # Behavioural parameters
    p.add_argument("--kbRadius", dest="kbRadius", type=int,
                   required=False,
                   help="""Optionally, specify the radius of the window in kb
                   surrounding an R/qtl discovery to mark as our window
                   region(s) (default==50)""",
                   default=50)
    p.add_argument("--pvalue", dest="pvalue", type=float,
                   required=False,
                   help="""Optionally, specify the P-value to enforce for
                   R/qtl discovery (default==0.05)""",
                   default=0.05)
    p.add_argument("--lod", dest="lod", type=float,
                   required=False,
                   help="""Optionally, specify the minimum LOD to accept from
                   a R/qtl result (default==10.0)""",
                   default=10)
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse GO.obo
    go = obo_parser.GODag(args.goOboFile)
    
    # Parse metadata file
    metadataDict = parse_pops_metadata(args.metadataFile)
    assert set(metadataDict.keys()) == {'bulk1', 'bulk2', 'parent'}, \
        "Metadata file should have 3 populations: bulk1, bulk2, and parent!"
    
    # Parse GFF3 with NCLS indexing
    gff3Obj = ZS_GFF3IO.GFF3(args.gffFile, strict_parse=False)
    gff3Obj.create_ncls_index(typeToIndex="gene")
    
    # Parse FASTA file
    genomeFASTA_obj = ZS_SeqIO.FASTA(args.fastaFile)
    
    # Parse SNP index-like file for outlier SNPs & delta SNP index values
    snpPositions = get_snp_positions_from_indexlike(args.snpFile, args.snpColumns[0:2])
    deltaSnpIndex = get_delta_snpindex_from_indexlike(args.snpFile, args.snpColumns)
    
    # Parse VCF data for outlier SNP genotypes
    snpGenotypes = haplotypes.get_genotypes_from_vcf(args.vcfFile, snpPositions, imputeMissing=False)
    # validate_snpGenotypes(snpGenotypes) # raise a flag if something looks amiss
    
    # Get the proximity reporting dict (including effects prediction)
    proximityDict = generate_bsa_proximity_dict(gff3Obj, genomeFASTA_obj, snpGenotypes)
    
    # Parse the scanone file (if relevant)
    if args.scanoneFile != None:
        scanoneDict = parse_scanone_file(args.scanoneFile, args.pvalue, args.lod)
        scanoneWindows = convert_scanone_to_windows(scanoneDict, args.kbRadius)
    else:
        scanoneWindows = {}
    
    # Parse the scantwo file (if relevant)
    if args.scantwoFile != None:
        scantwoDict = parse_scantwo_file(args.scantwoFile, args.pvalue, args.lod)
        scantwoWindows = convert_scantwo_to_windows(scantwoDict, args.kbRadius)
    else:
        scantwoWindows = {}
    
    # Parse annotation file for relevant entries
    mrnaDict = { proximityDict[key]["mRNA"]: key for key in proximityDict.keys() }
    annotDict = parse_annot_file(args.annotFile, mrnaDict)
    
    # Combine everything into a table output
    needsReplace = set()
    with open(args.outputFileName, "w") as fileOut:
        # Write header line
        fileOut.write("{0}\n".format("\t".join([
            "#contig", "gene_ID", "mRNA_ID", "strand",
            "coords", "gene_name", "GO_IDs", "GO_names",
            "within_pos", "within_location", "left_pos",
            "left_distance", "right_pos", "right_distance",
            "scan1Window", "scan2Window", "parentGT",
            "bulk1GT", "bulk2GT", "delta_SNPindex"#,
            #"desirable_delta"
        ])))
        
        # Write content lines
        for geneID, snpDict in proximityDict.items():
            annot = annotDict[geneID]
            
            # Check if this is in a window region
            try:
                for windowRange in scanoneWindows[snpDict["contig"]]:
                    if snpDict["start"] in windowRange or snpDict["end"] in windowRange:
                        isScanoneWindow = True
                    else:
                        isScanoneWindow = False
            except:
                isScanoneWindow = False
            
            try:
                for windowRange in scantwoWindows[snpDict["contig"]]:
                    if snpDict["start"] in windowRange or snpDict["end"] in windowRange:
                        isScantwoWindow = True
                    else:
                        isScantwoWindow = False
            except:
                isScantwoWindow = False
            
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
            
            # Format delta SNP index info
            deltaFormat = format_delta(deltaSnpIndex, snpDetails, snpDict["contig"])
            
            # Figure out what extreme of SNP index is desirable
            """Depending on whether the phenotype of the parent is high or low value,
            we'll want a +ve or -ve delta value.
            
            For example, let's say the parent has 0/0 genotype and is low value.
            In that instance, a 1/1 genotype would be desirable. Hence, we'd want to
            find a SNP with a +1 delta (high value bulk SNP index will be 1 since it
            will have all alternative alleles; low value bulk SNP index will be 0 since
            it will have all reference alleles).
            
            Alternatively, if the low value parent has 1/1 genotype, we'd want a 0/0
            genotype. In that case, we'd prefer a SNP with -1 delta (high value bulk SNP
            index == 0 since it has all reference alleles; low value bulk SNP index
            == 1 since it has all alternate alleles).
            
            This realisation makes me see that the current delta SNP index is non-intuitive
            and could be calculated better as:
                SNP-index (at a position) = Count of HIGH-VALUE base/ Count of reads aligned.
            That'd mean SNP-index == 1 will always be good. Oh well.
            """
            #desirableDelta = determine_desirable_delta(snpGtFormat, args.parentalPhenotype)
            
            # Provide a quick overview of whether a gene contains any desirable deltas
            ## Not done, won't be implemented unless I figure out how QTL-seq calculates deltas
            
            # Format output line
            outputLine = "{contig}\t{geneID}\t{mrnaID}\t{strand}\t{coords}\t{geneName}\
\t{gos}\t{goNames}\t{within_pos}\t{within_location}\t{left}\
\t{left_dist}\t{right}\t{right_dist}\t{scan1}\t{scan2}\
\t{gtFormat}\t{delta}\n".format(
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
                scan1 = "Y" if isScanoneWindow else ".",
                scan2 = "Y" if isScantwoWindow else ".",
                gtFormat = "\t".join(snpGtFormat),
                delta = deltaFormat#,
                #desirableDelta = desirableDelta
            )
            
            # Write output row as line
            fileOut.write(outputLine)
    
    if len(needsReplace) != 0:
        print("Some GOs need replacing for this to work correctly!")
        print(f"These are: {needsReplace}")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

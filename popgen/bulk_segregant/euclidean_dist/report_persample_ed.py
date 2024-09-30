#! python3
# report_persample_ed.py
# Script to combine several input files into a cohesive report file
# to understand the results of a SNP prediction project where the
# goals are to obtain information on segregating populations, but where
# SNP prediction has occurred per-sample. It builds upon the
# calculate_persample_ed.py script which may have been filtered.

# Load normal/pip packages
import os, argparse, sys
from goatools import obo_parser

# Load ZS_IO Class code
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))) # 4 dirs up is where we find Function_packages
from Function_packages import ZS_GFF3IO, ZS_GO

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))) # 3 dirs up is where we find windows
import windows # pulls relevant functions from snp_proximity_report.py
import haplotypes

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.edistFile):
        print(f'I am unable to locate the Euclideand distance file ({args.edistFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.annotFile):
        print(f'I am unable to locate the annotation file ({args.annotFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.gff3File):
        print(f'I am unable to locate the GFF3 file ({args.gff3File})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.goOboFile):
        print(f'I am unable to locate the GO .obo file ({args.goOboFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
    # Validate numeric arguments
    if 0 > args.euclideanDistance:
        print("--ed must be 0 (no filtering) or higher")
        quit()
    if args.radiusSize < 1:
        print("--radiusSize must be a positive integer")
        quit()
    # Validate output file location
    for outputSuffix in [".gene_report.tsv", ".variant_report.tsv"]:
        if os.path.isfile(args.outputPrefix + outputSuffix):
            print(f'File already exists at output location ({args.outputPrefix + outputSuffix})')
            print('Make sure you specify a unique file name and try again.')
            quit()

def parse_ed_file(tsvFile, bulkAlleles=[], bulkOccurrence=None,
                  HEADER_VALUES = ["CHROM", "POSI", "euclideanDist", "bulk1_alleles", "bulk2_alleles"]):
    '''
    Parameters:
        tsvFile -- a string pointing to the TSV file containing relevant statistics
        bulkAlleles -- OPTIONAL; a list of two integers indicating the maximum number of alleles
                       in each bulk OR an empty list OR None to indicate that bulk occurrence
                       should not be considered (default=[])
        bulkOccurrence -- OPTIONAL; a float value indicating the minimum fraction of occurrence
                          for one of the two bulks to be considered for plotting OR None to indicate
                          that bulk occurrence should not be considered (default=None)
    Returns:
        dotsX -- a dict pairing chromosome IDs (keys) to a list of integers indicating
                 the position for each SNP (values)
        dotsY -- a dict pairing chromosome IDs (keys) to a list of floats indicating the
                 statistic for each SNP (values)
    '''
    # Check that bulk values make sense
    if bulkAlleles != [] and bulkOccurrence != None:
        assert all([ isinstance(val, int) for val in bulkAlleles ]), "Bulk alleles must be integers!"
        assert all([ val >= 2 for val in bulkAlleles ]), "Bulk alleles must be >= 2!"
        assert isinstance(bulkOccurrence, float) or isinstance(bulkOccurrence, int), \
            "Bulk occurrence must be an int or float!"
        assert 0 < bulkOccurrence <= 1, "Bulk occurrence must be >0 and <=1!"
        shouldBulk = True
    else:
        shouldBulk = False
    
    # Iterate through file and grab statistical values
    edDict = {}
    with open(tsvFile, "r") as fileIn:
        firstLine = True
        for line in fileIn:
            sl = line.rstrip("\r\n ").split("\t")
            
            # Handle header lines
            if firstLine:
                assert all([ hv in sl for hv in HEADER_VALUES ]), "Header line doesn't contain expected values!"
                
                chromIndex = sl.index(HEADER_VALUES[0])
                posIndex = sl.index(HEADER_VALUES[1])
                distIndex = sl.index(HEADER_VALUES[2])
                bulk1Index = sl.index(HEADER_VALUES[3])
                bulk2Index = sl.index(HEADER_VALUES[4])
                
                firstLine = False
            
            # Handle content lines
            else:
                # Parse out relevant details from this line
                chrom, pos, edist, bulk1, bulk2 = sl[chromIndex], sl[posIndex], \
                    sl[distIndex], sl[bulk1Index], sl[bulk2Index]
                
                # Check bulk occurrence and skip if applicable
                if shouldBulk:
                    # Parse out the bulk allele counts
                    bulk1Alleles, bulk2Alleles = bulkAlleles
                    
                    # Check if the bulk occurrence is met
                    if ((int(bulk1) / bulk1Alleles) < bulkOccurrence) and ((int(bulk2) / bulk2Alleles) < bulkOccurrence):
                        continue
                
                # Store data
                edDict.setdefault(chrom, {})
                edDict[chrom][int(pos)] = {
                    "bulk1_alleles": int(bulk1),
                    "bulk2_alleles": int(bulk2),
                    "euclidean_distance": float(edist) if edist != "." else "."
                }
                
    return edDict

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

def determine_snp_proximity_from_matches(pos, matches):
    '''
    Parameters:
        pos -- an integer indicating the position of a SNP in the genome
        matches -- list containing ZS_GFF3IO.Feature objects that are
                   the result of a NCLS search for a SNP position; the
                   NCLS search is assumed to have been performed over
                   the "gene" feature type.
    Returns:
        proximalGene -- a string indicating the gene ID of the nearest gene
        location -- a string indicating where the SNP was located
                    i.e., "CDS", "UTR", "intron", or "intergenic".
    '''
    if matches == []:
        return "intergenic"
    else:
        bestMatch = "intron"
        for geneFeature in matches:
            for childFeature in geneFeature.retrieve_all_children():
                if int(pos) >= childFeature.start and int(pos) <= childFeature.end:
                    # Update best match according to feature type
                    if childFeature.type == "CDS":
                        bestMatch = "CDS"
                    elif childFeature.type == "exon":
                        if bestMatch != "CDS":
                            bestMatch = "UTR"
    return bestMatch

def generate_proximity_dicts(gff3Obj, snpsDict, radiusSize=50000):
    '''
    Parameters:
        gff3Obj -- a ZS_GFF3IO.GFF3 object with NCLS indexing performed on
                   gene features.
        snpsDict -- a dict with structure like:
                    {
                        'contig1': {
                            pos1: < contents don't matter >,
                            pos2: ...,
                            ...
                        },
                        'contig2': { ... },
                        ...
                    }
        radiusSize -- [OPTIONAL] an integer indicating the radius surrounding a gene that
                      you want to consider as being 'proximal' to a gene; default
                      == 50,000 bp.
    Returns:
        snpDict -- a dictionary with structure like:
                   {
                       'contig1': {
                           pos1: {
                               'left_of': [ 'geneID1', 'geneID2', ... ],
                               'left_dist': [ distance1, distance2, ... ],
                               'right_of': [ 'geneID1', 'geneID2', ... ],
                               'right_dist': [ distance1, distance2, ... ]
                           },
                           pos2: { ... },
                           ...
                       },
                       'contig2': { ... },
                   }
        geneDict -- a dictionary with structure like:
                    {
                        'contig1': {
                            'geneID1': {
                                'left_snp': [ pos1, pos2, ... ],
                                'left_dist': [ distance1, distance2, ... ],
                                'right_snp': [ pos1, pos2, ... ],
                                'right_dist': [ distance1, distance2, ... ]
                            },
                            'geneID2': { ... },
                            ...
                        },
                        'contig2': { ... },
                    }
    '''
    def setdefault_snp(snpDict, contig, pos):
        snpDict.setdefault(contig, {})
        snpDict[contig].setdefault(pos, {
            "within": [],
            "within_type": [],
            "left_of": [],
            "left_dist": [],
            "right_of": [],
            "right_dist": []
        })
    
    def setdefault_gene(geneDict, contig, geneID):
        geneDict.setdefault(contig, {})
        geneDict[contig].setdefault(geneID, {
            "within_snp": [],
            "within_type": [],
            "left_snp": [],
            "left_dist": [],
            "right_snp": [],
            "right_dist": []
        })
    
    snpDict = {}
    geneDict = {}
    for contig, posDict in snpsDict.items():
        for pos, _data in posDict.items():
            pos = int(pos) # just for testing
            
            # Locate the SNP in relation to gene features
            matches = gff3Obj.ncls_finder(pos, pos, "contig", contig)
            
            # Handle scenario where SNP is located within a gene
            if matches != []:
                location, geneID = haplotypes.determine_snp_location_from_matches(pos, matches)
                
                # Store data in SNP dict
                setdefault_snp(snpDict, contig, pos)
                snpDict[contig][pos]["within"].append(geneID)
                snpDict[contig][pos]["within_type"].append(location)
                
                # Store data in gene dict
                setdefault_gene(geneDict, contig, geneID)
                geneDict[contig][geneID]["within_snp"].append(pos)
                geneDict[contig][geneID]["within_type"].append(location)
                
            # Handle scenario where SNP is (potentially) located near a gene
            else:
                nearestLeft, nearestRight = windows.locate_genes_near_snp(gff3Obj, contig, pos,
                                                                          radiusSize=radiusSize,
                                                                          mrnaIndexing=False) # NCLS has gene indexing
                if nearestLeft[1] != None: # gene is left of snp; hence snp is right of gene
                    distance, geneID = nearestLeft
                    
                    # Store data in SNP dict
                    setdefault_snp(snpDict, contig, pos)
                    snpDict[contig][pos]["right_of"].append(geneID)
                    snpDict[contig][pos]["right_dist"].append(str(distance))
                    
                    # Store data in gene dict
                    setdefault_gene(geneDict, contig, geneID)
                    geneDict[contig][geneID]["right_snp"].append(pos)
                    geneDict[contig][geneID]["right_dist"].append(distance)
                
                if nearestRight[1] != None: # gene is right of snp; hence snp is left of gene
                    distance, geneID = nearestRight
                    
                    # Store data in SNP dict
                    setdefault_snp(snpDict, contig, pos)
                    snpDict[contig][pos]["left_of"].append(geneID)
                    snpDict[contig][pos]["left_dist"].append(str(distance))
                    
                    # Store data in gene dict
                    setdefault_gene(geneDict, contig, geneID)
                    geneDict[contig][geneID]["left_snp"].append(pos)
                    geneDict[contig][geneID]["left_dist"].append(distance)
    return snpDict, geneDict

def format_snp_types(snpDict):
    '''
    Helper function for formatting SNP data into within/left/right categories.
    
    Parameters:
        snpDict -- a dictionary with structure like:
                   {
                       'within_snp': [ pos1, pos2, ... ],
                       'within_type': [ "UTR", "CDS", "intron", ... ],
                       'left_snp': [ pos1, pos2, ... ],
                       'left_dist': [ distance1, distance2, ... ],
                       'right_snp': [ pos1, pos2, ... ],
                       'right_dist': [ distance1, distance2, ... ]
                   }
    Returns:
        snpDetails -- a dictionary with structure like:
                      {
                          < same keys as input, with string formatted outputs ... >
                      }
    '''
    WTYPE_ORDER = ["CDS", "UTR", "intron"]
    wTypes = set(snpDict["within_type"])
    return {
        "within_pos": len(snpDict["within_snp"]) if snpDict["within_snp"] != [] else ".",
        "within_type": ",".join(
            [ f'{wType}={snpDict["within_type"].count(wType)}'
            for wType in WTYPE_ORDER
            if wType in wTypes])
            if snpDict["within_type"] != [] else ".",
        "left": len(snpDict["left_snp"]) if snpDict["left_snp"] != [] else ".",
        "left_dist": sum(snpDict["left_dist"]) / len(snpDict["left_dist"]) if snpDict["left_dist"] != [] else ".",
        "right": len(snpDict["right_snp"]) if snpDict["right_snp"] != [] else ".",
        "right_dist": sum(snpDict["right_dist"]) / len(snpDict["right_dist"]) if snpDict["right_dist"] != [] else "."
    }

def format_euclidean_distances(edDict, contigID, snpDict):
    '''
    Helper function to format Euclidean distances for eventual output in tabular format.
    
    Parameters:
        edDict -- a dictionary resulting from parse_ed_file()
        contigID -- a string indicating the contig ID of the SNPs being considered
        snpDict -- a dictionary with structure like:
                   {
                       'within_snp': [ pos1, pos2, ... ],
                       'within_type': [ "UTR", "CDS", "intron", ... ],
                       'left_snp': [ pos1, pos2, ... ],
                       'left_dist': [ distance1, distance2, ... ],
                       'right_snp': [ pos1, pos2, ... ],
                       'right_dist': [ distance1, distance2, ... ]
                   }
    Returns:
        formatDict -- a dictionary with structure like:
                      {
                          'left': ( _min, _max, _mean ),
                          'right': ( _min, _max, _mean ),
                          'within': ( _min, _max, _mean ),
                          'total': mean of all SNPs
                      }
    '''
    leftDist = [
        edDict[contigID][pos]["euclidean_distance"]
        for pos in snpDict["left_snp"]
    ]
    rightDist = [
        edDict[contigID][pos]["euclidean_distance"]
        for pos in snpDict["right_snp"]
    ]
    withinDist = [
        edDict[contigID][pos]["euclidean_distance"]
        for pos in snpDict["within_snp"]
    ]
    totalAvg = (sum(leftDist) + sum(rightDist) + sum(withinDist)) / (len(leftDist) + len(rightDist) + len(withinDist))
    return {
        "left": get_edist_distribution(leftDist),
        "right": get_edist_distribution(rightDist),
        "within": get_edist_distribution(withinDist),
        "total": totalAvg
    }

def get_edist_distribution(edistList):
    '''
    Helper function to get the distribution of Euclidean distances in terms
    of their min, max, and median.
    
    Parameters:
        edistList -- a list of Euclidean distances
    Returns:
        _min -- the minimum Euclidean distance in the list.
        _max -- the maximum Euclidean distance in the list.
        _mean -- the mean Euclidean distance in the list.
    '''
    _min = min(edistList) if edistList != [] else "."
    _max = max(edistList) if edistList != [] else "."
    _mean = sum(edistList) / len(edistList) if edistList != [] else "."
    return _min, _max, _mean

def get_annot(annotDict, geneID, allowAbsence=False):
    '''
    Helper function to get annotation details for a gene ID with tolerance for absence.
    
    Parameters:
        annotDict -- a dictionary resulting from parse_annot_file()
        geneID -- a string indicating the gene ID for which we want annotation details
        allowAbsence -- [OPTIONAL] a boolean indicating whether we should allow the
                        absence of a gene in the annotation file; default == False
    Returns:
        annot -- the value of annotDict[geneID] if it exists, or None if it does not
    '''
    try:
        return annotDict[geneID]
    except Exception as e:
        if allowAbsence:
            return None
        else:
            raise e

def main():
    # User input
    usage = """%(prog)s receives several files associated with a per-sample SNP
    prediction project dealing with two segregating populations. It combines relevant
    information into two comprehensive report TSV files. The first file is gene-centric
    and provides one output row per gene with information on all the SNPs proximal to it.
    The second file is variant-centric and provides one output row per SNP with information
    on which gene(s) it is proximal to.
    
    Required input files are:
    1) Difference ratio TSV file from calculate_persample_ed.py.
    2) GOextended annotation table from annotation_table_extend_GOs.py or BINge's
    generate_annotation_table.py.
    3) A GFF3 file with gene annotations.
    4) A GO .obo file.
    """
    
    # Required
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-d", dest="edistFile",
                   required=True,
                   help="Specify the Euclidean distance TSV file name")
    p.add_argument("-a", dest="annotFile",
                   required=True,
                   help="Specify the annotation table file name")
    p.add_argument("-g", dest="gff3File",
                   required=True,
                   help="Specify the genome annotation GFF3 file name")
    p.add_argument("-obo", dest="goOboFile",
                   required=True,
                   help="Specify the GO .obo file name")
    p.add_argument("-out", dest="outputPrefix",
                   required=True,
                   help="Specify the prefix for the two output files")
    # Optional
    p.add_argument("--ed", dest="euclideanDistance",
                   type=float,
                   required=False,
                   help="""Optionally, only report a gene if it is proximal to
                   at least one SNP with >= the indicated Euclidean distance value, or
                   a SNP if it has >= the Euclidean distance; default
                   == 0 (i.e., no filtering occurs, all results are reported)""",
                   default=0)
    p.add_argument("--radius", dest="radiusSize",
                   type=int,
                   required=False,
                   help="""Optionally, specify the radius surrounding a gene that you
                   want to consider as being 'proximal' to a gene; default == 50000 bp""",
                   default=50000)
    p.add_argument("--allowAbsence", dest="allowAbsence",
                   required=False,
                   action="store_true",
                   help="""Optionally, provide this flag if you want to allow genes
                   present in the GFF3 to be absent from the annotation file; genes
                   not found in the annotation file will be skipped""",
                   default=False)
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse GO.obo
    go = obo_parser.GODag(args.goOboFile)
    
    # Parse Euclidean distance file
    edDict = parse_ed_file(args.edistFile)
    
    # Parse GFF3 with NCLS indexing
    gff3Obj = ZS_GFF3IO.GFF3(args.gff3File, strict_parse=False)
    gff3Obj.create_ncls_index(typeToIndex="gene")
    
    # Produce proximity dicts
    snpDict, geneDict = generate_proximity_dicts(gff3Obj, edDict, args.radiusSize)
    
    # Parse annotation file for relevant entries
    mrnaDict = {
        mrnaFeature.ID : geneID
        for geneData in geneDict.values()
        for geneID in geneData.keys()
        for mrnaFeature in gff3Obj[geneID].mRNA
    }
    annotDict = parse_annot_file(args.annotFile, mrnaDict, set(mrnaDict.values()))
    
    # Write gene-centric output
    queriedGOs = {}
    needsReplace = set()
    with open(args.outputPrefix + ".gene_report.tsv", "w") as fileOut:
        # Write header line
        fileOut.write("{0}\n".format("\t".join([
            "#contig", "gene_ID", "strand",
            "coords", "gene_name", "GO_IDs", "GO_names",
            "within_num_snps", "within_types", "left_num_snps",
            "left_avg_distance", "right_num_snps", "right_avg_distance",
            "left_edist_min", "left_edist_max", "left_edist_mean",
            "right_edist_min", "right_edist_max", "right_edist_mean",
            "within_edist_min", "within_edist_max", "within_edist_mean",
            "total_edist_mean"
        ])))
        
        # Write content lines
        for contigID, geneData in geneDict.items():
            for geneID, snpData in geneData.items():
                # Get annotation details with tolerance for absence (if applicable)
                annot = get_annot(annotDict, geneID, args.allowAbsence)
                if annot == None:
                    continue
                
                # Skip if this gene does not meet our --euclideanDistance threshold
                meetsThreshold = any([
                    edDict[contigID][pos]["euclidean_distance"] >= args.euclideanDistance
                    for snpKey in ["left_snp", "right_snp"]
                    for pos in snpData[snpKey]
                ])
                if not meetsThreshold:
                    continue
                
                # Get GO names from IDs
                if annot["gos"] == ".":
                    goNames = "."
                else:
                    # Fix obsoletions in GO terms
                    replacedGOs = ZS_GO.fix_obsoletions(annot["gos"].split("; "), go, queriedGOs)
                    
                    # Get the GO names now
                    goNames = []
                    for term in replacedGOs:
                        try:
                            goNames.append(go.get(term).name)
                        except:
                            needsReplace.add(term)
                    goNames = "; ".join(goNames)
                
                # Format SNP details
                snpDetails = format_snp_types(snpData)
                
                # Obtain ED values
                formattedRatios = format_euclidean_distances(edDict, contigID, snpData)
                
                # Format output line
                outputLine = "{contig}\t{geneID}\t{strand}\t{coords}\t{geneName}\
\t{gos}\t{goNames}\t{within_pos}\t{within_type}\t{left}\
\t{left_dist}\t{right}\t{right_dist}\
\t{left_edist_min}\t{left_edist_max}\t{left_edist_mean}\
\t{right_edist_min}\t{right_edist_max}\t{right_edist_mean}\
\t{within_edist_min}\t{within_edist_max}\t{within_edist_mean}\t{total_edist_mean}\n".format(
                    contig = contigID,
                    geneID = geneID,
                    strand = gff3Obj[geneID].strand,
                    coords = "{0}-{1}".format(*gff3Obj[geneID].coords),
                    geneName = annot["name"],
                    gos = annot["gos"],
                    goNames = goNames,
                    within_pos = snpDetails["within_pos"],
                    within_type = snpDetails["within_type"],
                    left = snpDetails["left"],
                    left_dist = snpDetails["left_dist"],
                    right = snpDetails["right"],
                    right_dist = snpDetails["right_dist"],
                    left_edist_min = formattedRatios["left"][0],
                    left_edist_max = formattedRatios["left"][1],
                    left_edist_mean = formattedRatios["left"][2],
                    right_edist_min = formattedRatios["right"][0],
                    right_edist_max = formattedRatios["right"][1],
                    right_edist_mean = formattedRatios["right"][2],
                    within_edist_min = formattedRatios["within"][0],
                    within_edist_max = formattedRatios["within"][1],
                    within_edist_mean = formattedRatios["within"][2],
                    total_edist_mean = formattedRatios["total"]
                )
                
                # Write output line
                fileOut.write(outputLine)
    
    if len(needsReplace) != 0:
        print("Some GOs were not found in the .obo file and need to be replaced for the gene-centric report.")
        print(f"These are: {needsReplace}")
    
    # Write variant-centric output
    with open(args.outputPrefix + ".variant_report.tsv", "w") as fileOut:
        # Write header line
        fileOut.write("{0}\n".format("\t".join([
            "#contig", "position", "left_of", "right_of",
            "within", "within_type", "bulk1_alleles", "bulk2_alleles",
            "euclidean_distance"
        ])))
        
        # Write content lines
        for contigID, snpData in snpDict.items():
            for pos, snpDetails in snpData.items():
                # Get Euclidean distance
                edData = edDict[contigID][pos]
                
                # Skip if this SNP does not meet our --euclideanDistance threshold
                meetsThreshold = edData["euclidean_distance"] >= args.euclideanDistance
                if not meetsThreshold:
                    continue
                
                # Format gene names
                withinName = "." # set default condition
                for geneID in snpDetails["within"]: # there's only ever 1 value in this list
                    annot = get_annot(annotDict, geneID, args.allowAbsence)
                    if annot == None:
                        continue
                    withinName = f'{geneID} ({annot["name"]})'
                leftName = "."
                for geneID in snpDetails["left_of"]: # iterating thru is just a convenience
                    annot = get_annot(annotDict, geneID, args.allowAbsence)
                    if annot == None:
                        continue
                    leftName = f'{geneID} ({annot["name"]})'
                rightName = "."
                for geneID in snpDetails["right_of"]:
                    annot = get_annot(annotDict, geneID, args.allowAbsence)
                    if annot == None:
                        continue
                    rightName = f'{geneID} ({annot["name"]})'
                
                # Format output line
                outputLine = "{contig}\t{position}\t{left}\t{right}\
\t{within}\t{withinTypes}\t{bulk1Alleles}\t{bulk2Alleles}\t{euclideanDistance}\n".format(
                    contig = contigID,
                    position = str(pos),
                    left = leftName,
                    right = rightName,
                    within = withinName,
                    withinTypes = snpDetails["within_type"][0] if snpDetails["within_type"] != [] else ".",
                    bulk1Alleles = edData["bulk1_alleles"],
                    bulk2Alleles = edData["bulk2_alleles"],
                    euclideanDistance = edData["euclidean_distance"]
                )
                
                # Write output line
                fileOut.write(outputLine)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

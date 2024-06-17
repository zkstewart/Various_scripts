#! python3
# report_persample_diffratio.py
# Script to combine several input files into a cohesive report file
# to understand the results of a SNP prediction project where the
# goals are to obtain information on segregating populations, but where
# SNP prediction has occurred per-sample. It builds upon the
# calculate_persample_diffratio.py script which may have been filtered.

# Load normal/pip packages
import os, argparse, sys, re
from goatools import obo_parser

# Load ZS_IO Class code
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))) # 4 dirs up is where we find Function_packages
from Function_packages import ZS_GFF3IO, ZS_SeqIO, ZS_GO, ZS_VCFIO

# Load functions from other scripts
from filter_persample_diffratio import parse_persample_diffratio_file

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))) # 3 dirs up is where we find windows
import windows # pulls relevant functions from snp_proximity_report.py
import haplotypes

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.diffratioFile):
        print(f'I am unable to locate the difference ratio file ({args.diffratioFile})')
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
    if 0 > args.differenceRatio or 1 < args.differenceRatio:
        print("--differenceRatio must be in the range 0 (no filtering) to 1 " + 
                "(retain only variants that completely segregate between the two populations)")
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

def format_difference_ratios(diffratioDict, contigID, snpDict):
    '''
    Helper function to format difference ratios for eventual output in tabular format.
    
    Parameters:
        diffratioDict -- a dictionary resulting from parse_persample_diffratio_file()
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
    leftDiff = [
        diffratioDict[contigID][pos]["differenceRatio"]
        for pos in snpDict["left_snp"]
    ]
    rightDiff = [
        diffratioDict[contigID][pos]["differenceRatio"]
        for pos in snpDict["right_snp"]
    ]
    withinDiff = [
        diffratioDict[contigID][pos]["differenceRatio"]
        for pos in snpDict["within_snp"]
    ]
    totalAvg = (sum(leftDiff) + sum(rightDiff) + sum(withinDiff)) / (len(leftDiff) + len(rightDiff) + len(withinDiff))
    return {
        "left": get_diffratio_distribution(leftDiff),
        "right": get_diffratio_distribution(rightDiff),
        "within": get_diffratio_distribution(withinDiff),
        "total": totalAvg
    }

def get_diffratio_distribution(diffRatioList):
    '''
    Helper function to get the distribution of difference ratios in terms
    of their min, max, and median.
    
    Parameters:
        diffRatioList -- a list of difference ratios.
    Returns:
        _min -- the minimum difference ratio in the list.
        _max -- the maximum difference ratio in the list.
        _mean -- the mean difference ratio in the list.
    '''
    _min = min(diffRatioList) if diffRatioList != [] else "."
    _max = max(diffRatioList) if diffRatioList != [] else "."
    _mean = sum(diffRatioList) / len(diffRatioList) if diffRatioList != [] else "."
    return _min, _max, _mean

def main():
    # User input
    usage = """%(prog)s receives several files associated with a per-sample SNP
    prediction project dealing with two segregating populations. It combines relevant
    information into two comprehensive report TSV files. The first file is gene-centric
    and provides one output row per gene with information on all the SNPs proximal to it.
    The second file is variant-centric and provides one output row per SNP with information
    on which gene(s) it is proximal to.
    
    Required input files are:
    1) Difference ratio TSV file from calculate_persample_diffratio.py.
    2) GOextended annotation table from annotation_table_extend_GOs.py or BINge's
    generate_annotation_table.py.
    3) A GFF3 file with gene annotations.
    4) A GO .obo file.
    """
    
    # Required
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-d", dest="diffratioFile",
                   required=True,
                   help="Specify the difference ratio TSV file name")
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
    p.add_argument("--differenceRatio", dest="differenceRatio",
                   type=float,
                   required=False,
                   help="""Optionally, only report a gene if it is proximal to
                   at least one SNP with >= the indicated difference ratio value, or
                   a SNP if it has >= the difference ratio; default
                   == 0 (i.e., no filtering occurs, all results are reported)""",
                   default=0)
    p.add_argument("--radius", dest="radiusSize",
                   type=int,
                   required=False,
                   help="""Optionally, specify the radius surrounding a gene that you
                   want to consider as being 'proximal' to a gene; default == 50000 bp""",
                   default=50000)
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse GO.obo
    go = obo_parser.GODag(args.goOboFile)
    
    # Parse difference ratio file
    diffratioDict = parse_persample_diffratio_file(args.diffratioFile)
    
    # Parse GFF3 with NCLS indexing
    gff3Obj = ZS_GFF3IO.GFF3(args.gff3File, strict_parse=False)
    gff3Obj.create_ncls_index(typeToIndex="gene")
    
    # Produce proximity dicts
    snpDict, geneDict = generate_proximity_dicts(gff3Obj, diffratioDict, args.radiusSize)
    
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
            "left_diff_min", "left_diff_max", "left_diff_mean",
            "right_diff_min", "right_diff_max", "right_diff_mean",
            "within_diff_min", "within_diff_max", "within_diff_mean",
            "total_diff_mean"
        ])))
        
        # Write content lines
        for contigID, geneData in geneDict.items():
            for geneID, snpData in geneData.items():
                annot = annotDict[geneID]
                
                # Skip if this gene does not meet our --differenceRatio threshold
                meetsThreshold = any([
                    diffratioDict[contigID][pos]["differenceRatio"] >= args.differenceRatio
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
                
                # Obtain difference ratio values
                formattedRatios = format_difference_ratios(diffratioDict, contigID, snpData)
                
                # Format output line
                outputLine = "{contig}\t{geneID}\t{strand}\t{coords}\t{geneName}\
\t{gos}\t{goNames}\t{within_pos}\t{within_type}\t{left}\
\t{left_dist}\t{right}\t{right_dist}\
\t{left_diff_min}\t{left_diff_max}\t{left_diff_mean}\
\t{right_diff_min}\t{right_diff_max}\t{right_diff_mean}\
\t{within_diff_min}\t{within_diff_max}\t{within_diff_mean}\t{total_diff_mean}\n".format(
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
                    left_diff_min = formattedRatios["left"][0],
                    left_diff_max = formattedRatios["left"][1],
                    left_diff_mean = formattedRatios["left"][2],
                    right_diff_min = formattedRatios["right"][0],
                    right_diff_max = formattedRatios["right"][1],
                    right_diff_mean = formattedRatios["right"][2],
                    within_diff_min = formattedRatios["within"][0],
                    within_diff_max = formattedRatios["within"][1],
                    within_diff_mean = formattedRatios["within"][2],
                    total_diff_mean = formattedRatios["total"]
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
            "#contig", "position", "variant", "left_of", "right_of",
            "within", "within_type", "bulk1_alleles", "bulk2_alleles",
            "difference_ratio"
        ])))
        
        # Write content lines
        for contigID, snpData in snpDict.items():
            for pos, snpDetails in snpData.items():
                # Get difference ratio
                diffData = diffratioDict[contigID][pos]
                
                # Skip if this SNP does not meet our --differenceRatio threshold
                meetsThreshold = diffData["differenceRatio"] >= args.differenceRatio
                if not meetsThreshold:
                    continue
                
                # Format gene names
                withinName = "." # set default condition
                for geneID in snpDetails["within"]: # there's only ever 1 value in this list
                    annot = annotDict[geneID]
                    withinName = f'{geneID} ({annot["name"]})'
                leftName = "."
                for geneID in snpDetails["left_of"]: # iterating thru is just a convenience
                    annot = annotDict[geneID]
                    leftName = f'{geneID} ({annot["name"]})'
                rightName = "."
                for geneID in snpDetails["right_of"]:
                    annot = annotDict[geneID]
                    rightName = f'{geneID} ({annot["name"]})'
                
                # Format output line
                outputLine = "{contig}\t{position}\t{variant}\t{left}\t{right}\
\t{within}\t{withinTypes}\t{bulk1Alleles}\t{bulk2Alleles}\t{differenceRatio}\n".format(
                    contig = contigID,
                    position = str(pos),
                    variant = diffData["variant"],
                    left = leftName,
                    right = rightName,
                    within = withinName,
                    withinTypes = snpDetails["within_type"][0] if snpDetails["within_type"] != [] else ".",
                    bulk1Alleles = diffData["bulk1_alleles"],
                    bulk2Alleles = diffData["bulk2_alleles"],
                    differenceRatio = diffData["differenceRatio"]
                )
                
                # Write output line
                fileOut.write(outputLine)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

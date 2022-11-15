#! python3
# snp_proximity_report.py
# Script to take in a GFF3 file and a VCF and
# identify genes near to or containing SNPs/variants.

import os, argparse, math, sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))) # 3 dirs up is where we find dependencies
from Function_packages import ZS_GFF3IO

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
    if args.windowsTable != None:
        if not os.path.isfile(args.windowsTable):
            print('I am unable to locate the windows table file (' + args.windowsTable + ')')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    if args.annotationTable != None:
        if not os.path.isfile(args.annotationTable):
            print('I am unable to locate the annotation table file (' + args.annotationTable + ')')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    # Validate output file location
    if os.path.isfile(args.outputFilePrefix + ".snp_proximity.tsv"):
        print(f'File already exists at output location ({args.outputFilePrefix + ".snp_proximity.tsv"})')
        print('Make sure you specify a unique file name and try again.')
        quit()
    elif os.path.isfile(args.outputFilePrefix + ".gene_proximity.tsv"):
        print(f'File already exists at output location ({args.outputFilePrefix + ".gene_proximity.tsv"})')
        print('Make sure you specify a unique file name and try again.')
        quit()
    # Validate other parameters
    if args.windowsColumns != None:
        if len(args.windowsColumns) != 3:
            print("windowsColumns, if provided, must consist of 3 separate column names (no spaces please!)")
            quit()

def get_snp_positions_from_vcf(vcfFile):
    '''
    This function is designed to parse VCFs where comments (#) can be safely ignored.
    Any other line is expected to have column 1 be the chromosome/contig name, and column
    2 is an integer of the position of the variant.
    
    However, there can be other file formats which fit this description e.g., the output
    from QTL-seq. Hence, this script can also handle file formats where the first line
    lacks comment (#) heading but still should be ignored.
    '''
    snpPositions = {}
    with open(vcfFile, "r") as fileIn:
        firstLine = True
        for line in fileIn:
            l = line.rstrip("\r\n").split("\t")
            
            # Handle header lines
            if line.startswith("#"):
                firstLine = False
                continue
            elif firstLine is True:
                firstLine = False
                if not l[1].isdigit():
                    continue
            
            # Extract relevant details
            chrom = l[0]
            pos = int(l[1])
            
            # Store in dictionary
            if chrom not in snpPositions:
                snpPositions[chrom] = []
            snpPositions[chrom].append(pos)
    
    return snpPositions

def generate_snp_proximity_dict(gff3Obj, snpPositions, windows=None):
    '''
    Reads in a GFF3 and locates SNPs in relation to these features.
    
    Parameters:
        gff3Obj -- a ZS_GFF3IO.GFF3 object with NCLS indexing of the gene features
        snpPositions -- a dictionary as produced by get_snp_positions_from_vcf()
                        with structure like:
                        {
                            'contig1': [pos1, pos2, ...],
                            'contig2': ...,
                        }
        windows -- a dictionary with structure like:
                   {
                       chromID1: [ [start1, end1], [start2, end2], ... ],
                       chromID2: [ ... ],
                       ...
                   }
                   OR, None if this is irrelevant.
    Returns:
        snpProximityDict -- a dictionary with structure like:
                            {
                                'contig1': {
                                    'pos1': [
                                        insideGeneID,
                                        insideGeneLocation,
                                        insideGeneCoords,
                                        leftGeneID,
                                        leftGeneDistance,
                                        leftGeneStrand,
                                        rightGeneID,
                                        rightGeneDistance,
                                        rightGeneStrand
                                    ]
                                    'pos2': ...
                                },
                                'contig2': ...
                            }
        geneProximityDict -- a dictionary with structure like:
                             {
                                 'gene1': [contigPosObject1, contigPosObject2]
                             }
                             ... where contigPosObject is a pointer to the same list
                             indexed by snpProximityDict.
    '''
    snpProximityDict = {}
    geneProximityDict = {}
    
    for contig, positions in snpPositions.items():
        snpProximityDict[contig] = {}
        for pos in positions:
            # Establish this SNP's results storage structure
            snpResultDict = {
                "insideGeneID": ".",
                "insideGeneLocation": ".",
                "insideGeneCoords": ".",
                "leftGeneID": ".",
                "leftGeneDistance": ".",
                "leftGeneStrand": ".",
                "rightGeneID": ".",
                "rightGeneDistance": ".",
                "rightGeneStrand": "."
            }
            if windows != None:
                snpResultDict["isInsideWindow"] = "."
            
            # Establish whether this SNP is inside a window (if relevant)
            if windows != None:
                snpResultDict["isInsideWindow"] = "Y" if locate_snp_within_window(contig, pos, windows) else "N"
            
            # Scenario 1: SNP located within a gene
            matches = gff3Obj.ncls_finder(pos, pos, "contig", contig)
            if matches != []: # i.e., if the SNP is located within a gene
                snpLocation, geneID = locate_snp_within_gene(matches, pos)
                snpResultDict["insideGeneID"] = geneID
                snpResultDict["insideGeneLocation"] = snpLocation
                snpResultDict["insideGeneCoords"] = gff3Obj[geneID].coords
                
                geneProximityDict.setdefault(geneID, {
                        "snpWithin": [],
                        "snpRight": [],
                        "snpLeft": []
                    }
                )
                geneProximityDict[geneID]["snpWithin"].append(snpLocation)
                if windows != None:
                    geneProximityDict[geneID]["isInsideWindow"] = snpResultDict["isInsideWindow"][:]
                
            
            # Scenario 2: SNP located not inside, but near a gene
            else:
                nearestLeft, nearestRight = locate_genes_near_snp(gff3Obj, contig, pos)
                if nearestLeft[1] != None:
                    snpResultDict["leftGeneID"] = nearestLeft[1]
                    snpResultDict["leftGeneDistance"] = nearestLeft[0]
                    snpResultDict["leftGeneStrand"] = gff3Obj[nearestLeft[1]].strand
                    
                    geneProximityDict.setdefault(nearestLeft[1], {
                        "snpWithin": [],
                        "snpRight": [],
                        "snpLeft": []
                        }
                    )
                    geneProximityDict[nearestLeft[1]]["snpRight"].append(pos) # if gene is left of SNP, SNP is right of gene!
                    if windows != None:
                        geneProximityDict[nearestLeft[1]]["isInsideWindow"] = snpResultDict["isInsideWindow"][:]
                    
                if nearestRight[1] != None:
                    snpResultDict["rightGeneID"] = nearestRight[1]
                    snpResultDict["rightGeneDistance"] = nearestRight[0]
                    snpResultDict["rightGeneStrand"] = gff3Obj[nearestRight[1]].strand
                    
                    geneProximityDict.setdefault(nearestRight[1], {
                            "snpWithin": [],
                            "snpRight": [],
                            "snpLeft": []
                        }
                    )
                    geneProximityDict[nearestRight[1]]["snpLeft"].append(pos)
                    if windows != None:
                        geneProximityDict[nearestRight[1]]["isInsideWindow"] = snpResultDict["isInsideWindow"][:]
            
            # Store whatever results we have found if relevant
            if any([value != "." for value in snpResultDict.values()]):
                snpProximityDict[contig][pos] = snpResultDict
    
    return snpProximityDict, geneProximityDict

def locate_snp_within_gene(matches, snpPos):
    bestLocation = None
    bestID = None
    for feature in matches:
        # Look for CDS overlap
        try:
            cds = False
            for cdsFeature in feature.CDS:
                if cdsFeature.start <= snpPos and snpPos <= cdsFeature.end: # if it overlaps...
                    cds = True
                    break
        except:
            cds = False
        # Look for exon overlap (if CDS overlap wasn't found)
        if cds is False:
            try:
                exon = False
                for exonFeature in feature.exon:
                    if exonFeature.start <= snpPos and snpPos <= exonFeature.end: # if it overlaps...
                        exon = True
                        break
            except:
                exon = False
        # Get the ID to index this under
        if feature.type == "mRNA":
            geneID = feature.Parent
        else:
            geneID = feature.ID
        # Get the position of the SNP relative to the gene
        if cds:
            bestLocation = "CDS"
            bestID = geneID
            return bestLocation, bestID # immediately return since it can't get better
        elif exon:
            bestLocation = "UTR" if bestLocation != "CDS" else bestLocation
            bestID = geneID if bestLocation != "CDS" else bestID
        else:
            bestLocation = "intron" if (bestLocation != "CDS" and bestLocation != "UTR") else bestLocation
            bestID = geneID if (bestLocation != "CDS" and bestLocation != "UTR") else bestID
    return bestLocation, bestID

def locate_genes_near_snp(gff3Obj, contig, snpPos, radiusSize=50000):
    '''
    This function assumes SNPs have been filtered so that SNPs overlapping a gene are
    absent. It won't cause any bugs, but the notion of finding a gene to the left
    and right of a SNP that's inside of a gene loses meaning since we'll find the
    same gene in both directions.
    
    Parameters:
        radiusSize -- an integer indicating the distance to check surrounding
                      either side of a SNP for a gene
    '''
    # Get genes within radius length of the SNP position
    start = (snpPos - radiusSize) if (snpPos - radiusSize) >= 0 else 0
    end = snpPos + radiusSize
    mrnaMatches = gff3Obj.ncls_finder(start, end, "contig", contig)
    geneMatches = [gff3Obj[match.Parent] for match in mrnaMatches]
    
    nearestLeft = [math.inf, None]
    nearestRight = [math.inf, None]
    checkedIDs = set()
    for geneFeature in geneMatches:
        if geneFeature.ID in checkedIDs:
            continue
        
        isLeft = geneFeature.end < snpPos
        isRight = geneFeature.start > snpPos
        if isLeft:
            dist = snpPos - geneFeature.end
            if dist < nearestLeft[0]:
                nearestLeft = [dist, geneFeature.ID]
        elif isRight:
            dist = geneFeature.start - snpPos
            if dist < nearestRight[0]:
                nearestRight = [dist, geneFeature.ID]
        checkedIDs.add(geneFeature.ID)
    
    return nearestLeft, nearestRight

def find_genes_nearest_snp(chromosomeGeneRegions, snpPositions, numNearest):
    '''
    This function is deprecated, but I'll leave it in here just in case I ever
    want it to jog my memory of how it used to behave.
    '''
    genes = []
    for chrom, positions in snpPositions.items():
        for pos in positions:
            if chrom not in chromosomeGeneRegions:
                continue

            nNearest = []
            for gStart, gEnd, geneID in chromosomeGeneRegions[chrom]:
                isLeft = gEnd < pos
                isRight = gStart > pos
                if isLeft:
                    dist = pos - gEnd
                    nNearest.append([dist, geneID])
                elif isRight:
                    dist = gStart - pos
                    nNearest.append([dist, geneID])
            nNearest.sort(key = lambda x: x[0])
            
            for n in range(numNearest):
                if n < len(nNearest): # only append as many entries as exist if numNearest is large
                    genes.append(nNearest[n][1])
            
    return list(set(genes))

def parse_windows_table(windowsTableFile, windowsColumns):
    '''
    Intended to receive a CSV or TSV file containing at least 3 columns
    which indicate the chromosome/contig, start, and end positions of
    predicted QTL windows.
    
    Parameters:
        windowsTableFile -- a string indicating the location of the windows
                            CSV / TSV file.
        windowsColumns -- a list containing 3 string values which correspond to
                          values in the header of the table file.
    Returns:
        windows -- a dictionary with structure like:
                   {
                       chromID1: [ [start1, end1], [start2, end2], ... ],
                       chromID2: [ ... ],
                       ...
                   }
    '''
    windows = {}
    
    skipLine = True
    with open(windowsTableFile, "r") as fileIn:
        for line in fileIn:
            l = line.rstrip("\r\n ").replace("'", "").replace('"', "")
            sl = l.split("\t") if "\t" in l else l.split(",")
            # Skip any lines prior to our header line; handle the header when we find it
            if skipLine == True:
                headerLine = all([c in sl for c in windowsColumns])
                if headerLine == True:
                    skipLine = False
                    colIndices = [sl.index(c) for c in windowsColumns]
            # Handle content lines
            else:
                chrom, start, end = [sl[i] for i in colIndices]
                windows.setdefault(chrom, [])
                windows[chrom].append([int(start), int(end)])
    return windows

def locate_snp_within_window(contig, pos, windows):
    '''
    Return
    
    Parameters:
        contig -- a string indicating which contig/chromosome we want to check
        pos -- a string or integer indicating the position in the chromosome to check
        windows -- a dictionary with structure like:
                   {
                       chromID1: [ [start1, end1], [start2, end2], ... ],
                       chromID2: [ ... ],
                       ...
                   }
    Returns:
        isInsideWindow -- a boolean indicating whether the identified positon is inside
                          of a window (True) or not (False)
    '''
    assert str(pos).isdigit(), \
        "locate_snp_within_window needs pos to be an integer or a str representation of an int"
    
    if contig not in windows:
        return False
    else:
        for start, end in windows[contig]:
            if int(pos) >= start and int(pos) <= end:
                return True
    return False

def parse_annotation_table_geneids(annotationTableFile, geneSet):
    '''
    Specifically parses an annotation table file and returns annotations for the
    genes specified.
    
    Parameters:
        annotationTableFile -- a string indicating the location of a gene annotation
                               table produced using my (ZKS)'s pipeline.
        geneSet -- a set containing string values which correspond to one or more
                    rows in the annotation table.
    Returns:
        geneAnnotationDict -- a dictionary with structure like:
                              {
                                  'geneID1': ['geneName1', ["GO_ID1", "GO_ID2"]],
                                  'geneID2': ['geneName2', ["GO_ID3", "GO_ID4"]],
                                  ...
                              }
    '''
    geneAnnotationDict = {}
    with open(annotationTableFile, "r") as fileIn:
        for line in fileIn:
            l = line.rstrip("\r\n ")
            if l.startswith("#") or l == "":
                continue
            
            sl = l.split("\t")
            geneID = sl[0]
            if geneID not in geneSet:
                geneID = geneID.rsplit(".", maxsplit=1)[0] # attempt a rescue
                if geneID not in geneSet:
                    continue
            geneName = sl[3].split(" [")[0]
            GOs = sl[18]
            geneAnnotationDict[geneID] = [geneName, GOs]
    
    return geneAnnotationDict

## Main
def main():
    # User input
    usage = """%(prog)s reads in GFF3 and VCF files to obtain genes that are surrounding
    or contain SNPs in the VCF file. The output file is a TSV indicating genes and
    which SNP they contain/are near to.
    
    NOTE: This script can also handle the QTL-seq snp_index.p##.tsv files. In this format
    column 1 == contig ID and column 2 == contig position as an integer. The first line
    can be a header which will be ignored.
    """
    p = argparse.ArgumentParser(description=usage)
    ## Required
    p.add_argument("-v", dest="vcf",
                   required=True,
                   help="Input VCF containing only SNPs of interest (filter it beforehand)")
    p.add_argument("-g", dest="gff3",
                   required=True,
                   help="Input GFF3 file")
    p.add_argument("-o", dest="outputFilePrefix",
                   required=True,
                   help="Output prefix for .snp_proximity.tsv and .gene_proximity.tsv files")
    ## Optional
    p.add_argument("--annotationTable", dest="annotationTable",
                   required=False,
                   help="Optionally, specify the gene annotation table file (produced by ZKS pipeline scripts)")
    p.add_argument("--windowsTable", dest="windowsTable",
                   required=False,
                   help="Optionally, specify the QTL windows table file (must be TSV or CSV file)")
    p.add_argument("--windowsColumns", dest="windowsColumns", nargs="+",
                   required=False,
                   help="""Optionally, indicate which columns identify relevant window details;
                   must be ordered to specify the contig/chromosome, its start, and then its end
                   position (default=='CHROM start end')""",
                   default=["CHROM", "start", "end"])
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse GFF3 with NCLS indexing
    gff3Obj = ZS_GFF3IO.GFF3(args.gff3, strict_parse=False)
    gff3Obj.create_ncls_index(typeToIndex="mRNA")
    
    # Parse VCF as dictionary indexing contig: [positions]
    snpPositions = get_snp_positions_from_vcf(args.vcf)
    
    # Parse windows table file (if relevant)
    if args.windowsTable != None:
        windows = parse_windows_table(args.windowsTable, args.windowsColumns)
    else:
        windows = None
    
    # Get the proximity reporting dicts
    snpProximityDict, geneProximityDict = generate_snp_proximity_dict(gff3Obj, snpPositions, windows)
    
    # Parse the annotation table (if relevant)
    if args.annotationTable != None:
        annotations = parse_annotation_table_geneids(args.annotationTable, set(geneProximityDict.keys()))
    else:
        annotations = None
    
    # Produce output file 1: SNP-centric report
    with open(args.outputFilePrefix + ".snp_proximity.tsv", "w") as fileOut:
        # Handle header
        header = [
            "#contig", "pos", "inside",
            "at location", "left of", "at distance",
            "right of", "at distance"
        ]
        if windows != None:
            header.append("inside window")
        fileOut.write("\t".join(header) + "\n")
        
        # Handle content
        for contig in snpProximityDict.keys():
            orderedPositions = list(snpProximityDict[contig].keys())
            orderedPositions.sort()
            
            for pos in orderedPositions:
                snpResultsDict = snpProximityDict[contig][pos]
                contents = list(map(str, 
                    [contig, pos,
                    snpResultsDict["insideGeneID"], snpResultsDict["insideGeneLocation"],
                    
                    snpResultsDict["rightGeneID"] + f" ({snpResultsDict['rightGeneStrand']})" \
                        if snpResultsDict["rightGeneID"] != "." else ".",
                    snpResultsDict["rightGeneDistance"], 
                    
                    snpResultsDict["leftGeneID"] + f" ({snpResultsDict['leftGeneStrand']})" \
                        if snpResultsDict["leftGeneID"] != "." else ".",
                    snpResultsDict["leftGeneDistance"]]
                ))
                if windows != None:
                    contents.append(snpResultsDict["isInsideWindow"])
                
                fileOut.write("{0}\n".format("\t".join(contents)))
    
    # Produce output file 2: Gene-centric report
    with open(args.outputFilePrefix + ".gene_proximity.tsv", "w") as fileOut:
        # Handle header
        header = [
            "#gene", "coords", "strand", "contig",
            "snp within", "left snp at pos", "right snp at pos"
        ]
        if annotations != None:
            header.insert(1, "GO terms")
            header.insert(1, "gene name")
        if windows != None:
            header.append("inside window")
        fileOut.write("\t".join(header) + "\n")
        
        # Handle content
        orderedGenes = list(geneProximityDict.keys())
        orderedGenes.sort()
        for geneID in orderedGenes:
            geneResultsDict = geneProximityDict[geneID]
            geneCoords = f"{gff3Obj[geneID].start}-{gff3Obj[geneID].end}"
            geneStrand = gff3Obj[geneID].strand
            geneContig = gff3Obj[geneID].contig
            
            contents = [
                geneID, geneCoords, geneStrand, geneContig,
                "; ".join(geneResultsDict["snpWithin"]) \
                    if geneResultsDict["snpWithin"] != [] else ".",
                "; ".join(map(str, geneResultsDict["snpLeft"])) \
                    if geneResultsDict["snpLeft"] != [] else ".",
                "; ".join(map(str, geneResultsDict["snpRight"])) \
                    if geneResultsDict["snpRight"] != [] else "." # {5}
            ]
            if annotations != None:
                contents.insert(1, annotations[geneID][1])
                contents.insert(1, annotations[geneID][0])
            if windows != None:
                contents.append(geneResultsDict["isInsideWindow"])
            
            fileOut.write("{0}\n".format("\t".join(contents)))
    
    # Let user know everything went swimmingly
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

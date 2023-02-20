#! python3
# bsa_report_pipeline.py
# Script to combine several input files into a cohesive report file
# to understand the results of a bulked segregant analysis.

# Load normal/pip packages
import os, argparse, sys
from copy import deepcopy

# Load ZS_IO Class code
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))) # 4 dirs up is where we find Function_packages
from Function_packages import ZS_GFF3IO, ZS_SeqIO

# Load functions from other scripts
from add_qtl_info_to_proximity_reports import parse_scanone_file, parse_scantwo_file

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))) # 3 dirs up is where we find windows
import windows # pulls relevant functions from snp_proximity_report.py
import haplotypes

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.goFile):
        print(f'I am unable to locate the GOextended annotation file ({args.goFile})')
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
    if len(args.snpColumns) != 2:
        print("snpColumns should provide two values corresponding to header labels")
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
    
    Optional input files are:
    1) R/qtl scanone file.
    2) R/qtl scantwo file.
    
    Alongside these, you should specify several parameters to control the behaviour
    of this script.
    """
    
    # Required
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-go", dest="goFile",
                   required=True,
                   help="Specify the output file name")
    p.add_argument("-vcf", dest="vcfFile",
                   required=True,
                   help="Specify the output file name")
    ## TBD - second/multiple VCF combined??
    p.add_argument("-gff", dest="gffFile",
                   required=True,
                   help="Specify the output file name")
    p.add_argument("-fasta", dest="fastaFile",
                   required=True,
                   help="Specify the output file name")
    p.add_argument("-snp", dest="snpFile",
                   required=True,
                   help="Specify the output file name")
    p.add_argument("--snpColumns", dest="snpColumns", nargs="+",
                   required=False,
                   help="""Optionally, specify the header values of the two columns
                   to parse from the SNP index-like file
                   (default == 'CHROM POSI' """,
                   default=["CHROM", "POSI"])
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
    
    # Parse GFF3 with NCLS indexing
    gff3Obj = ZS_GFF3IO.GFF3(args.gffFile, strict_parse=False)
    gff3Obj.create_ncls_index(typeToIndex="gene")
    
    # Parse FASTA file
    genomeFASTA_obj = ZS_SeqIO.FASTA(args.fastaFile)
    
    # Parse SNP index-like file for outlier SNPs
    snpPositions = get_snp_positions_from_indexlike(args.snpFile, args.snpColumns)
    
    # Parse VCF data for outlier SNP genotypes
    snpGenotypes = haplotypes.get_genotypes_from_vcf(args.vcfFile, snpPositions)
    
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
    annotDict = parse_annot_file(args.goFile, mrnaDict)
    
    # Combine everything into a table output
    with open(args.outputFileName, "w") as fileOut:
        # Write header line
        ## TBD...
        
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
            goNames = "..." # placefiller, TBD
            
            # Format SNP details
            snpDetails = format_snp_types(snpDict)
            
            # Format output line
            outputLine = "{geneID}\t{mrnaID}\t{strand}\t{coords}\t{geneName}\
            \t{gos}\t{goNames}\t{within_pos}\t{within_location}\t{left}\
            \t{left_dist}\t{right}\t{right_dist}\t{scan1}\t{scan2}\n".format(
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
            )
            
            # Write output row as line
            fileOut.write(outputLine)
        
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

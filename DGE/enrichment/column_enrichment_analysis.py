#!/usr/bin/env python3

import os
import sys
import pandas as pd

from scipy.stats import fisher_exact, boschloo_exact

def validate_args(args):
    # Validate input file location depending on type of input
    args.tableFileName = os.path.abspath(args.tableFileName)
    if not os.path.isfile(args.tableFileName):
        raise FileNotFoundError(f"-i '{args.tableFileName}' is not a file!")
    
    ACCEPTED_SUFFIXES = [".tsv", ".csv", ".xlsx"]
    if not any([ args.tableFileName.endswith(suffix) for suffix in ACCEPTED_SUFFIXES ]):
        raise ValueError("-i file must end with " + " or ".join(ACCEPTED_SUFFIXES))
    
    # Validate output file
    if os.path.exists(args.outputFileName):
        raise FileExistsError(f"-o output file '{args.outputFileName}' already exists!")
    args.outputFileName = os.path.abspath(args.outputFileName)
    
    if not os.path.isdir(os.path.dirname(args.outputFileName)):
        raise FileNotFoundError(f"-o parent directory '{os.path.dirname(args.outputFileName)}' does not exist!")

def parse_table_column(fileName, columnHeaders):
    if fileName.endswith(".xlsx"):
        return _parse_excel_column(excelFileName, columnHeaders)
    elif fileName.endswith(".tsv"):
        return _parse_text_column(excelFileName, columnHeaders, "\t")
    elif fileName.endswith(".csv"):
        return _parse_text_column(excelFileName, columnHeaders, ",")
    else:
        raise NotImplementedError(f"'{fileName}' ends with a file suffix that is not handled by parse_table_column()")
    degs = set()
    
    return degs

def _parse_text_column(textFileName, columnHeaders, delimiter):
    columns = []
    firstLine = True
    with open(textFileName, "r") as fileIn:
        for line in fileIn:
            sl = line.rstrip().split(delimiter)
            
            if firstLine:
                columnIndices = [None] * len(columnHeaders)
                for i, header in enumerate(columnHeaders):
                    try:
                        columnIndices[i] = sl.index(header)
                    except ValueError:
                        raise KeyError(f"'{header}' is not a column header in '{textFileName}'")
                firstLine = False
            else:
                columns.append([ sl[colIndex] for colIndex in columnIndices ])
    return columns

def _parse_excel_column(excelFileName, columnHeaders):
    columns = []
    df = pd.read_excel(excelFileName)
    try:
        for header in columnHeaders:
            try:
                columns.append(df[columneaders].to_list())
            except KeyError:
                raise KeyError(f"'{header}' is not a column header in '{excelFileName}'")
    return 

def make_contingency_table(annotationDict, degIDs):
    # Obtain our ID sets
    fullSet = set(annotationDict.keys())
    degSet = set([ x for x in degIDs if x in fullSet ])
    
    # Partition counts
    topleft = sum([ 1 if annotationDict[x] == 1 else 0 for x in degSet ]) # ncRNA DEG
    bottomleft = sum([ 1 if annotationDict[x] == 0 else 0 for x in degSet ]) # protein DEG
    topright = sum([ 1 if annotationDict[x] == 1 else 0 for x in fullSet if not x in degSet ]) # ncRNA non-DEG
    bottomright = sum([ 1 if annotationDict[x] == 0 else 0 for x in fullSet if not x in degSet ]) # protein non-DEG
    
    # Format the 2x2 table
    contingencyTable = [
            #  DEG     non-DEG
            [topleft, topright],      # ncRNA
            [bottomleft, bottomright] # protein
        ]
    
    return contingencyTable

def run_enrichment_test(contingencyTable, test="fisher"):
    '''
    Assumes a table where the left column are the DEGs and the first row
    are proteins. Or in more general terms, the top left cell is the
    "most important" value and the bottom right cell is the "least important"
    value.
    '''
    (topleft, topright), (bottomleft, bottomright) = contingencyTable
    
    # Run the exact test
    if test == "fisher":
        statistic, p = fisher_exact(contingencyTable)
    elif test == "boschloo":
        statistic, p = boschloo_exact(contingencyTable)
    else:
        raise ValueError(f"run_enrichment_test doesn't know what a '{test}' is")
    
    # Figure out if it's over or underrepresented
    leftRatio = topleft / (topleft + bottomleft)
    rightRatio = topright / (topright + bottomright)
    
    if leftRatio > rightRatio:
        representation = "over"
    else:
        representation = "under"
    
    return statistic, p, representation

def print_contingency_table(contingencyTable, leftColLabel, rightColLabel, firstRowLabel, secondRowLabel):
    (topleft, topright), (bottomleft, bottomright) = contingencyTable
    
    leftColWidth = max(len(leftColLabel), len(str(topleft)), len(str(bottomleft)))
    rightColWidth = max(len(rightColLabel), len(str(topright)), len(str(bottomright)))
    
    border = "+-" + ("-"*leftColWidth) + "-+-" + ("-"*rightColWidth) + "-+"
    header = "".join(["| ", # buffer left
                    " "*(leftColWidth-len(leftColLabel)),
                    leftColLabel,
                    " | ",
                    " "*(rightColWidth-len(rightColLabel)),
                    rightColLabel,
                    " |" # buffer right
                    ])
    row1 = "".join(["| ", # buffer left
                    " "*(leftColWidth-len(str(topleft))),
                    str(topleft),
                    " | ",
                    " "*(rightColWidth-len(str(topright))),
                    str(topright),
                    " |" # buffer right
                    ])
    row2 = "".join(["| ", # buffer left
                    " "*(leftColWidth-len(str(bottomleft))),
                    str(bottomleft),
                    " | ",
                    " "*(rightColWidth-len(str(bottomright))),
                    str(bottomright),
                    " |" # buffer right
                    ])
    
    maxRowWidth = max(len(firstRowLabel), len(secondRowLabel))
    
    yspace = " "*(2 + maxRowWidth + 1)
    yborder = "+-" + ("-"*maxRowWidth) + "-"
    ylabel1 = "| " + " "*(maxRowWidth-len(firstRowLabel)) + firstRowLabel + " "
    ylabel2 = "| " + " "*(maxRowWidth-len(secondRowLabel)) + secondRowLabel + " "
    
    formattedTable = [
        yspace + border,
        yspace + header,
        yborder + border,
        ylabel1 + row1,
        yborder + border,
        ylabel2 + row2,
        yborder + border
    ]
    for line in formattedTable:
        print(line)
    return formattedTable

def interpret_enrichment_test(table, p, representation, experiment,
                              leftColLabel, rightColLabel, firstRowLabel, secondRowLabel,
                              testApplied="Fisher's exact", pCutoff=0.05):
    print(f"# Contigency table for {experiment}")
    print_contingency_table(table, leftColLabel, rightColLabel, firstRowLabel, secondRowLabel)
    if p > pCutoff:
        print(f"> No statistical difference in '{firstRowLabel}' vs. '{secondRowLabel}' proportions " + 
              f"between the '{leftColLabel}' and '{rightColLabel}' groups")
    else:
        print(f"> '{firstRowLabel}' is {representation}represented relative to '{secondRowLabel}' " + 
              f"when comparing the '{leftColLabel}' and '{rightColLabel}' groups")
        print(f"> {testApplied} P-value = {p}")

def main():
    ##### USER INPUT SECTION
    usage = """%(prog)s concatenates FASTA files together, ensuring
    that proper newline spacing is maintained. It offers some filtering
    ability on sequence length, and can also multiline the output.
    
    NOTE TO SELF: This is a nightmare and almost impossible to generalise like I've been doing.
    It needs specifically formatted inputs to know the proper class membership / sizes and whatnot.
    """
    # Required arguments
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="tableFileName",
                   required=True,
                   help="Location of table file; expects suffix of .txt or .csv or .xlsx")
    p.add_argument("-c", dest="columnHeaders",
                   required=True,
                   nargs=3,
                   help="""Specify two headers: the first identifies the sequence IDs,
                   the second relates the sequences to a binary classifier (e.g., is DEG/is
                   not DEG), and the third assigns the sequence to a group""")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Location to write enrichment results")
    # Optional arguments
    p.add_argument("--labels", dest="binaryLabels",
                   required=False,
                   nargs=2,
                   help="""Optionally, specify the labels for use when
                   presenting results; default == 'deg non-de'""",
                   default=["deg", "non-de"])
    p.add_argument("--binary", dest="binaryValue",
                   required=False,
                   help="""Optionally, specify the value which identifies one of
                   the binary classes (the second column e.g., the non-DE group); default == '.'""",
                   default=".")
    p.add_argument("--group", dest="groupValue",
                   required=False,
                   help="""Optionally, specify the value which identifies one of
                   the group classes (the third column e.g., the non-DE group); default == '.'""",
                   default=None)
    
    args = p.parse_args()
    validate_args(args)
    
    ####
    bingeFasta = "/mnt/f/plant_group/ted/binge/BINge_clustering_representatives.aa"

    protFile = "/mnt/f/plant_group/ted/annotation/MMseqs2_results.filtered.tsv"
    ncrnaFile1 = "poncirus_trifoliata.filtered.outfmt6"
    ncrnaFile2 = "arabidopsis_thaliana.filtered.outfmt6"

    budpc1 = "/mnt/f/plant_group/ted/DGE/deseq2/latent_pc1.bud.tsv"
    leafpc1 = "/mnt/f/plant_group/ted/DGE/deseq2/latent_pc1.leaf.tsv"
    budtc = "/mnt/f/plant_group/ted/DGE/deseq2/timecourse.bud.tsv"
    leaftc = "/mnt/f/plant_group/ted/DGE/deseq2/timecourse.leaf.tsv"
    ####

    # Parse table columns
    binaryColumn, groupColumn = parse_table_column(args.tableFileName, args.columnHeaders)
    
    
    
    # Parse DESeq2 outputs
    budpc1IDs = parse_deseq2_ids(budpc1)
    leafpc1IDs = parse_deseq2_ids(leafpc1)

    budtcIDs = parse_deseq2_ids(budtc)
    leaftcIDs = parse_deseq2_ids(leaftc)

    # Parse BLAST outputs
    protParsed = Outfmt6(protFile, evalue=1e-3, numHits=1)
    protParsed.parse_to_dict()
    protDict = protParsed.resultsDict

    ncrnaParsed1 = Outfmt6(ncrnaFile1, evalue=1e-3, numHits=1)
    ncrnaParsed1.parse_to_dict()
    ncrnaDict1 = ncrnaParsed1.resultsDict

    ncrnaParsed2 = Outfmt6(ncrnaFile2, evalue=1e-3, numHits=1)
    ncrnaParsed2.parse_to_dict()
    ncrnaDict2 = ncrnaParsed2.resultsDict

    # Get a full list of sequence IDs
    clusterSequences = {}
    with open(bingeFasta, "r") as fileIn:
        for line in fileIn:
            if line.startswith(">"):
                clusterID = line[1:].rstrip().split(" ")[0]
                seqID = line.rstrip().split("=")[1]
                clusterSequences[seqID] = clusterID

    # Annotate each sequence as protein or ncRNA
    "I think there might be a genuinely interesting result w/r/t citrus genome architecture and ncRNA occurrence at the chromosome terminals???"
    annotationDict = {}
    for seqID, clusterID in clusterSequences.items():
        protHit = protDict[seqID][0] if seqID in protDict else None
        ncrnaHit1 = ncrnaDict1[clusterID][0] if clusterID in ncrnaDict1 else None
        ncrnaHit2 = ncrnaDict2[clusterID][0] if clusterID in ncrnaDict2 else None
        
        # Select the best ncRNA hit
        if ncrnaHit1 == None and ncrnaHit2 == None:
            ncrnaHit = None
        elif ncrnaHit1 == None and ncrnaHit2 != None:
            ncrnaHit = ncrnaHit2
        elif ncrnaHit1 != None and ncrnaHit2 == None:
            ncrnaHit = ncrnaHit1
        elif ncrnaHit1.score > ncrnaHit2.score:
            ncrnaHit = ncrnaHit1
        else:
            ncrnaHit = ncrnaHit2
        
        # Skip if no hit found
        if protHit == None and ncrnaHit == None:
            pass
        # Handle only one hit
        elif protHit != None and ncrnaHit == None:
            annotationDict[clusterID] = 0 # 0 == protein
        elif protHit == None and ncrnaHit != None:
            annotationDict[clusterID] = 1 # 1 == ncRNA
        # Compare score for protein and ncRNA
        elif protHit.score >= ncrnaHit.score: # prefer an equivalent protein hit over a ncRNA hit
            annotationDict[clusterID] = 0
        else:
            annotationDict[clusterID] = 1

    ####

    # Format 2x2 tables for each DGE experiment
    budpc1Table = make_contingency_table(annotationDict, budpc1IDs)
    budpc1_stat, budpc1_p, budpc1_representation = run_enrichment_test(budpc1Table, test="fisher")

    leafpc1Table = make_contingency_table(annotationDict, leafpc1IDs)
    leafpc1_stat, leafpc1_p, leafpc1_representation = run_enrichment_test(leafpc1Table, test="fisher")

    budtcTable = make_contingency_table(annotationDict, budtcIDs)
    budtc_stat, budtc_p, budtc_representation = run_enrichment_test(budtcTable, test="fisher")

    leaftcTable = make_contingency_table(annotationDict, leaftcIDs)
    leaftc_stat, leaftc_p, leaftc_representation = run_enrichment_test(leaftcTable, test="fisher")

    # Interpret the result
    leftColLabel = "DEG"
    rightColLabel = "non-DE"

    firstRowLabel = "ncRNA"
    secondRowLabel = "protein"

    experiment = "bud PC1"
    interpret_enrichment_test(budpc1Table, budpc1_p, budpc1_representation,
                            experiment, leftColLabel, rightColLabel,
                            firstRowLabel, secondRowLabel)
    print("")

    experiment = "leaf PC1"
    interpret_enrichment_test(leafpc1Table, leafpc1_p, leafpc1_representation,
                            experiment, leftColLabel, rightColLabel,
                            firstRowLabel, secondRowLabel)
    print("")

    experiment = "bud timecourse"
    interpret_enrichment_test(budtcTable, budtc_p, budtc_representation,
                            experiment, leftColLabel, rightColLabel,
                            firstRowLabel, secondRowLabel)
    print("")

    experiment = "leaf timecourse"
    interpret_enrichment_test(leaftcTable, leaftc_p, leaftc_representation,
                            experiment, leftColLabel, rightColLabel,
                            firstRowLabel, secondRowLabel)

#! python 3
# contaminant_curate.py
# A script/program designed to curate a transcriptome assembly
# by removing sequences that get significant hits to taxonomic ranks
# that are indicative of contamination e.g., bacteria, archaea,
# or viruses. This code is designed with animals in mind,
# and won't be useful otherwise.

import os, argparse, sys
from Bio import SeqIO
from ete3 import NCBITaxa
from collections import Counter

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))) # 3 dirs up is where we find dependencies
from Function_packages import ZS_BlastIO


def validate_args(args):
    # Validate input data locations
    if not os.path.isfile(args.transcriptomeFile):
        print(f'I am unable to locate the input transcriptome file ({args.transcriptomeFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.blastFile):
        print(f'I am unable to locate the input BLAST outfmt6 file ({args.blastFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate numeric arguments
    if args.contamEvalue < 0:
        print("contamEvalue must be number >= 0")
        quit()
    if args.rescueEvalue < 0:
        print("rescueEvalue must be number >= 0")
        quit()
    if not (0 <= args.ovlCutoff <= 1):
        print("ovlCutoff must be float in the range of 0 to 1 (inclusive)")
        quit()
    # Handle file output
    if os.path.isfile(args.outputFileName):
        print('The specified output file already exists. This program will not allowing overwriting.')
        print('Specify a different file name, or move the existing one and try again.')
        quit()

def parse_annots(annotationTable, evalueCutoff):
    '''
    Parses an annotation table for data specifically of interest to this script (contaminant_curate.py)
    
    Parameters:
        annotationTable -- a string indicating the location of an annotation table produced by
                           the Genome_analysis_scripts/annotation_table scripts.
        evalueCutoff -- a float specifying the maximum allowed E-value for indexing purposes
    Returns:
        annotDict -- a dictionary with structure like:
                     {
                         'seqID1': {
                             'taxons': [taxonInt1, taxonInt2, ...],
                             'alignment_length': [alignLenInt1, alignLenInt2, ...],
                             'evalues': [evalueFloat1, evalueFloat2, ...]
                         },
                         'seqID2' { ... },
                         ...
                     }
    '''
    annotDict = {}
    with open(annotationTable, "r") as fileIn:
        for line in fileIn:
            l = line.rstrip("\r\n ")
            sl = l.split("\t")
            
            # Skip irrelevant lines
            if sl == []:
                continue
            # Handle header line
            elif l.startswith("#"):
                taxonIndex = sl.index("NCBI_taxonomy_of_hits")
                alignIndex = sl.index("Alignment_length")
                evalueIndex = sl.index("Expect_value")
            # Handle content lines
            else:
                seqID, taxons, alignLens, evalues = sl[0], sl[taxonIndex], sl[alignIndex], sl[evalueIndex]
                
                # Skip if we have no hits
                if evalues == ".": # any variable will do here other than seqID
                    continue
                
                # Extract information into lists
                taxons = taxons.replace(" ", "").replace("]", "").split("[")
                alignLens = list(map(int, alignLens.replace(" ", "").replace("]", "").split("[")))
                evalues = list(map(float, evalues.replace(" ", "").replace("]", "").split("[")))
                
                # Drop any Evalues that don't meet cutoff
                evalues = [e for e in evalues if e <= evalueCutoff]
                taxons = taxons[0:len(evalues)] # can do it this way since it's sorted
                alignLens = alignLens[0:len(evalues)]
                
                # Store hit in dictionary (if relevant)
                if evalues != []:
                    annotDict[seqID] = {
                        "taxons": taxons,
                        "alignment_length": alignLens,
                        "evalues": evalues
                    }
    
    return annotDict

def modify_annot_taxons(annotDict):
    '''
    Receives the dictionary output by parse_annots() and adds another key in
    for sequences to indicate whether it belongs to a desired kingdom or a putative
    contaminant kingdom
    '''
    ncbi = NCBITaxa()
    ANIMAL_TAXID = 33208 # this is the one we want to see
    
    for annots in annotDict.values():
        kingdoms = []
        for taxon in annots["taxons"]:
            lineage = ncbi.get_lineage(taxon)
            if ANIMAL_TAXID in lineage:
                kingdoms.append("a") # a for animal
            else:
                kingdoms.append("c") # c for contaminant

        annots["kingdoms"] = kingdoms

def contaminant_curate_transcriptome(transcriptomeRecords, annotDict,
                             contamEvalue, rescueEvalue, ovlCutoff,
                             outputFileName):
    '''
    Parameters:
        transcriptomeRecords -- a Bio.SeqIO.parse() generator object
        annotDict -- a dictionary object obtained by parsing an annotation table
                     produced by the Genome_analysis_scripts/annotation_table scripts.
                     It will have a structure like:
                     {
                         'seqID1': {
                             'taxons': [taxonInt1, taxonInt2, ...],
                             'alignment_length': [alignLenInt1, alignLenInt2, ...],
                             'evalues': [evalueFloat1, evalueFloat2, ...],
                             'kingdoms': [a, c, ...] # a for animal, c for contaminant
                         },
                         'seqID2' { ... },
                         ...
                     }
        contamEvalue -- a float value indicating the threshold to enforce when identifying hits
                        to putative contaminants
        rescueValue -- a float value indicating the threshold to enforce when identifying hits
                       to putative non-contaminants (that will result in rescue of the sequence)
        ovlCutoff -- a float value in the range of 0 -> 1 that will set a threshold for how much
                     of the query sequence must be aligned against the contaminant for the hit to
                     be considered genuine
        outputFileName -- a string indicating the location to write the filtered
                          transcriptome file to
    Returns:
        contamCounts -- a dictionary with structure like:
                        {
                            taxonInt1: numContaminantsWithThisTaxon,
                            taxonInt2: ... ,
                            ...
                        }
    '''
    ncbi = NCBITaxa()
    contamCounts = {"num_contaminants": 0}
    
    with open(outputFileName, "w") as fileOut:
        for record in transcriptomeRecords:
            # Curate if we have BLAST hit(s)
            isContaminant = False
            if record.id in annotDict:
                annots = annotDict[record.id]
                seqLen = len(record)
                
                # Figure out whether we have any contaminant or animal hits that pass cutoffs
                isContaminant = False
                contaminants = []
                for i in range(len(annots["kingdoms"])):
                    # Get details for this iteration
                    thisTaxon = annots["taxons"][i]
                    thisKingdom = annots["kingdoms"][i]
                    thisEvalue = annots["evalues"][i]
                    thisAlignLen = annots["alignment_length"][i]
                    thisAlignOvl = thisAlignLen / seqLen
                    
                    # Handle animals
                    if thisKingdom == "a" and thisEvalue <= rescueEvalue and thisAlignOvl <= ovlCutoff:
                        isContaminant = False
                        break # we're going to rescue this sequence!
                    
                    # Handle contaminants
                    elif thisKingdom == "c" and thisEvalue <= contamEvalue and thisAlignOvl <= ovlCutoff:
                        isContaminant = True
                        contaminants.append(thisTaxon)
            
            # If we found no evidence of it being a contaminant, write to file
            if not isContaminant:
                fileOut.write(">{0}\n{1}\n".format(
                    record.description,
                    str(record.seq)
                ))
            # If we think it is a contaminant, get some info from the sequence
            else:
                contamCounts["num_contaminants"] += 1 # log the absolute number of sequences dropped
                
                # Store information from the best contaminant hit
                bestContamTaxon = contaminants[0]
                lineage = ncbi.get_lineage(bestContamTaxon)
                for lineageTaxonID in lineage:
                    contamCounts.setdefault(lineageTaxonID, 0)
                    contamCounts[lineageTaxonID] += 1
    return contamCounts

def main():
    usage = """%(prog)s receives a transcriptome and attempts to curate it of
    contaminant contigs using an annotation table produced by ZSK' scripts.
    
    Notably, this program assumes the transcriptome is of animal origin, and hence
    any microorganism hits will be considered to be contamination. Statistical
    breakdown of contaminant sequences will be presented to help identify
    whether there's any sort of pattern to the contaminants being identified.
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="transcriptomeFile",
                   required=True,
                   help="Specify the location of the input EVG transcriptome file")
    p.add_argument("-a", dest="annotationTable",
                   required=True,
                   help="Specify an annotation table file")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Specify where the curated transcriptome should be written")
    # Opts
    p.add_argument("--contam-evalue", dest="contamEvalue", type=float,
                   required=False,
                   help="""Optionally, specify an E-value to use for identifying
                   contaminant BLAST hits (default == 1e-5)""",
                   default=1e-5)
    p.add_argument("--rescue-evalue", dest="rescueEvalue", type=float,
                   required=False,
                   help="""Optionally, specify an E-value to use for identifying
                   non-contaminant BLAST hits that, if detected, will result in the
                   rescue of the sequence (default == 1e-5)""",
                   default=1e-5)
    p.add_argument("--overlap", dest="ovlCutoff", type=float,
                   required=False,
                   help="""Optionally, the proportion of the contig that must be
                   matched for a hit to be considered significant
                   (value from 0>1; default == 0.4)""",
                   default=0.4)
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse annotation file
    annotDict = parse_annots(args.annotationTable, max(args.contamEvalue, args.rescueEvalue))
    modify_annot_taxons(annotDict)
    
    # Load transcriptome records
    records = SeqIO.parse(open(args.transcriptomeFile, "r"), "fasta")
    
    # Curate contaminants from transcriptome
    contamCounts = contaminant_curate_transcriptome(records, annotDict,
        args.contamEvalue, args.rescueEvalue,
        args.ovlCutoff, args.outputFileName
    )
    
    # Tabulate statistical information regarding contamination in this transcriptome
    ncbi = NCBITaxa()
    uniqueTaxons = set([taxonID for taxonID in contamCounts.keys() if taxonID != "num_contaminants"])
    taxid2name = ncbi.get_taxid_translator(uniqueTaxons)
    
    namedContamCounts = [
        [taxid2name[taxonID], count]
            for taxonID, count in contamCounts.items()
                if taxonID != "num_contaminants"
    ]
    namedContamCounts.sort(key = lambda x: -x[1])
    
    with open(args.outputFileName + ".stats", "w") as fileOut:
        fileOut.write("#rank\tnum_contaminants\n") # header
        for taxonName, count in namedContamCounts:
            fileOut.write(f"{taxonName}\t{count}\n")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

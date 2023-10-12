#! python3
# gff3_rename_by_blast.py
# Uses the ZS_GFF3IO.GFF3 Class to parse a GFF3 file alongside a
# BLAST outfmt6 file to rename the gene models based on their best match.

import os, argparse, sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Function_packages import ZS_GFF3IO, ZS_BlastIO

PROBLEM_SEQ_BITS = ["PREDICTED: ", "LOW QUALITY PROTEIN: "]

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.gff3File):
        print('I am unable to locate the GFF3 file (' + args.gff3File + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.blastFile):
        print('I am unable to locate the BLAST outfmt6 file (' + args.blastFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if args.idMapFile != None:
        if not os.path.isfile(args.idMapFile):
            print('I am unable to locate the ID mapping file (' + args.idMapFile + ')')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print('File already exists at output location (' + args.outputFileName + ')')
        print('Make sure you specify a unique file name and try again.')
        quit()

def tsv_to_dict(textFile):
    outDict = {}
    with open(textFile, 'r') as fileIn:
        for line in fileIn:
            if line.startswith("#"):
                continue
            
            sl = line.rstrip("\r\n ").split("\t")
            assert len(sl) == 2, \
                f"{textFile} does not have two columns on line '{line}'"
            outDict[sl[0]] = sl[1]
    return outDict

def main():
    ##### USER INPUT SECTION
    usage = """%(prog)s reads in a GFF3 file and an outfmt6 file to rename the
    gene features in a GFF3 according to the best hit. An additional layer of
    functionality comes from the idMapFile which can be used to relate the best
    hit to a more meaningful name.
    """
    # Required
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-g", dest="gff3File",
                   required=True,
                   help="Input GFF3 file")
    p.add_argument("-b", dest="blastFile",
                   required=True,
                   help="Input BLAST outfmt6 file")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file name.")
    # Optional
    p.add_argument("--idmap", dest="idMapFile",
                   required=False,
                   help="""Optionally, specify a TSV pairing GFF3 mRNA IDs to
                   the name you'd rather associate with that hit""")
    
    args = p.parse_args()
    validate_args(args)
    
    # Load GFF3 file
    gff3Obj = ZS_GFF3IO.GFF3(args.gff3File, strict_parse=False)
    
    # Load BLAST file
    blastResultsObj = ZS_BlastIO.BLAST_Results(args.blastFile)
    blastResultsObj.num_hits = 1
    blastResultsObj.parse_blast_hit_coords()
    blastResultsDict = blastResultsObj.results
    
    # Load ID mapping if relevant
    if args.idMapFile != None:
        idMapDict = tsv_to_dict(args.idMapFile)
    else:
        idMapDict = {}
    
    # Update GFF3 based on BLAST hits
    handledIDs = set()
    idsCount = {}
    for mRNA_ID, blastHitList in blastResultsDict.items():
        bestBlast = blastHitList[0]
        try:
            geneFeature = gff3Obj[gff3Obj[mRNA_ID].Parent]
        except KeyError:
            try:
                assert gff3Obj[mRNA_ID].Parent.rsplit("_", maxsplit=1)[0] in idsCount
            except AssertionError:
                print(gff3Obj[mRNA_ID].Parent)
                print(idsCount)
                print("##")
                print("Major issue encountered during GFF3 update")
                print(f"{mRNA_ID} couldn't be found in your GFF3")
                print("I expected this based on your BLAST input. You should try to figure out what's going wrong...")
                print("Sorry. Program exit now.")
                quit()
        
        # Skip if we've looked at this gene feature already
        if geneFeature.ID in handledIDs:
            continue
        handledIDs.add(geneFeature.ID) # now store it to prevent future runs
        
        # Get best BLAST for all mRNAs of this gene
        for mrnaFeature in geneFeature.mRNA:
            if mrnaFeature.ID in blastResultsDict:
                if blastResultsDict[mrnaFeature.ID][0][-1] < bestBlast[-1]:
                    bestBlast = blastResultsDict[mrnaFeature.ID][0]
        bestID = bestBlast[0]
        
        # Get mapped ID if relevant
        bestID = bestID if not bestID in idMapDict else idMapDict[bestID]
        
        # Fix it up to be nice and all
        for problemBit in PROBLEM_SEQ_BITS:
            if problemBit in bestID:
                bestID = bestID.replace(problemBit, "")
        
        bestID = bestID.replace(" ", "_") # no empty spaces in gene ID!
        
        # Append a count to the gene ID to prevent identical IDs
        idsCount.setdefault(bestID, 0)
        idsCount[bestID] += 1
        bestID += f"_{idsCount[bestID]}"
        
        # Update this gene feature and all child features
        geneFeature.update_id(geneFeature.ID, bestID)
    
    # Write output file
    with open(args.outputFileName, "w") as fileOut:
        for geneFeature in gff3Obj.types["gene"]:
            fileOut.write(geneFeature.format_as_gff3() + "\n")
            for mrnaFeature in geneFeature.mRNA:
                fileOut.write(mrnaFeature.format_as_gff3() + "\n")
                for childFeature in mrnaFeature.retrieve_all_children():
                    fileOut.write(childFeature.format_as_gff3() + "\n")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

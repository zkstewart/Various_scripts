#! python3
# prep_edaseq_norm.py
# Script to generate a file that is necessary for EDAseq
# normalisation to commence.

# Import necessary packages
import os, argparse, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))) # 3 dirs up is where we find dependencies
from Function_packages import ZS_SeqIO, ZS_GFF3IO

# Define script-level functions
def validate_args(args):  
    # Validate input file location
    if not os.path.isfile(args.gff3File):
        print("I am unable to locate the GFF3 file ({0})".format(args.gff3File))
        print("Make sure you\'ve typed the file name or location correctly and try again.")
        quit()
    if not os.path.isfile(args.fastaFile):
        print("I am unable to locate the FASTA file ({0})".format(args.fastaFile))
        print("Make sure you\'ve typed the file name or location correctly and try again.")
        quit()
    if args.mappingFile != None:
        if not os.path.isfile(args.mappingFile):
            print("I am unable to locate the mapping file ({0})".format(args.mappingFile))
            print("Make sure you\'ve typed the file name or location correctly and try again.")
            quit()
    
    # Validate output file
    if os.path.isfile(args.outputFileName):
        print("Specifed output file already exists. This program does not allow file overwriting.")
        print("Specify a new name or move the existing file and try again.")
        quit()

def parse_mapping_file(mappingFile):
    mappingDict = {}
    with open(mappingFile, "r") as fileIn:
        for line in fileIn:
            line = line.strip("\r\n ")
            if line.startswith("#") or line == "":
                continue
            
            sl = line.split("\t")
            mappingDict[sl[0]] = sl[1]
    return mappingDict

def create_edaseq_gene_details(GFF3_obj, FASTA_obj, mappingDict=None):
    edaseqTableDict = {} # as a dict so we get one entry per mappingDict value (if applicable)
    for FastASeq_obj in FASTA_obj:
        # Get our sequence ID
        if mappingDict != None:
            if FastASeq_obj.id not in mappingDict:
                continue
            else:
                seqID = mappingDict[FastASeq_obj.id]
        else:
            seqID = FastASeq_obj.id
        
        # Get our GFF3 entry
        assert seqID in GFF3_obj, \
            "'{0}' sequence not found in GFF3; files aren't matched, can't proceed".format(seqID)
        GFF3_feature = GFF3_obj[seqID]
        assert GFF3_feature.type in ["gene", "mRNA", "lnc_RNA"], \
            "'{0}' sequence is not a gene/mRNA/lnc_RNA, can't proceed".format(seqID)
        
        # Figure out if our sequence is exon or CDS
        # exonLengths = []
        # cdsLengths = []
        # RNA_features = [GFF3_feature] if GFF3_feature.type == "mrna" else GFF3_feature.types["mRNA"]
        
        # for RNA_feature in RNA_features:
        #     assert "CDS" in RNA_feature.types and "exon" in RNA_feature.types, \
        #         "'{0}' mRNA feature does not contain CDS and exon types, can't proceed".format(RNA_feature.ID)
        #     exonLengths.append(0)
        #     cdsLengths.append(0)
            
        #     try: # this try-except clause will catch lnc_RNA features
        #         for cds_feature in RNA_feature.types["CDS"]:
        #             cdsLengths[-1] += cds_feature.end - cds_feature.start + 1
        #     except:
        #         pass
        #     for exon_feature in RNA_feature.types["exon"]:
        #         exonLengths[-1] += exon_feature.end - exon_feature.start + 1
        
        # exonDistances = [max(len(FastASeq_obj.seq), l) - min(len(FastASeq_obj.seq), l) for l in exonLengths]
        # cdsDistances = [max(len(FastASeq_obj.seq), l) - min(len(FastASeq_obj.seq), l) for l in cdsLengths]
        # isCDS = True if min(cdsDistances) <= min(exonDistances) else False
        
        # Calculate the GC content for this feature
        gcCount = [1 if base in ["C", "G"] else 0 for base in FastASeq_obj.seq.upper()]
        gcPct = sum(gcCount) / len(gcCount)
        
        # Figure out mRNA feature to use (if relevant)
        # mrnaIndex = exonLengths.index(max(exonLengths))
        
        # Format details for this gene/mRNA
        "metaSeqR2 wants a dumb format so we just have to give it to them"
        detailsTable = []
        # iterType = "CDS" if isCDS else "exon"
        # for feature in mRNA_features[mrnaIndex].types[iterType]:
        #     detailsTable.append([
        #         feature.contig, feature.start, feature.end,
        #         seqID, gcPct
        #     ])
        start = 1 # fake start since it's only used for length calculations
        end = len(FastASeq_obj.seq)
        detailsTable.append([
            GFF3_feature.contig, start, end,
            seqID, gcPct
        ])
        
        # Store in dictionary
        edaseqTableDict[seqID] = detailsTable
    return edaseqTableDict

def main():
    usage = """%(prog)s will produce a table that can be used by EDAseq to perform
    read normalisation. It requires a GFF3 input and the corresponding FASTA file
    """
    # Required
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-g", dest="gff3File", required=True,
        help="Specify GFF3 file")
    p.add_argument("-f", dest="fastaFile", required=True,
        help="Specify FASTA file")
    p.add_argument("-o", dest="outputFileName", required=True,
        help="output file name")
    # Optional
    p.add_argument("-m", dest="mappingFile", required=False,
        help="Optionally, provide a mapping file relating FASTA IDs (left) to GFF3 IDS (right)")
    args = p.parse_args()
    validate_args(args)
    
    # Parse GFF3
    GFF3_obj = ZS_GFF3IO.GFF3(args.gff3File, strict_parse = False)
    
    # Parse FASTA
    FASTA_obj = ZS_SeqIO.FASTA(args.fastaFile)
    
    # Parse mapping file (if applicable)
    if args.mappingFile != None:
        mappingDict = parse_mapping_file(args.mappingFile)
    else:
        mappingDict = None
    
    # Generate EDAseq gene details table
    edaseqTableDict = create_edaseq_gene_details(GFF3_obj, FASTA_obj, mappingDict)
    
    # Write output
    with open(args.outputFileName, "w") as fileOut:
        fileOut.write("{0}\n".format("\t".join([
            "chromosome", "start", "end",
            "gene_id", "gc_content"
        ])))
        for tableRows in edaseqTableDict.values():
            for row in tableRows:
                fileOut.write("{0}\n".format("\t".join(map(str, row))))
    
    print('Program completed successfully!')

if __name__ == "__main__":
    main()

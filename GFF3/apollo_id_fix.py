#! python3
# apollo_id_fix.py
# Script to receive a GFF3 obtained from Apollo and update
# IDs to not look terrible. It will pretty up the GFF3
# format too, akin to how PASA does things.

import os, sys, argparse
from Bio import SeqIO

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) # 2 dirs up is where we find dependencies
from Function_packages import ZS_GFF3IO

def validate_args(args):
    # Validate input data locations
    if not os.path.isfile(args.gff3File):
        print('I am unable to locate the input GFF3 file (' + args.gff3File + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.fastaFile):
        print('I am unable to locate the input genome FASTA file (' + args.fastaFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Handle file output
    if os.path.exists(args.outputGff3):
        print('The specified output GFF3 already exists. This program will not allowing overwriting.')
        print('Move the existing file, or specify a new file name and try again.')

def main():
    usage = """%(prog)s receives a GFF3 produced by Apollo, alongside a FASTA
    containing it's protein translations for mRNA sequences, and produces an updated
    GFF3 with better gene IDs and PASA-like pretty formatting.
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-g", dest="gff3File",
                   required=True,
                   help="Specify the location of the input GFF3 file")
    p.add_argument("-f", dest="fastaFile",
                   required=True,
                   help="Specify the location of the input FASTA file")
    p.add_argument("-o", dest="outputGff3",
                   required=True,
                   help="Output location for modified GFF3")
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse GFF3
    gff3Object = ZS_GFF3IO.GFF3(args.gff3File)
    
    # Parse FASTA
    records = SeqIO.to_dict(SeqIO.parse(args.fastaFile, "fasta"))
    
    # Create updated apollo annotation
    outputLines = []
    contigCount = {}
    for geneFeature in gff3Object.types["gene"]:
        contigCount.setdefault(geneFeature.contig, 1)
        
        # Write comment line demarcating gene
        geneID = f"{geneFeature.contig}_manual_{contigCount[geneFeature.contig]}"
        outputLines.append(f"# APOLLO_MANUAL: {geneID}; num_mRNA={len(geneFeature.mRNA)}")
        
        # Create line for gene-level feature
        geneLine = "\t".join(map(str, [
            geneFeature.contig, geneFeature.source, geneFeature.type,
            geneFeature.start, geneFeature.end, geneFeature.score,
            geneFeature.strand, geneFeature.frame,
            f"ID={geneID};Name=apollo_{geneID};date_creation={geneFeature.date_creation}"
        ]))
        outputLines.append(geneLine)
        contigCount[geneFeature.contig] += 1
        
        # Create lines for mRNA-level feature
        mrnaCount = 1
        for mrnaFeature in geneFeature.mRNA:
            mrnaID = f"{geneID}.mrna{mrnaCount}" # geneID carries over
            mrnaLine = "\t".join(map(str, [
                mrnaFeature.contig, mrnaFeature.source, mrnaFeature.type,
                mrnaFeature.start, mrnaFeature.end, mrnaFeature.score,
                mrnaFeature.strand, mrnaFeature.frame,
                f"ID={mrnaID};Parent={geneID};Name=apollo_{mrnaID};date_creation={mrnaFeature.date_creation}"
            ]))
            outputLines.append(mrnaLine)
            mrnaCount += 1
            
            # Sort child features appropriately
            if geneFeature.strand == "+":
                exonFeatures = sorted(mrnaFeature.exon, key = lambda x: (x.start, x.end))
                cdsFeatures = sorted(mrnaFeature.CDS, key = lambda x: (x.start, x.end))
            else:
                exonFeatures = sorted(mrnaFeature.exon, key = lambda x: (-x.end, -x.start))
                cdsFeatures = sorted(mrnaFeature.CDS, key = lambda x: (-x.end, -x.start))
            
            # Create lines for exon-level feature
            exonCount = 1
            for exonFeature in exonFeatures:
                exonID = f"{mrnaID}.exon{exonCount}"
                exonLine = "\t".join(map(str, [
                    exonFeature.contig, exonFeature.source, exonFeature.type,
                    exonFeature.start, exonFeature.end, exonFeature.score, 
                    exonFeature.strand, exonFeature.frame,
                    f"ID={exonID};Parent={mrnaID};Name=apollo_{exonID}"
                ]))
                outputLines.append(exonLine)
                exonCount += 1
            
            # Create lines for CDS-level feature
            cdsCount = 1
            for cdsFeature in cdsFeatures:
                cdsID = f"{mrnaID}.cds{cdsCount}"
                cdsLine = "\t".join(map(str, [
                    cdsFeature.contig, cdsFeature.source, cdsFeature.type,
                    cdsFeature.start, cdsFeature.end, cdsFeature.score,
                    cdsFeature.strand, cdsFeature.frame,
                    f"ID={cdsID};Parent={mrnaID};Name=apollo_{cdsID}"
                ]))
                outputLines.append(cdsLine)
                cdsCount += 1
        
        # Write comment line containing sequence
        mrnaCount = 1
        for mrnaFeature in geneFeature.mRNA:
            mrnaID = f"{geneID}.mrna{mrnaCount}"
            
            seq = str(records[mrnaFeature.ID].seq)
            outputLines.append(f"# PROT {mrnaID}\t{seq}")
            mrnaCount += 1
    
    # Write artifact lines to file
    artifactFeatures = []
    for typeKey in gff3Object.types.keys():
        if "artifact" in typeKey:
            for artifactFeature in gff3Object.types[typeKey]:
                artifactFeatures.append(artifactFeature)
    
    if len(artifactFeatures) > 0:
        artifactFeatures.sort(key = lambda x: (x.contig, x.start, x.end))
        
        outputLines.append(f"# INDEL ARTIFACTS; num={len(artifactFeatures)}")
        
        contigCount = {}
        for artifactFeature in artifactFeatures:
            contigCount.setdefault(geneFeature.contig, 1)
            artifactID = f"{artifactFeature.contig}_artifact_{contigCount[artifactFeature.contig]}"
            
            attributes = f"ID={artifactID};Name=apollo_{artifactID};" + \
                    f"justification={artifactFeature.justification};date_creation={artifactFeature.date_creation}"
            if hasattr(artifactFeature, "residues"):
                attributes += f";residues={artifactFeature.residues}"
            
            artifactLine = "\t".join(map(str, [
                artifactFeature.contig, artifactFeature.source, artifactFeature.type,
                artifactFeature.start, artifactFeature.end, artifactFeature.score,
                artifactFeature.strand, artifactFeature.frame,
                attributes
            ]))
            outputLines.append(artifactLine)
            contigCount[geneFeature.contig] += 1
    
    # Write file
    with open(args.outputGff3, "w") as fileOut:
        fileOut.write("\n".join(outputLines))
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

#! python3

import sys
from Bio import SeqIO

sys.path.append(r"C:\git\Various_scripts")
from Function_packages import ZS_GFF3IO

gff3 = "Annotations-F1_Ma8D_418_trimmed_to_Markers.final.gff3"
fasta = "manual_418_models.aa"
out = "manual_418_models.gff3"

# Parse GFF3
gff3Object = ZS_GFF3IO.GFF3(gff3)

# Parse FASTA
records = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))

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
        geneFeature.start, geneFeature.end, geneFeature.strand,
        geneFeature.score, geneFeature.frame,
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
            mrnaFeature.start, mrnaFeature.end, mrnaFeature.strand,
            mrnaFeature.score, mrnaFeature.frame,
            f"ID={mrnaID};Parent={geneID};Name=apollo_{mrnaID};date_creation={mrnaFeature.date_creation}"
        ]))
        outputLines.append(mrnaLine)
        mrnaCount += 1
        
        # Create lines for exon-level feature
        exonCount = 1
        for exonFeature in mrnaFeature.exon:
            exonID = f"{mrnaID}.exon{exonCount}"
            exonLine = "\t".join(map(str, [
                exonFeature.contig, exonFeature.source, exonFeature.type,
                exonFeature.start, exonFeature.end, exonFeature.strand,
                exonFeature.score, exonFeature.frame,
                f"ID={exonID};Parent={mrnaID};Name=apollo_{exonID}"
            ]))
            outputLines.append(exonLine)
            exonCount += 1
        
        # Create lines for CDS-level feature
        cdsCount = 1
        for cdsFeature in mrnaFeature.CDS:
            cdsID = f"{mrnaID}.cds{cdsCount}"
            cdsLine = "\t".join(map(str, [
                cdsFeature.contig, cdsFeature.source, cdsFeature.type,
                cdsFeature.start, cdsFeature.end, cdsFeature.strand,
                cdsFeature.score, cdsFeature.frame,
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
        if artifactFeature.type == "insertion_artifact":
            attributes += f";residues={artifactFeature.residues}"
        
        artifactLine = "\t".join(map(str, [
            artifactFeature.contig, artifactFeature.source, artifactFeature.type,
            artifactFeature.start, artifactFeature.end, artifactFeature.strand,
            artifactFeature.score, artifactFeature.frame,
            attributes
        ]))
        outputLines.append(artifactLine)
        contigCount[geneFeature.contig] += 1

# Write file
with open(out, "w") as fileOut:
    fileOut.write("\n".join(outputLines))

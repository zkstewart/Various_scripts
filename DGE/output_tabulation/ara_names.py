#! python3

import os

def blasttab_best_hits(blastTab, evalue, numHits, SUFFIX_TRIM=False):
        from itertools import groupby
        # Preliminary declaration of grouper function & main dictionary to hold onto the alignment details of all selected hits
        grouper = lambda x: x.split('\t')[0]
        outDict = {}
        # Iterate through the blastTab file looking at each group of hits to a single query sequence
        with open(blastTab, 'r') as fileIn:
                for key, value in groupby(fileIn, grouper):
                        # Make the group value amenable to multiple iterations, and sort it to present most significant first, least significant last
                        value = list(value)
                        for i in range(len(value)):
                                value[i] = value[i].rstrip('\n').rstrip('\r').split('\t')
                        value.sort(key = lambda x: (float(x[10]),-float(x[11])))        # Sorts by E-value (lowest first) and bitscore (highest first)
                        # Pull out the [X] best hits
                        bestHits = []
                        for val in value:
                                if float(val[10]) <= evalue and len(bestHits) < numHits:
                                        # Alter the target name if it has UniRef prefix
                                        if val[1].startswith('UniRef'):
                                                val[1] = val[1].split('_')[1]           # This handles normal scenarios ("UniRef100_UPI0000") as well as MMseqs2 weird ID handling ("UniRef100_UPI0000_0")
                                                bestHits.append(val)
                                        else:
                                                bestHits.append(val)
                        # Skip this hit if we found no matches which pass E-value cut-off
                        if bestHits == []:
                                continue
                        # Dig into the hits and retrieve our best hit which has a mapping in the idmapping file
                        bestMapped = '.'
                        for val in value:
                                if val[1].startswith('UniRef'):
                                        val[1] = val[1].split('_', maxsplit=1)[1]
                                if float(val[10]) <= evalue:
                                        bestMapped = val[1] + ' (' + val[10] + ')'
                                        break
                        # Process line to format it for output
                        formattedList = []
                        for i in range(len(bestHits)):
                                if i == 0:
                                        for x in range(1, 12):
                                                formattedList.append([bestHits[i][x]])
                                else:
                                        for x in range(1, 12):
                                                formattedList[x-1].append('[' + bestHits[i][x] + ']')
                        for i in range(len(formattedList)):
                                formattedList[i] = ''.join([formattedList[i][0], ' ', *formattedList[i][1:]])
                        # Add our best mapped value onto the end & save to our outDict
                        if SUFFIX_TRIM is True:
                            index = bestHits[0][0].rsplit(".", maxsplit=1)[0]
                            if index == "chr4.comp19742_c0_seq1":
                                index = "chr4.comp19742_c0_seq1.path1"
                            outDict[index] = formattedList + [bestMapped]
                        else:
                            outDict[bestHits[0][0]] = formattedList + [bestMapped]
        return outDict

def parse_araport_gff(gff3FileName):
    araportDict = {}
    with open(gff3FileName, "r") as fileIn:
        for line in fileIn:
            # Skip comments
            if line.startswith("#"):
                continue
            sl = line.rstrip("\r\n").split("\t")
            # Skip irrelevant lines
            if sl[2] != "mRNA" and 'ATMG' not in sl[8]: # ATMG accounts for mitochondrial genome genes
                continue
            # Get the relevant attributes
            note = [attribute[5:] for attribute in sl[8].split(";") if attribute.startswith("Note")]
            if note == []:
                if "product=" in sl[8]:
                    note = [attribute[8:] for attribute in sl[8].split(";") if attribute.startswith("product")]
                else:
                    continue
            else:
                note = note[0]
            mrnaID = [attribute[3:] for attribute in sl[8].split(";") if attribute.startswith("ID")][0]
            # Save to dict
            araportDict[mrnaID] = note
    return araportDict

## HARD-CODED WORKING
idsFolders = ["timecourse_DE", "pairwise_DE", "venn_timecourse_DE"]
tairBlast = r"F:\plant_group\plant_rnaseq\annotations\mango\mango_tair10_mms2SEARCH_sorted.m8"
tairGFF3 = r"F:\plant_group\plant_rnaseq\annotations\databases\Araport11_GFF3_genes_transposons.Mar92021.gff"
outputFileNames = ["tc_dge_tair_columns.txt", "pw_dge_tair_columns.txt", "venn_dge_tair_columns.txt"]
evalue = 0.001

# Parse the IDs files to get a parent list of IDs to get TAIR names for
geneIDsList = []
for folder in idsFolders:
    geneIDsList.append(set())
    headerLines = 1
    
    for file in os.listdir(folder):
        if file.endswith("_results.tsv"):
            file = os.path.join(folder, file)
            with open(file, "r") as fileIn:
                for line in fileIn:
                    if headerLines > 0:
                        headerLines -= 1
                        continue
                    sl = line.rstrip("\r\n").split("\t")
                    geneIDsList[-1].add(sl[0])
        elif file.endswith("_results.ids"):
            file = os.path.join(folder, file)
            with open(file, "r") as fileIn:
                for line in fileIn:
                    l = line.rstrip(" \r\n")
                    geneIDsList[-1].add(l)

# Parse the blast-tab file
outDict = blasttab_best_hits(tairBlast, evalue, 1, SUFFIX_TRIM=True)

# Parse the TAIR10 annotation GFF3
araportDict = parse_araport_gff(tairGFF3)

# Hard-coded dict for handling obsoleted genes
obsoleteDict = {
    "AT2G25050.2": "Actin-binding FH2 (Formin Homology) protein",
    "AT2G45770.2": "chloroplast SRP receptor homolog",
    "ATMG00860.1": "Uncharacterized mitochondrial protein",
    "ATMG00580.1": "NADH dehydrogenase subunit 4",
    "AT5G50011.1": "Conserved upstream opening reading frame relative to major ORF AT5G50010.1",
    "AT1G27590.1": "Unknown protein",
    "AT2G33440.1": "RNA-binding (RRM/RBD/RNP motifs) family protein",
    "AT3G09470.1": "Major facilitator superfamily protein",
    "ATMG01250.1": "RNA-directed DNA polymerase (reverse transcriptase)",
    "ATMG00060.1": "Mitochondrial NADH dehydrogenase subunit 5",
    "ATMG01280.1": "Encodes a cytochrome c oxidase subunit II",
    "AT4G16360.3": "5'-AMP-activated protein kinase beta-2 subunit protein",
    "ATMG00300.1": "Gag-Pol-related retrotransposon family protein",
    "ATMG00285.1": "Nadh-ubiquinone oxidoreductase chain 2",
    "ATMG00660.1": "Hypothetical protein",
    "AT2G36480.3": "Pre-mRNA cleavage complex 2 Pcf11-like protein",
    "AT5G25752.1": "Arabidopsis rhomboid-like protein 11",
    "AT2G44270.2": "Cytoplasmic tRNA 2-thiolation protein 1",
    "AT2G33430.1": "Multiple organellar rna editing factor 2",
    "ATMG00860.1": "Uncharacterized mitochondrial protein",
    "AT2G25050.2": "Actin-binding FH2 (Formin Homology) protein",
    "ATMG00810.1": "Uncharacterized mitochondrial protein",
    "AT4G16360.3": "5'-AMP-activated protein kinase beta-2 subunit protein",
    "ATMG00710.1": "Polynucleotidyl transferase, ribonuclease H-like superfamily protein",
    "AT2G44270.2": "Encodes ROL5, a repressor of lrx1 mutants that develop aberrant root hairs",
    "AT3G09470.1": "Major facilitator superfamily protein",
    "AT5G18280.2": "Apyrase 2",
    "ATMG00910.1": "cytochrome B/B6 protein",
    "AT2G33440.1": "RNA-binding (RRM/RBD/RNP motifs) family protein",
    "ATMG00120.1": "Uncharacterized mitochondrial protein",
    "ATMG00750.1": "Uncharacterized mitochondrial protein",
    "AT1G27590.1": "CONTAINS InterPro DOMAIN/s: Protein of unknown function DUF3453",
    "ATMG00580.1": "NADH dehydrogenase subunit 4",
    "ATMG00060.1": "Mitochondrial NADH dehydrogenase subunit 5",
    "ATMG00285.1": "Encodes a subunit of mitochondrial NAD(P)H dehydrogenase that is trans-spliced from two precursors",
    "ATMG00410.1": "ATPase subunit 6",
    "AT2G42240.2": "RNA-binding (RRM/RBD/RNP motifs) family protein",
    "AT2G36480.3": "ENTH/VHS family protein",
    "AT2G33430.1": "Multiple organellar rna editing factor 2",
    "AT3G54670.3": "Structural maintenance of chromosomes (SMC) family protein",
    "ATMG00850.1": "Uncharacterized mitochondrial protein",
    "AT2G45770.2": "chloroplast SRP receptor homolog, alpha subunit CPFTSY",
    "AT5G25752.1": "Chloroplast-localized rhomboid-like protein"
}

# Loop through ID file to rename the genes (if applicable), order the output appropriately, and identify gaps in the BLAST-tab file
for i in range(len(outputFileNames)):
    geneIDs = geneIDsList[i]
    outputFileName = outputFileNames[i]
    with open(outputFileName, 'w') as fileOut:
        fileOut.write("#geneID\ttair10_id\ttair10_name\ttair10_link\n")
        for geneID in geneIDs:
            # Obtain relevant details from BLAST-tab dictionary
            if geneID in outDict:
                values = outDict[geneID]
                tair10ID = values[0].rstrip(" ")
                tair10Evalue = float(values[9].rstrip(" "))
                if tair10Evalue > evalue:
                    tair10ID = "."
                    tair10Name = "."
                    tair10Evalue = "."
                    tair10Link = "."
                else:
                    try:
                        tair10Name = araportDict[tair10ID]
                    except:
                        tair10Name = obsoleteDict[tair10ID]
                    tair10Link = r"https://www.arabidopsis.org/servlets/TairObject?type=gene&name={0}".format(tair10ID)
            else:
                tair10ID = "."
                tair10Name = "."
                tair10Evalue = "."
                tair10Link = "."
            # Fix list values??
            if isinstance(tair10Name, list):
                tair10Name = tair10Name[0]
            # Replace URL characters
            tair10Name = tair10Name.replace("%2C", ",").replace("%3B", ";")
            # Write results to file
            fileOut.write("{0}\t{1}\t{2}\t{3}\n".format(
                geneID, tair10ID, tair10Name, tair10Link
            ))
# Done!
print('Program completed successfully!')

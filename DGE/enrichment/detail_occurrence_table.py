#! python3

from goatools import obo_parser

species = ["act", "exa", "nem", "par", "tel"]
goSuffix = ".GOs.chisq_results.tsv"
pfamSuffix = ".PFAMs.chisq_results.tsv"

newGoSuffix = ".GOs.chisq_results.detailed.tsv"
newPfamSuffix = ".PFAMs.chisq_results.detailed.tsv"

# Parse GO obo file
oboFile = r"F:\anemone_expansions\go-basic.obo"
goOboDict = obo_parser.GODag(oboFile)

# Parse PFAM clans file
clansFile = r"F:\anemone_expansions\Pfam-A.clans.tsv"
pfamDict = {}
with open(clansFile, "r") as fileIn:
    for line in fileIn:
        sl = line.rstrip("\r\n").split("\t")
        pfamID, clanID, _, shortName, longName = sl
        pfamDict[shortName] = [longName, pfamID]

# Parse and update GO tables
PROBLEM_GOS = {
    "GO:0000187": "obsolete activation of MAPK activity",
    'GO:0042766': "obsolete nucleosome mobilization",
    'GO:0007257': "obsolete activation of JUN kinase activity",
    'GO:0031936': "obsolete negative regulation of chromatin silencing",
}

for sp in species:
    firstLine = True
    with open(sp + goSuffix, "r") as fileIn, open(sp + newGoSuffix, "w") as fileOut:
        for line in fileIn:
            sl = line.rstrip("\r\n").split("\t")
            if firstLine == True:
                sl.insert(1, "name")
                firstLine = False
            else:
                goTerm = sl[0]
                try:
                    goName = goOboDict[goTerm].name
                except:
                    goName = PROBLEM_GOS[goTerm]
                sl.insert(1, goName)
            fileOut.write("{0}\n".format("\t".join(sl)))

# Parse and update PFAM tables
for sp in species:
    firstLine = True
    with open(sp + pfamSuffix, "r") as fileIn, open(sp + newPfamSuffix, "w") as fileOut:
        for line in fileIn:
            sl = line.rstrip("\r\n").split("\t")
            if firstLine == True:
                sl = ["PFAM_ID", "PFAM_name", "PFAM_description", "P_value", "represented"]
                firstLine = False
            else:
                shortName = sl[0]
                longName, pfamID = pfamDict[shortName]
                sl = [pfamID, shortName, longName, sl[1], sl[2]]
            fileOut.write("{0}\n".format("\t".join(sl)))

print("Program completed successfully!")
#! python3

from scipy.stats import chi2_contingency
from multipy.fdr import lsu
import numpy as np

species = ["act", "exa", "nem", "par", "tel"]
goSuffix = ".goseq.gos"
pfamSuffix = ".goseq.domains"

gosDict = {"all": {}}
pfamsDict = {"all": {}}
geneNumber = {}

P_VAL_CUTOFF = 0.05

# Get domain/GO occurrence data
for sp in species:
    gosDict[sp] = {}
    geneNumber[sp] = 0
    with open(sp + goSuffix, "r") as fileIn:
        for line in fileIn:
            model, gos = line.rstrip("\r\n").split("\t")
            geneNumber[sp] += 1
            if gos == "0":
                continue
            for go in gos.split("; "):
                gosDict[sp].setdefault(go, 0) # for tracking species occurrence
                gosDict["all"].setdefault(go) # for keeping a list of all GOs annotated
                gosDict[sp][go] += 1
    
    pfamsDict[sp] = {}
    with open(sp + pfamSuffix, "r") as fileIn:
        for line in fileIn:
            model, domains = line.rstrip("\r\n").split("\t")
            if domains == "0":
                continue
            for dom in domains.split("; "):
                pfamsDict[sp].setdefault(dom, 0) # for tracking species occurrence
                pfamsDict["all"].setdefault(dom) # for keeping a list of all domains annotated
                pfamsDict[sp][dom] += 1

# Perform statistical testing
goResults = []
pfamResults = []
for sp in species:
    spGoResults = []
    for go in gosDict["all"]:
        # Get this species' count data
        spOccurring = 0 if go not in gosDict[sp] else gosDict[sp][go]
        spNotOccurring = geneNumber[sp] - spOccurring
        # Get all other species' count data
        otherOccurring, otherNotOccurring = 0, 0
        for otherSp in species:
            if otherSp == sp:
                continue
            else:
                _occurring = 0 if go not in gosDict[otherSp] else gosDict[otherSp][go]
                _notOccurring = geneNumber[otherSp] - _occurring
                otherOccurring += _occurring
                otherNotOccurring += _notOccurring
        # Make a 2x2 contingency table and run chi2 analysis
        contingencyTable = [
            [spOccurring, spNotOccurring],
            [otherOccurring, otherNotOccurring]
        ]
        chi2, p, dof, expected = chi2_contingency(contingencyTable)
        # Figure out if it's over or underrepresented
        representation = "over" if expected[0][0] < spOccurring else "under" if expected[0][0] > spOccurring else "neither"
        # Store result
        spGoResults.append([go, p, representation])
    goResults.append(spGoResults)
    
    spPfamResults = []
    for dom in pfamsDict["all"]:
        # Get this species' count data
        spOccurring = 0 if dom not in pfamsDict[sp] else pfamsDict[sp][dom]
        spNotOccurring = geneNumber[sp] - spOccurring
        # Get all other species' count data
        otherOccurring, otherNotOccurring = 0, 0
        for otherSp in species:
            if otherSp == sp:
                continue
            else:
                _occurring = 0 if dom not in pfamsDict[otherSp] else pfamsDict[otherSp][dom]
                _notOccurring = geneNumber[otherSp] - _occurring
                otherOccurring += _occurring
                otherNotOccurring += _notOccurring
        # Make a 2x2 contingency table and run chi2 analysis
        contingencyTable = [
            [spOccurring, spNotOccurring],
            [otherOccurring, otherNotOccurring]
        ]
        chi2, p, dof, expected = chi2_contingency(contingencyTable)
        # Figure out if it's over or underrepresented
        representation = "over" if expected[0][0] < spOccurring else "under" if expected[0][0] > spOccurring else "neither"
        # Store result
        spPfamResults.append([dom, p, representation])
    pfamResults.append(spPfamResults)

# FDR correction
for i in range(len(goResults)):
    spGoResults = goResults[i]
    significant_pvals = lsu(np.array([p for _, p, _ in spGoResults]), q=P_VAL_CUTOFF)
    # Retain only significant results
    updatedResults = []
    for x in range(len(spGoResults)):
        if significant_pvals[x] == True:
            updatedResults.append(spGoResults[x])
    # Update the list now
    goResults[i] = updatedResults

for i in range(len(pfamResults)):
    spPfamResults = pfamResults[i]
    significant_pvals = lsu(np.array([p for _, p, _ in spPfamResults]), q=P_VAL_CUTOFF)
    # Retain only significant results
    updatedResults = []
    for x in range(len(spPfamResults)):
        if significant_pvals[x] == True:
            updatedResults.append(spPfamResults[x])
    # Update the list now
    pfamResults[i] = updatedResults

# Write nice output tables
for i in range(len(species)):
    sp = species[i]
    spGoResults = goResults[i]
    spPfamResults = pfamResults[i]
    # Sort results for better interpretation
    spGoResults.sort(key = lambda x: (x[2], x[1]))
    spPfamResults.sort(key = lambda x: (x[2], x[1]))
    # Write GO results
    with open(sp + ".GOs.chisq_results.tsv", "w") as fileOut:
        fileOut.write("GO_term\tP_value\trepresented\n")
        for result in spGoResults:
            fileOut.write("{0}\n".format("\t".join(map(str, result))))
    # Write PFAM results
    with open(sp + ".PFAMs.chisq_results.tsv", "w") as fileOut:
        fileOut.write("PFAM_domain\tP_value\trepresented\n")
        for result in spPfamResults:
            fileOut.write("{0}\n".format("\t".join(map(str, result))))

print("Program completed successfully!")
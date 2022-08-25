#! python3

import os, time

directories = ["pairwise_DE", "timecourse_DE", "venn_timecourse_DE"]

clusterGODir = r"timecourse_DE\clustering_goseq"

clusterResultsDir = "timecourse_DE\clustering"

clustersFileSuffix = "_patterns_and_clusters.tsv"

###

samplePrefixes = ["12.052.037", "Thai.wild", "Ampalam", "Keitt", "1243", "Laurina.Lombok"]

allVennIDsFile = r"venn_timecourse_DE\all_tc_results.ids"


################
#              #
# PARSING DATA #
#              #
################

# Parse results files to get a parent table of all DE genes & results for individual samples
annotTableDict = {}
annotSamplesDict = {}
for dir in directories:
    for file in os.listdir(dir):
        if file.endswith("_results.tsv"):
            fileName = os.path.join(dir, file)
            with open(fileName, "r") as fileIn:
                firstLine = True
                for line in fileIn:
                    if firstLine:
                        firstLine = False
                        continue
                    else:
                        sl = line.rstrip("\r\n").split("\t")
                        # Extract relevant details
                        geneID = sl[0]
                        isUpRegulated = True if float(sl[1]) > 0 else False
                        geneName = sl[4]
                        GOs, KEGGs = sl[7:9]
                        URL = sl[9]
                        sampleName = os.path.basename(fileName).split("_results")[0]
                        # Store details in relevant dicts
                        annotTableDict[geneID] = [geneName, GOs, KEGGs, URL]
                        annotSamplesDict.setdefault(sampleName, {})
                        annotSamplesDict[sampleName][geneID] = isUpRegulated

# Parse GO files to get a parent table of all DE genes & results for individual samples
goTableDict = {}
goSamplesDict = {}
for dir in directories + [clusterGODir]:
    for file in os.listdir(dir):
        if file.endswith("_GOseq.tsv"):
            fileName = os.path.join(dir, file)
            with open(fileName, "r") as fileIn:
                firstLine = True
                for line in fileIn:
                    if firstLine:
                        firstLine = False
                        continue
                    else:
                        sl = line.rstrip("\r\n").split("\t")
                        # Extract relevant details
                        goTerm = sl[0]
                        isOverRepresented = True if float(sl[1]) < float(sl[2]) else False
                        description = sl[5]
                        sampleName = os.path.basename(fileName).split("_GOseq")[0]
                        # Store details in relevant dicts
                        goTableDict[goTerm] = description
                        goSamplesDict.setdefault(sampleName, {})
                        goSamplesDict[sampleName][goTerm] = isOverRepresented

# Parse KEGG files to get a parent table of all DE genes & results for individual samples
keggTableDict = {}
keggSamplesDict = {}
for dir in directories  + [clusterGODir]:
    for file in os.listdir(dir):
        if file.endswith("_KEGGseq.tsv"):
            fileName = os.path.join(dir, file)
            with open(fileName, "r") as fileIn:
                firstLine = True
                sampleName = os.path.basename(fileName).split("_KEGGseq")[0]
                keggSamplesDict.setdefault(sampleName, {})
                for line in fileIn:
                    if firstLine:
                        firstLine = False
                        continue
                    else:
                        sl = line.rstrip("\r\n").split("\t")
                        # Extract relevant details
                        keggTerm = sl[0]
                        isOverRepresented = True if float(sl[1]) < float(sl[2]) else False
                        description = sl[5]
                        # Store details in relevant dicts
                        keggTableDict[keggTerm] = description
                        keggSamplesDict[sampleName][keggTerm] = isOverRepresented

# Parse clustering files for patterns and clusters
clusterIDsDict = {}
for file in os.listdir(clusterResultsDir):
    if file.endswith(clustersFileSuffix):
        fileName = os.path.join(clusterResultsDir, file)
        with open(fileName, "r") as fileIn:
            firstLine = True
            sampleName = os.path.basename(fileName).split(clustersFileSuffix)[0]
            clusterIDsDict.setdefault(sampleName, {})
            
            for line in fileIn:
                if firstLine is True:
                    firstLine = False
                    continue
                
                l = line.rstrip("\r\n ")
                if l == "":
                    continue
                else:
                    geneID, expressionValue, timeNum, pattern, clusterNum = l.split("\t")
                    if clusterNum == "NA":
                        clusterNum = None
                    clusterIDsDict[sampleName].setdefault(geneID, [pattern, clusterNum, {}])
                    clusterIDsDict[sampleName][geneID][2][int(timeNum)] = float(expressionValue)

# Parse the all_tc Venn grouping IDs
allVennIDs = []
with open(allVennIDsFile, "r") as fileIn:
    for line in fileIn:
        l = line.rstrip("\r\n ")
        if l == "": continue
        allVennIDs.append(l)


#################
#               #
# ORDERING DATA #
#               #
#################

##### SAMPLE ORDERS

tcSampleOrder = [
    'time1',
    'time2',
    'time3'
]

pairwiseSampleOrder = [
    '12.052.037_vs_Thai.wild_time1',
    '12.052.037_vs_Ampalam_time1',
    '12.052.037_vs_Keitt_time1',
    '12.052.037_vs_1243_time1',
    '12.052.037_vs_Laurina.Lombok_time1',
    'Thai.wild_vs_Ampalam_time1',
    'Thai.wild_vs_Keitt_time1',
    'Thai.wild_vs_1243_time1',
    'Thai.wild_vs_Laurina.Lombok_time1',
    'Ampalam_vs_Keitt_time1',
    'Ampalam_vs_1243_time1',
    'Ampalam_vs_Laurina.Lombok_time1',
    'Keitt_vs_1243_time1',
    'Keitt_vs_Laurina.Lombok_time1',
    '1243_vs_Laurina.Lombok_time1',
    '12.052.037_vs_Thai.wild_time2',
    '12.052.037_vs_Ampalam_time2',
    '12.052.037_vs_Keitt_time2',
    '12.052.037_vs_1243_time2',
    '12.052.037_vs_Laurina.Lombok_time2',
    'Thai.wild_vs_Ampalam_time2',
    'Thai.wild_vs_Keitt_time2',
    'Thai.wild_vs_1243_time2',
    'Thai.wild_vs_Laurina.Lombok_time2',
    'Ampalam_vs_Keitt_time2',
    'Ampalam_vs_1243_time2',
    'Ampalam_vs_Laurina.Lombok_time2',
    'Keitt_vs_1243_time2',
    'Keitt_vs_Laurina.Lombok_time2',
    '1243_vs_Laurina.Lombok_time2',
    '12.052.037_vs_Thai.wild_time3',
    '12.052.037_vs_Ampalam_time3',
    '12.052.037_vs_Keitt_time3',
    '12.052.037_vs_1243_time3',
    '12.052.037_vs_Laurina.Lombok_time3',
    'Thai.wild_vs_Ampalam_time3',
    'Thai.wild_vs_Keitt_time3',
    'Thai.wild_vs_1243_time3',
    'Thai.wild_vs_Laurina.Lombok_time3',
    'Ampalam_vs_Keitt_time3',
    'Ampalam_vs_1243_time3',
    'Ampalam_vs_Laurina.Lombok_time3',
    'Keitt_vs_1243_time3',
    'Keitt_vs_Laurina.Lombok_time3',
    '1243_vs_Laurina.Lombok_time3'
]

clusterSuffixOrder = [
    '_cluster_1',
    '_cluster_2',
    '_cluster_3',
    '_cluster_4',
    '_cluster_5',
    '_cluster_6',
    '_time_1_direction_D',
    '_time_1_direction_U',
    '_time_2_direction_D',
    '_time_2_direction_U',
    '_time_3_direction_D',
    '_time_3_direction_U'
]

vennOrder = [
    '12.052.037_tc_only',
    'Thai.wild_tc_only',
    'Ampalam_tc_only',
    'Keitt_tc_only',
    '1243_tc_only',
    'Laurina.Lombok_tc_only',
    'all_tc_only'
]

##### GENE ORDERS

pwGeneOrder = set()
for key, value in annotSamplesDict.items():
    if "tc" in key:
        continue
    else:
        pwGeneOrder.update(list(value.keys()))
pwGeneOrder = list(pwGeneOrder)
pwGeneOrder.sort()

vennGeneOrder = set()
for key, value in annotSamplesDict.items():
    if "only" not in key:
        continue
    else:
        vennGeneOrder.update(list(value.keys()))
vennGeneOrder.update(allVennIDs)
vennGeneOrder = list(vennGeneOrder)
vennGeneOrder.sort()

##### TERM ORDERS

tcGoOrders = []
for sample in samplePrefixes:
    tcGoOrders.append(set())
    tcGoOrders[-1].update(goSamplesDict[sample + "_tc"].keys())
    for suffix in clusterSuffixOrder:
        sampleGroup = sample + suffix
        tcGoOrders[-1].update(goSamplesDict[sampleGroup].keys())
for i in range(len(tcGoOrders)):
    tcGoOrders[i] = list(tcGoOrders[i])
    tcGoOrders[i].sort()

tcKeggOrders = []
for sample in samplePrefixes:
    tcKeggOrders.append(set())
    tcKeggOrders[-1].update(keggSamplesDict[sample + "_tc"].keys())
    for suffix in clusterSuffixOrder:
        sampleGroup = sample + suffix
        tcKeggOrders[-1].update(keggSamplesDict[sampleGroup].keys())
for i in range(len(tcKeggOrders)):
    tcKeggOrders[i] = list(tcKeggOrders[i])
    tcKeggOrders[i].sort()

pwGoOrder = set()
for pwKey in pairwiseSampleOrder:
    pwGoOrder.update(list(goSamplesDict[pwKey].keys()))
pwGoOrder = list(pwGoOrder)
pwGoOrder.sort()

pwKeggOrder = set()
for pwKey in pairwiseSampleOrder:
    pwKeggOrder.update(list(keggSamplesDict[pwKey].keys()))
pwKeggOrder = list(pwKeggOrder)
pwKeggOrder.sort()

vennGoOrder = set()
for vennKey in vennOrder:
    vennGoOrder.update(list(goSamplesDict[vennKey].keys()))
vennGoOrder = list(vennGoOrder)
vennGoOrder.sort()

vennKeggOrder = set()
for vennKey in vennOrder:
    vennKeggOrder.update(list(keggSamplesDict[vennKey].keys()))
vennKeggOrder = list(vennKeggOrder)
vennKeggOrder.sort()

###################
#                 #
# FORMATTING DATA #
#                 #
###################

# Update order values to be more suited for writing to file output
clusterOrderHeader = []
for value in clusterSuffixOrder:
    if "_cluster_" in value:
        newValue = value.split("_")[-1]
    else:
        newValue = value.split("_time_")[1].replace("_direction", "")
    clusterOrderHeader.append(newValue)

pairwiseSampleOrderHeader = []
for value in pairwiseSampleOrder:
    newValue = value.rsplit("_", maxsplit=1)[0]
    pairwiseSampleOrderHeader.append(newValue)


# Generate header lines for all tables
tcOutputAnnotTable = [
    ["gene details", "expression at time"],
    ["gene ID", "predicted gene name", "GO terms", "KEGG terms", "URL", "cluster", "pattern", *tcSampleOrder]
]

pwOutputAnnotTable = [
    ["gene details", "time1", "time2", "time3"],
    ["gene ID", "predicted gene name", "GO terms", "KEGG terms", "URL", *pairwiseSampleOrderHeader]
]

vennOutputAnnotTable = [
    ["gene details", "expression in group", "expression across all"],
    ["gene ID", "predicted gene name", "GO terms", "KEGG terms", "URL", *vennOrder]
]

####

tcOutputTermTable = [
    ["term details", "abundance in cluster", "abundance at time"],
    ["term ID", "term name", "tc", *clusterOrderHeader]
]

pwOutputTermTable = [
    ["term details", "time1", "time2", "time3"],
    ["term ID", "term name", *pairwiseSampleOrderHeader]
]

vennOutputTermTable = [
    ["term details", "abundance in group"],
    ["term ID", "term name", *vennOrder]
]

# Generate rows for DEG tables
## tc tables for each sample
TIMES=range(1, 4)
tcOutputTables = []
for sample in samplePrefixes:
    tcOutputTables.append(tcOutputAnnotTable[:])
    for geneID in clusterIDsDict[sample].keys():
        name, GOs, KEGGs, URL = annotTableDict[geneID]
        pattern, cluster, expressionDict = clusterIDsDict[sample][geneID]
        if cluster == None:
            cluster = "."
        row = [geneID, name, GOs, KEGGs, URL, cluster, pattern]
        for time in TIMES:
            expressionValue = expressionDict[time]
            row.append(str(expressionValue))
        tcOutputTables[-1].append(row)

## pw table
for geneID in pwGeneOrder:
    name, GOs, KEGGs, URL = annotTableDict[geneID]
    row = [geneID, name, GOs, KEGGs, URL]
    for sample in pairwiseSampleOrder:
        try:
            isUpregulated = annotSamplesDict[sample][geneID]
            row.append("up" if isUpregulated is True else "down")
        except:
            row.append(".")
    pwOutputAnnotTable.append(row)

## venn table
for geneID in vennGeneOrder:
    name, GOs, KEGGs, URL = annotTableDict[geneID]
    row = [geneID, name, GOs, KEGGs, URL]
    if geneID in allVennIDs:
        regulationCount = [0, 0] # up, down
        for prefix in samplePrefixes:
            isUpregulated = annotSamplesDict[prefix + "_tc"][geneID]
            if isUpregulated:
                regulationCount[0] += 1
            else:
                regulationCount[1] += 1
        
        row += ["."]*(len(vennOrder)-1)
        row.append(f"{regulationCount[0]}_{regulationCount[1]}")
    else:
        for sample in vennOrder:
            try:
                isUpregulated = annotSamplesDict[sample][geneID]
                row.append("up" if isUpregulated is True else "down")
            except:
                row.append(".")
        row.append(".")
    vennOutputAnnotTable.append(row)

# Generate rows for annotation tables
## tc tables for each sample
tcTermOutputTables = []
for i in range(len(samplePrefixes)):
    sample = samplePrefixes[i]
    tcTermOutputTables.append(tcOutputTermTable[:])
    
    # GO terms
    for term in tcGoOrders[i]:
        row = [term, goTableDict[term]]
        # Check if abundant in overall TC
        try:
            isAbundant = goSamplesDict[sample + "_tc"][term]
            row.append("over" if isAbundant is True else "under")
        except:
            row.append(".")
        # Check if abundant in cluster groups
        for suffix in clusterSuffixOrder:
            sampleGroup = sample + suffix
            try:
                isAbundant = goSamplesDict[sampleGroup][term]
                row.append("over" if isAbundant is True else "under")
            except:
                row.append(".")
        tcTermOutputTables[-1].append(row)
    # KEGG terms
    for term in tcKeggOrders[i]:
        row = [term, keggTableDict[term]]
        # Check if abundant in overall TC
        try:
            isAbundant = keggSamplesDict[sample + "_tc"][term]
            row.append("over" if isAbundant is True else "under")
        except:
            row.append(".")
        # Check if abundant in cluster groups
        for suffix in clusterSuffixOrder:
            sampleGroup = sample + suffix
            try:
                isAbundant = keggSamplesDict[sampleGroup][term]
                row.append("over" if isAbundant is True else "under")
            except:
                row.append(".")
        tcTermOutputTables[-1].append(row)

## pw table
for term in pwGoOrder:
    row = [term, goTableDict[term]]
    for sample in pairwiseSampleOrder:
        try:
            isAbundant = goSamplesDict[sample][term]
            row.append("over" if isAbundant is True else "under")
        except:
            row.append(".")
    pwOutputTermTable.append(row)
for term in pwKeggOrder:
    row = [term, keggTableDict[term]]
    for sample in pairwiseSampleOrder:
        
        try:
            isAbundant = keggSamplesDict[sample][term]
            row.append("over" if isAbundant is True else "under")
        except:
            row.append(".")
        
    pwOutputTermTable.append(row)

## venn table
for term in vennGoOrder:
    row = [term, goTableDict[term]]
    for sample in vennOrder:
        try:
            isAbundant = goSamplesDict[sample][term]
            row.append("over" if isAbundant is True else "under")
        except:
            row.append(".")
    vennOutputTermTable.append(row)
for term in vennKeggOrder:
    row = [term, keggTableDict[term]]
    for sample in vennOrder:
        try:
            isAbundant = keggSamplesDict[sample][term]
            row.append("over" if isAbundant is True else "under")
        except:
            row.append(".")
    vennOutputTermTable.append(row)

################
#              #
# WRITING DATA #
#              #
################

################## DEGS Files

## tc file
for i in range(len(samplePrefixes)):
    sample = samplePrefixes[i]
    outputAnnotFileName = f"{sample}_dge_results_table.tsv"
    with open(outputAnnotFileName, "w") as fileOut:
        fileOut.write("\n".join(
            ["\t".join(row) for row in tcOutputTables[i]]
        ))

## pw file
pwOutputAnnotFileName = "pw_dge_results_table.tsv"
with open(pwOutputAnnotFileName, "w") as fileOut:
    fileOut.write("\n".join(
        ["\t".join(row) for row in pwOutputAnnotTable]
    ))

## venn file
vennOutputFileName = "venn_dge_results_table.tsv"
with open(vennOutputFileName, "w") as fileOut:
    fileOut.write("\n".join(
        ["\t".join(row) for row in vennOutputAnnotTable]
    ))

################## Terms Files

## tc file
for i in range(len(samplePrefixes)):
    sample = samplePrefixes[i]
    outputAnnotFileName = f"{sample}_goseq_results_table.tsv"
    with open(outputAnnotFileName, "w") as fileOut:
        fileOut.write("\n".join(
            ["\t".join(row) for row in tcTermOutputTables[i]]
        ))

## pw file
pwOutputTermFileName = "pw_goseq_results_table.tsv"
with open(pwOutputTermFileName, "w") as fileOut:
    fileOut.write("\n".join(
        ["\t".join(row) for row in pwOutputTermTable]
    ))

## venn file
vennOutputTermFileName = "venn_goseq_results_table.tsv"
with open(vennOutputTermFileName, "w") as fileOut:
    fileOut.write("\n".join(
        ["\t".join(row) for row in vennOutputTermTable]
    ))

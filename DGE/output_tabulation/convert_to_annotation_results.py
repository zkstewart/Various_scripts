#! python3

import shutil, os

annotationTable = r"F:\plant_group\mango_cultivars\annotations\mango_GOextended_table.tsv"

keggMap = r"F:\plant_group\mango_cultivars\annotations\manindi_flc.aa_keggmap.tsv"
keggNames = r"F:\plant_group\mango_cultivars\annotations\manindi_flc.aa_keggnames.tsv"

directories = ["pairwise_DE", "timecourse_DE", "venn_timecourse_DE"]

tcFileSuffix = "tc_results.tsv"

# Parse KEGG information
keggNamesDict = {}
with open(keggNames, "r") as fileIn:
    for line in fileIn:
        if line == "": continue
        
        sl = id, name = line.rstrip("\r\n ").split("\t")
        keggNamesDict[id] = name

keggMapDict = {}
with open(keggMap, "r") as fileIn:
    for line in fileIn:
        if line == "": continue
        
        sl = geneID, keggAnnot = line.rstrip("\r\n ").split("\t")
        keggMapDict[geneID] = keggAnnot

# Parse annotation information
class Annotation:
    def __init__(self, sl):
        # Extract best value from all sl values
        for i in range(len(sl)):
            sl[i] = sl[i].split(" [")[0]
        
        # Store values
        self.query, self.source, self.target, self.name, self.taxon, \
        self.length, self.percent, self.alignLength, self.mismatch, \
        self.gaps, self.qstart, self.qend, self.tstart, self.tend, self.evalue, \
        self.bitscore, self.bestAccess, self.bestGO, self.bestGOParents = sl
        
annotTableDict = {}
with open(annotationTable, "r") as fileIn:
    for line in fileIn:
        if line == "" or line.startswith("#"): continue
        
        sl = line.rstrip("\r\n ").split("\t")
        annotObj = Annotation(sl)
        annotTableDict[annotObj.query] = annotObj

# Get full path to directories
for i in range(len(directories)):
    directories[i] = os.path.abspath(directories[i])

# Go through directories and change files
def handle_dir_contents(dir):
    currentDir = dir
    os.chdir(dir)
    
    # Loop through all files
    for file in os.listdir():
        # Handle GOseq files
        if file.endswith("_GOseq.tsv"):
            "No changes needed"
            pass
        
        # Handle KEGGseg files
        elif file.endswith("_KEGGseq.tsv"):
            # Create new file
            with open(file, "r") as fileIn, open(file + ".tmp", "w") as fileOut:
                for line in fileIn:
                    sl = line.rstrip("\r\n ").replace('"', '').split("\t") # fix any previous changes
                    sl = sl[0:5] # drop any concatenated columns during development
                    if sl[0].startswith("map"):
                        sl[0] = sl[0][3:] # drop the extra bit added to the start
                    
                    # Handle header line
                    if line.startswith("category\t"):
                        sl.append("term")
                    
                    # Handle body lines
                    else:
                        term = keggNamesDict[sl[0]]
                        sl[0] = 'map{0}'.format(sl[0]) # make Excel show things better
                        sl.append(term)
                    
                    # Write output line
                    fileOut.write("{0}\n".format("\t".join(sl)))
            
            # Replace original with new file
            shutil.move(file + ".tmp", file)
        
        # Handle DESeq2 files
        elif file.endswith("_results.tsv") and not "_cluster_" in file:
            # Create new file
            with open(file, "r") as fileIn, open(file + ".tmp", "w") as fileOut:
                for line in fileIn:
                    sl = line.rstrip("\r\n ").replace('"', '').split("\t") # fix any previous changes
                    
                    # Drop any columns we don't care about here / were concatenated during development
                    if sl[-1] == "URL" or "http" in sl[-1] or sl[-1] == ".":
                        sl = [sl[0], sl[1], sl[2]]
                    else:
                        sl = [sl[0], sl[2], sl[6]]
                    
                    # Handle header line
                    if sl[1] == "log2FoldChange":
                        sl = [
                            "gene_ID", "log2FoldChange", "pvalue_adj", "target_accession",
                            "gene_name", "percentage_identity", "evalue", "GOs", "KEGGs",
                            "URL"
                        ]
                    
                    # Handle body lines
                    else:
                        # Retrieve relevant details
                        keggTerm = keggMapDict[sl[0]]
                        if keggTerm != "0":
                            splitKeggTerm = keggTerm.split("; ")
                            for i in range(len(splitKeggTerm)):
                                splitKeggTerm[i] = "map" + splitKeggTerm[i]
                            keggTerm = "; ".join(splitKeggTerm)
                        else:
                            keggTerm = "."
                        annotObj = annotTableDict[sl[0]]
                        
                        # Format URL for redirection to database
                        if annotObj.target.startswith("UPI"): # handle UniParc sequences
                            url = "https://www.uniprot.org/uniparc/{0}".format(annotObj.target)
                        elif annotObj.target == ".": # handle blanks?
                            url = "."
                        else: # handle everything else?
                            url = "https://www.uniprot.org/uniprot/{0}".format(annotObj.target)
                        
                        # Update body line
                        sl.append(annotObj.target)
                        sl.append(annotObj.name)
                        sl.append(annotObj.percent)
                        sl.append(annotObj.evalue)
                        sl.append(annotObj.bestGOParents)
                        sl.append(keggTerm)
                        sl.append(url)
                    # Write output line
                    fileOut.write("{0}\n".format("\t".join(sl)))
            
            # Replace original with new file
            shutil.move(file + ".tmp", file)
            
        # Recursively handle subdirectories
        elif os.path.isdir(file):
            handle_dir_contents(file)
            os.chdir(currentDir) # come back here after our excursion
        
        # Note any other files
        else:
            print("{0} file was skipped".format(file))

for dir in directories:
    # Loop through all files recursively and process
    handle_dir_contents(dir)

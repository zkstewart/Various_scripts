#! python 3
# compare_tree_against_taxa_labels.py
# Script to assess the distance between a provided phylogeny
# and the species it makes up to see how much agreement there
# is between the tree and the known taxonomic groupings.

import os, argparse, statistics
from ete3 import Tree
from itertools import combinations

def validate_args(args):
    # Validate input data locations
    if not os.path.isfile(args.tsvFile):
        print('I am unable to locate the input TSV file (' + args.tsvFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.nwkFile):
        print('I am unable to locate the input NWK file (' + args.nwkFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Handle file output
    if os.path.isfile(args.outputFileName):
        print('The specified output file already exists. This program will not allow overwriting.')
        print('Specify a different output file name or move/rename the existing file and try again.')

def is_a_tsv(fileName):
    numRows = None
    with open(fileName, "r") as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n").split("\t")
            if sl == []:
                continue
            
            if numRows != None and len(sl) != numRows:
                return False
            
            numRows = len(sl)
    return True

def get_col_indices(tsvFileName, idColumn, labels):
    with open(tsvFileName, "r") as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n").split("\t")
            assert idColumn in sl, \
                f"{idColumn} not found within TSV header, can't proceed"
            assert all([l in sl for l in labels]), \
                "Not all labels found within TSV header, can't proceed"
            
            idColIndex = sl.index(idColumn)
            labColIndices = [sl.index(l) for l in labels]
            
            return idColIndex, labColIndices

class TaxonomyTSV:
    def __init__(self, tsvFileName, idColIndex, labelColIndices):
        self.parse_tsv(tsvFileName, idColIndex, labelColIndices)
    
    def parse_tsv(self, tsvFileName, idColIndex, labelColIndices):
        # Set values to check what params define this object
        self.fileName = tsvFileName
        self.idColIndex = idColIndex
        self.labelColIndices = labelColIndices
        
        # Parse TSV into dictionary groupings
        self.taxaDict = {}
        with open(tsvFileName, "r") as fileIn:
            for line in fileIn:
                if line.startswith("#"):
                    continue
                
                sl = line.rstrip("\r\n").split("\t")
                sampleID = sl[idColIndex]
                labels = [sl[i] for i in labelColIndices]
                for label in labels:
                    if label == ".":
                        continue
                    self.taxaDict.setdefault(label, [])
                    self.taxaDict[label].append(sampleID)
        
        # Retrieve non-redundant pairwise comparisons from groupings (where possible)
        self.taxaPairs = {}
        for taxon, sampleIDs in self.taxaDict.items():
            if len(sampleIDs) < 2:
                continue
            
            self.taxaPairs[taxon] = list(combinations(sampleIDs, 2))

def average_tree_distance_from_pairwise_comparisons(ete3Tree_obj, pairwiseDict):
    '''
    Parameters:
        ete3Tree_obj -- an ete3.Tree() object.
        pairwiseDict -- a dict with structure like:
                        {
                            taxonRank: [
                                (sampleID_1, sampleID_2),
                                (sampleID_1, sampleID_3),
                                ...
                            ],
                            ...
                        }
    '''
    # Get the distance for ALL pairwise comparisons
    "Yeah, this'll be slow"
    distsDict = {}
    for taxonRank, pairedSampleIDs in pairwiseDict.items():
        distsDict[taxonRank] = []
        for sampleID_A, sampleID_B in pairedSampleIDs:
            node_A = ete3Tree_obj&sampleID_A
            dist = node_A.get_distance(sampleID_B)
            distsDict[taxonRank].append(dist)
    
    # Get the average for each taxon rank
    averageDistsDict = {}
    for taxonRank, dists in distsDict.items():
        averageDistsDict[taxonRank] = statistics.mean(dists)
    
    return averageDistsDict

def main():
    usage = """%(prog)s receives a .tsv formatted file and a tree in .nwk format.
    Using the taxonomic labels in the .tsv file alongside the species IDs that can
    be found in the .nwk file, it will derive the amount of agreement between the tree
    and known taxonomic groupings based on labels provided in the .tsv.
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-t", dest="tsvFile", required=True,
                help="Specify the location of the input TSV file")
    p.add_argument("-n", dest="nwkFile", required=True,
                help="Specify the location of the input NWK file")
    p.add_argument("-o", dest="outputFileName", required=True,
                help="Output file name to be created")
    p.add_argument("--idColumn", dest="idColumn", required=True,
                help="Specify the header value for the column containing IDs as seen in the .nwk file")
    p.add_argument("--labels", dest="labels", required=True, nargs="+",
                help="Specify one or more of the header values for columns containing relevant labels")
    # Opts
    p.add_argument("--root", dest="root", required=False,
                help="Optionally, specify which taxon ID should be the root of the tree (if unrooted currently)")
    
    args = p.parse_args()
    validate_args(args)
    
    # Validate that the input file is a TSV
    isTsv = is_a_tsv(args.tsvFile)
    assert isTsv, \
        "TSV file format doesn't check out; some rows have fewer columns than others"
    
    # Validate and retrieve indices for ID and label columns
    idColIndex, labelColIndices = get_col_indices(args.tsvFile, args.idColumn, args.labels)
    
    # Parse nwk as an ETE3 Tree
    tree = Tree(args.nwkFile)
    
    # If root was specified, validate that the taxon exists and root the tree
    try:
        rootNode = tree&args.root
        tree.set_outgroup(rootNode)
    except:
        print(f"{args.root} not found within the tree file; program closing now")
        quit()
    
    # Normalise branch lengths 
    tree.convert_to_ultrametric()
    for node in tree.iter_descendants("postorder"):
        node.dist = 1.0
    
    # Parse input TSV into a helpful structure for getting pairwise comparisons
    taxaTSV = TaxonomyTSV(args.tsvFile, idColIndex, labelColIndices)
    
    # Solve distances
    averageDistsDict = average_tree_distance_from_pairwise_comparisons(tree, taxaTSV.taxaPairs)
    
    # Write sorted output
    taxonOrder = list(averageDistsDict.keys())
    taxonOrder.sort()
    with open(args.outputFileName, "w") as fileOut:
        fileOut.write("#rank\tdist\n")
        for taxonID in taxonOrder:
            fileOut.write("{0}\t{1}\n".format(taxonID, averageDistsDict[taxonID]))
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

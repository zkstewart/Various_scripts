#! python3
# parse_iqtree.py
# This script takes an IQTREE file and obtains a Newick file out of it

import os, argparse
from Bio import Phylo
from io import StringIO

# Define functions
def validate_args(args):
    # Validate input file
    args.iqtreeFile = os.path.abspath(args.iqtreeFile)
    if not os.path.isfile(args.iqtreeFile):
        raise FileNotFoundError(f"-i '{args.iqtreeFile}' does not exist or is not a file!")
    
    if args.renameFile != None:
        args.renameFile = os.path.abspath(args.renameFile)
        if not os.path.isfile(args.renameFile):
            raise FileNotFoundError(f"--rename '{args.renameFile}' does not exist or is not a file!")
    
    # Validate output file
    args.outputFileName = os.path.abspath(args.outputFileName)
    if os.path.exists(args.outputFileName):
        raise FileExistsError(f"-o '{args.outputFileName}' already exists and will not be rewritten.")
    parentDir = os.path.dirname(args.outputFileName)
    if not os.path.isdir(parentDir):
        raise ValueError(f"Output cannot be written to '{args.outputFileName}' as its parent directory ({parentDir}) does not exist!")

def parse_iqtree_nwk(fileName):
    with open(fileName, "r") as fileIn:
        for line in fileIn:
            if line.startswith("("):
                return line.rstrip()

def parse_fam_file(famFileName):
    famDict = {}
    numInfam = {}
    with open(famFileName, "r") as fileIn:
        for line in fileIn:
            fid, iid, p1, p2, sex, pheno = line.rstrip().split("\t")
            
            numInfam.setdefault(fid, 0)
            numInfam[fid] += 1
            famDict[iid] = f"{fid}.{numInfam[fid]}"
    return famDict

def write_tree_to_file(outputFileName, tree):
    Phylo.write(tree, outputFileName, "newick")

## Main
def main():
    # User input
    usage = """%(prog)s will receive a .iqtree file and extract a .nwk file with
    some limited options for reformatting the file.
    """
    p = argparse.ArgumentParser(description=usage)
    # Reqs
    p.add_argument("-i", dest="iqtreeFile",
                   required=True,
                   help="Input .iqtree file")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Specify location to write rerooted tree")
    # Opts
    p.add_argument("--rename", dest="renameFile",
                   required=False,
                   help="""Optionally provide a .fam PLINK2 metadata file
                   to automatically rename samples""",
                   default=None)
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse the tree out from .iqtree file
    nwkString = parse_iqtree_nwk(args.iqtreeFile)
    tree = Phylo.read(StringIO(nwkString), "newick")
    
    # Rename samples as needed
    if args.renameFile != None:
        renameDict = parse_fam_file(args.renameFile)
        
        for terminal in tree.get_terminals():
            if terminal.name in renameDict:
                newName = renameDict[terminal.name]
                print(f"{terminal.name} -> {newName}")
                terminal.name = newName
    
    # Write output
    write_tree_to_file(args.outputFileName, tree)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

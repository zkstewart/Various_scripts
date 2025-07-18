#! python3
# reroot_tree.py
# This script reroots a phylogenetic tree with a specified outgroup.

import os, argparse
from Bio import Phylo

# Define functions
def validate_args(args):
    # Validate input file
    if not os.path.isfile(args.inputTree):
        raise FileNotFoundError(f"Input tree file '{args.inputTree}' does not exist.")
    # Validate output file
    if os.path.exists(args.outputTree):
        raise FileExistsError(f"Output tree file '{args.outputTree}' already exists.")

## Main
def main():
    # User input
    usage = """%(prog)s will receive a phylogenetic tree in Newick format and reroot
    it with a specified outgroup species.
    """
    p = argparse.ArgumentParser(description=usage)
    # Reqs
    p.add_argument("-i", dest="inputTree",
                   required=True,
                   help="Input newick tree file to reroot")
    p.add_argument("-r", dest="rootSpecies",
                   required=True,
                   help="Species to reroot the tree with")
    p.add_argument("-o", dest="outputTree",
                   required=True,
                   help="Specify location to write rerooted tree")
    args = p.parse_args()
    validate_args(args)
    
    # Load the tree from the Newick file
    tree = Phylo.read(args.inputTree, "newick")
    
    # Reroot the tree at the new root node
    tree.root_with_outgroup(args.rootSpecies)
    
    # Save the rerooted tree to a new Newick file
    Phylo.write(tree, args.outputTree, "newick")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

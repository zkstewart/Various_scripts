#! python3
# prune_tree.py
# This script will prune a phylogenetic tree by removing specified samples.

import os, argparse, shutil, subprocess, platform
from Bio import Phylo

# Define functions
def validate_args(args):
    def _not_found_error(program):
        raise FileNotFoundError(f"{program} not discoverable in your system PATH and was not specified as an argument.")
    def _specified_wrong_error(program, path):
        raise FileNotFoundError(f"{program} was not found at the indicated location '{path}'")
    
    # Validate program discoverability
    if args.nwkit == None:
        args.nwkit = shutil.which("nwkit")
        if args.nwkit == None:
            _not_found_error("nwkit")
    else:
        if not os.path.isfile(args.nwkit):
            _specified_wrong_error("nwkit", args.nwkit)
    
    # Validate input file
    if not os.path.isfile(args.inputTree):
        raise FileNotFoundError(f"Input tree file '{args.inputTree}' does not exist.")
    # Validate output file
    if os.path.exists(args.outputTree):
        raise FileExistsError(f"Output tree file '{args.outputTree}' already exists.")

def parse_leaf_names(tree):
    '''
    Parses the leaf names from a Phylo tree object.
    
    Parameters:
        tree -- a Phylo tree object
    
    Returns:
        A list of leaf names in the tree.
    '''
    return [leaf.name for leaf in tree.get_terminals()]

def run_nwkit(treeFile, samples, outputFile, nwkitPath="nwkit"):
    '''
    Runs nwkit to prune the samples from the input tree file.
    
    Parameters:
        treeFile -- a string indicating the location of the tree file to be pruned
        samples -- a list of strings indicating the samples to prune from the tree
        outputFile -- a string indicating the location to write the pruned tree
        nwkitPath -- (OPTIONAL) a string indicating the location of the nwkit executable
    '''
    # Format the samples into a nwkit pattern
    pattern = "|".join(samples)
    
    # Format nwkit command
    cmd = [
        nwkitPath, "prune", "--infile", treeFile, "--outfile", outputFile,
        "--pattern", f'"{pattern}"'
    ]
    
    # Run nwkit
    if platform.system() != "Windows":
        run_nwkit = subprocess.Popen(" ".join(cmd), shell = True,
                                      stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
    else:
        run_nwkit = subprocess.Popen(cmd, shell = True,
                                      stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
    nwkitout, nwkiterr = run_nwkit.communicate()
    if not "Number of leaves in output tree" in nwkiterr.decode("utf-8"):
        raise Exception('nwkit error text below\n' + nwkiterr.decode("utf-8"))
    
    # Print information about the pruning
    print(nwkiterr.decode("utf-8").strip())

## Main
def main():
    # User input
    usage = """%(prog)s will receive a phylogenetic tree in Newick format and remove
    specified samples from the tree.
    """
    p = argparse.ArgumentParser(description=usage)
    # Reqs
    p.add_argument("-i", dest="inputTree",
                   required=True,
                   help="Input newick tree file to reroot")
    p.add_argument("-p", dest="prunedSamples",
                   required=True,
                   nargs="+",
                   help="Specify one or more samples to prune from the tree")
    p.add_argument("-o", dest="outputTree",
                   required=True,
                   help="Specify location to write rerooted tree")
    # Optional
    p.add_argument("--nwkit", dest="nwkit",
                   required=False,
                   help="""Optionally, specify the nwkit executable file location
                   if it is not discoverable in the path""",
                   default=None)
    p.add_argument("--keep", dest="keep",
                   required=False,
                   action="store_true",
                   help="""Optionally, provide this flag to indicate that the samples
                   provided through -p should be kept in the tree, rather than pruned.""",
                   default=False)
    args = p.parse_args()
    validate_args(args)
    
    # If --keep, obtain the samples to prune
    if args.keep:
        # Parse leaf names from the tree
        tree = Phylo.read(args.inputTree, "newick")
        allSamples = parse_leaf_names(tree)
        
        # Determine which samples to prune
        toPrune = [ s for s in allSamples if s not in args.prunedSamples ]
        args.prunedSamples = toPrune
    
    run_nwkit(args.inputTree, args.prunedSamples, args.outputTree, nwkitPath=args.nwkit)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

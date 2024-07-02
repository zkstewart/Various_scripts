#! python3
# changed_families_at_nodes.py
# Script to identify gene families that have expanded or
# contracted at specific nodes.

import os, argparse
from CAFE5 import parse_change_tab, parse_probabilities_tab

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.changeFile):
        print('I am unable to locate the *_change.tab file (' + args.changeFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.probabilitiesFile):
        print('I am unable to locate the *_branch_probabilities.tab file (' + args.probabilitiesFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate numeric arguments
    if not all([x > 0 for x in args.nodes]):
        print("All node numbers must be greater than 0.")
        quit()
    # Handle file overwrites
    if os.path.exists(args.outputFileName):
        print(f"{args.outputFileName} already exists, and I will not overwrite it.")
        print("Delete/move whatever exists here, or specify a different output name when you try again.")
        quit()

def main():
    # User input
    usage = """%(prog)s reads in the *(Base or Gamma)_change.tab file and the *_family_results.txt file
    from CAFE in order to identify gene families that have expanded or contracted at nodes specified
    by you. The output will be a TSV listing gene families as well as whether they, cumulatively,
    expanded or contracted. For -n you should specify node numbers e.g., if you have a column in the
    change file labelled "species<7>" you should specify 7 as the node number, and for internal nodes
    labelled as "<13>" you should specify 13 as the node number. NOTE: for a family to be
    assessed at all, it must have changed significantly at least once in the given nodes.
    """
    p = argparse.ArgumentParser(description=usage)
    # Reqs
    p.add_argument("-c", dest="changeFile",
                   required=True,
                   help="Input *_change.tab file name")
    p.add_argument("-p", dest="probabilitiesFile",
                   required=True,
                   help="Input *_branch_probabilities.tab file name")
    p.add_argument("-n", dest="nodes",
                   required=True,
                   type=int,
                   nargs="+",
                   help="Specify one or more node numbers to assess")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file prefix for the chi-square results")
    # Opts
    p.add_argument("--pvalue", dest="pvalue",
                   required=False,
                   type=float,
                   help="""Optionally specify the p-value threshold for significance;
                   default == 0.05""",
                   default=0.05)
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse *_change.tab file
    changeDict = parse_change_tab(args.changeFile)
    
    # Parse *_branch_probabilities.txt file
    significantDict = parse_probabilities_tab(args.probabilitiesFile,
                                              args.pvalue)
    
    # Find families that have a significant change at the specified nodes
    sigFamilies = set()
    for nodeNum in args.nodes:
        sigFamilies.update(significantDict[str(nodeNum)])
    
    # Sum the changes of significantly changing families at the specified nodes
    changeSum = {}
    for family in sigFamilies:
        changeSum[family] = 0
        for nodeNum in args.nodes:
            changeSum[family] += changeDict[str(nodeNum)][family]
    
    # Write output file
    with open(args.outputFileName, "w") as outFile:
        outFile.write("FamilyID\tChange\n")
        for family, change in changeSum.items():
            if change > 0:
                change = "expansion"
            elif change < 0:
                change = "contraction"
            else:
                continue
            outFile.write(f"{family}\t{change}\n")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

#! python3
# orthofinder_group_containing_sequences
# Script to produce an abbreviated Orthogroups.csv file which contains only
# orthogroups which contain any of a list of sequence IDs

import os, argparse

def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.orthogroupsFileName):
        raise FileNotFoundError(f"-i file '{args.orthogroupsFileName}' could not be found or is not a file!")
    if not os.path.isfile(args.textFile):
        raise FileNotFoundError(f"-t file '{args.textFile}' could not be found or is not a file!")
    # Handle file overwrites
    if os.path.exists(args.outputFileName):
        raise FileExistsError(f"-o file '{args.outputFileName}' already exists!")

def parse_list_file(fileName):
    seqIDs = []
    with open(fileName, "r") as fileIn:
        for line in fileIn:
            l = line.strip()
            if l != "":
                seqIDs.append(l)
    return set(seqIDs)

def main():
    ### USER INPUT
    usage = """%(prog)s will parse the Orthogroups.tsv file output by Orthofinder and 
    produce a modified output file with only orthogroups which contain sequences from
    the specified species of interest (SOI). Specification of SOI is in the form of
    providing the name of the column(s) you want; if specifying more than one column,
    separate each value with a space e.g., -s species1 species2.
    """

    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="orthogroupsFileName",
                required=True,
                help="Specify the orthogroup file")
    p.add_argument("-t", dest="textFile",
                required=True,
                help="Text file listing sequence IDs to look for")
    p.add_argument("-o", dest="outputFileName",
                required=True,
                help="Output file name")

    args = p.parse_args()
    validate_args(args)
    
    # Parse text file
    seqIDs = parse_list_file(args.textFile)
    
    # Parse Orthofinder file
    isHeader = True
    with open(args.orthogroupsFileName, "r") as fileIn, open(args.outputFileName, "w") as fileOut:
        for line in fileIn:
            if isHeader:
                fileOut.write(line)
                isHeader = False
            else:
                sl = line.strip().split("\t")
                
                wroteLine = False
                for memberIDs in sl[1:]:
                    if wroteLine:
                        break
                    for member in memberIDs.split(", "):
                        if member in seqIDs:
                            fileOut.write(line)
                            wroteLine = True
                            break
    
    print('Program completed successfully!')

if __name__ == "__main__":
    main()

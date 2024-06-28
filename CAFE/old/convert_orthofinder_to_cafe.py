#! python3
# convert_orthofinder_to_cafe.py
# Script to take the OrthoFinder output file 
# (Orthogroups.GeneCount.tsv) and make it amenable
# to CAFE (or other gene family history reconstruction) operations.

import os, argparse

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.tsvFileName):
        print('I am unable to locate the Orthogroups.GeneCount.tsv file (' + args.tsvFileName + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate output file location
    if os.path.isfile(args.outputFileName):
        print('File already exists at output location (' + args.outputFileName + ')')
        print('Make sure you specify a unique file name and try again.')
        quit()

def main():
    # User input
    usage = """%(prog)s reads in the Orthogroups.GeneCount.tsv file producted by OrthoFinder
    and reformats it for use by CAFE5 or DupliPHY-ML.
    """
    ## Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-t", dest="tsvFileName",
                   required=True,
                   help="Input Orthogroups.GeneCount.tsv file location")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file name for the modified tsv file")
    ## Opts
    p.add_argument("--dupliphy", dest="dupliphyFormat",
                   required=False,
                   action="store_true",
                   help="Optionally specify if you'd like the output to be in DupliPHY-ML format",
                   default=False)
    args = p.parse_args()
    validate_args(args)
    
    # Produce output file
    with open(args.tsvFileName, "r") as fileIn, open(args.outputFileName, "w") as fileOut:
        firstLine = True
        for line in fileIn:
            sl = line.rstrip("\r\n").split("\t")
            
            # Handle header line
            if firstLine == True:
                if args.dupliphyFormat:
                    fileOut.write("FAMILY\t{0}".format("\t".join(sl[1:-1]))) # drop the Total column at end
                else:
                    fileOut.write("Desc\tFamily ID\t{0}".format("\t".join(sl[1:-1])))
                firstLine = False
            # Handle content lines
            else:
                familyID, counts = sl[0], sl[1:-1] # drop the Total column at end
                if args.dupliphyFormat:
                    fileOut.write("{0}\t{1}".format(familyID, "\t".join(counts)))
                else:
                    fileOut.write("(null)\t{0}\t{1}".format(familyID, "\t".join(counts)))
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

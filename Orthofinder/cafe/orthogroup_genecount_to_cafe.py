#! python3
# orthogroup_genecount_to_cafe.py
# Script intended to take the Orthogroups.GeneCount.tsv file
# from OrthoFinder and modify its format to be suitable for
# input to CAFE5 / DupliPHY-ML.

import argparse, os

def validate_args(args):
    # Validate input data location
    if not os.path.isfile(args.geneCountFile):
        print('I am unable to locate the directory where the Orthogroup.GeneCount file is (' + args.geneCountFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Handle file output
    if os.path.isfile(args.outputFileName):
        print(args.outputFileName + ' already exists. Specify a new file name or move/rename the existing file.')
        quit()

if __name__ == "__main__":
    #### USER INPUT SECTION
    usage = """%(prog)s will take the Orthogroups.GeneCount.tsv file produced
    by OrthoFinder and modify it for use by CAFE5 or DupliPHY. This involves light reformatting
    of the original file including appending a new 'Desc' column which will be
    populated with '(null)' values.
    """
    p = argparse.ArgumentParser(description=usage)
    # Required
    p.add_argument("-i", dest="geneCountFile",
                   required=True,
                   help="""Specify the location of the directory containing
                   single copy gene FASTA files""")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="""Specify the location to write output files to; this
                   location will be populated with orthogroup files""")
    # Optional
    p.add_argument("--new_ids", dest="newIDsList",
                   required=False,
                   nargs="+",
                   help="Optionally, specify multiple space-separated IDs to rename headers to",
                   default=None)
    p.add_argument("--dupliphy", dest="dupliphyFormat",
                   required=False,
                   action="store_true",
                   help="Optionally specify if you'd like the output to be in DupliPHY-ML format",
                   default=False)
    args = p.parse_args()
    validate_args(args)
    
    # Parse through input file, generating output
    with open(args.geneCountFile, "r") as fileIn, open(args.outputFileName, "w") as fileOut:
        for line in fileIn:
            sl = line.rstrip("\r\n ").split("\t")
            # Handle header line
            if sl[0] == "Orthogroup":
                sampleIDs = sl[1:-1] # omit first ('Orthogroup') and last ('Total')
                
                # Update sample IDs if relevant
                if args.newIDsList != None:
                    assert len(sampleIDs) == len(args.newIDsList), \
                        f"--new_ids number of values is not equal to the number of values in {sampleIDs}"
                    sampleIDs = args.newIDsList
                
                # Write header line to file
                if args.dupliphyFormat:
                    fileOut.write("FAMILY\t{0}\n".format(
                        "\t".join(sampleIDs)
                    ))
                else:
                    fileOut.write("Desc\tFamily ID\t{0}\n".format(
                        "\t".join(sampleIDs)
                    ))
            
            # Handle body lines
            else:
                if args.dupliphyFormat:
                    fileOut.write("{0}\n".format(
                        "\t".join(sl[:-1]) # omit last ('Total') value
                    ))
                else:
                    fileOut.write("(null)\t{0}\n".format(
                        "\t".join(sl[:-1])
                    ))
    
    # Done!
    print('Program completed successfully!')

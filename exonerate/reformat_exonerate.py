#! python3
# reformat_exonerate.py
# Reformats exonerate's output into a standard GFF3 format.

import os, argparse

def validate_args(args):
    # Validate input data location
    if not os.path.isfile(inputFile):
        raise FileNotFoundError((f"I am unable to locate the file at '{args.inputFile}'. " + 
                                "Make sure you've typed the name or location correctly and try again."))
    # Handle file output
    if os.path.exists(args.outputLocation):
        raise FileExistsError((f"File already exists at output location ({args.outputLocation}). " + 
                                "Make sure you specify a unique file name and try again."))
    elif not os.path.isdir(os.path.dirname(os.path.abspath(args.outputLocation))):
        FileNotFoundError((f"Output file '{args.outputLocation}' would be written to a non-existent directory" + 
                            "If you provide a full path, make sure its parent directories exist; " +
                            "otherwise, provide a file name only."))

def convert_exonerate_gff(exonerateFile, outputFileName):
    '''
    This function will take an exonerate output file and reformat
    it into a GFF3-style format.
    
    Function copy pasted from old code, so it's written poorly.
    
    Parameters:
        exonerateFile -- a string pointing to the location of the exonerate
                         output file.
        outputFileName -- a string pointing to the location where the output
                          file should be written.
    '''
    # Setup
    geneIDDict = {}
    # Main function
    with open(exonerateFile, "r") as fileIn, open(outputFileName, "w") as fileOut:
        for line in fileIn:
            # Skip non-GFF2 lines
            sl = line.split('\t')
            if len(sl) != 9: # GFF lines have 9 entries
                continue
            try:
                int(sl[3]) # This corresponds to start coordinate, and should be numeric
                int(sl[4]) # This corresponds to end coordinate, and should be numeric
                assert sl[6] in ['+', '-'] # This corresponds to orientation, and should be specified as + or -
            except:
                continue  # If any of the above tests fail we know that this isn't a GFF2 line
            
            # Skip similarity line since it's irrelevant
            if sl[2] == 'similarity':
                continue
            
            # Parse details
            details = sl[8].rstrip('\n').split(' ; ')
            detailDict = {}
            for i in range(len(details)):
                if details[i] == '':
                    continue
                splitDetail = details[i].split(' ', maxsplit=1)
                detailDict[splitDetail[0]] = splitDetail[1]
            
            # Reformat detail lines to be GFF3-style
            if sl[2] == 'gene':
                geneID, geneIDDict = exonerate_geneid_produce(sl[0], detailDict['sequence'], geneIDDict)    # This will carry over into CDS/exon lines below this
                name = 'exonerate_' + geneID
                sl[8] = 'ID=' + geneID + ';Name=' + name + ';Sequence=' + detailDict['sequence'] + ';identity=' + detailDict['identity'] + ';similarity=' + detailDict['similarity']
                exonCount = 1
                cdsCount = 1
            elif sl[2] == 'cds':
                sl[2] = 'CDS'
                sl[8] = 'ID=' + geneID + '.mrna1.CDS' + str(cdsCount) + ';Parent=' + geneID + '.mrna1'
                cdsCount += 1
            elif sl[2] == 'exon':
                sl[8] = 'ID=' + geneID + '.mrna1.exon' + str(exonCount) + ';Parent=' + geneID + '.mrna1'
                exonCount += 1
            else:
                continue    # Skip any other lines; these include 'intron' and 'splice5'/'splice3' which are implied based on other GFF details
            
            # Write line to file
            fileOut.write('\t'.join(sl) + '\n')
            
            # Format an additional mRNA line under gene lines
            if sl[2] == 'gene':
                sl[2] = 'mRNA'
                sl[8] = 'ID=' + geneID + '.mrna1;Name=' + name + ';Sequence=' + detailDict['sequence'] + ';Parent=' + geneID + ';identity=' + detailDict['identity'] + ';similarity=' + detailDict['similarity']
                fileOut.write('\t'.join(sl) + '\n')

def main():
    #### USER INPUT SECTION
    usage = """%(prog)s will take in a raw exonerate output and reformat it to
    a GFF3-style format more suitable for modern analyses.
    """
    p = argparse.ArgumentParser(description=usage)
    # Required
    p.add_argument("-i", dest="inputFile",
                   required=True,
                   help="Specify locations of exonerate output file")
    p.add_argument("-o", dest="outputLocation",
                   required=True,
                   help="Output file name for reformatted exonerate output")
    
    args = p.parse_args()
    validate_args(args)
    
    # Run the main function
    convert_exonerate_gff(args.inputFile, args.outputLocation)
    
    # Done!
    print('Program completed successfully!')

if __name__ == "__main__":
    main()

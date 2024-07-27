#! python3
# geneious_to_fasta.py
# Simple script to open a .geneious file and extract
# sequences. Geneious files are actually just zipped
# XMLs, so we can open them programmatically. This
# script assumes a format that might not be universally
# true though.

import os, argparse, xmltodict
from zipfile import ZipFile

def validate_args(args):
    # Validate input location
    if not os.path.isfile(args.geneiousFile):
        raise FileNotFoundError((f"I am unable to locate the input file ({args.geneiousFile}). " + 
                                "Make sure you've typed the name or location correctly and try again."))
    # Handle file output
    if os.path.exists(args.outputFileName):
        raise FileExistsError((f"File already exists at output location ({args.outputFileName}). " + 
                                "Make sure you specify a unique file name and try again."))
    elif not os.path.isdir(os.path.dirname(os.path.abspath(args.outputFileName))):
            FileNotFoundError((f"Output file '{args.outputFileName}' would be written to a non-existent directory" + 
                               "If you provide a full path, make sure its parent directories exist; " +
                               "otherwise, provide a file name only."))

def main():
    #### USER INPUT SECTION
    usage = """%(prog)s will receive a .geneious file and attempt
    to extract the sequences contained within to a FASTA file.
    There's a very good chance this script only works in a narrow range
    of scenarios when it comes to GENEIOUS outputs.
    """
    p = argparse.ArgumentParser(description=usage)
    # Required
    p.add_argument("-i", dest="geneiousFile",
                   required=True,
                   help="Specify input .geneious file location")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="""Specify location to write output FASTA file""")
    
    args = p.parse_args()
    validate_args(args)
    
    # Load GENEIOUS .xml into dict
    with ZipFile(args.geneiousFile) as zf:
        if not len(zf.namelist()) == 1:
            raise ValueError(("Geneious input file has multiple contents; I " +
                            "only know how to handle 1"))
        contents = zf.namelist()[0]
        with zf.open(contents) as fileIn:
            geneiousDict = xmltodict.parse(fileIn)
    
    # Parse out sequences as FASTAs
    with open(args.outputFileName, "w") as fileOut:
        for documentDict in geneiousDict["geneious"]["geneiousDocument"][1:]: # first value is irrelevant
            geneID = documentDict["originalElement"]["XMLSerialisableRootElement"]["name"]
            seq = documentDict["originalElement"]["XMLSerialisableRootElement"]["charSequence"]
            fileOut.write(f">{geneID}\n{seq}\n")

    print("Program completed successfully!")

if __name__ == "__main__":
    main()

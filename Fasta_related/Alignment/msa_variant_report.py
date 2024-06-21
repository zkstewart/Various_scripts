#! python3
# msa_variant_report.py
# Script for working with one or more MSA files
# to generate a report on any variants present.

import sys, argparse, os

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))) # 3 dirs up is where we find dependencies
from Function_packages import ZS_SeqIO, ZS_AlignIO

def validate_args(args):
    # Validate input data location
    if not os.path.isdir(args.msaDir):
        print('I am unable to locate the directory where the MSA FASTA files are (' + args.msaDir + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate logic of argument combination
    if args.asCodons and not args.isNucleotide:
        print("You've specified --asCodons without --isNucleotide; provide them together or not at all")
        quit()
    if args.translationTable != 1 and not args.asCodons:
        print("You've specified --translationTable without --asCodons; provide them together or not at all")
        quit()
    # Validate numeric arguments
    # if args.translationTable < 1 or args.translationTable > 31:
    #     print('The translation table must be an integer between 1 and 31')
    #     quit()
    # Handle file output
    if os.path.isfile(args.outputFileName):
        print(f'File already exists at output location ({args.outputFileName})')
        print('Make sure you specify a unique file name and try again.')
        quit()

def main():
    #### USER INPUT SECTION
    usage = """%(prog)s will receive one or more input MSA FASTA files located within a specified
    directory and identify any variants present in the alignment. In this context, a variant is
    a nucleotide or amino acid residue that deviates from the consensus sequence. An output
    TSV will be generated listing all variants found, along with their positions in the alignment.
    Note that input can be provided as nucleotide or amino acid sequences; if providing a nucleotide
    file, you can specify whether you'd like to consider codons (which get translated to amino acids)
    or individual positions.
    """
    p = argparse.ArgumentParser(description=usage)
    # Required
    p.add_argument("-i", dest="msaDir",
                   required=True,
                   help="Specify the location of the directory containing one or more MSA FASTA files")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output report file name")
    # Optional
    p.add_argument("--isNucleotide", dest="isNucleotide",
                   required=False,
                   action="store_true",
                   help="""Optionally provide this argument if the MSAs contain nucleotide
                   sequences; otherwise, we assume they are protein alignments""",
                   default=False)
    p.add_argument("--asCodons", dest="asCodons",
                   required=False,
                   action="store_true",
                   help="""Optionally provide this argument if you're providing a nucleotide
                   file and would like to consider codons as the unit of comparison""",
                   default=False)
    # p.add_argument("--translationTable", dest="translationTable",
    #                require=False,
    #                type=int,
    #                help="""Optionally specify the NCBI numeric genetic code to utilise for
    #                CDS translation if --isNucleotide and --asCodons; this should be an integer
    #                from 1 to 31 (default == 1 i.e., Standard Code)""",
    #                default=1)
    
    args = p.parse_args()
    validate_args(args)
    
    # Locate all FASTA files
    files = [os.path.join(args.msaDir, file) for file in os.listdir(args.msaDir)]
    
    # Iterate through each file and store variants
    variantDict = {}
    for file in files:
        filePrefix = os.path.splitext(os.path.basename(file))[0]
        
        # Load in the aligned FASTA file
        FASTA_obj = ZS_SeqIO.FASTA(file, isAligned=True)
        FASTA_obj.make_uppercase() # make sure comparison isn't case-sensitive
        
        # Locate variants in this MSA file
        variantDict[filePrefix] = ZS_AlignIO.MSA.locate_variants_from_msa(
            FASTA_obj, args.isNucleotide, args.asCodons
        )
    
    # Exit if no variants found
    if all([len(variantDict[file]) == 0 for file in variantDict]):
        print("No variants found in any of the provided MSA files; program will exit without writing output.")
        quit()
    
    # Otherwise, write output report
    else:
        with open(args.outputFileName, "w") as fileOut:
            # Write header
            fileOut.write("#gene\tposition_number\tconsensus_residue\tvariant_residue\tseqs_with_variant\n")
            
            # Write content lines
            for geneID, positionDict in variantDict.items():
                for position, detailsDict in positionDict.items():
                    for variantResidue, seqList in detailsDict["variants"].items():
                        fileOut.write(f"{geneID}\t{position+1}\t{detailsDict['consensus']}\t" +
                                      f"{variantResidue}\t{', '.join(seqList)}\n")
    
    # Done!
    print('Program completed successfully!')

if __name__ == "__main__":
    main()

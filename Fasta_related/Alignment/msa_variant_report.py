#! python3
# msa_variant_report.py
# Script for working with one or more MSA files
# to generate a report on any variants present.

import sys, argparse, os

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))) # 3 dirs up is where we find dependencies
from Function_packages import ZS_SeqIO, ZS_AlignIO

def validate_args(args):
    # Validate input data location
    args.files = []
    for location in args.msas:
        if os.path.isfile(location):
            args.files.append(location)
        elif os.path.isdir(location):
            args.files.extend([
                os.path.join(location, x)
                for x in os.listdir(location)
                if x.lower().endswith(tuple(args.suffixes))
            ])
        else:
            raise FileNotFoundError((f"I am unable to locate the MSA file or directory at '{location}'. " + 
                                    "Make sure you've typed the name or location correctly and try again."))
    # Validate logic of argument combination
    if args.asCodons and not args.isNucleotide:
        raise ValueError("You've specified --asCodons without --isNucleotide; provide them together or not at all")
    # Handle file output
    if args.reportFormat == "per_variant":
        if os.path.exists(args.outputLocation):
            raise FileExistsError((f"File already exists at output location ({args.outputLocation}). " + 
                                   "Make sure you specify a unique file name and try again."))
        elif not os.path.isdir(os.path.dirname(os.path.abspath(args.outputLocation))):
            FileNotFoundError((f"Output file '{args.outputLocation}' would be written to a non-existent directory" + 
                               "If you provide a full path, make sure its parent directories exist; " +
                               "otherwise, provide a file name only."))
    elif args.reportFormat == "per_sequence":
        if not os.path.isdir(os.path.dirname(os.path.abspath(args.outputLocation))):
            FileNotFoundError((f"Output file '{args.outputLocation}' would be written to a non-existent directory" + 
                               "If you provide a full path, make sure its parent directories exist; " +
                               "otherwise, provide a directory name only."))
        if os.path.isdir(args.outputLocation):
            print((f"# Output location '{args.outputLocation}' already exists. Be warned that if you're trying to " +
                   "resume a previous run, you should have deleted any flawed or incomplete files before " + 
                   "proceeding since files will NOT be overwritten."))
            print("# msa_variant_report.py execution will continue...")
        else:
            try:
                os.mkdir(args.outputLocation)
                print(f"# Output directory '{args.outputLocation}' created as part of argument validation.")
            except:
                raise FileNotFoundError((f"Unable to create output directory '{args.outputLocation}' for " + 
                                         "unknown reasons. Please check the path and try again."))

def report_per_variant(outputFileName, variantDict):
    # Write report file
    with open(outputFileName, "w") as fileOut:
        # Write header
        fileOut.write("#gene\tposition_number\tconsensus_residue\tvariant_residue\tseqs_with_variant\n")
        
        # Write content lines
        for fileName, pairedValue in variantDict.items():
            geneID = fileName.rsplit(".", maxsplit=1)[0]
            positionDict = pairedValue[0] # pairedValue = [positionDict, fullFileName]
            for position, detailsDict in positionDict.items():
                for variantResidue, seqList in detailsDict["variants"].items():
                    fileOut.write(f"{geneID}\t{position+1}\t{detailsDict['consensus']}\t" +
                                f"{variantResidue}\t{', '.join(seqList)}\n")

def report_per_sequence(outputDirectory, variantDict):
    for fileName, pairedValue in variantDict.items():
        geneID = fileName.rsplit(".", maxsplit=1)[0]
        positionDict, fullFileName = pairedValue
        
        # Establish dictionary for storing reformatted data
        FASTA_obj = ZS_SeqIO.FASTA(fullFileName, isAligned=True)
        reformattedDict = { x.id : [] for x in FASTA_obj.seqs }
        
        # Reformat variants for per-sequence reporting
        for position, detailsDict in positionDict.items():
            for variantResidue, seqList in detailsDict["variants"].items():
                for seqID in seqList:
                    reformattedDict[seqID].append(f"{position+1}:{detailsDict['consensus']}:{variantResidue}")
        
        # Derive report file name and skip if it exists
        reportFileName = os.path.join(outputDirectory, f"{geneID}.report.tsv")
        if os.path.exists(reportFileName):
            print(f"# WARNING: File already exists at '{reportFileName}' and will not be overwritten; skipping...")
            continue
        
        # Write report file
        with open(reportFileName, "w") as fileOut:
            # Write header
            fileOut.write("#sequence_id\tvariants\n")
            
            # Write content lines
            for seqID, variantList in reformattedDict.items():
                variants = ", ".join(variantList) if variantList != [] else "."
                fileOut.write(f"{seqID}\t{variants}\n")

def main():
    #### USER INPUT SECTION
    usage = """%(prog)s will receive one or more input MSA FASTA files specified by the -i argument.
    Inputs can be specified as individual files or as a directory containing multiple files.
    This script will identify any variants present in these aligned files. In this context, a variant is
    a nucleotide or amino acid residue that deviates from the consensus sequence. An output
    TSV will be generated listing all variants found, along with their positions in the alignment.
    Note that input can be provided as nucleotide or amino acid sequences; if providing a nucleotide
    file, you can specify whether you'd like to consider codons (which get translated to amino acids)
    or individual positions.
    """
    p = argparse.ArgumentParser(description=usage)
    # Required
    p.add_argument("-i", dest="msas",
                   required=True,
                   nargs="+",
                   help="Specify locations of MSA files and/or director(y/ies) containing MSA files")
    p.add_argument("-f", dest="reportFormat",
                   required=True,
                   choices=["per_variant", "per_sequence"],
                   help="Specify the format of the output report file(s)")
    p.add_argument("-o", dest="outputLocation",
                   required=True,
                   help="""Output file name (if -f == 'per_variant')
                   or directory (if 'per_sequence')""")
    # Optional
    p.add_argument("--suffixes", dest="suffixes",
                   required=False,
                   nargs="+",
                   help="""Optionally indicate the suffixes of the MSA files you'd like to
                   obtain from any directories specified in -i; default == 
                   '.fasta .fas .fa .faa .fna .aa .cds .trans .nucl'""",
                   default=[".fasta", ".fas", ".fa", ".faa", ".fna", ".aa", ".cds", ".trans", ".nucl"])
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
    p.add_argument("--reportUntilStop", dest="reportUntilStop",
                   required=False,
                   action="store_true",
                   help="""Optionally provide this argument if you do not want to see variants
                   reported for a sequence after a stop codon is encountered in that sequence
                   (variants will be reported in other sequences if they don't encounter that
                   stop codon)""",
                   default=False)
    
    args = p.parse_args()
    validate_args(args)
    
    # Iterate through each file and store variants
    variantDict = {}
    for file in args.files:
        # Load in the aligned FASTA file
        FASTA_obj = ZS_SeqIO.FASTA(file, isAligned=True)
        FASTA_obj.make_uppercase() # make sure comparison isn't case-sensitive
        
        # Locate variants in this MSA file
        variantDict[os.path.basename(file)] = [ZS_AlignIO.MSA.locate_variants_from_msa(
            FASTA_obj, args.isNucleotide, args.asCodons, reportUntilStop=args.reportUntilStop
        ), file] # store file location for report_per_sequence
    
    # Exit if no variants found
    if all([len(variantDict[file]) == 0 for file in variantDict]):
        raise ValueError("No variants found in any of the provided MSA files; program will exit without writing output.")
    
    # Otherwise, write output report
    else:
        if args.reportFormat == "per_variant":
            report_per_variant(args.outputLocation, variantDict)
        elif args.reportFormat == "per_sequence":
            report_per_sequence(args.outputLocation, variantDict)
    
    # Done!
    print('Program completed successfully!')

if __name__ == "__main__":
    main()

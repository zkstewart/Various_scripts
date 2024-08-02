#! python3
# filter_fasta_by_snp_num.py
# Script to reduce a FASTA down to just the contigs that
# have >N SNPs predicted in the VCF file

import os, argparse
from Bio import SeqIO

def validate_args(args):
    # Validate input data locations
    if not os.path.isfile(args.fastaFile):
        raise FileNotFoundError((f"I am unable to locate the FASTA file at '{args.fastaFile}'. " + 
                                "Make sure you've typed the name or location correctly and try again."))
    if not os.path.isfile(args.vcfFile):
        raise FileNotFoundError((f"I am unable to locate the VCF file at '{args.vcfFile}'. " + 
                                "Make sure you've typed the name or location correctly and try again."))
    # Handle numeric parameters
    if args.numSNPs < 0:
        raise ValueError("numSNPs must be a positive integer")
    # Handle file output
    if os.path.exists(args.outputLocation):
        raise FileExistsError((f"File already exists at output location ({args.outputLocation}). " + 
                                "Make sure you specify a unique file name and try again."))
    elif not os.path.isdir(os.path.dirname(os.path.abspath(args.outputLocation))):
        FileNotFoundError((f"Output file '{args.outputLocation}' would be written to a non-existent directory" + 
                            "If you provide a full path, make sure its parent directories exist; " +
                            "otherwise, provide a file name only."))

def main():
    usage = """%(prog)s will filter a FASTA file down to just the contigs that have >=--numSNPs
    predicted in the VCF file. The output will be in BED format. This is a necessary step prior
    to running ReLERNN.
    """
    # Establish main parser
    p = argparse.ArgumentParser(description=usage)
    
    # Set arguments shared by subparsers
    ## Required arguments
    p.add_argument("-f", dest="fastaFile",
                    required=True,
                    help="Specify the location of the input Euclidean distance file")
    p.add_argument("-v", dest="vcfFile",
                    required=True,
                    help="Specify the location of the genome FASTA file")
    p.add_argument("-o", dest="outputLocation",
                    required=True,
                    help="Output BED file")
    ## Opts
    p.add_argument("--numSNPs", dest="numSNPs",
                    required=False,
                    type=int,
                    help="""Optionally, indicate the minimum number of SNPs a contig
                    must have for it to be included in the output (default=250; recommended
                    for ReLERNN)""",
                    default=250)
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse VCF and count SNPs on each contig
    countsDict = {}
    with open(args.vcfFile, "r") as fileIn:
        for line in fileIn:
            if line.startswith("#") or line.startswith("CHROM"):
                continue
            line = line.strip().split("\t")
            contig = line[0]
            countsDict.setdefault(contig, 0)
            countsDict[contig] += 1
    
    # Check that at least 1 contig meets the threshold
    if not any([v >= args.numSNPs for v in countsDict.values()]):
        raise ValueError("No contigs in the FASTA file have >=--numSNPs SNPs predicted in the VCF file")
    
    # Write out the contigs that meet the threshold
    with open(args.fastaFile, "r") as fileIn, open(args.outputLocation, "w") as fileOut:
        for record in SeqIO.parse(fileIn, "fasta"):
            if countsDict.get(record.id, 0) >= args.numSNPs:
                fileOut.write(f"{record.id}\t0\t{len(record)}\n")
    
    # Give details on the output
    print("# filter_fasta_by_snp_num.py results:")
    print(f"# Num. contigs input: {len(countsDict)}")
    print(f"# Num. contigs output: {sum([1 for v in countsDict.values() if v >= args.numSNPs])}")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

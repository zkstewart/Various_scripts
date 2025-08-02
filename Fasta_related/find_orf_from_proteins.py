#! python3
# find_orf_from_proteins.py
# Takes matching FASTA files with mRNA and protein sequences,
# and derives the CDS using this information. Sequence IDs can
# differ between the mRNA and protein files, but they are otherwise
# expected to be matching and equivalently ordered.

import os, argparse
from Bio import SeqIO

STOP_CODONS = ["TAA", "TAG", "TGA"]

def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.mrnaFile):
        raise FileNotFoundError(f"-m '{args.mrnaFile} does not exist!")
    if not os.path.isfile(args.proteinFile):
        raise FileNotFoundError(f"-m '{args.proteinFile} does not exist!")
    # Validate that input files match
    for record1, record2 in iterate_paired_fastas(args.mrnaFile, args.proteinFile):
        pass # will naturally raise an error if the files don't have an equal number of sequences
    # Validate output file locations
    if os.path.exists(args.cdsFile):
        raise FileExistsError(f"-o '{args.cdsFile}' already exists! Move/rename/delete this file or set a different -o value.")

def iterate_paired_fastas(fasta1, fasta2):
    '''
    Parameters:
        fasta1 / fasta2 -- a string indicating the FASTA file location
    Yields:
        record1 / record2 -- Bio.SeqIO.Seq records for the two input FASTA files
    '''
    with open(fasta1, "r") as file1, open(fasta2, "r") as file2:
        records1 = SeqIO.parse(file1, "fasta")
        records2 = SeqIO.parse(file2, "fasta")
        for record1, record2 in zip(records1, records2):
            yield record1, record2

def main():
    usage = """%(prog)s will receive matching mRNA and protein FASTA files, and
    use this to derive the CDS from the mRNA sequences. The CDS will be output
    as its own FASTA file.
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-m", dest="mrnaFile",
                   required=True,
                   help="Specify the input mRNA FASTA file")
    p.add_argument("-p", dest="proteinFile",
                   required=True,
                   help="Specify the input protein FASTA file")
    p.add_argument("-o", dest="cdsFile",
                   required=True,
                   help="Specify the output CDS FASTA file name")
    # Optional
    p.add_argument("--proteinID", dest="proteinID",
                   required=False,
                   action="store_true",
                   help="""Set this flag to take the ID from the protein,
                   not the mRNA""",
                   default=False)
    
    args = p.parse_args() # sets .conversionFile; .fastaFile
    validate_args(args)
    
    with open(args.cdsFile, "w") as fileOut:
        for mrnaRecord, proteinRecord in iterate_paired_fastas(args.mrnaFile, args.proteinFile):
            proteinSeq = str(proteinRecord.seq)
            
            # Find strand and frame where translation was obtained
            found = False
            for strand, nuc in [(+1, mrnaRecord.seq), (-1, mrnaRecord.seq.reverse_complement())]:
                for frame in range(3):
                    length = 3 * ((len(mrnaRecord)-frame) // 3)
                    frameNuc = str(nuc[frame:frame+length])
                    frameProt = str(nuc[frame:frame+length].translate(table=1))
                    if proteinSeq in frameProt:
                        found = True
                        
                        # Extract just the region of the protein
                        startPosition = frameProt.find(proteinSeq) * 3
                        endPosition = startPosition + (len(proteinSeq) * 3)
                        cdsSeq = frameNuc[startPosition:endPosition]
                        
                        # See if we can extend for a stop codon
                        if proteinSeq[-1] != "*":
                            extraCodon = frameNuc[endPosition:endPosition+3].upper()
                            if extraCodon in STOP_CODONS:
                                cdsSeq += extraCodon
                        break
                if found == True:
                    break
            if found == False:
                raise ValueError(f"Could not find '{proteinSeq}' within the mRNA '{str(mrnaRecord.seq)}'")
            
            # Write to file
            outputID = proteinRecord.description if args.proteinID else mrnaRecord.description
            fileOut.write(f">{outputID}\n{cdsSeq}\n")
            
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

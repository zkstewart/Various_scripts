#! python3
# crate_single_copy_MSA.py
# Script intended to work with the Single_Copy_Orthologue_Sequences
# output directory from OrthoFinder and create an aligned MSA
# for later phylogenetic analysis. It can also work with things
# other than OrthoFinder output though.

import sys, argparse, os, platform
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) # 2 dirs up is where we find dependencies
from Function_packages import ZS_SeqIO, ZS_AlignIO

def validate_args(args):
    # Validate input data location
    if not os.path.isdir(args.singleCopyOrthDir):
        print('I am unable to locate the directory where the single copy FASTA files are (' + args.singleCopyOrthDir + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isdir(args.mafftDir):
        print('I am unable to locate the directory where the MAFFT executables are (' + args.alignmentsDir + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    else:
        if platform.system() == "Windows":
            if not os.path.isfile(os.path.join(args.mafftDir, "mafft.bat")):
                raise Exception("{0} does not exist".format(os.path.join(args.mafftDir, "mafft.bat")))
        else:
            if not os.path.isfile(os.path.join(args.mafftDir, "mafft")) and not os.path.isfile(os.path.join(args.mafftDir, "mafft.exe")):
                raise Exception("mafft or mafft.exe does not exist at {0}".format(args.mafftDir))
    # Validate numeric inputs
    if args.threads < 1:
        print("threads arg needs to be greater than zero")
        quit()
    if not (0 <= args.noninfo_pct <= 1):
        print("noninfo_pct must be in the range of 0 to 1 (inclusive)")
        quit()
    # Handle file output
    if os.path.isfile(args.outputFileName):
        print(args.outputFileName + ' already exists. Specify a new file name or move/rename the existing file.')
        quit()

def trim_noninformative_flanks(FASTA_obj, isNucleotide=True, ALLOWED_NONINFO_PCT=0.25):
    '''
    Trims a ZS_SeqIO.FASTA object to remove non-informative flanks i.e.,
    sequence that is only gap ("-") or unknown ("N"/"X") across
    entire columns on the left and right sides of the MSA.
    
    INFO_DROP_PCT is a little non-intuitive, so I'll try to explain it.
    In short, this function will calculate the percentage of positions
    are non-informative. As soon as we reach a position that has a non-informative
    percentage LESS THAN the provided value, we'll stop trimming. If all
    positions are informative, the percentage calculated will be 0. If all positions
    are non-informative, the percentage will instead be 1. So in essence, if you want
    to stop trimming when 75% of sequences have informative positions, you should
    set the value to be 0.25.
    
    Params:
        FASTA_obj -- a ZS_SeqIO.FASTA instance
        isNucleotide -- a boolean indicating whether the MSA contains nucleotide sequences
                        or, if False, it is assumed to be amino acid sequence.
        ALLOWED_NONINFO_PCT -- a float value indicating what percentage of the MSA members
                               are allowed to be non-informative without being trimmed. If left
                               at the default 0.25, then 75% or more of the sequences must
                               be informative to avoid having the head or tail position trimmed.
    Returns:
        pctTrimmed -- a float value indicating what proportion of the alignment was trimmed.
    '''
    # Validate input types and value ranges
    assert isinstance(ALLOWED_NONINFO_PCT, float) or isinstance(ALLOWED_NONINFO_PCT, int)
    assert 0 <= ALLOWED_NONINFO_PCT <= 1, "ALLOWED_NONINFO_PCT must be between 0 or 1 (inclusive of 0 and 1)"
    assert isinstance(isNucleotide, bool)
    
    # Specify ambiguous character depending on sequence type
    if isNucleotide:
        ambiguous = "n"
    else:
        ambiguous = "x"
    
    # Trim the start
    startTrim = 0
    for i in range(len(FASTA_obj[0].gap_seq)):
        isNonInformative = [FastASeq_obj.gap_seq[i].lower() in [ambiguous, "-"] for FastASeq_obj in FASTA_obj if FastASeq_obj.id != "Codons"]
        nonInfoPct = sum(isNonInformative) / len(isNonInformative)
        
        if nonInfoPct > ALLOWED_NONINFO_PCT:
            startTrim += 1
        else:
            break
    
    # Trim the end
    endTrim = 0
    for i in range(len(FASTA_obj[0].gap_seq)-1, -1, -1):
        isNonInformative = [FastASeq_obj.gap_seq[i].lower() in [ambiguous, "-"] for FastASeq_obj in FASTA_obj if FastASeq_obj.id != "Codons"]
        nonInfoPct = sum(isNonInformative) / len(isNonInformative)
        
        if nonInfoPct > ALLOWED_NONINFO_PCT:
            endTrim += 1
        else:
            break
    
    # Calculate statistics
    pctTrimmed = (startTrim + endTrim) / len(FASTA_obj[0].gap_seq)
    
    # Trim to boundaries
    FASTA_obj.trim_left(startTrim, asAligned=True)
    FASTA_obj.trim_right(endTrim, asAligned=True)
    
    return pctTrimmed

def replace_nonstandard_aminoacids(FASTA_obj):
    '''
    Replaces amino acid characters that aren't within the standard list of amino acids.
    Prevents issues with MAFFT alignment and later phylogenetic analysis.
    
    Note that it will make all sequences upper case.
    '''
    STANDARD_AMINO_ACIDS = [
        "A", "C", "D", "E", "F",
        "G", "H", "I", "K", "L",
        "M", "N", "P", "Q", "R",
        "S", "T", "V", "W", "Y",
    ]
    for FastASeq_obj in FASTA_obj:
        if FastASeq_obj.gap_seq != None:
            seq = FastASeq_obj.gap_seq.upper()
        else:
            seq = FastASeq_obj.seq.upper()
        
        for i in range(len(seq)):
            letter = seq[i]
            if letter != "-" and letter not in STANDARD_AMINO_ACIDS:
                seq = seq[0:i] + "X" + seq[i+1:]
        
        if FastASeq_obj.gap_seq != None:
            FastASeq_obj.gap_seq = seq
            FastASeq_obj.seq = seq.replace("-", "")
        else:
            FastASeq_obj.seq = seq

if __name__ == "__main__":
    #### USER INPUT SECTION
    usage = """%(prog)s will take all the files inside the Single_Copy_Orthologue_Sequences
    directory produced by OrthoFinder and create a concatenated, aligned MSA suitable for
    phylogenetic analysis. Realistically, this will also work for any directory (OrthoFinder
    -produced or not) that contains a single file for every gene to be concenated into a single
    copy gene sequences MSA.
    
    You can set new IDs for each sequences through the --new_ids argument, but note that this
    needs to be ordered the same as the FASTA files themselves. For that matter, the FASTA files
    should ALL be ordered IDENTICALLY. It's kind of important.
    """

    p = argparse.ArgumentParser(description=usage)
    # Required
    p.add_argument("-i", dest="singleCopyOrthDir", required=True,
                help="Specify the location of the directory containing single copy gene FASTA files")
    p.add_argument("-m", dest="mafftDir", required=True,
                help="Specify the directory where MAFFT executables are located")
    p.add_argument("-o", dest="outputFileName", required=True,
                help="Specify the location to write output files to; this location will be populated with orthogroup files")
    # Optional
    p.add_argument("-t", dest="threads", type=int, default=1,
                help="Optionally specify the number of threads to run MAFFT alignment with (default == 1).")
    p.add_argument("--new_ids", dest="newIDsList", nargs="+", default=None,
                help="Optionally, specify multiple space-separated IDs to rename sequences to")
    p.add_argument("--noninfo_pct", dest="noninfo_pct", type=float, default=0.25,
                help="""Optionally specify where trimming should occur on the basis of the
                allowed amount of sequence that can be non-informative i.e., be gap
                or ambiguous position (default == 0.25).""")
    p.add_argument("--is_nucleotide", dest="is_nucleotide", action="store_true", default=False,
                help="""Optionally provide this argument if the MSAs contain nucleotide
                sequences; otherwise, we assume they are protein alignments""")
    args = p.parse_args()
    validate_args(args)
    
    # Locate all FASTA files
    files = [os.path.join(args.singleCopyOrthDir, file) for file in os.listdir(args.singleCopyOrthDir)]
    
    # Load FASTA files
    fastaObjs = []
    for file in files:
        f = ZS_SeqIO.FASTA(file)
        f.make_uppercase()
        fastaObjs.append(f)
    
    # Rename if applicable
    if args.newIDsList != None:
        for FASTA_obj in fastaObjs:
            FASTA_obj.set_ids_via_list(args.newIDsList)
    
    # Replace characters that don't conform to standard amino acid coding
    if not args.is_nucleotide:
        for FASTA_obj in fastaObjs:
            replace_nonstandard_aminoacids(FASTA_obj)
    
    # Align FASTA objects
    mafftAligner = ZS_AlignIO.MAFFT(args.mafftDir) # set up here for use later
    mafftAligner.set_threads(args.threads)
    for FASTA_obj in fastaObjs:
        mafftAligner.run(FASTA_obj)
    
    # Perform quality trimming of MSAs
    for FASTA_obj in fastaObjs:
        trim_noninformative_flanks(FASTA_obj, isNucleotide=args.is_nucleotide, ALLOWED_NONINFO_PCT=args.noninfo_pct)
    
    # Concatenate FASTA objects
    for x in range(1, len(fastaObjs)):
        fastaObjs[0].concat(fastaObjs[x]) # it all gets concatenated into the first FASTA object
    
    # Write output
    fastaObjs[0].write(args.outputFileName, asAligned=True)
    
    # Done!
    print('Program completed successfully!')

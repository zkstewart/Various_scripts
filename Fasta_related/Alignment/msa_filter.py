#! python3
# msa_filter.py
# Script for working with one or more MSA files
# to filter rows and columns based on provided criteria

import sys, argparse, os
import concurrent.futures
from itertools import repeat

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))) # 3 dirs up is where we find dependencies
from Function_packages import ZS_SeqIO

def validate_args(args):
    # Validate input data location
    if not os.path.isdir(args.msaDir):
        print('I am unable to locate the directory where the MSA FASTA files are (' + args.msaDir + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate numeric inputs
    if not (0 <= args.noninfo_pct <= 1):
        print("noninfo_pct must be in the range of 0 to 1 (inclusive)")
        quit()
    if args.threads < 1:
        print("threads argument must be a positive integer greater than 0")
        quit()
    # Ensure paired args are provided
    if args.filterID != None and len(args.filterValues) == 0:
        print("If filterID is specified, you must also specify one or more filterValues")
        quit()
    # Handle file output
    if not os.path.isdir(args.outputDir):
        try:
            os.mkdir(args.outputDir)
            print("Created '{0}' directory as part of argument validation".format(args.outputDir))
        except:
            print("Wasn't able to create '{0}' directory; does '{1}' actually exist?".format(args.outputDir, os.path.dirname(args.outputDir)))

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

def drop_FASTA_rows_by_id(FASTA_obj, rowsToDrop):
    '''
    Function to take a ZS_SeqIO.FASTA object and, based on the IDs in the list
    rowsToDrop, remove any rows that match those IDs. This will modify
    the FASTA object in place.
    
    Parameters:
        FASTA_obj -- a ZS_SeqIO.FASTA object.
        rowsToDrop -- a list containing strings corresponding to the .id values
                      within FASTA_obj.
    '''
    rowIndices = [i for i in range(len(FASTA_obj)) if FASTA_obj[i].id in rowsToDrop]
    assert len(rowIndices) == len(rowsToDrop), \
        "Not all values in rowsToDrop were discovered in the FASTA object!"
    
    for index in rowIndices[::-1]:
        del FASTA_obj.seqs[index]

def drop_aligned_FASTA_columns_by_value(FASTA_obj, filterRowID, filterValues):
    '''
    Receives an aligned ZS_SeqIO.FASTA object and eliminates any columns containing
    any values in the filterValues list that are found within the sequence identified
    by filterRowID.
    
    A practical example of its use is when a sequence is in the MSA called "Codons", and
    this line contains values like 4 or 5 to indicate poor quality sequence, and 1/2/3 for
    the codon numbers. In this case, filterRowID == "Codons" and filterValues == ["4", "5"].
    
    Parameters:
        FASTA_obj -- a ZS_SeqIO.FASTA object.
        filterRowID -- a string corresponding to a sequence in FASTA_obj by its
                       .id attribute.
        filterValues -- a list containing strings corresponding to values within the
                        sequence (identified by filterRowID)'s .gap_seq attribute.
    Returns:
        FASTA_obj -- a modified deepcopy of the original input FASTA_obj.
    '''
    assert filterRowID in FASTA_obj.ids, \
        "{0} could not be found within the provided FASTA object!".format(filterRowID)
    assert FASTA_obj.isAligned is True, \
        "FASTA_obj must be aligned to be used by this function!"
    filterFastASeq_obj = FASTA_obj[filterRowID]
    
    # Figure out the coordinate ranges to drop using an open:shut algorithm
    open = False
    coords = []
    for i in range(len(filterFastASeq_obj.gap_seq)):
        character = filterFastASeq_obj.gap_seq[i]
        
        if character in filterValues:
            if not open:
                coords.append([i, None])
                open = True
        else:
            if open:
                coords[-1][1] = i
                open = False
    if len(coords) != 0 and coords[-1][1] == None:
        coords[-1][1] = i
    
    # Drop columns in reverse
    for start, end in coords[::-1]:
        FASTA_obj = FASTA_obj.cut_cols(start, end)
    
    return FASTA_obj

def process_handler(args, fileName, FASTA_obj):
    # Filter columns if relevant
    if args.filterID != None:
        FASTA_obj = drop_aligned_FASTA_columns_by_value(FASTA_obj, args.filterID, args.filterValues)
    
    # Drop any rows if relevant
    if len(args.rowsToDrop) != 0:
        drop_FASTA_rows_by_id(FASTA_obj, args.rowsToDrop)
    
    # Perform quality trimming of MSA
    trim_noninformative_flanks(FASTA_obj, isNucleotide=args.is_nucleotide, ALLOWED_NONINFO_PCT=args.noninfo_pct)
    
    # Replace characters that don't conform to standard amino acid coding
    if not args.is_nucleotide:
        replace_nonstandard_aminoacids(FASTA_obj)
    
    # Write FASTA to file
    outputFileName = os.path.join(args.outputDir, os.path.basename(fileName))
    FASTA_obj.write(outputFileName, asAligned=True)

def main():
    #### USER INPUT SECTION
    usage = """%(prog)s offers two main capabilities. Firstly, it can drop rows based on their ID. 
    It requires an exact match, and multiple rows can be specified (or none). Secondly, it can drop
    columns based on their value. Specifically, it requires a single row's ID to be specified, and
    one or more values contained in that row which will cause the corresponding column to be removed.
    Row IDs and filter values provided are CASE-SENSITIVE! Outputs will be written to the provided 
    directory with names identical to the input. File overwrites will NOT be permitted and this program
    end if this is detected
    """
    p = argparse.ArgumentParser(description=usage)
    # Required
    p.add_argument("-i", dest="msaDir", required=True,
                help="Specify the location of the directory containing one or more MSA FASTA files")
    p.add_argument("-o", dest="outputDir", required=True,
                help="Output directory location")
    # Optional
    p.add_argument("--threads", dest="threads", required=False, type=int, default=1,
                help="Optionally specify how many threads to run with")
    p.add_argument("--rows_to_drop", dest="rowsToDrop", nargs="+", default=[],
                help="Optionally, specify one or more space-separated IDs to drop these sequences from the MSA")
    p.add_argument("--filter_id", dest="filterID", default=None,
                help="""Optionally specify a single row ID which contains values that will mark columns to be dropped""")
    p.add_argument("--filter_values", dest="filterValues", nargs="+", default=[],
                help="Optionally, specify one or more space-separated values that mark columns to be dropped")
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
    files = [os.path.join(args.msaDir, file) for file in os.listdir(args.msaDir)]
    
    # Stop processing if file overwriting would be necessary
    for file in files:
        if os.path.isfile(os.path.join(args.outputDir, os.path.basename(file))):
            print("'{0}' file already exists in the output directory; program will end now".format(os.path.basename(file)))
            quit()
    
    # Load FASTA files
    fastaObjs = []
    for file in files:
        f = ZS_SeqIO.FASTA(file, isAligned=True)
        f.make_uppercase()
        fastaObjs.append(f)
    
    # Ensure that any provided row IDs actually exist
    if len(args.rowsToDrop) != 0:
        dropsCount = sum([FastASeq_obj.id in args.rowsToDrop for FASTA_obj in fastaObjs for FastASeq_obj in FASTA_obj])
        if dropsCount != (len(args.rowsToDrop) * len(fastaObjs)):
            print("Not all FASTA files contain the rowsToDrop values; we can't process this")
    
    if args.filterID != None:
        filterCount = sum([FastASeq_obj.id == args.filterID for FASTA_obj in fastaObjs for FastASeq_obj in FASTA_obj])
        if filterCount != len(fastaObjs):
            print("Not all FASTA files contain the rowsToDrop values; we can't process this")
    
    # Filtering
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as executor:
        executor.map(process_handler, repeat(args, len(files)), files, fastaObjs)
    
    # Done!
    print('Program completed successfully!')

if __name__ == "__main__":
    main()

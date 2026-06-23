#!/usr/bin/env python3

import argparse
import codecs
import gzip
import os
import pandas as pd

from Bio.SeqIO.QualityIO import FastqGeneralIterator
from collections import Counter
from copy import deepcopy
from contextlib import contextmanager
from skbio.alignment import pair_align
from xlsxwriter.utility import xl_col_to_name

COMPLEMENT = {
    "A": "T", "a": "t",
    "T": "A", "t": "a",
    "C": "G", "c": "g",
    "G": "C", "g": "c"
}

def align(queryString, targetString, mode):
    '''
    Use the skbio pairwise alignment function to align two sequences
    to each other.
    
    Parameters:
        queryString / targetString -- a string to be used as the query or target
        mode -- a string equal to "global" or "local" to influence the type
                of alignment to perform
    Returns:
        score -- a float indicating the alignment score
        (queryAlign, targetAlign) -- multiple sequence alignment formatted strings
    '''
    score, paths, matrices = pair_align(queryString, targetString, mode=mode)
    return score, paths[0].to_aligned([queryString, targetString])

def alignment_check(queryAlign, targetAlign, mismatchTolerance=1):
    '''
    Check whether the alignment occurs without indels (- gaps) and with
    tolerance for 1bp mismatch.
    
    Parameters:
        queryAlign / targetAlign -- ultiple sequence alignment formatted strings
        mismatchTolerance -- an integer indicating how many sequence mismatches
                             can be tolerated before discarding the alignment as
                             "not good enough"
    Returns:
        isGood -- a boolean indicating whether the query>target alignment
                  meets the criteria of being indel-free and with 1bp mismatch
                  allowance
        (qStart, qEnd) -- integers indicating the position in the query where
                          the alignment starts and ends
        tSegment -- a string giving the segment of the target sequence where
                    the alignment occurred
    '''    
    qUngapped = queryAlign.replace("-", "")
    qStart = queryAlign.find(qUngapped[0])
    qEnd = len(queryAlign) - queryAlign[::-1].find(qUngapped[-1])
    
    qSegment = queryAlign[qStart:qEnd]
    tSegment = targetAlign[qStart:qEnd]
    
    if "-" in qSegment or "-" in tSegment:
        return False, (None, None), None
    
    differences = len(qSegment) - sum([ x==y for x,y in zip(qSegment,tSegment)])
    if differences > mismatchTolerance:
        return False, (None, None), None
    
    return True, (qStart, qEnd), tSegment

def find_index(leftFlank, rightFlank, readSeq, mismatchTolerance=1, indexLength=8):
    '''
    Pipeline function for identifying an index sequence between
    a left and right flanking sequence.
    
    The read sequence can be provided in its standard (as sequenced)
    or reverse complement orientation.
    
    Note that this function returns the lSegment and rSegment values
    which is a holdover from deprecated program behaviour. The values
    are still emitted just in case it enables future behaviours.
    
    Parameters:
        leftFlank -- the left flanking sequence as a string OR None
        rightFlank -- the right flanking sequence as a string OR None
        readSeq -- the sequenced read as a string
        mismatchTolerance -- an integer indicating how many basepairs can differ
                             (default == 1)
        indexLength -- an integer indicating how many basepairs to obtain as the index
                       IF only a single flank is given; if both left and right flanks exist,
                       the index is taken to be the gap sequence regardless of length
    Returns:
        indexSeq -- a string of the i5 or i7 index, OR None if no
                    valid match to the flanks occurred
        lSegment -- a string of the left sequence as found within the
                    sequenced read, OR None if no valid match occurred
        rSegment -- a string of the right sequence as found within the
                    sequenced read, OR None if no valid match occurred
    '''
    # Set defaults
    "If the flank is None, we should default the 'is it good' statement to 'yes'"
    lIsGood = True
    rIsGood = True
    
    # Check left and right flank alignment correctness
    if not leftFlank is None:
        lScore, lAlign = align(leftFlank, readSeq, mode="global")
        lIsGood, (lStart, lEnd), lSegment = alignment_check(*lAlign, mismatchTolerance=mismatchTolerance)
    if not rightFlank is None:
        rScore, rAlign = align(rightFlank, readSeq, mode="global")
        rIsGood, (rStart, rEnd), rSegment = alignment_check(*rAlign, mismatchTolerance=mismatchTolerance)
    if not (lIsGood and rIsGood):
        return None, None, None
    
    # Obtain results contigent on flanks provided
    if not (leftFlank is None or rightFlank is None): # both flanks are provided
        # Check that the arrangement of sequence segments is as expected
        lSide = "left" if lStart < rStart else "right" # is the "left flank" actually to the left of the gap?
        if lSide != "left":
            return None, None, None
        
        isGap = lEnd < rStart if lSide == "left" else lStart > rEnd # is there a gap between the left and right flank alignments?
        if not isGap:
            return None, None, None

        # If both flanks are given, obtain the gap as the putative index
        gapStart = readSeq.find(lSegment) + len(lSegment)
        gapEnd = readSeq.find(rSegment)
        indexSeq = readSeq[gapStart:gapEnd]
        return indexSeq, lSegment, rSegment
    else: # only one flank was provided
        # Handle index finding from the left flank
        if not leftFlank is None:
            gapStart = readSeq.find(lSegment) + len(lSegment)
            gapEnd = gapStart + indexLength
            indexSeq = readSeq[gapStart:gapEnd]
            return indexSeq, lSegment, None # no right segment
        # Handle index finding from the right flank
        else:
            gapEnd = readSeq.find(rSegment)
            gapStart = gapEnd - indexLength
            indexSeq = readSeq[gapStart:gapEnd]
            return indexSeq, rSegment, None # no right segment

def reverse_complement(sequence):
    '''
    Parameters:
        sequence -- a nucleotide sequence as a string
    Returns:
        sequencerc -- the reverse complement of the input nucleotide
    '''
    return "".join(COMPLEMENT.get(nuc, nuc) for nuc in reversed(sequence))

def fillup_blanks(_counter, NUM_TO_COUNT):
    '''
    Makes sure the list returned by a Counter object
    has at least NUM_TO_COUNT values within it, to enable
    pandas DataFrame creation. Blank values will be imputed
    as "NA" with a count of 0.
    
    Parameters:
        _counter -- a collections.Counter object
        NUM_TO_COUNT -- an integer giving the number of
                        values needed to be listed.
    Returns:
        _counts -- a list of tuples with format like:
                   [
                       ("sequence1": count1),
                       ("sequence2": count2),
                       ...
                   ]
    '''
    _counts = _counter.most_common(NUM_TO_COUNT)
    if len(_counts) < NUM_TO_COUNT:
        for _ in range(NUM_TO_COUNT - len(_counts)):
            _counts.append(("NA", 0))
    return _counts

def validate_args(args):
    '''
    Check validity of argparse values.
    '''
    # Validate input file locations depending on type of input
    args.files = []
    for inputLocation in args.inputLocations:
        inputLocation = os.path.abspath(inputLocation)
        
        # Handle directories
        if os.path.isdir(inputLocation):
            foundFiles = False
            for f in os.listdir(inputLocation):
                file = os.path.join(inputLocation, f)
                if os.path.isfile(file) and file.endswith(args.fileSuffix):
                    args.files.append(file)
                    foundFiles = True
            if not foundFiles:
                raise FileNotFoundError(f"'{inputLocation}' does not contain any files.")
        
        # Handle files
        elif os.path.isfile(inputLocation) and file.endswith(args.fileSuffix):
            args.files.append(inputLocation)
        else:
            raise FileNotFoundError(f"'{inputLocation}' is not a valid file or directory, ")
    
    # Validate flanking sequences
    if len(args.leftFlank) == 0 and len(args.rightFlank) == 0: # no flank is given
        raise ValueError(f"At least one value must be given to --left and/or --right")
    
    for flankList in [args.leftFlank, args.rightFlank]:
        for flank in flankList:
            if len(flank) < 2:
                raise ValueError(f"Flanking sequence '{flank}' is too short to make sense!")
    
    if len(args.leftFlank) > 0 and len(args.rightFlank) > 0: # both --left and --right are given
        if len(args.leftFlank) != len(args.rightFlank):
            raise ValueError(f"Number of values given to -l does not match the number given to -r !")
        if len(args.leftFlank) == 0:
            raise ValueError("At least one value must be given to -l !")
    else: # only --left or --right is given
        if len(args.leftFlank) > 0:
            args.rightFlank = [None]*len(args.leftFlank) # pad out with Nones
        else:
            args.leftFlank = [None]*len(args.rightFlank) # pad out with Nones
    
    # Validate numeric arguments
    if args.indexLength < 1:
        raise ValueError(f"--length must be given a value of 1bp or longer to make any sense!")
    if args.numMismatches < 0:
        raise ValueError(f"--mismatches must be given a value >= 0")
    
    # Validate output file location
    args.outputFileName = os.path.abspath(args.outputFileName)
    if not args.outputFileName.endswith(".xlsx"):
        print("-o file name did not end in the expected .xlsx suffix; this has been automatically appended")
        args.outputFileName += ".xlsx"
    
    if os.path.exists(args.outputFileName):
        raise FileExistsError(f"File already exists at output location ({args.outputFileName}). " +
                              "Make sure you specify a unique file name and try again.")
    elif not os.path.isdir(os.path.dirname(args.outputFileName)):
        raise FileNotFoundError(f"Output directory does not exist ({os.path.dirname(args.outputFileName)}). " +
                                "Make sure you specify a parent directory for your output file which already exists.")

def get_paired_files(files, fileSuffix):
    '''
    Parameters:
        files -- a list of strings for all files to be paired up
        fileSuffix -- a string allowing us to locate the value '1' or '2'
                      immediately preceding the file ending
    Returns:
        forwardReads / reverseReads -- a list of strings of file names ordered to
                                       pair up with the corresponding list
    '''
    # Locate files from the directory
    forwardReads = []
    reverseReads = []
    for file in files:
        if file.endswith(f"1{fileSuffix}"):
            forwardReads.append(file)
        elif file.endswith(f"2{fileSuffix}"):
            reverseReads.append(file)
        else:
            raise ValueError(f"{file} ends with the expected suffix '{fileSuffix}' but is not preceeded by a 1 or 2!")
    forwardReads.sort()
    reverseReads.sort()
    
    # Validate that paired files match
    assert len(forwardReads) == len(reverseReads), \
        f"Number of reads don't match for forward ({len(forwardReads)}) and reverse ({len(reverseReads)}) files"
    for i in range(len(forwardReads)):
        prefix = os.path.commonprefix([forwardReads[i], reverseReads[i]])
        assert prefix != "", \
            "forward and reverse read pairs don't have a common prefix?"
        assert forwardReads[i].startswith(f"{prefix}1") and reverseReads[i].startswith(f"{prefix}2"), \
            f"forward and reverse reads don't start with a recognised prefix ({prefix} should preceed a 1 or 2)"
    
    # Return files
    return forwardReads, reverseReads if reverseReads != [] else None

def get_codec(fileName):
    '''
    Checks the input file for its UTF encoding.
    
    Returns:
        codecValue -- a string of either "utf-8" or "utf-16"
    '''
    try:
        f = codecs.open(fileName, encoding='utf-8', errors='strict')
        for line in f:
            break
        f.close()
        return "utf-8"
    except:
        try:
            f = codecs.open(fileName, encoding='utf-16', errors='strict')
            for line in f:
                break
            f.close()
            return "utf-16"
        except UnicodeDecodeError:
            print(f"Can't tell what codec '{fileName}' is!!")

@contextmanager
def read_gz_file(filename):
    '''
    Flexibly open a file in "with <...> as variable" context with handling
    of either .gz or standard (uncompressed) formatting.
    '''
    if filename.endswith(".gz"):
        with gzip.open(filename, "rt") as f:
            yield f
    else:
        with open(filename, "r", encoding=get_codec(filename)) as f:
            yield f

def main():
    usage = """%(prog)s will take one or more pairs of FASTQ files and locate
    sequences (such as the Illumina P5/7 and TruSeq reads 1/2) expected to flank an index
    sequence. Returns tabulated results in Excel format.
    """
    # Required arguments
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="inputLocations",
                   required=True,
                   nargs="+",
                   help="Location of file(s) or directories containing FASTQ files to parse through")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file name for the resulting Excel document")
    p.add_argument("--left", dest="leftFlank",
                   required=False,
                   nargs="+",
                   help="One or more sequences expected to occur to the left of the index",
                   default=[])
    p.add_argument("--right", dest="rightFlank",
                   required=False,
                   nargs="+",
                   help="""One or more sequences expected to occur to the right of the index; input
                   should be given as pairs to the -l values""",
                   default=[])
    # Optional arguments
    p.add_argument("--top", dest="numTopIndices",
                   required=False,
                   type=int,
                   help="""Optionally, specify the number of top indices to report
                   (default: 20)""",
                   default=20)
    p.add_argument("--mismatches", dest="numMismatches",
                   required=False,
                   type=int,
                   help="""Optionally, specify the maximum number of mismatches allowed
                   (default: 1)""",
                   default=1)
    p.add_argument("--length", dest="indexLength",
                   required=False,
                   type=int,
                   help="""Optionally, specify the length of the index sequence that borders
                   a single flanking sequence
                   (default: 8)""",
                   default=8)
    p.add_argument("--suffix", dest="fileSuffix",
                   required=False,
                   help="""Optionally, specify the suffix which uniquely identifies all relevant
                   read files; this must immediately precede the number '1' or '2' for the
                   forward and reverse files, respectively (default: '.fq.gz')""",
                   default=".fq.gz")
    
    args = p.parse_args()
    validate_args(args)
    
    # Pair up files
    forwardReads, reverseReads = get_paired_files(args.files, args.fileSuffix)
    print(f"# Detected {len(forwardReads)} file pairs for processing")
    
    # Iterate over sample files to locate indices
    DEFAULT_DICT = {
        1: { "standard": { i: {} for i in range(len(args.leftFlank)) },
             "revcomp": { i: {} for i in range(len(args.leftFlank)) }
        },
        2: { "standard": { i: {} for i in range(len(args.leftFlank)) },
             "revcomp": { i: {} for i in range(len(args.leftFlank)) }
        }
    }
    
    sampleIndices = {}
    #sampleFlanks = {}
    for forwardFile, reverseFile in zip(forwardReads, reverseReads):
        sampleID = os.path.commonprefix([os.path.basename(forwardFile), os.path.basename(reverseFile)]).rstrip("._ ")
        
        sampleIndices[sampleID] = deepcopy(DEFAULT_DICT)
        #sampleFlanks[sampleID] = deepcopy(DEFAULT_DICT)
        
        for pairIndex, readFile in zip([1, 2], [forwardFile, reverseFile]):
            with read_gz_file(readFile) as fileIn:
                for title, seq, qual in FastqGeneralIterator(fileIn):
                    rcseq = reverse_complement(seq)
                    
                    for flankIndex, (leftFlank, rightFlank) in enumerate(zip(args.leftFlank, args.rightFlank)):
                        # Align the flanks to the sequence
                        index, foundLeft, foundRight = find_index(leftFlank, rightFlank, seq,
                                                                  mismatchTolerance=args.numMismatches,
                                                                  indexLength=args.indexLength)
                        foundSeq, foundOrientation = seq, "standard" # set default value
                        if index is None:
                            index, foundLeft, foundRight = find_index(leftFlank, rightFlank, rcseq,
                                                                      mismatchTolerance=args.numMismatches,
                                                                      indexLength=args.indexLength)
                            foundSeq, foundOrientation = rcseq, "revcomp" # reset as revcomp values
                        
                        # Store result if a valid match occurred
                        if not index is None:
                            sampleIndices[sampleID][pairIndex][foundOrientation][flankIndex].setdefault(
                                index, 0
                            )
                            sampleIndices[sampleID][pairIndex][foundOrientation][flankIndex][index] += 1
                            
                            # sampleFlanks[sampleID][pairIndex][foundOrientation][flankIndex].setdefault(
                            #     arrangement, 0
                            # )
                            # sampleFlanks[sampleID][pairIndex][foundOrientation][flankIndex][arrangement] += 1
    
    # Format results into a table
    writer = pd.ExcelWriter(args.outputFileName, engine = "xlsxwriter")
    for flankIndex in range(len(args.leftFlank)):
        # Format results for this flank into a list
        dfList = []
        for sampleID, pairIndexDict in sampleIndices.items():
            for pairIndex, orientationDict in pairIndexDict.items():
                for orientation, flankDict in orientationDict.items():
                    values = flankDict[flankIndex]
                    count = fillup_blanks(Counter(values), args.numTopIndices)
                    thisRow = [sampleID, pairIndex, orientation]
                    for pair in count:
                        thisRow.extend(pair)
                    dfList.append(thisRow)
        
        # Convert list into pandas DataFrame
        columns = ["sample", "file_pair", "orientation"]
        for i in range(args.numTopIndices):
            columns.extend([ f"index{i}", f"count{i}" ])
        df = pd.DataFrame(dfList, columns=columns)
        
        # Store DataFrame as a sheet in the Excel document
        sheetName = f"flank{flankIndex+1}"
        df.to_excel(writer, sheet_name=sheetName, index=True)
        workbook = writer.book
        worksheet = writer.sheets[sheetName]
        
        # Set column width based on contents
        columnLength = max([len(str(x)) for x in df.index])
        worksheet.set_column(0, 0, columnLength+1)
        
        for i, header in enumerate(df.columns):
            columnLength = max(df[header].astype(str).map(len).max(), len(header))
            worksheet.set_column(i+1, i+1, columnLength+1)
    
    writer.close()
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

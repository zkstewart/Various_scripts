#! python3
# exome_benchmark.py
# Program to enable discovery of exon sequences
# from genome sequences on the basis of exome
# sequencing alignments

import sys, argparse, os, statistics
sys.path.append(os.path.dirname(os.path.dirname(__file__))) # 2 dirs up is where we find dependencies
from Function_packages import ZS_SeqIO
from exome_liftover import ssw_parasail

def validate_args(args):
    # Validate input data location
    if not os.path.isdir(args.testDir):
        print('I am unable to locate the test directory where FASTA files are (' + args.testDir + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isdir(args.expectationDir):
        print('I am unable to locate the expectation directory where FASTA files are (' + args.expectationDir + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Handle file output
    if os.path.isfile(args.outputFileName):
        print('Output file name already exists; we will not allow file overwriting to prevent data loss')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    elif os.path.isdir(args.outputFileName):
        print('You\'ve specified an existing directory as outputFileName; you need to specify a file name that does not exist')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()

def calculate_sequence_similarity(FastASeq_obj1, FastASeq_obj2):
    '''
    This function takes two FastASeq objects and calculates several
    metrics for determining their similarity to one another. It works
    by aligning the two sequences and looking for gap positions and
    positions where the two sequences are the same.
    
    Params:
        FastASeq_obj1 -- a ZS_SeqIO.FastASeq object.
        FastASeq_obj2 -- a ZS_SeqIO.FastASeq object.
        
    Return:
        gaps -- the number of gaps obtained through alignment.
        identical -- the number of positions shared in common.
        substitutions -- the number of non-gapped positions that are different.
        percentageIdentical -- an overall ratio of how identical the two sequences are.
                               It's a simple calculation of identical / len(total)
        percentageNongappy -- an overall ratio of how "similar" the two sequences are.
                              It's the same as the other percentage but allowing for 
                              substitutions.
    '''
    
    # Easy scenario -- two sequences are identical
    if FastASeq_obj1.seq == FastASeq_obj2.seq:
        return [0, len(FastASeq_obj1.seq), 0, 1.0, 1.0]
    # Other scenarios -- they differ
    else:
        # Run alignment
        queryAlign, targetAlign, startIndex, _ = ssw_parasail(FastASeq_obj1.seq, FastASeq_obj2.seq)
        
        # Check how many gaps are present
        gaps = startIndex + (len(FastASeq_obj1.seq) - len(targetAlign)) # Unaligned residues are considered gaps
        gaps += sum([1 for i in range(len(queryAlign)) if queryAlign[i] == "-" or targetAlign[i] == "-"])

        # Check how many identical positions are present
        identical = sum([1 for i in range(len(queryAlign)) if queryAlign[i] == targetAlign[i]])
        
        # Check how many substitutions are present
        substitutions = sum([1 for i in range(len(queryAlign)) if queryAlign[i] != targetAlign[i] and queryAlign[i] != "-" and targetAlign[i] != "-"])
        
        # Calculate overall identical percentage
        percentageIdentical = identical / max(len(FastASeq_obj1.seq), len(FastASeq_obj2.seq))
        
        # Calculate nongappy percentage
        percentageNongappy = (identical+substitutions) / max(len(FastASeq_obj1.seq), len(FastASeq_obj2.seq))
        
        return [gaps, identical, substitutions, percentageIdentical, percentageNongappy]
        
def main():
    usage = """%(prog)s receives two directories of FASTA files with the same names. The
    "expectation" directory contains our gold-standard sequences we want to match. The
    "test" directory contains our new sequences we want to check. The output is a tabular
    file indicating per-sequence details, as well as overall statistics of how well our
    test files match our expectation files.
    
    Note that files can be missing from "test", but not from "expectation".
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-t", dest="testDir", required=True,
                help="Specify the directory where aligned FASTA files are located")
    p.add_argument("-e", dest="expectationDir", required=True,
                help="Specify the directory where HMMER executables are located")
    p.add_argument("-o", dest="outputFileName", required=True,
                help="Output name for tabular file")
    
    args = p.parse_args()
    validate_args(args)
    
    # Locate all files
    testFilesSuffix = os.listdir(args.testDir)[0].rsplit(".", maxsplit=1)[1] # Handle .fa or .fasta
    testFilesBase = [file.rsplit(".", maxsplit=1)[0] for file in os.listdir(args.testDir)]
    
    expectationFilesSuffix = os.listdir(args.expectationDir)[0].rsplit(".", maxsplit=1)[1]
    expectationFilesBase = [file.rsplit(".", maxsplit=1)[0] for file in os.listdir(args.expectationDir)]
    
    # Find files that can be used for benchmarking
    benchmarkingFilesBase = [file for file in testFilesBase if file in expectationFilesBase]
    
    # Compare like-for-like sequences and save data for each comparison
    results = []
    for fileBase in benchmarkingFilesBase:
        # Derive file names for comparison
        testFileName = os.path.join(args.testDir, "{0}.{1}".format(fileBase, testFilesSuffix))
        expectationFileName = os.path.join(args.expectationDir,"{0}.{1}".format(fileBase, expectationFilesSuffix))
    
        # Load in files as FASTA objects
        testFASTA = ZS_SeqIO.FASTA(testFileName)
        expectationFASTA = ZS_SeqIO.FASTA(expectationFileName)
        
        # Determine sequence similarity
        r = gaps, identical, substitutions, percentageIdentical, percentageNongappy = calculate_sequence_similarity(testFASTA[0], expectationFASTA[0])

        # Store result
        results.append(r)
    
    # Compute some high-level statistics
    ## Gaps
    gapsMedian = statistics.median([r[0] for r in results])
    gapsMax = max([r[0] for r in results])
    print("Gaps median = {0}".format(gapsMedian))
    print("Gaps max = {0}".format(gapsMax))

    ## Percentages
    pIdenticalMedian = statistics.median([r[3] for r in results])
    pNongappyMedian = statistics.median([r[4] for r in results])
    print("Percentage identical median = {0}".format(pIdenticalMedian))
    print("Percentage nongappy median = {0}".format(pNongappyMedian))
    
    ## Locate problem sequences
    PROBLEMATIC = 0.70
    problems = [testFilesBase[i] for i in range(len(results)) if results[i][3] < PROBLEMATIC]
    
    ## Problems to follow up on
    # problems[3] 'ENSSHAP00000000194-mx' has a huge gap - what's the deal?
    ## ANS: It's also part of many other sequences, so it's probably good =D
    import pyperclip
    for p in problems:
        # Derive file names for comparison
        testFileName = os.path.join(args.testDir, "{0}.{1}".format(p, testFilesSuffix))
        expectationFileName = os.path.join(args.expectationDir,"{0}.{1}".format(p, expectationFilesSuffix))
    
        # Load in files as FASTA objects
        testFASTA = ZS_SeqIO.FASTA(testFileName)
        expectationFASTA = ZS_SeqIO.FASTA(expectationFileName)
        
        # Store in clipboard for analysis
        pyperclip.copy("\n".join([testFASTA[0].get_str(), expectationFASTA[0].get_str()]))
        
        print("Okay, next?")
        button = input()
        
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

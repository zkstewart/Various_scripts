#! python3
# copynum_prediction_pipeline.py
# Script to automatically run Control-FREEC analysis for the
# prediction of gene copy number variation in a set of BAM files
# using a GFF3 file and a reference genome FASTA file.

import os, argparse, sys
from Bio import SeqIO

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))) # 3 dirs up is where we find Function_packages
from Function_packages import ZS_Utility, ZS_BAMIO

# Define functions
def validate_args(args):
    def _not_specified_error(program):
        print(f"ERROR: {program} not discoverable in your system PATH and was not specified as an argument.")
        quit()
    def _not_found_error(program, path):
        print(f"ERROR: {program} was not found at the location indicated ('{path}')")
        quit()
    
    # Validate input file locations
    if not os.path.isfile(args.gff3File):
        print('I am unable to locate the GFF3 file (' + args.gff3File + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.fastaFile):
        print('I am unable to locate the genome FASTA file (' + args.fastaFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isdir(args.bamDirectory):
        print('I am unable to locate the BAM directory (' + args.bamDirectory + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate BAM suffix
    foundABAM = False
    for file in os.listdir(args.bamDirectory):
        if file.endswith(args.bamSuffix):
            foundABAM = True
            break
    if not foundABAM:
        print(f'No BAM files with suffix "{args.bamSuffix}" found in directory "{args.bamDirectory}"')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate program discoverability
    if args.bamParser is None:
        # Check for the presence of samtools or sambamba
        samtools = ZS_Utility.wsl_which("samtools")
        sambamba = ZS_Utility.wsl_which("sambamba")
        if samtools is None and sambamba is None:
            _not_specified_error("samtools or sambamba are")
        
        # Preferentially pick sambamba, or fallback to samtools
        elif sambamba is None:
            args.bamParser = samtools
        else:
            args.bamParser = sambamba
    else:
        assert os.path.basename(args.bamParser) in ["samtools", "sambamba"], \
            "The --bamParser argument must point to either 'samtools' or 'sambamba'"
        
        if not os.path.isfile(args.bamParser):
            _not_found_error(os.path.basename(args.bamParser), args.bamParser)
    
    if args.samtools is None:
        args.samtools = ZS_Utility.wsl_which("samtools")
        if args.samtools is None:
            _not_specified_error("samtools")
    else:
        if not os.path.isfile(args.samtools):
            _not_found_error("samtools", args.samtools)
    
    if args.picard is None:
        args.picard = ZS_Utility.wsl_which("picard")
        if args.picard is None:
            _not_specified_error("picard")
    else:
        if not os.path.isfile(args.picard):
            _not_found_error("picard", args.picard)
    
    if args.freec is None:
        args.freec = ZS_Utility.wsl_which("freec")
        if args.freec is None:
            _not_specified_error("freec")
    else:
        if not os.path.isfile(args.freec):
            _not_found_error("freec", args.freec)
    # Validate numeric arguments
    if args.cpus < 1:
        print("cpus must be a positive integer; fix this and try again.")
        quit()
    if len(args.ploidy) != len(set(args.ploidy)):
        print("ploidy values must not be duplicated; fix this and try again.")
        quit()
    for ploidyNum in args.ploidy:
        if ploidyNum < 1:
            print(f"ploidy values must be positive integers (unlike '{ploidyNum}'); fix this and try again.")
            quit()
    # Validate output file location
    if os.path.isdir(args.outputDirectory) and os.listdir(args.outputDirectory) != []:
        print(f"Output directory '{args.outputDirectory}' already exists; I'll write output files here.")
        print("But, I won't overwrite any existing files, so beware that if a previous run had issues, " +
              "you may need to delete/move files first.")
    if not os.path.isdir(args.outputDirectory):
        os.makedirs(args.outputDirectory)
        print(f"Output directory '{args.outputDirectory}' has been created as part of argument validation.")

def explode_and_rename_for_freec(fastaFile, outputRenameTSV, outputDirectory):
    '''
    Parameters:
        fastaFile -- a string indicating the location of a genome to make
                     Control-FREEC compliant
        outputRenameTSV -- a string indicating the location to write a TSV
                           file that maps the original sequence IDs to the
                           new sequence IDs.
        outputDirectory -- a string indicating the location to write the
                           exploded FASTA files to.
    '''
    
    with open(fastaFile, "r") as fastaIn, open(outputRenameTSV, "w") as renameOut:
        records = SeqIO.parse(fastaIn, 'fasta')
        ongoingCount = 1
        for record in records:
            # Write the exploded FASTA file
            newSeqID = f"chr{ongoingCount}"
            outputFileName = os.path.join(outputDirectory, newSeqID + ".fasta")
            if not os.path.exists(outputFileName):
                with open(outputFileName, 'w') as fastaOut:
                    fastaOut.write(f">{newSeqID}\n{str(record.seq)}\n")
            ongoingCount += 1
            
            # Add to the renaming TSV
            renameOut.write(f"{record.id}\t{newSeqID}\n")

def parse_rename_tsv(renameTSVFileName):
    '''
    Parameters:
        renameTSVFileName -- a TSV file with two columns, the first being the
                             original sequence ID and the second being the new
                             sequence ID; this file should not have a header row.
    Returns:
        renameDict -- a dictionary with the original sequence ID as the key and
                      the new sequence ID as the value.
    '''
    
    renameDict = {}
    with open(renameTSVFileName, "r") as renameIn:
        for line in renameIn:
            sl = line.strip().split("\t")
            renameDict[sl[0]] = sl[1]
    
    return renameDict

def freec_bam_reheader(renameTSVFileName, bamFile, outputFile, samtoolsPath, picardPath):
    '''
    Parameters:
        renameTSVFileName -- a TSV file with two columns, the first being the
                             original sequence ID and the second being the new
                             sequence ID; this file should not have a header row.
        bamFile -- a string indicating the location to write of the BAM file to
                   reheader.
        outputFile -- a string indicating the location to write the reheadered
                      BAM to.
        samtoolsPath -- a string indicating the location of the samtools executable.
        picardPath -- a string indicating the location of the picard executable.
    '''
    
    # Parse the renaming TSV
    renameDict = parse_rename_tsv(renameTSVFileName)
    
    # Obtain the current header for the BAM file
    tmpHeaderFile = ZS_Utility.tmp_file_name_gen("freec_bam_reheader_tmp", ".raw.sam")
    ZS_BAMIO.StandardProgramRunners.samtools_view_header(bamFile, tmpHeaderFile, samtoolsPath)
    
    # Modify the header to reflect the new sequence IDs
    tmpRenamedFile = ZS_Utility.tmp_file_name_gen("freec_bam_reheader_tmp", ".renamed.sam")
    with open(tmpHeaderFile, "r") as headerIn, open(tmpRenamedFile, "w") as headerOut:
        for line in headerIn:
            if line.startswith("@SQ"):
                sl = line.strip().split("\t")
                seqID = sl[1].split(":")[1]
                sl[1] = f"SN:{renameDict[seqID]}"
                headerOut.write("\t".join(sl) + "\n")
            else:
                headerOut.write(line)
    
    # Clean up first temporary file
    os.remove(tmpHeaderFile)
    
    # Reheader the BAM file
    ZS_BAMIO.StandardProgramRunners.picard_ReplaceSamHeader(bamFile, tmpRenamedFile, outputFile, picardPath)
    
    # Clean up second temporary file
    os.remove(tmpRenamedFile)

## Main
def main():
    # User input
    usage = """%(prog)s reads in GFF3, a genome FASTA, and a directory containing one or more
    BAM files. It will rename the FASTA in Control-FREEC compliant fashion, update the BAM files
    to have headers reflecting the change in sequence IDs, and then run Control-FREEC for each
    sample/BAM file.
    """
    p = argparse.ArgumentParser(description=usage)
    # Reqs
    p.add_argument("-g", dest="gff3File",
                   required=True,
                   help="Input GFF3 file")
    p.add_argument("-f", dest="fastaFile",
                   required=True,
                   help="Input genome FASTA file")
    p.add_argument("-b", dest="bamDirectory",
                   required=True,
                   help="Input directory containing BAM file(s)")
    p.add_argument("-o", dest="outputDirectory",
                   required=True,
                   help="Specify location to write output files to")
    # Opts (programs)
    p.add_argument("--samtools", dest="samtools",
                   required=False,
                   help="""Optionally, specify the samtools executable file
                   if it is not discoverable in the path""",
                   default=None)
    p.add_argument("--picard", dest="picard",
                   required=False,
                   help="""Optionally, specify the picard executable file
                   if it is not discoverable in the path""",
                   default=None)
    p.add_argument("--freec", dest="freec",
                   required=False,
                   help="""Optionally, specify the freec executable file
                   if it is not discoverable in the path""",
                   default=None)
    # Opts (Control-FREEC)
    p.add_argument("--freecBamParser", dest="bamParser",
                   required=False,
                   help="""Optionally, specify the sambamba executable OR
                   samtools executable file if neither are discoverable in the path;
                   if both are discoverable, sambamba will be preferred""",
                   default=None)
    p.add_argument("--cpus", dest="cpus",
                   required=False,
                   type=int,
                   help="""Optionally, specify the number of CPUs to use for
                   multi-threadable operations as in Control-FREEC; default == 1""",
                   default=1)
    p.add_argument("--ploidy", dest="ploidy",
                   required=False,
                   type=int,
                   nargs="+",
                   help="""Optionally, specify the one or more ploidy numbers to test for
                   in each sample; default == 2""",
                   default=[2])
    p.add_argument("--mateOrientation", dest="mateOrientation",
                   required=False,
                   choices=["0", "RF", "FR", "FF"],
                   help="""Specify the mate orientation for Control-FREEC;
                   0 for single end or sorted BAM, RF for mate pairs, FR for paired end,
                   FF for SOLiD; default == '0'""",
                   default="0")
    # Opts (behavioural)
    p.add_argument("--bamSuffix", dest="bamSuffix",
                   required=False,
                   help="""Indicate the suffix of the BAM files to look for;
                   default == '.sorted.bam'""",
                   default=".sorted.bam")
    p.add_argument("--relaxed", dest="relaxedParsing",
                   required=False,
                   action='store_true',
                   help="""Optionally specify whether we should use relaxed GFF3 parsing.""",
                   default=False)
    
    args = p.parse_args()
    validate_args(args)
    
    # Explode and rename FASTA for Control-FREEC
    explosionDir = os.path.join(args.outputDirectory, "explosion")
    renameTSVFileName = os.path.join(args.outputDirectory, "id_renaming_list.tsv")
    os.makedirs(explosionDir, exist_ok=True)
    
    if not os.path.exists(renameTSVFileName) or not \
        os.path.exists(os.path.join(args.outputDirectory, "explosion_was_successful.flag")):
            explode_and_rename_for_freec(args.fastaFile, renameTSVFileName, explosionDir)
            open(os.path.join(args.outputDirectory, "explosion_was_successful.flag"), "w").close()
    else:
        print(f"sequence explosion has already been performed; skipping.")
    
    # Reheader BAM files with new sequence IDs
    bamReheaderDir = os.path.join(args.outputDirectory, "bam_reheader")
    os.makedirs(explosionDir, exist_ok=True)
    if not os.path.exists(os.path.join(args.outputDirectory, "reheader_was_successful.flag")):
        # Locate BAM files
        bamFiles = [
            os.path.join(args.bamDirectory, file)
            for file in os.listdir(args.bamDirectory)
            if file.endswith(args.bamSuffix)
        ]
        
        # Reheader each BAM file
        for bamFile in bamFiles:
            outputBamFile = os.path.join(bamReheaderDir, os.path.basename(bamFile))
            if not os.path.exists(outputBamFile):
                freec_bam_reheader(renameTSVFileName, bamFile, outputBamFile, args.samtools, args.picard)
        open(os.path.join(args.outputDirectory, "reheader_was_successful.flag"), "w").close()
    else:
        print(f"bam reheadering has already been performing; skipping.")
    
    # Let user know everything went swimmingly
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

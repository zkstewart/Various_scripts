#! python3
# copynum_prediction_pipeline.py
# Script to automatically run Control-FREEC analysis for the
# prediction of gene copy number variation in a set of BAM files
# using a GFF3 file and a reference genome FASTA file.

import os, argparse, sys, platform, subprocess
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
    tmpHeaderFile = ZS_Utility.tmp_file_name_gen(
        os.path.join(os.getcwd(), "freec_bam_reheader_tmp"), "raw.sam")
    ZS_BAMIO.StandardProgramRunners.samtools_view_header(bamFile, tmpHeaderFile, samtoolsPath)
    
    # Modify the header to reflect the new sequence IDs
    tmpRenamedFile = ZS_Utility.tmp_file_name_gen(
        os.path.join(os.getcwd(), "freec_bam_reheader_tmp"), "renamed.sam")
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

def generate_chr_len_file(explosionDir, outputFileName):
    '''
    Parameters:
        explosionDir -- a string indicating the directory containing exploded contigs
        outputFileName -- a string indicating the location to write the chrLen file to
    '''
    # Parse the FASTA files for sequence lengths
    chrLenList = []
    for fastaFile in os.listdir(explosionDir):
        # Parse in a single record
        record = SeqIO.read(os.path.join(explosionDir, fastaFile), "fasta")
        
        # Extract components for the chrLen file
        seqID = fastaFile.rsplit(".", maxsplit=1)[0]
        seqNum = int(seqID.split("chr")[1])
        seqLen = len(record.seq)
        chrLenList.append([seqNum, seqID, seqLen])
    
    # Sort the list by sequence number
    chrLenList.sort(key=lambda x: x[0])
    
    # Write the chromosome length file
    with open(outputFileName, "w") as chrLenOut:
        for seqNum, seqID, seqLen in chrLenList:
            chrLenOut.write(f"{seqNum}\t{seqID}\t{seqLen}\n")

def generate_freec_conf_file(bamFile, explosionDir, chrLenFile, ploidyNums,
                             mateOrientation, outputDir, outputConfFile, bamSuffix,
                             cpus=1, controlBamFile=None, sambambaPath=None):
    '''
    Parameters:
        bamFile -- a string indicating the location of the BAM file to use
        explosionDir -- a string indicating the location of the exploded FASTA files
        chrLenFile -- a string indicating the location of the chromosome length file to use
        ploidyNums -- a list of integers indicating the ploidy numbers to test for
        mateOrientation -- a string indicating the mate orientation to use for Control-FREEC
        outputDir -- a string indicating the location to write the Control-FREEC output to
        outputConfFile -- a string indicating the location to write the Control-FREEC
                          configuration file to
        bamSuffix -- a string indicating the suffix of the BAM files
        cpus -- OPTIONAL; an integer indicating the number of CPUs to use for Control-FREEC
                computations; default == 1
        controlBamFile -- OPTIONAL; a string indicating the location of a control BAM file
                          to use for Control-FREEC; default == None
        sambambaPath -- OPTIONAL; a string indicating the location of the sambamba executable
                        to use for Control-FREEC in lieu of samtools; default == None
    '''
    # Upgrade file locations to absolute paths with WSL formatting
    bamFile = ZS_Utility.convert_to_wsl_if_not_unix(os.path.abspath(bamFile))
    explosionDir = ZS_Utility.convert_to_wsl_if_not_unix(os.path.abspath(explosionDir))
    chrLenFile = ZS_Utility.convert_to_wsl_if_not_unix(os.path.abspath(chrLenFile))
    outputDir = ZS_Utility.convert_to_wsl_if_not_unix(os.path.abspath(outputDir))
    
    # Check if the control BAM (if applicable) is the same as the sample BAM
    sampleIsControl = False
    if controlBamFile != None:
        controlBamPrefix = os.path.basename(controlBamFile).replace(bamSuffix, "")
        sampleBamPrefix = os.path.basename(bamFile).replace(bamSuffix, "")
        sampleIsControl = controlBamPrefix == sampleBamPrefix
    
    # Generate the Control-FREEC configuration file
    with open(outputConfFile, "w") as confOut:
        # Write the general settings
        confOut.write(f"[general]\n\n")
        confOut.write(f"chrLenFile={chrLenFile}\n")
        confOut.write("ploidy={0}\n".format(",".join(map(str, ploidyNums))))
        confOut.write(f"chrFiles={explosionDir}\n")
        confOut.write(f"outputDir={outputDir}\n\n") # linebreak for behavioural params
        confOut.write(f"breakPointThreshold=.8\n")
        confOut.write(f"window=50000\n\n") # linebreak for computational params
        confOut.write(f"maxThreads={cpus}\n")
        confOut.write(f"numberOfProcesses={cpus}\n")
        if sambambaPath is not None:
            confOut.write(f"sambamba={ZS_Utility.convert_to_wsl_if_not_unix(sambambaPath)}\n")
        
        # Write general settings contingent on control presence
        if (controlBamFile is not None) and (sampleIsControl is False):
            confOut.write(f"\nintercept=0\n")
            confOut.write(f"degree=1\n")
            confOut.write(f"forceGCcontentNormalization=0\n")
        
        # Write sample settings
        confOut.write(f"\n[sample]\n\n")
        confOut.write(f"mateFile={bamFile}\n")
        confOut.write(f"inputFormat=BAM\n")
        confOut.write(f"mateOrientation={mateOrientation}\n\n")
        
        # Write control settings (if relevant)
        if (controlBamFile is not None) and (sampleIsControl is False):
            controlBamFile = ZS_Utility.convert_to_wsl_if_not_unix(os.path.abspath(controlBamFile))
            confOut.write(f"[control]\n\n")
            confOut.write(f"mateFile={controlBamFile}\n")
            confOut.write(f"inputFormat=BAM\n")
            confOut.write(f"mateOrientation={mateOrientation}\n\n")

def freec_worked(outputDir, EXPECTED_SUFFIXES):
    return len(
        [ f for f in os.listdir(outputDir) for s in EXPECTED_SUFFIXES if f.endswith(s) ]
    ) == len(EXPECTED_SUFFIXES)

def run_freec(confFile, outputDir, freecPath,
              EXPECTED_SUFFIXES=["_CNVs", "_info.txt", "_ratio.txt", "_sample.cpn"]):
    '''
    Parameters:
        confFile -- a string indicating the location of the freec .conf file
        outputDir -- a string indicating the location where freec will write
                     outputs to; needed for validating if the program run successfully
        freecPath -- a string indicating the location of the freec executable
        EXPECTED_SUFFIXES -- OPTIONAL; a list of strings indicating the suffixes
                             of files that freec should output; used to check for
                             successful program completion
    '''
    # Construct the cmd for subprocess
    cmd = ZS_Utility.base_subprocess_cmd(freecPath)
    cmd += [
        "-conf", ZS_Utility.convert_to_wsl_if_not_unix(confFile)
    ]
    
    if platform.system() != "Windows":
        cmd = " ".join(cmd)
    
    # Run the command
    run_freec = subprocess.Popen(cmd, shell = True,
                                 stdout = subprocess.PIPE,
                                 stderr = subprocess.PIPE)
    freecout, freecerr = run_freec.communicate()
    
    # Check file outputs to see if there was an error
    if not freec_worked(outputDir, EXPECTED_SUFFIXES):
        raise Exception(("ERROR: run_freec encountered an error; have a look " +
                        f'at the stdout ({freecout.decode("utf-8")}) and stderr ' + 
                        f'({freecerr.decode("utf-8")}) to make sense of this.'))

def parse_freec_info_file(infoFileName):
    '''
    Parameters:
        infoFileName -- a string indicating the location of the Control-FREEC
                        info file to parse.
    Returns:
        infoDict -- a dictionary with the parsed information from the info file
                    where keys correspond to column 1, and values to column 2.
    '''
    infoDict = {}
    with open(infoFileName, "r") as fileIn:
        for line in fileIn:
            sl = line.strip().split("\t")
            infoDict[sl[0]] = sl[1]
    return infoDict

def tabulate_ploidy_estimates(freecBaseDir, outputFileName):
    '''
    Parameters:
        freecBaseDir -- a string indicating the location where freec folders are
                        located
        outputFileName -- a string indicating the location to write the ploidy
                          estimates to
    '''
    with open(outputFileName, "w") as fileOut:
        # Write header line
        fileOut.write("sample_id\tploidy\n")
        
        # Iterate through each sample folder
        for freecOutputDir in os.listdir(freecBaseDir):
            if not freecOutputDir.endswith(".conf"):
                # Format the *_info.txt file path
                infoFileName = os.path.join(freecBaseDir, freecOutputDir, freecOutputDir + "_info.txt")
                
                # Parse the Control-FREEC info file
                infoDict = parse_freec_info_file(infoFileName)
                
                # Write to file
                fileOut.write(f"{freecOutputDir}\t{infoDict['Output_Ploidy']}\n")

def tabulate_copynum_estimates(freecBaseDir, outputFileName):
    '''
    Parameters:
        freecBaseDir -- a string indicating the location where freec folders are
                        located
        outputFileName -- a string indicating the location to write the copy number
                          estimates to
    '''
    raise NotImplementedError()
    
    with open(outputFileName, "w") as fileOut:
        # Write header line
        fileOut.write("sample_id\tploidy\n")
        
        # Iterate through each sample folder
        for freecOutputDir in os.listdir(freecBaseDir):
            if not freecOutputDir.endswith(".conf"):
                # Format the *_info.txt file path
                infoFileName = os.path.join(freecBaseDir, freecOutputDir, freecOutputDir + "_info.txt")
                
                # Parse the Control-FREEC info file
                infoDict = parse_freec_info_file(infoFileName)
                
                # Write to file
                fileOut.write(f"{freecOutputDir}\t{infoDict['Output_Ploidy']}\n")

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
    p.add_argument("--controlSample", dest="controlSample",
                   required=False,
                   help="""Optionally, if you have a control BAM file e.g., reads
                   from the reference genome species/variety, specify it here so as
                   to configure Control-FREEC to use it as control; specify the file
                   prefix prior to the --bamSuffix value!""",
                   default=None)
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
    os.makedirs(bamReheaderDir, exist_ok=True)
    
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
    
    # Generate chrLenFile for use with Control-FREEC
    chrLenFile = os.path.join(args.outputDirectory, "chrLenFile.tsv")
    if not os.path.exists(os.path.join(args.outputDirectory, "chrLen_was_successful.flag")) \
        or not os.path.exists(chrLenFile):
            generate_chr_len_file(explosionDir, chrLenFile)
            open(os.path.join(args.outputDirectory, "chrLen_was_successful.flag"), "w").close()
    else:
        print(f"chrLenFile has already been generated; skipping.")
    
    # Run Control-FREEC for each BAM file
    EXPECTED_SUFFIXES = ["_CNVs", "_info.txt", "_ratio.txt", "_sample.cpn"] # used for program resumption
    freecBaseDir = os.path.join(args.outputDirectory, "freec_output")
    os.makedirs(freecBaseDir, exist_ok=True)
    
    if not os.path.exists(os.path.join(args.outputDirectory, "freec_was_successful.flag")):
        # Locate reheadered BAM files
        bamFiles = [ os.path.join(bamReheaderDir, file) for file in os.listdir(bamReheaderDir) ]
        
        # Ensure the control BAM exists (if applicable)
        if args.controlSample != None:
            controlBamFile = os.path.join(bamReheaderDir, args.controlSample + args.bamSuffix)
            if not os.path.exists(controlBamFile):
                print(f"Control BAM file '{controlBamFile}' not found; program must end now.")
                quit()
        else:
            controlBamFile = None
        
        # Iterate through each BAM file
        for bamFile in bamFiles:
            bamBase = os.path.basename(bamFile).replace(args.bamSuffix, "")
            workingDir = os.path.join(freecBaseDir, bamBase)
            
            # Generate .conf file
            confFileName = os.path.join(freecBaseDir, bamBase + ".conf")
            generate_freec_conf_file(bamFile, explosionDir, chrLenFile, args.ploidy,
                                     args.mateOrientation, workingDir, confFileName, args.bamSuffix,
                                     cpus=args.cpus, controlBamFile=controlBamFile,
                                     sambambaPath=None if "samtools" in args.bamParser else args.bamParser)
            
            # Run Control-FREEC
            if (not os.path.exists(workingDir)) or (not freec_worked(workingDir, EXPECTED_SUFFIXES)):
                os.makedirs(workingDir, exist_ok=True)
                run_freec(confFileName, workingDir, args.freec, EXPECTED_SUFFIXES)
            else:
                print(f"Control-FREEC has already been run for '{bamBase}'; skipping.")
        
        open(os.path.join(args.outputDirectory, "freec_was_successful.flag"), "w").close()
    else:
        print(f"Control-FREEC has already been performed; skipping.")
    
    # Tabulate sample ploidy results
    ploidyTableFile = os.path.join(args.outputDirectory, "ploidy_estimates.tsv")
    if not os.path.exists(ploidyTableFile) or not \
        os.path.exists(os.path.join(args.outputDirectory, "ploidy_tabulation_was_successful.flag")):            
            tabulate_ploidy_estimates(freecBaseDir, ploidyTableFile)
            open(os.path.join(args.outputDirectory, "ploidy_tabulation_was_successful.flag"), "w").close()
    else:
        print(f"ploidy tabulation file has already been generated; skipping.")
    
    # Tabulate gene copy number results
    copynumTableFile = os.path.join(args.outputDirectory, "gene_copy_numbers.tsv")
    if not os.path.exists(copynumTableFile) or not \
        os.path.exists(os.path.join(args.outputDirectory, "copynum_tabulation_was_successful.flag")):            
            tabulate_copynum_estimates(freecBaseDir, copynumTableFile) # not implemented yet
            open(os.path.join(args.outputDirectory, "copynum_tabulation_was_successful.flag"), "w").close()
    else:
        print(f"copy number table file has already been generated; skipping.")
    
    # Let user know everything went swimmingly
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

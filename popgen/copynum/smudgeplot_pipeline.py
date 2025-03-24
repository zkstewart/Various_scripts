#! python3
# smudgeplot_pipeline.py
# Script to automatically run smudgeplot analysis for the
# prediction of genome ploidy values in a set of FASTA/Q files.

import os, argparse, sys, platform, subprocess

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))) # 3 dirs up is where we find Function_packages
from Function_packages import ZS_Utility, ZS_SeqIO

# Define functions
def validate_args(args):
    def _not_specified_error(program):
        raise FileNotFoundError(f"{program} not discoverable in your system PATH and was not specified as an argument.")
    def _not_found_error(program, path):
        raise FileNotFoundError(f"{program} was not found at the location indicated ('{path}')")
    
    # Validate input file locations
    if not os.path.isdir(args.fileDirectory):
        raise FileNotFoundError('I am unable to locate the FASTA/Q directory (' + args.fileDirectory + ')')
    
    # Impute file suffix if not specified
    if args.fileSuffix is None:
        if args.fileType == "fasta" or args.fileType == "fasta_m":
            args.fileSuffix = ".fasta"
        elif args.fileType == "fastq":
            args.fileSuffix = ".fq.gz"
        else:
            raise NotImplementedError(f"File type '{args.fileType}' not recognized; please specify a file suffix.")
    
    # Validate file suffix
    foundAFile = False
    for file in os.listdir(args.fileDirectory):
        if file.endswith(args.fileSuffix):
            foundAFile = True
            break
    if not foundAFile:
        raise FileNotFoundError(f'No FASTA/Q files with suffix "{args.fileSuffix}" found in directory "{args.fileDirectory}"')
    
    # Validate program discoverability
    if args.smudgeplot is None:
        args.smudgeplot = ZS_Utility.wsl_which("smudgeplot.py")
        if args.smudgeplot is None:
            _not_specified_error("smudgeplot.py")
    else:
        if not os.path.isfile(args.smudgeplot):
            _not_found_error("smudgeplot.py", args.smudgeplot)
    
    if args.fastk is None:
        args.fastk = ZS_Utility.wsl_which("FastK")
        if args.fastk is None:
            _not_specified_error("FastK")
    else:
        if not os.path.isfile(args.fastk):
            _not_found_error("FastK", args.fastk)
    
    if args.kmc is None:
        args.kmc = ZS_Utility.wsl_which("kmc")
        if args.kmc is None:
            _not_specified_error("kmc")
    else:
        if not os.path.isfile(args.kmc):
            _not_found_error("kmc", args.kmc)
    
    if args.kmc_tools is None:
        args.kmc_tools = ZS_Utility.wsl_which("kmc_tools")
        if args.kmc_tools is None:
            _not_specified_error("kmc_tools")
    else:
        if not os.path.isfile(args.kmc_tools):
            _not_found_error("kmc_tools", args.kmc_tools)
    
    # Validate output file location
    if os.path.isdir(args.outputDirectory) and os.listdir(args.outputDirectory) != []:
        print(f"Output directory '{args.outputDirectory}' already exists; I'll write output files here.")
        print("But, I won't overwrite any existing files, so beware that if a previous run had issues, " +
              "you may need to delete/move files first.")
    if not os.path.isdir(args.outputDirectory):
        os.makedirs(args.outputDirectory)
        print(f"Output directory '{args.outputDirectory}' has been created as part of argument validation.")

# Functions for newer FastK-based smudgeplot
def run_fastk(reads, fastkTableFile, cpus, mem, fastkPath):
    '''
    Parameters:
        reads -- a string indicating the location of the FASTA/Q file(s) to process; paired
                 reads should be formatted as a string with square brackets giving the
                 alternative characters for the two files e.g., "sample_R[12].fastq.gz"
        fastkTableFile -- a string indicating the location to write the FastK table file to
        cpus -- an integer indicating how many threads to run FastK with
        mem -- an integer indicating how many gigabytes of memory to run FastK with
        fastkPath -- a string indicating the location of the FastK executable
    '''
    # Construct the cmd for subprocess
    cmd = ZS_Utility.base_subprocess_cmd(fastkPath)
    cmd += [
        f"-t4", "-k31", f"-M{mem}", f"-T{cpus}",
        ZS_Utility.convert_to_wsl_if_not_unix(reads),
        f"-N{ZS_Utility.convert_to_wsl_if_not_unix(fastkTableFile)}"
    ]
    
    if platform.system() != "Windows":
        cmd = " ".join(cmd)
    
    # Run the command
    run_fastk = subprocess.Popen(cmd, shell = True,
                                 stdout = subprocess.PIPE,
                                 stderr = subprocess.PIPE)
    fastkout, fastkerr = run_fastk.communicate()
    
    # Check to see if there was an error
    if not "Total Resources:" in fastkerr.decode("utf-8"):
        raise Exception(("ERROR: run_fastk encountered an error; have a look " +
                        f'at the stdout ({fastkout.decode("utf-8")}) and stderr ' + 
                        f'({fastkerr.decode("utf-8")}) to make sense of this.'))

def run_smudgeplot_hetmers(fastkTableFile, kmerPairsFile, L, cpus, smudgeplotPath):
    '''
    Not to be confused with run_smudgeplot_hetkmers from the older KMC-based smudgeplot.
    
    Parameters:
        fastkTableFile -- a string indicating the location of the FastK table file
        kmerPairsFile -- a string indicating the output file name for the hetmers result
        L -- the 
        cpus -- an integer indicating how many threads to run smudgeplot with
        smudgeplotPath -- a string indicating the location of the smudgeplot.py executable
    '''
    # Construct the cmd for subprocess
    cmd = ZS_Utility.base_subprocess_cmd(smudgeplotPath)
    cmd += [
        "hetmers", "-L", str(L), "-t", str(cpus),
        "-o", ZS_Utility.convert_to_wsl_if_not_unix(kmerPairsFile),
        ZS_Utility.convert_to_wsl_if_not_unix(fastkTableFile)
    ]
    
    if platform.system() != "Windows":
        cmd = " ".join(cmd)
    
    # Run the command
    run_smudge_hetmers = subprocess.Popen(cmd, shell = True,
                                          stdout = subprocess.PIPE,
                                          stderr = subprocess.PIPE)
    smudgeout, smudgeerr = run_smudge_hekmers.communicate()
    
    # Check to see if there was an error
    if (smudgeout.decode("utf-8") != "") or (not "Done!" in smudgeerr.decode("utf-8")):
        raise Exception(("ERROR: run_smudgeplot_hetmers encountered an error; have a look " +
                        f'at the stdout ({smudgeout.decode("utf-8")}) and stderr ' + 
                        f'({smudgeerr.decode("utf-8")}) to make sense of this.'))

def run_smudgeplot_all(kmerPairsFile, outputPrefix, smudgeplotPath):
    '''    
    Parameters:
        kmerPairsFile -- a string indicating the output file name for the hetmers result
        outputPrefix -- a string indicating the prefix for output file names
        smudgeplotPath -- a string indicating the location of the smudgeplot.py executable
    '''
    # Construct the cmd for subprocess
    cmd = ZS_Utility.base_subprocess_cmd(smudgeplotPath)
    cmd += [
        "all", "-o", ZS_Utility.convert_to_wsl_if_not_unix(outputPrefix),
        f"{ZS_Utility.convert_to_wsl_if_not_unix(kmerPairsFile)}_text.smu"
    ]
    
    if platform.system() != "Windows":
        cmd = " ".join(cmd)
    
    # Run the command
    run_smudge_hetmers = subprocess.Popen(cmd, shell = True,
                                          stdout = subprocess.PIPE,
                                          stderr = subprocess.PIPE)
    smudgeout, smudgeerr = run_smudge_hekmers.communicate()
    
    # Check to see if there was an error
    if (smudgeout.decode("utf-8") != "") or (not "Done!" in smudgeerr.decode("utf-8")):
        raise Exception(("ERROR: run_smudgeplot_hetmers encountered an error; have a look " +
                        f'at the stdout ({smudgeout.decode("utf-8")}) and stderr ' + 
                        f'({smudgeerr.decode("utf-8")}) to make sense of this.'))

def fastk_pipeline(smudgeplotDir, kmerDir, fastaqsDir, pair, samplePrefix, args):
    # Format pairs as string with square bracketed alternatives
    if len(pair) == 2:
        commonPrefix = os.path.commonprefix(pair)
        suffixes = [ p.split(commonPrefix)[1] for p in pair ]
        commonSuffix = os.path.commonprefix([ s[::-1] for s in suffixes ])[::-1]
        alternatives = [ suffixes[i].rsplit(commonSuffix)[0] for i in range(len(suffixes)) ]
        
        reads = os.path.join(
            fastaqsDir,  
            f"{commonPrefix}[{alternatives[0]}{alternatives[1]}]{commonSuffix}"
        )
    else:
        reads = os.path.join(fastaqsDir, pair[0])
    
    # Run FastK
    fastkTableFile = os.path.join(kmerDir, samplePrefix + "_table")
    if not os.path.exists(fastkTableFile):
        run_fastk(reads, fastkTableFile, args.cpus, args.mem, args.fastk)
    else:
        print(f"FastK has already been run for '{samplePrefix}'; skipping.")
    
    # Run smudgeplot hetmers
    kmerPairsFile = os.path.join(smudgeplotDir, samplePrefix + "_kmerpairs")
    if not os.path.exists(kmerPairsFile + "_text.smu"):
        run_smudgeplot_hetmers(fastkTableFile, kmerPairsFile, args.smudgeErroneousLower,
                               args.cpus, args.smudgeplot)
    else:
        print(f"smudgeplot hetmers has already been run for '{samplePrefix}'; skipping.")
    
    # Run smudgeplot all
    resultsFile = os.path.join(smudgeplotDir, samplePrefix + "_fastkout")
    if not os.path.exists(resultsFile):
        run_smudgeplot_all(kmerPairsFile, resultsFile, args.smudgeplot)
    else:
        print(f"smudgeplot all has already been run for '{samplePrefix}'; skipping.")

# Functions for older KMC-based smudgeplot
def generate_atfile_file(pair, fastaqDir, atFilesName):
    '''
    Parameters:
        pair -- a list containing one or two strings indicating the location of the FASTA/Q file(s)
                making up a pair
        fastaqDir -- a string indicating the location of the FASTA/Q files in the pair
        atFilesName -- a string indicating the location to write the @FILES file to
    '''
    with open(atFilesName, "w") as fileOut:
        for fastaqFile in pair:
            fileOut.write(os.path.join(fastaqDir, fastaqFile) + "\n")

def run_kmc(atFilesName, fileType, kmcdbPrefix, cpus, mem, tmpDir, kmcPath):
    '''
    Parameters:
        atFilesName -- a string indicating the location of the @FILES file containing
                       input file names
        fileType -- a string indicating the type of input file(s) (fasta, fasta_m, or fastq)
        kmcdbPrefix -- a string indicating the basename for kmc to write the k-mer
                       database to (note: kmc will add .kmc_pre and .kmc_suf extension)
        cpus -- an integer indicating how many threads to run KMC with
        mem -- an integer indicating how many gigabytes of memory to run KMC with
        tmpDir -- a string indicating the location to write KMC temporary files to
        kmcPath -- a string indicating the location of the kmc executable
    '''
    # Construct the cmd for subprocess
    cmd = ZS_Utility.base_subprocess_cmd(kmcPath)
    cmd += [
        "-fa" if fileType == "fasta" else "-fm" if fileType == "fasta_m" else "-fq",
        "-k21", f"-t{cpus}", f"-m{mem}", "-ci1", "-cs10000",
        "@" + ZS_Utility.convert_to_wsl_if_not_unix(atFilesName),
        ZS_Utility.convert_to_wsl_if_not_unix(kmcdbPrefix),
        ZS_Utility.convert_to_wsl_if_not_unix(tmpDir)
    ]
    
    if platform.system() != "Windows":
        cmd = " ".join(cmd)
    
    # Run the command
    run_kmc = subprocess.Popen(cmd, shell = True,
                               stdout = subprocess.PIPE,
                               stderr = subprocess.PIPE)
    kmcout, kmcerr = run_kmc.communicate()
    
    # Check to see if there was an error
    if (not "Stats:" in kmcout.decode("utf-8")) and (not "100%" in kmcerr.decode("utf-8")):
        raise Exception(("ERROR: run_kmc encountered an error; have a look " +
                        f'at the stdout ({kmcout.decode("utf-8")}) and stderr ' + 
                        f'({kmcerr.decode("utf-8")}) to make sense of this.'))

def run_kmctools_histogram(kmcdbPrefix, outFileName, kmc_toolsPath):
    '''
    Parameters:
        bamFile -- a string indicating the location of the BAM file to process
        outFileName -- a string indicating the output file name for the histogram
        kmc_toolsPath -- a string indicating the location of the kmc_tools executable
    '''
    # Construct the cmd for subprocess
    cmd = ZS_Utility.base_subprocess_cmd(kmc_toolsPath)
    cmd += [
        "transform", ZS_Utility.convert_to_wsl_if_not_unix(kmcdbPrefix),
        "histogram", ZS_Utility.convert_to_wsl_if_not_unix(outFileName),
        "-cx10000"
    ]
    
    if platform.system() != "Windows":
        cmd = " ".join(cmd)
    
    # Run the command
    run_kmc = subprocess.Popen(cmd, shell = True,
                               stdout = subprocess.PIPE,
                               stderr = subprocess.PIPE)
    kmcout, kmcerr = run_kmc.communicate()
    
    # Check to see if there was an error
    if (not kmcout.decode("utf-8") == "") and (not "100%" in kmcerr.decode("utf-8")):
        raise Exception(("ERROR: run_kmctools_histogram encountered an error; have a look " +
                        f'at the stdout ({kmcout.decode("utf-8")}) and stderr ' + 
                        f'({kmcerr.decode("utf-8")}) to make sense of this.'))

def run_smudgeplot_cutoff(kmcHistogramFile, cutoff, smudgeplotPath):
    '''
    Parameters:
        kmcHistogramFile -- a string indicating the location of the KMC histogram file
        cutoff -- a string of "L" or "U" for lower or upper cutoffs, respectively
        smudgeplotPath -- a string indicating the location of the smudgeplot.py executable
    Returns:
        stdout -- a string indicating the stdout from the smudgeplot cutoff command
    '''
    # Construct the cmd for subprocess
    cmd = ZS_Utility.base_subprocess_cmd(smudgeplotPath)
    cmd += [
        "cutoff", ZS_Utility.convert_to_wsl_if_not_unix(kmcHistogramFile),
        cutoff
    ]
    
    if platform.system() != "Windows":
        cmd = " ".join(cmd)
    
    # Run the command
    run_smudge_cutoff = subprocess.Popen(cmd, shell = True,
                                         stdout = subprocess.PIPE,
                                         stderr = subprocess.PIPE)
    smudgeout, smudgerr = run_smudge_cutoff.communicate()
    
    # Check to see if there was an error
    if (smudgeout.decode("utf-8") == "") or (not "Done!" in smudgerr.decode("utf-8")):
        raise Exception(("ERROR: run_smudgeplot_cutoff encountered an error; have a look " +
                        f'at the stdout ({smudgeout.decode("utf-8")}) and stderr ' + 
                        f'({smudgerr.decode("utf-8")}) to make sense of this.'))
    
    # Return cutoff value as printed to stdout
    return smudgeout.decode("utf-8").rstrip("\r\n ")

def run_kmctools_dump(kmcdbPrefix, lCutoff, uCutoff, outFileName, kmc_toolsPath):
    '''
    Parameters:
        bamFile -- a string indicating the location of the BAM file to process
        lCutoff -- an integer/string indicating the lower cutoff value
        uCutoff -- an integer/string indicating the upper cutoff value
        outFileName -- a string indicating the output file name for the histogram
        kmc_toolsPath -- a string indicating the location of the kmc_tools executable
    '''
    # Construct the cmd for subprocess
    cmd = ZS_Utility.base_subprocess_cmd(kmc_toolsPath)
    cmd += [
        "transform", ZS_Utility.convert_to_wsl_if_not_unix(kmcdbPrefix),
        f"-ci{lCutoff}", f"-cx{uCutoff}",
        "dump", "-s", ZS_Utility.convert_to_wsl_if_not_unix(outFileName)
    ]
    
    if platform.system() != "Windows":
        cmd = " ".join(cmd)
    
    # Run the command
    run_kmc = subprocess.Popen(cmd, shell = True,
                               stdout = subprocess.PIPE,
                               stderr = subprocess.PIPE)
    kmcout, kmcerr = run_kmc.communicate()
    
    # Check to see if there was an error
    if (not kmcout.decode("utf-8") == "") and (not "100%" in kmcerr.decode("utf-8")):
        raise Exception(("ERROR: run_kmctools_dump encountered an error; have a look " +
                        f'at the stdout ({kmcout.decode("utf-8")}) and stderr ' + 
                        f'({kmcerr.decode("utf-8")}) to make sense of this.'))

def run_smudgeplot_hetkmers(kmcDumpFile, outputFileName, smudgeplotPath):
    '''
    Parameters:
        kmcDumpFile -- a string indicating the location of the KMC dump file
        outputFileName -- a string indicating the output file name for the hetkmers result
        smudgeplotPath -- a string indicating the location of the smudgeplot.py executable
    '''
    # Construct the cmd for subprocess
    cmd = ZS_Utility.base_subprocess_cmd(smudgeplotPath)
    cmd += [
        "hetkmers", "-o", ZS_Utility.convert_to_wsl_if_not_unix(outputFileName),
        ZS_Utility.convert_to_wsl_if_not_unix(kmcDumpFile)
    ]
    
    if platform.system() != "Windows":
        cmd = " ".join(cmd)
    
    # Run the command
    run_smudge_hetkmers = subprocess.Popen(cmd, shell = True,
                                         stdout = subprocess.PIPE,
                                         stderr = subprocess.PIPE)
    smudgeout, smudgeerr = run_smudge_hetkmers.communicate()
    
    # Check to see if there was an error
    if (smudgeout.decode("utf-8") != "") or (not "Done!" in smudgeerr.decode("utf-8")):
        raise Exception(("ERROR: run_smudgeplot_hetkmers encountered an error; have a look " +
                        f'at the stdout ({smudgeout.decode("utf-8")}) and stderr ' + 
                        f'({smudgeerr.decode("utf-8")}) to make sense of this.'))

def run_smudgeplot_plot(smudgeCoveragesFile, outputPrefix, samplePrefix, smudgeplotPath):
    '''
    Parameters:
        smudgeCoveragesFile -- a string indicating the location of the smudgeplot
                               _coverages.tsv file
        outputPrefix -- a string indicating the prefix for output file names
        samplePrefix -- a string to use for the plot title; typically the sample name
        smudgeplotPath -- a string indicating the location of the smudgeplot.py executable
    '''
    # Construct the cmd for subprocess
    cmd = ZS_Utility.base_subprocess_cmd(smudgeplotPath)
    cmd += [
        "plot", "-t", samplePrefix,
        "-o", ZS_Utility.convert_to_wsl_if_not_unix(outputPrefix),
        ZS_Utility.convert_to_wsl_if_not_unix(smudgeCoveragesFile)
        
    ]
    
    if platform.system() != "Windows":
        cmd = " ".join(cmd)
    
    # Run the command
    run_smudge_plot = subprocess.Popen(cmd, shell = True,
                                         stdout = subprocess.PIPE,
                                         stderr = subprocess.PIPE)
    smudgeout, smudgeerr = run_smudge_plot.communicate()
    
    # Check to see if there was an error
    if not "Done!" in smudgeerr.decode("utf-8"):
        raise Exception(("ERROR: run_smudge_plot encountered an error; have a look " +
                        f'at the stdout ({smudgeout.decode("utf-8")}) and stderr ' + 
                        f'({smudgeerr.decode("utf-8")}) to make sense of this.'))

def parse_smudge_ploidy(summaryFileName):
    '''
    Parameters:
        summaryFileName -- a string indicating the location of the smudgepot verbose summary file
    Returns:
        ploidyNum -- string value of the most likely ploidy number according to smudgeplot
    '''
    with open(summaryFileName, "r") as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n ").split("\t")
            if sl[0] == "* Proposed ploidy:":
                assert len(sl) == 2, \
                    f"ERROR: smudgeplot summary file '{summaryFileName}' has an unexpected format."
                return sl[1]

def kmc_pipeline(smudgeplotDir, kmerDir, fastaqsDir, pair, samplePrefix, kmerTmpDir, args):
    # Write @FILES file for kmc
    atFilesName = os.path.join(kmerDir, samplePrefix + ".atfile")
    if not os.path.exists(atFilesName):
        generate_atfile_file(pair, fastaqsDir, atFilesName)
    else:
        print(f"@FILE generation has already been run for '{samplePrefix}'; skipping.")
    
    # Run kmc
    EXPECTED_SUFFIXES = [".kmc_pre", ".kmc_suf"] # used for program resumption
    
    kmcdbPrefix = os.path.join(kmerDir, samplePrefix + "_db")
    if not all([ os.path.exists(kmcdbPrefix + suf) for suf in EXPECTED_SUFFIXES ]):
        run_kmc(atFilesName, args.fileType, kmcdbPrefix, args.cpus, args.mem, kmerTmpDir, args.kmc)
    else:
        print(f"kmc db has already been generated for '{samplePrefix}'; skipping.")
    
    # Run kmc_tools histogram
    histFileName = kmcdbPrefix + ".hist"
    if not os.path.exists(histFileName):
        run_kmctools_histogram(kmcdbPrefix, histFileName, args.kmc_tools)
    else:
        print(f"kmc_tools histogram has already been run for '{samplePrefix}'; skipping.")
    
    # Run kmc_tools dump
    dumpFileName = kmcdbPrefix + ".dump"
    if not os.path.exists(dumpFileName) or os.path.getsize(dumpFileName) == 0:
        # Get lower and upper cutoffs from smudgeplot
        lCutoff = run_smudgeplot_cutoff(histFileName, "L", args.smudgeplot)
        uCutoff = run_smudgeplot_cutoff(histFileName, "U", args.smudgeplot)
        
        # Now drop a dump (file)
        run_kmctools_dump(kmcdbPrefix, lCutoff, uCutoff, dumpFileName, args.kmc_tools)
    else:
        print(f"kmc_tools has already taken a dump at '{samplePrefix}'; skipping.")
    
    # Run smudgeplot hetkmers function
    hetkmersPrefix = os.path.join(smudgeplotDir, samplePrefix + "_hetkmers")
    hetkmersFileName = hetkmersPrefix + "_coverages.tsv"
    if not os.path.exists(hetkmersFileName):
        run_smudgeplot_hetkmers(dumpFileName, hetkmersPrefix, args.smudgeplot)
    else:
        print(f"smudgeplot hetkmers has already been run for '{samplePrefix}'; skipping.")
    
    # Run smudgeplot plot function
    plotPrefix = os.path.join(args.outputDirectory, samplePrefix + "_plot")
    plotFileName = plotPrefix + "_smudgeplot.png"
    if not os.path.exists(plotFileName):
        run_smudgeplot_plot(hetkmersFileName, plotPrefix, samplePrefix, args.smudgeplot)
    else:
        print(f"smudgeplot plot has already been run for '{samplePrefix}'; skipping.")

## Main
def main():
    # User input
    usage = """%(prog)s accepts a directory containing one or more FASTA/Q files in single or paired end.
    It will run smudgeplot for each sample and tabulate their most likely ploidy. Using -k kmc assumes
    an older version of smudgeplot whereas newer versions only support fastk.
    """
    p = argparse.ArgumentParser(description=usage)
    # Reqs
    p.add_argument("-i", dest="fileDirectory",
                   required=True,
                   help="Input directory containing FASTA/Q file(s)")
    p.add_argument("-f", dest="fileType",
                   required=True,
                   choices=["fasta", "fasta_m", "fastq"],
                   help="""Specify whether the input read files are FASTA (single line),
                   FASTA (multi-line) or FASTQ""")
    p.add_argument("-k", dest="kmerProgram",
                   required=True,
                   choices=["kmc", "fastk"],
                   help="""Specify whether to use kmc
                   or fastk for k-mer counting""")
    p.add_argument("-o", dest="outputDirectory",
                   required=True,
                   help="Specify location to write output files to")
    # Opts (programs)
    p.add_argument("--smudgeplot", dest="smudgeplot",
                   required=False,
                   help="""Optionally, specify the smudgeplot.py file
                   if it is not discoverable in the path""",
                   default=None)
    p.add_argument("--fastk", dest="fastk",
                   required=False,
                   help="""Optionally, specify the FastK file
                   if it is not discoverable in the path""",
                   default=None)
    p.add_argument("--kmc", dest="kmc",
                   required=False,
                   help="""Optionally, specify the kmc file
                   if it is not discoverable in the path""",
                   default=None)
    p.add_argument("--kmc_tools", dest="kmc_tools",
                   required=False,
                   help="""Optionally, specify the kmc_tools file
                   if it is not discoverable in the path""",
                   default=None)
    # Opts (computational)
    p.add_argument("--cpus", dest="cpus",
                   required=False,
                   type=int,
                   help="""Indicate how many threads to run KMC with;
                   default == 1""",
                   default=1)
    p.add_argument("--mem", dest="mem",
                   required=False,
                   type=int,
                   help="""Indicate how many gigabytes of memory to run KMC with;
                   default == 64""",
                   default=64)
    # Opts (behavioural)
    p.add_argument("--L", dest="smudgeErroneousLower",
                   required=False,
                   type=int,
                   help="""If using FastK, specify the 'count threshold below which
                   k-mers are considered erroneous'; this should try to accommodate
                   sequencing depth divided by genome size * expected ploidy; default == 6""",
                   default=6)
    p.add_argument("--fileSuffix", dest="fileSuffix",
                   required=False,
                   help="""Indicate the suffix of the FASTA/Q files to look for;
                   default == '.fasta' if -f fasta else '.fq.gz'""",
                   default=None)
    p.add_argument("--isSingleEnd", dest="isSingleEnd",
                   required=False,
                   action="store_true",
                   help="Specify if files are single-ended",
                   default=False)
    
    args = p.parse_args()
    validate_args(args)
    
    # Locate and validate FASTA/Q files
    forwardReads, reverseReads = ZS_SeqIO.FASTA.get_rnaseq_files(args.fileDirectory, args.fileSuffix, args.isSingleEnd)
    
    # Reformat and tie paired reads together
    forwardReads = [ os.path.basename(f) for f in forwardReads ]
    reverseReads = [ os.path.basename(r) for r in reverseReads ] if reverseReads is not None else None
    
    if reverseReads is not None:
        pairedReads = [ [f, r] for f, r in zip(forwardReads, reverseReads) ]
    else:
        pairedReads = [ [f] for f in forwardReads ]
    
    # Symlink FASTA/Q files into working directory
    fastaqsDir = os.path.abspath(os.path.join(args.outputDirectory, "fastqs" if args.fileType == "fastq" else "fastas"))
    os.makedirs(fastaqsDir, exist_ok=True)
    
    if not os.path.exists(os.path.join(args.outputDirectory, "fastaq_symlink_was_successful.flag")):
        # Symlink each FASTA/Q file pair
        for pair in pairedReads:
            for fastaqFile in pair:
                outputFastaqFile = os.path.join(fastaqsDir, fastaqFile)
                if not os.path.exists(outputFastaqFile):
                    os.symlink(os.path.join(args.fileDirectory, fastaqFile), outputFastaqFile)
        open(os.path.join(args.outputDirectory, "fastaq_symlink_was_successful.flag"), "w").close()
    else:
        print(f"fastaq symlinking has already been performing; skipping.")
    
    # Set up the working directory structure
    if args.kmerProgram == "kmc":
        kmerDir = os.path.abspath(os.path.join(args.outputDirectory, "kmc_files"))
        os.makedirs(kmerDir, exist_ok=True)
        
        kmerTmpDir = os.path.join(kmerDir, "tmp")
        os.makedirs(kmerTmpDir, exist_ok=True)
    elif args.kmerProgram == "fastk":
        kmerDir = os.path.abspath(os.path.join(args.outputDirectory, "fastk_files"))
        os.makedirs(kmerDir, exist_ok=True)
    
    smudgeplotDir = os.path.abspath(os.path.join(args.outputDirectory, "smudgeplot_files"))
    os.makedirs(smudgeplotDir, exist_ok=True)
    
    # Run kmer->smudgeplot for each FASTA/Q file pairing
    if not os.path.exists(os.path.join(args.outputDirectory, "pipeline_was_successful.flag")):
        # Iterate through each FASTA/Q file pair
        for pair in pairedReads:
            samplePrefix = pair[0].replace(args.fileSuffix, "")
            
            # Run kmc or fastk
            if args.kmerProgram == "kmc":
                kmc_pipeline(smudgeplotDir, kmerDir, fastaqsDir, pair, samplePrefix, kmerTmpDir, args)
            elif args.kmerProgram == "fastk":
                fastk_pipeline(smudgeplotDir, kmerDir, fastaqsDir, pair, samplePrefix, args)
            
            
        
        open(os.path.join(args.outputDirectory, "pipeline_was_successful.flag"), "w").close()
    else:
        print(f"kmer->smudgeplot pipeline has already been performed; skipping.")
    
    # Tabulate ploidy number results
    ploidyTableFile = os.path.join(args.outputDirectory, "ploidy_numbers.tsv")
    if not os.path.exists(os.path.join(args.outputDirectory, "tabulation_was_successful.flag")):
        # Begin generating result tabulation
        with open(ploidyTableFile, "w") as fileOut:
            # Write header line
            fileOut.write("sampleid\tploidy\n")
            
            # Iterate through samples
            for pair in pairedReads:
                samplePrefix = pair[0].replace(args.fileSuffix, "")
                
                # Parse the output summary file
                summaryFileName = os.path.join(args.outputDirectory, samplePrefix + "_plot_verbose_summary.txt")
                ploidyNum = parse_smudge_ploidy(summaryFileName)
                
                # Write result to file
                fileOut.write(f"{samplePrefix}\t{ploidyNum}\n")
        
        open(os.path.join(args.outputDirectory, "tabulation_was_successful.flag"), "w").close()
    else:
        print(f"tabulation has already been performed; skipping.")
    
    # Let user know everything went swimmingly
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

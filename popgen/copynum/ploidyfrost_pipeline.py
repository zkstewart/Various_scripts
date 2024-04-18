#! python3
# ploidyfrost_pipeline.py
# Script to automatically run ploidyfrost analysis for the
# prediction of genome ploidy values in a set of FASTQ files.

import os, argparse, sys, platform, subprocess

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))) # 3 dirs up is where we find Function_packages
from Function_packages import ZS_Utility, ZS_SeqIO

# Define functions
def validate_args(args):
    def _not_specified_error(program):
        print(f"ERROR: {program} not discoverable in your system PATH and was not specified as an argument.")
        quit()
    def _not_found_error(program, path):
        print(f"ERROR: {program} was not found at the location indicated ('{path}')")
        quit()
    
    # Validate input file locations
    if not os.path.isdir(args.fastqDirectory):
        print('I am unable to locate the FASTQ directory (' + args.fastqDirectory + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isdir(args.pfScriptDir):
        print('I am unable to locate the PloidyFrost scripts directory (' + args.pfScriptDir + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    
    # Impute and validate the location of the Filter.R and Drawfreq.R scripts
    args.pfFilter = os.path.join(args.pfScriptDir, "Filter.R")
    args.pfDrawfreq = os.path.join(args.pfScriptDir, "Drawfreq.R")
    if not os.path.isfile(args.pfFilter):
        print('I am unable to locate the PloidyFrost Filter.R script where expected (' + args.pfScriptDir + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.pfDrawfreq):
        print('I am unable to locate the PloidyFrost Drawfreq.R script where expected (' + args.pfScriptDir + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    
    # Validate FASTQ suffix
    foundAFASTQ = False
    for file in os.listdir(args.fastqDirectory):
        if file.endswith(args.fastqSuffix):
            foundAFASTQ = True
            break
    if not foundAFASTQ:
        print(f'No FASTQ files with suffix "{args.fastqSuffix}" found in directory "{args.fastqDirectory}"')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    
    # Validate program discoverability
    if args.ploidyfrost is None:
        args.ploidyfrost = ZS_Utility.wsl_which("PloidyFrost")
        if args.ploidyfrost is None:
            _not_specified_error("PloidyFrost")
    else:
        if not os.path.isfile(args.ploidyfrost):
            _not_found_error("PloidyFrost", args.ploidyfrost)
    
    if args.bifrost is None:
        args.bifrost = ZS_Utility.wsl_which("Bifrost")
        if args.bifrost is None:
            _not_specified_error("Bifrost")
    else:
        if not os.path.isfile(args.bifrost):
            _not_found_error("Bifrost", args.bifrost)
    
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
    
    if args.rscript is None:
        args.rscript = ZS_Utility.wsl_which("Rscript")
        if args.rscript is None:
            _not_specified_error("Rscript")
    else:
        if not os.path.isfile(args.rscript):
            _not_found_error("Rscript", args.rscript)
    
    # Validate output file location
    if os.path.isdir(args.outputDirectory) and os.listdir(args.outputDirectory) != []:
        print(f"Output directory '{args.outputDirectory}' already exists; I'll write output files here.")
        print("But, I won't overwrite any existing files, so beware that if a previous run had issues, " +
              "you may need to delete/move files first.")
    if not os.path.isdir(args.outputDirectory):
        os.makedirs(args.outputDirectory)
        print(f"Output directory '{args.outputDirectory}' has been created as part of argument validation.")

def generate_atfile_file(pair, fastqDir, atFilesName):
    '''
    Parameters:
        pair -- a list containing one or two strings indicating the location of the FASTQ file(s)
                making up a pair
        fastqDir -- a string indicating the location of the FASTQ files in the pair
        atFilesName -- a string indicating the location to write the @FILES file to
    '''
    with open(atFilesName, "w") as fileOut:
        for fastqFile in pair:
            fileOut.write(os.path.join(fastqDir, fastqFile) + "\n")

def run_kmc(atFilesName, kmcdbPrefix, cpus, mem, tmpDir, kmcPath):
    '''
    Parameters:
        bamFile -- a string indicating the location of the BAM file to process
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
        "-k25", f"-t{cpus}", f"-m{mem}", "-ci1", "-cs10000",
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
    
    # Check file outputs to see if there was an error
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
        "histogram", ZS_Utility.convert_to_wsl_if_not_unix(outFileName)
    ]
    
    if platform.system() != "Windows":
        cmd = " ".join(cmd)
    
    # Run the command
    run_kmc = subprocess.Popen(cmd, shell = True,
                               stdout = subprocess.PIPE,
                               stderr = subprocess.PIPE)
    kmcout, kmcerr = run_kmc.communicate()
    
    # Check file outputs to see if there was an error
    if (not "Stats:" in kmcout.decode("utf-8")) and (not "100%" in kmcerr.decode("utf-8")):
        raise Exception(("ERROR: run_kmc encountered an error; have a look " +
                        f'at the stdout ({kmcout.decode("utf-8")}) and stderr ' + 
                        f'({kmcerr.decode("utf-8")}) to make sense of this.'))

def run_ploidyfrost_cutoff(kmcHistogramFile, cutoff, ploidyfrostPath):
    '''
    Parameters:
        kmcHistogramFile -- a string indicating the location of the KMC histogram file
        cutoff -- a string of "L" or "U" for lower or upper cutoffs, respectively
        ploidyfrostPath -- a string indicating the location of the PloidyFrost executable
    Returns:
        stdout -- a string indicating the stdout from the ploidyfrost cutoff command
    '''
    # Construct the cmd for subprocess
    cmd = ZS_Utility.base_subprocess_cmd(ploidyfrostPath)
    cmd += [
        f"cutoff{cutoff}", ZS_Utility.convert_to_wsl_if_not_unix(kmcHistogramFile)
    ]
    
    if cutoff == "U":
        cmd.append("0.998")
    
    if platform.system() != "Windows":
        cmd = " ".join(cmd)
    
    # Run the command
    run_frost_cutoff = subprocess.Popen(cmd, shell = True,
                                        stdout = subprocess.PIPE,
                                        stderr = subprocess.PIPE)
    frostout, frosterr = run_frost_cutoff.communicate()
    
    # Check file outputs to see if there was an error
    if (frostout.decode("utf-8") == "") or (not frosterr.decode("utf-8") == ""):
        raise Exception(("ERROR: run_frost_cutoff encountered an error; have a look " +
                        f'at the stdout ({frostout.decode("utf-8")}) and stderr ' + 
                        f'({frosterr.decode("utf-8")}) to make sense of this.'))
    
    # Return cutoff value as printed to stdout
    return frostout.decode("utf-8").rstrip("\r\n ")

def run_kmctools_filter(kmcdbPrefix, atFilesName, lCutoff, cpus, outFileName, kmc_toolsPath):
    '''
    Parameters:
        kmcdbPrefix -- a string indicating the basename that KMC wrote the k-mer database to
        atFilesName -- a string indicating the location of the @FILES file containing
                       input file names
        lCutoff -- an integer/string indicating the lower cutoff value
        cpus -- an integer/string indicating the number of CPUs to run KMC with
        outFileName -- a string indicating the output file name for the histogram
        kmc_toolsPath -- a string indicating the location of the kmc_tools executable
    '''
    # Construct the cmd for subprocess
    cmd = ZS_Utility.base_subprocess_cmd(kmc_toolsPath)
    cmd += [
        f"-t{cpus}", "filter", "-hm", ZS_Utility.convert_to_wsl_if_not_unix(kmcdbPrefix),
        "@" + ZS_Utility.convert_to_wsl_if_not_unix(atFilesName),
        f"-ci{lCutoff}", ZS_Utility.convert_to_wsl_if_not_unix(outFileName)
    ]
    
    if platform.system() != "Windows":
        cmd = " ".join(cmd)
    
    # Run the command
    run_kmc = subprocess.Popen(cmd, shell = True,
                               stdout = subprocess.PIPE,
                               stderr = subprocess.PIPE)
    kmcout, kmcerr = run_kmc.communicate()
    
    # Check file outputs to see if there was an error
    if (not kmcout.decode("utf-8") == "") and (not "100%" in kmcerr.decode("utf-8")):
        raise Exception(("ERROR: run_kmctools_filter encountered an error; have a look " +
                        f'at the stdout ({kmcout.decode("utf-8")}) and stderr ' + 
                        f'({kmcerr.decode("utf-8")}) to make sense of this.'))

def run_bifrost_build(filteredFastqFile, cpus, outputFileName, bifrostPath):
    '''
    Parameters:
        kmcDumpFile -- a string indicating the location of the KMC dump file
        cpus -- an integer/string indicating the number of CPUs to run KMC with
        outputFileName -- a string indicating the output file name for the hetkmers result
        bifrostPath -- a string indicating the location of the Bifrost executable
    '''
    # Construct the cmd for subprocess
    cmd = ZS_Utility.base_subprocess_cmd(bifrostPath)
    cmd += [
        "build", "-i", "-d", "-k", "25", "-v",
        "-r", ZS_Utility.convert_to_wsl_if_not_unix(filteredFastqFile),
        "-o", ZS_Utility.convert_to_wsl_if_not_unix(outputFileName),
        "-t", str(cpus)
    ]
    
    if platform.system() != "Windows":
        cmd = " ".join(cmd)
    
    # Run the command
    run_bifrost_build = subprocess.Popen(cmd, shell = True,
                                         stdout = subprocess.PIPE,
                                         stderr = subprocess.PIPE)
    bifrostout, bifrosterr = run_bifrost_build.communicate()
    
    # Check file outputs to see if there was an error
    if (bifrostout.decode("utf-8") == "") or (not bifrosterr.decode("utf-8") == ""):
        raise Exception(("ERROR: run_bifrost_build encountered an error; have a look " +
                        f'at the stdout ({bifrostout.decode("utf-8")}) and stderr ' + 
                        f'({bifrosterr.decode("utf-8")}) to make sense of this.'))

def gunzip(gzFileName):
    '''
    Parameters:
        gzFileName -- a string indicating the location of a .gz file
    '''
    # Construct the cmd for subprocess
    cmd = ZS_Utility.base_subprocess_cmd("gunzip")
    cmd += [
        ZS_Utility.convert_to_wsl_if_not_unix(gzFileName)
    ]
    
    if platform.system() != "Windows":
        cmd = " ".join(cmd)
    
    # Run the command
    run_gunzip = subprocess.Popen(cmd, shell = True,
                                  stdout = subprocess.PIPE,
                                  stderr = subprocess.PIPE)
    gunzipout, gunziperr = run_gunzip.communicate()
    
    # Check file outputs to see if there was an error
    if (gunzipout.decode("utf-8") == "") or (not gunziperr.decode("utf-8") == ""):
        raise Exception(("ERROR: gunzip encountered an error; have a look " +
                        f'at the stdout ({gunzipout.decode("utf-8")}) and stderr ' + 
                        f'({gunziperr.decode("utf-8")}) to make sense of this.'))

def run_ploidyfrost_estimate(bifrostGfaFile, kmcdbPrefix, ploidyfrostFileName,
                             lCutoff, uCutoff, cpus, ploidyfrostPath):
    '''
    Parameters:
        bifrostGfaFile -- a string indicating the prefix of the Bifrost .gfa
        kmcdbPrefix -- a string indicating the location of the KMC db file
        ploidyfrostFileName -- a string indicating the output file name for the ploidyfrost result
        lCutoff -- an integer/string indicating the lower cutoff value
        uCutoff -- an integer/string indicating the upper cutoff value
        ploidyfrostPath -- a string indicating the location of the PloidyFrost executable
    '''
    # Construct the cmd for subprocess
    cmd = ZS_Utility.base_subprocess_cmd(ploidyfrostPath)
    cmd += [
        "-g", ZS_Utility.convert_to_wsl_if_not_unix(bifrostGfaFile),
        "-d", ZS_Utility.convert_to_wsl_if_not_unix(kmcdbPrefix),
        "-t", str(cpus), "-v", "-o", ZS_Utility.convert_to_wsl_if_not_unix(ploidyfrostFileName),
        "-l", str(lCutoff), "-u", str(uCutoff)
    ]
    
    if platform.system() != "Windows":
        cmd = " ".join(cmd)
    
    # Run the command
    run_frost_estimate = subprocess.Popen(cmd, shell = True,
                                          stdout = subprocess.PIPE,
                                          stderr = subprocess.PIPE)
    frostout, frosterr = run_frost_estimate.communicate()
    
    # Check file outputs to see if there was an error
    if (frostout.decode("utf-8") == "") or (not frosterr.decode("utf-8") == ""):
        raise Exception(("ERROR: run_ploidyfrost_estimate encountered an error; have a look " +
                        f'at the stdout ({frostout.decode("utf-8")}) and stderr ' + 
                        f'({frosterr.decode("utf-8")}) to make sense of this.'))

def run_ploidyfrost_filter(ploidyfrostFileName, outputFileName, pfFilterPath, rscriptPath):
    '''
    Parameters:
        ploidyfrostFileName -- a string indicating the prefix of the PloidyFrost
                               estimate files
        outputFileName -- a string indicating the output file name for filtered result
        pfFilterPath -- a string indicating the location of the ploidyfrost Filter.R file
        rscriptPath -- a string indicating the location of the Rscript executable
    '''
    # Construct the cmd for subprocess
    cmd = ZS_Utility.base_subprocess_cmd(rscriptPath)
    cmd += [
        ZS_Utility.convert_to_wsl_if_not_unix(pfFilterPath),
        "-i", ZS_Utility.convert_to_wsl_if_not_unix(ploidyfrostFileName),
        "-o", ZS_Utility.convert_to_wsl_if_not_unix(outputFileName),
        "-n", "6", "-s", "11", "-q", "0.05"
    ]
    
    if platform.system() != "Windows":
        cmd = " ".join(cmd)
    
    # Run the command
    run_frost_filter = subprocess.Popen(cmd, shell = True,
                                        stdout = subprocess.PIPE,
                                        stderr = subprocess.PIPE)
    frostout, frosterr = run_frost_filter.communicate()
    
    # Check file outputs to see if there was an error
    if (frostout.decode("utf-8") == "") or (not frosterr.decode("utf-8") == ""):
        raise Exception(("ERROR: run_ploidyfrost_filter encountered an error; have a look " +
                        f'at the stdout ({frostout.decode("utf-8")}) and stderr ' + 
                        f'({frosterr.decode("utf-8")}) to make sense of this.'))

def run_ploidyfrost_gmm(pfFilterFileName, outputFileName, ploidyfrostPath):
    '''
    Parameters:
        ploidyfrostFileName -- a string indicating the prefix of the PloidyFrost
                               estimate files
        outputFileName -- a string indicating the output file name for filtered result
        pfFilterPath -- a string indicating the location of the ploidyfrost Filter.R file
        ploidyfrostPath -- a string indicating the location of the PloidyFrost executable
    '''
    # Construct the cmd for subprocess
    cmd = ZS_Utility.base_subprocess_cmd(ploidyfrostPath)
    cmd += [
        "model", "-g", ZS_Utility.convert_to_wsl_if_not_unix(pfFilterFileName),
        "-l", "2", "-u", "10", "-q", "0.05",
        "-o", ZS_Utility.convert_to_wsl_if_not_unix(outputFileName)
    ]
    
    if platform.system() != "Windows":
        cmd = " ".join(cmd)
    
    # Run the command
    run_frost_model = subprocess.Popen(cmd, shell = True,
                                       stdout = subprocess.PIPE,
                                       stderr = subprocess.PIPE)
    frostout, frosterr = run_frost_model.communicate()
    
    # Check file outputs to see if there was an error
    if (frostout.decode("utf-8") == "") or (not frosterr.decode("utf-8") == ""):
        raise Exception(("ERROR: run_ploidyfrost_gmm encountered an error; have a look " +
                        f'at the stdout ({frostout.decode("utf-8")}) and stderr ' + 
                        f'({frosterr.decode("utf-8")}) to make sense of this.'))

def run_ploidyfrost_drawfreq(ploidyfrostFileName, outputFileName, samplePrefix, pfDrawPath, rscriptPath):
    '''
    Parameters:
        ploidyfrostFileName -- a string indicating the prefix of the PloidyFrost
                               estimate files
        outputFileName -- a string indicating the output file name for filtered result
        samplePrefix -- a string to use for the plot title; typically the sample name
        pfDrawPath -- a string indicating the location of the ploidyfrost Drawfreq.R file
        rscriptPath -- a string indicating the location of the Rscript executable
    '''
    # Construct the cmd for subprocess
    cmd = ZS_Utility.base_subprocess_cmd(rscriptPath)
    cmd += [
        ZS_Utility.convert_to_wsl_if_not_unix(pfDrawPath),
        "-f", ZS_Utility.convert_to_wsl_if_not_unix(ploidyfrostFileName),
        "-t", samplePrefix, "-p", "2",
        "-o", ZS_Utility.convert_to_wsl_if_not_unix(outputFileName)
    ]
    
    if platform.system() != "Windows":
        cmd = " ".join(cmd)
    
    # Run the command
    run_frost_drawfreq = subprocess.Popen(cmd, shell = True,
                                          stdout = subprocess.PIPE,
                                          stderr = subprocess.PIPE)
    frostout, frosterr = run_frost_drawfreq.communicate()
    
    # Check file outputs to see if there was an error
    if (frostout.decode("utf-8") == "") or (not frosterr.decode("utf-8") == ""):
        raise Exception(("ERROR: run_ploidyfrost_drawfreq encountered an error; have a look " +
                        f'at the stdout ({frostout.decode("utf-8")}) and stderr ' + 
                        f'({frosterr.decode("utf-8")}) to make sense of this.'))

## Main
def main():
    # User input
    usage = """%(prog)s accepts a directory containing one or more FASTQ files in single or paired end.
    It will run ploidyfrost for each sample and tabulate their most likely ploidy.
    """
    p = argparse.ArgumentParser(description=usage)
    # Reqs
    p.add_argument("-q", dest="fastqDirectory",
                   required=True,
                   help="Input directory containing FASTQ file(s)")
    p.add_argument("-o", dest="outputDirectory",
                   required=True,
                   help="Specify location to write output files to")
    p.add_argument("-p", dest="pfScriptDir",
                   required=True,
                   help="""Specify the location of the PloidyFrost script directory;
                   this should contain the Filter.R and Drawfreq.R scripts""")
    # Opts (programs)
    p.add_argument("--ploidyfrost", dest="ploidyfrost",
                   required=False,
                   help="""Optionally, specify the PloidyFrost executable file
                   if it is not discoverable in the path""",
                   default=None)
    p.add_argument("--bifrost", dest="bifrost",
                   required=False,
                   help="""Optionally, specify the Bifrost executable file
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
    p.add_argument("--rscript", dest="rscript",
                   required=False,
                   help="""Optionally, specify the Rscript file
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
    p.add_argument("--fastqSuffix", dest="fastqSuffix",
                   required=False,
                   help="""Indicate the suffix of the FASTQ files to look for;
                   default == '.fq.gz'""",
                   default=".fq.gz")
    p.add_argument("--isSingleEnd", dest="isSingleEnd",
                   required=False,
                   action="store_true",
                   help="Specify if files are single-ended",
                   default=False)
    
    args = p.parse_args()
    validate_args(args)
    
    # Locate and validate FASTQ files
    forwardReads, reverseReads = ZS_SeqIO.FASTA.get_rnaseq_files(args.fastqDirectory, args.fastqSuffix, args.isSingleEnd)
    
    # Reformat and tie paired reads together
    forwardReads = [ os.path.basename(f) for f in forwardReads ]
    reverseReads = [ os.path.basename(r) for r in reverseReads ] if reverseReads is not None else None
    
    if reverseReads is not None:
        pairedReads = [ [f, r] for f, r in zip(forwardReads, reverseReads) ]
    else:
        pairedReads = [ [f] for f in forwardReads ]
    
    # Symlink FASTQ files into working directory
    fastqsDir = os.path.abspath(os.path.join(args.outputDirectory, "fastqs"))
    os.makedirs(fastqsDir, exist_ok=True)
    
    if not os.path.exists(os.path.join(args.outputDirectory, "fastq_symlink_was_successful.flag")):
        # Symlink each FASTQ file pair
        for pair in pairedReads:
            for fastqFile in pair:
                outputFastqFile = os.path.join(fastqsDir, fastqFile)
                if not os.path.exists(outputFastqFile):
                    os.symlink(os.path.join(args.fastqDirectory, fastqFile), outputFastqFile)
        open(os.path.join(args.outputDirectory, "fastq_symlink_was_successful.flag"), "w").close()
    else:
        print(f"fastq symlinking has already been performing; skipping.")
    
    # Set up the working directory structure
    kmcDir = os.path.abspath(os.path.join(args.outputDirectory, "kmc_files"))
    os.makedirs(kmcDir, exist_ok=True)
    
    kmcTmpDir = os.path.join(kmcDir, "tmp")
    os.makedirs(kmcTmpDir, exist_ok=True)
    
    ploidyfrostDir = os.path.abspath(os.path.join(args.outputDirectory, "ploidyfrost_files"))
    os.makedirs(ploidyfrostDir, exist_ok=True)
    
    # Run kmc->ploidyfrost for each FASTQ file pairing
    if not os.path.exists(os.path.join(args.outputDirectory, "pipeline_was_successful.flag")):
        # Iterate through each FASTQ file pair
        for pair in pairedReads:
            samplePrefix = pair[0].replace(args.fastqSuffix, "")
            
            # Write @FILES file for kmc
            atFilesName = os.path.join(kmcDir, samplePrefix + ".atfile")
            if not os.path.exists(atFilesName):
                generate_atfile_file(pair, fastqsDir, atFilesName)
            else:
                print(f"@FILE generation has already been run for '{samplePrefix}'; skipping.")
            
            # Run kmc
            EXPECTED_SUFFIXES = [".kmc_pre", ".kmc_suf"] # used for program resumption
            
            kmcdbPrefix = os.path.join(kmcDir, samplePrefix + "_db")
            if not all([ os.path.exists(kmcdbPrefix + suf) for suf in EXPECTED_SUFFIXES ]):
                run_kmc(atFilesName, kmcdbPrefix, args.cpus, args.mem, kmcTmpDir, args.kmc)
            else:
                print(f"kmc db has already been generated for '{samplePrefix}'; skipping.")
            
            # Run kmc_tools histogram
            histFileName = kmcdbPrefix + ".hist"
            if not os.path.exists(histFileName):
                run_kmctools_histogram(kmcdbPrefix, histFileName, args.kmc_tools)
            else:
                print(f"kmc_tools histogram has already been run for '{samplePrefix}'; skipping.")
            
            # Run kmc_tools filter
            lCutoff = None # used for later to see if we need to re-compute this value
            
            filterFileName = os.path.join(kmcDir, samplePrefix + "_filtered.fq")
            if not os.path.exists(filterFileName):
                lCutoff = run_ploidyfrost_cutoff(histFileName, "L", args.ploidyfrost)
                run_kmctools_filter(kmcdbPrefix, atFilesName, lCutoff, args.cpus, filterFileName, args.kmc_tools)
            else:
                print(f"kmc_tools filter has already been run for '{samplePrefix}'; skipping.")
            
            # Run bifrost for CBDG graph construction            
            bifrostPrefix = os.path.join(kmcDir, samplePrefix + "_dbg")
            bifrostFileName = bifrostPrefix + ".gfa"
            bifrostGzFileName = bifrostFileName + ".gz"
            
            if (not os.path.exists(bifrostFileName)) and (not os.path.exists(bifrostGzFileName)):
                run_bifrost_build(filterFileName, args.cpus, bifrostPrefix, args.bifrost)
            else:
                print(f"bifrost build has already been run for '{samplePrefix}'; skipping.")
            
            # Gunzip the GFA file
            if os.path.exists(bifrostGzFileName):
                gunzip(bifrostGzFileName)
            else:
                print(f"bifrost gfa has already been gunzipped for '{samplePrefix}'; skipping.")
            
            # Run ploidyfrost for ploidy prediction
            ploidyfrostFileName = os.path.join(ploidyfrostDir, samplePrefix + "_pf")
            if not os.path.exists(ploidyfrostFileName):
                lCutoff = run_ploidyfrost_cutoff(histFileName, "L", args.ploidyfrost) if lCutoff is None else lCutoff
                uCutoff = run_ploidyfrost_cutoff(histFileName, "U", args.ploidyfrost)
                run_ploidyfrost_estimate(bifrostFileName, kmcdbPrefix, ploidyfrostFileName,
                                         lCutoff, uCutoff, args.cpus, args.ploidyfrost)
            else:
                print(f"ploidyfrost estimate has already been run for '{samplePrefix}'; skipping.")
            
            # Run ploidyfrost filter
            pfFilterPrefix = ploidyfrostFileName + "filt"
            pfFilterFileName = pfFilterPrefix + "-filtered_allele_frequency.txt"
            if not os.path.exists(pfFilterFileName):
                run_ploidyfrost_filter(ploidyfrostFileName, pfFilterPrefix, args.pfFilter, args.rscript)
            else:
                print(f"ploidyfrost filter has already been run for '{samplePrefix}'; skipping.")
            
            # Run ploidyfrost gmm
            gmmFileName = os.path.join(args.outputDirectory, samplePrefix + "_pf-gmm")
            if not os.path.exists(gmmFileName):
                run_ploidyfrost_gmm(pfFilterFileName, gmmFileName, args.ploidyfrost)
            else:
                print(f"ploidyfrost gmm has already been run for '{samplePrefix}'; skipping.")
            
            # Run ploidyfrost drawfreq
            drawfreqFileName = os.path.join(args.outputDirectory, samplePrefix + "_pf-drawfreq")
            if not os.path.exists(drawfreqFileName):
                run_ploidyfrost_drawfreq(pfFilterFileName, drawfreqFileName, samplePrefix, args.pfDrawFreq, args.rscript)
            else:
                print(f"ploidyfrost drawfreq has already been run for '{samplePrefix}'; skipping.")
        
        open(os.path.join(args.outputDirectory, "pipeline_was_successful.flag"), "w").close()
    else:
        print(f"kmc->ploidyfrost pipeline has already been performed; skipping.")
    
    # Let user know everything went swimmingly
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

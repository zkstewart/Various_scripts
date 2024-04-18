#! python3
# smudgeplot_pipeline.py
# Script to automatically run smudgeplot analysis for the
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
    if args.smudgeplot is None:
        args.smudgeplot = ZS_Utility.wsl_which("smudgeplot.py")
        if args.smudgeplot is None:
            _not_specified_error("smudgeplot.py")
    else:
        if not os.path.isfile(args.smudgeplot):
            _not_found_error("smudgeplot.py", args.smudgeplot)
    
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
        atFilesName -- a string indicating the location of the @FILES file containing
                       input file names
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
    
    # Check file outputs to see if there was an error
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
    
    # Check file outputs to see if there was an error
    if (smudgeout.decode("utf-8") == "") or (not "Done!" in smudgerr.decode("utf-8")):
        raise Exception(("ERROR: run_smudgeplot_cutoff encountered an error; have a look " +
                        f'at the stdout ({smudgeout.decode("utf-8")}) and stderr ' + 
                        f'({smudgerr.decode("utf-8")}) to make sense of this.'))
    
    # Return cutoff value as printed to stdout
    return smudgeout.decode("utf-8")

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
    
    # Check file outputs to see if there was an error
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
    
    # Check file outputs to see if there was an error
    if (smudgeout.decode("utf-8") != "") or (not "Done!" in smudgeerr.decode("utf-8")):
        raise Exception(("ERROR: run_smudgeplot_hetkmers encountered an error; have a look " +
                        f'at the stdout ({smudgeout.decode("utf-8")}) and stderr ' + 
                        f'({smudgeerr.decode("utf-8")}) to make sense of this.'))

def run_smudgeplot_plot(smudgeCoveragesFile, outputFileName, samplePrefix, smudgeplotPath):
    '''
    Parameters:
        smudgeCoveragesFile -- a string indicating the location of the smudgeplot
                               _coverages.tsv file
        samplePrefix -- a string to use for the plot title; typically the sample name
        outputFileName -- a string indicating the output file name for the plot result
        smudgeplotPath -- a string indicating the location of the smudgeplot.py executable
    '''
    # Construct the cmd for subprocess
    cmd = ZS_Utility.base_subprocess_cmd(smudgeplotPath)
    cmd += [
        "plot", "-t", samplePrefix,
        "-o", ZS_Utility.convert_to_wsl_if_not_unix(outputFileName),
        ZS_Utility.convert_to_wsl_if_not_unix(smudgeCoveragesFile)
        
    ]
    
    if platform.system() != "Windows":
        cmd = " ".join(cmd)
    
    # Run the command
    run_smudge_plot = subprocess.Popen(cmd, shell = True,
                                         stdout = subprocess.PIPE,
                                         stderr = subprocess.PIPE)
    smudgeout, smudgeerr = run_smudge_plot.communicate()
    
    # Check file outputs to see if there was an error
    if (smudgeout.decode("utf-8") != "") or (not "Done!" in smudgeerr.decode("utf-8")):
        raise Exception(("ERROR: run_smudge_plot encountered an error; have a look " +
                        f'at the stdout ({smudgeout.decode("utf-8")}) and stderr ' + 
                        f'({smudgeerr.decode("utf-8")}) to make sense of this.'))

## Main
def main():
    # User input
    usage = """%(prog)s accepts a directory containing one or more FASTQ files in single or paired end.
    It will run smudgeplot for each sample and tabulate their most likely ploidy.
    """
    p = argparse.ArgumentParser(description=usage)
    # Reqs
    p.add_argument("-q", dest="fastqDirectory",
                   required=True,
                   help="Input directory containing FASTQ file(s)")
    p.add_argument("-o", dest="outputDirectory",
                   required=True,
                   help="Specify location to write output files to")
    # Opts (programs)
    p.add_argument("--smudgeplot", dest="smudgeplot",
                   required=False,
                   help="""Optionally, specify the smudgeplot.py file
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
    
    smudgeplotDir = os.path.abspath(os.path.join(args.outputDirectory, "smudgeplot_files"))
    os.makedirs(smudgeplotDir, exist_ok=True)
    
    # Run kmc->smudgeplot for each FASTQ file pairing
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
            
            # Run kmc_tools dump
            dumpFileName = kmcdbPrefix + ".dump"
            if not os.path.exists(dumpFileName):
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
            plotFileName = os.path.join(args.outputDirectory, samplePrefix + "_plot.txt")
            if not os.path.exists(hetkmersFileName):
                run_smudgeplot_plot(hetkmersFileName, plotFileName, samplePrefix, args.smudgeplot)
            else:
                print(f"smudgeplot plot has already been run for '{samplePrefix}'; skipping.")
        
        open(os.path.join(args.outputDirectory, "pipeline_was_successful.flag"), "w").close()
    else:
        print(f"kmc->smudgeplot pipeline has already been performed; skipping.")
    
    # Let user know everything went swimmingly
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

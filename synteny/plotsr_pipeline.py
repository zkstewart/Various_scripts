#! python3
# plotsr_pipeline.py
# Script to take in a directory of FASTA files and perform pairwise minimap2 alignments
# with subsequent syri structural variant calling and plotting of results with plotsr.

import os, argparse, sys, shutil, platform, subprocess

# Define functions
def validate_args(args):
    def _not_found_error(program):
        raise FileNotFoundError(f"{program} not discoverable in your system PATH and was not specified as an argument.")
    def _specified_wrong_error(program, path):
        raise FileNotFoundError(f"{program} was not found at the indicated location '{path}'")
    
    # Locate and validate input files
    args.fastaFiles = []
    fastaPrefixes = set()
    for location in args.inputLocations:
        location = os.path.abspath(location)
        # Handle files
        if os.path.isfile(location):
            args.fastaFiles.append(location)
        # Handle directories
        elif os.path.isdir(location):
            foundAny = False
            for f in os.listdir(location):
                if any([ f.endswith(x) for x in args.fastaSuffix ]):
                    fastaPrefix = os.path.basename(location).rsplit(".", maxsplit=1)[0]
                    if fastaPrefix in fastaPrefixes:
                        raise ValueError(f"Duplicate FASTA prefix found: '{fastaPrefix}'")
                    
                    args.fastaFiles.append(os.path.join(location, f))
                    foundAny = True
            if not foundAny:
                raise FileNotFoundError(f"No FASTA files found in directory '{location}' ending with a suffix in '{args.fastaSuffix}'")
        # Handle other cases
        else:
            raise FileNotFoundError(f"Input FASTA file or directory '{location}' not found!")
    if len(args.fastaFiles) < 2:
        raise ValueError(f"At least two FASTA files are required for pairwise alignment; only {len(args.fastaFiles)} found.")
    
    # Validate numeric arguments
    if args.threads < 1:
        raise ValueError(f"--threads must be greater than or equal to 1.")
    
    # Validate program discoverability
    if args.minimap2 == None:
        args.minimap2 = shutil.which("minimap2")
        if args.minimap2 == None:
            _not_found_error("minimap2")
    else:
        if not os.path.isfile(args.minimap2):
            _specified_wrong_error("minimap2", args.minimap2)
    
    if args.samtools == None:
        args.samtools = shutil.which("samtools")
        if args.samtools == None:
            _not_found_error("samtools")
    else:
        if not os.path.isfile(args.samtools):
            _specified_wrong_error("samtools", args.samtools)
    
    if args.syri == None:
        args.syri = shutil.which("syri")
        if args.syri == None:
            _not_found_error("syri")
    else:
        if not os.path.isfile(args.syri):
            _specified_wrong_error("syri", args.syri)
    
    if args.plotsr == None:
        args.plotsr = shutil.which("plotsr")
        if args.plotsr == None:
            _not_found_error("plotsr")
    else:
        if not os.path.isfile(args.plotsr):
            _specified_wrong_error("plotsr", args.plotsr)
    
    # Validate output file location
    args.outputDirectory = os.path.abspath(args.outputDirectory)
    if os.path.isdir(args.outputDirectory) and os.listdir(args.outputDirectory) != []:
        print(f"Output directory '{args.outputDirectory}' already exists; I'll write output files here.")
        print("But, I won't overwrite any existing files, so beware that if a previous run had issues, " +
              "you may need to delete/move files first.")
    if not os.path.isdir(args.outputDirectory):
        os.makedirs(args.outputDirectory)
        print(f"Output directory '{args.outputDirectory}' has been created as part of argument validation.")

def run_sorted_minimap2(fasta1, fasta2, outputFileName, preset,
                        minimap2Path, samtoolsPath, threads=1):
    '''
    Runs minimap2 to align two FASTA files and generate a sorted BAM file.
    
    Parameters:
        fasta1 -- a string indicating the location of the first FASTA file
        fasta2 -- a string indicating the location of the second FASTA file
        outputFileName -- a string indicating the location of the output BAM file
        preset -- a string indicating the preset to use for minimap2
        minimap2Path -- a string indicating the location of the minimap2 executable file
        samtoolsPath -- a string indicating the location of the samtools executable file
        threads -- (OPTIONAL) an integer indicating the number of threads to use for minimap2
    '''
    # Format syri command
    cmd = [
        minimap2Path, "-ax", preset, "-t", str(threads), "--eqx", fasta1, fasta2,
        "|",
        samtoolsPath, "sort", "-O", "BAM",  "-o", outputFileName, "-"
    ]
    
    # Run syri
    if platform.system() != "Windows":
        runMinimap2 = subprocess.Popen(" ".join(cmd), shell = True,
                                       stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
    else:
        runMinimap2 = subprocess.Popen(cmd, shell = True,
                                       stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
    mm2Out, mm2Err = runMinimap2.communicate()
    if not "Real time:" in mm2Err.decode("utf-8"):
            raise Exception('minimap2 error text below\n' + mm2Err.decode("utf-8"))

def run_syri(bam, fasta1, fasta2, prefix, outputDir, syriPath):
    '''
    Runs syri to call structural variants between the two FASTA files.
    
    Parameters:
        bam -- a string indicating the location of the BAM file, generated by minimap2,
               to be used as input for syri
        fasta1 -- a string indicating the location of the first FASTA file
        fasta2 -- a string indicating the location of the second FASTA file
        prefix -- a string indicating the prefix to use for the output files
        outputDir -- a string indicating the location to write syri outputs
        syriPath -- a string indicating the location of the syri executable file
    '''
    # Format syri command
    cmd = [
        syriPath, "-c", bam, "-r", fasta1, "-q", fasta2,
        "-F", "B", "--prefix", prefix, "--dir", outputDir
    ]
    
    # Run syri
    if platform.system() != "Windows":
        runSyri = subprocess.Popen(" ".join(cmd), shell = True,
                                   stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
    else:
        runSyri = subprocess.Popen(cmd, shell = True,
                                   stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
    syriOut, syriErr = runSyri.communicate()
    print(f"stdout: {syriOut.decode('utf-8')}") # for debugging during script development
    print(f"stderr: {syriErr.decode('utf-8')}")
    #if not "INFO: Finished running" in syriErr.decode("utf-8"):
    #    raise Exception('syri error text below\n' + syriErr.decode("utf-8"))

def run_plotsr(syriFiles, genomeTxtFile, outputFileName, plotsrPath):
    '''
    Runs plotsr to visualise structural variants between syri output files.
    
    Parameters:
        syriFiles -- a list of strings indicating the locations of the syri output files
        genomeTxtFile -- a string indicating the location of the genomes.txt file
        outputFileName -- a string indicating the location to write the plotsr output
        plotsrPath -- a string indicating the location of the plotsr executable file
    '''
    # Format syri command
    cmd = [plotsrPath]
    for syriFile in syriFiles:
        cmd += ["--sr", syriFile]
    cmd += ["--genomes", genomeTxtFile, "-o", outputFileName]
    
    # Run syri
    if platform.system() != "Windows":
        runPlotsr = subprocess.Popen(" ".join(cmd), shell = True,
                                     stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
    else:
        runPlotsr = subprocess.Popen(cmd, shell = True,
                                     stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
    plotsrOut, plotsrErr = runPlotsr.communicate()
    print(f"stdout: {plotsrOut.decode('utf-8')}") # for debugging during script development
    print(f"stderr: {plotsrErr.decode('utf-8')}")
    #if not "INFO: Finished running" in plotsrErr.decode("utf-8"):
    #    raise Exception('plotsr error text below\n' + plotsrErr.decode("utf-8"))

## Main
def main():
    # User input
    usage = """%(prog)s ...
    """
    p = argparse.ArgumentParser(description=usage)
    # Reqs
    p.add_argument("-i", dest="inputLocations",
                   required=True,
                   nargs="+",
                   help="Input one or more file names and/or directories containing FASTA files")
    p.add_argument("-o", dest="outputDirectory",
                   required=True,
                   help="Specify location to write output files to")
    # Opts (minimap2)
    p.add_argument("--preset", dest="preset",
                   required=False,
                   choices=["asm5", "asm10", "asm20"],
                   help="""Optionally, specify the preset to use for minimap2;
                   default == 'asm20'""",
                   default="asm20")
    p.add_argument("--threads", dest="threads",
                   required=False,
                   type=int,
                   help="""Optionally, specify the number of threads to use for minimap2;
                   default == 1""",
                   default=1)
    # Opts (programs)
    p.add_argument("--fastaSuffix", dest="fastaSuffix",
                   required=False,
                   nargs="+",
                   help="""Optionally, specify one or more suffixes to use when scanning for FASTA
                   files in input directories; default == '.fasta'""",
                   default=[".fasta"])
    p.add_argument("--minimap2", dest="minimap2",
                   required=False,
                   help="""Optionally, specify the minimap2 executable file
                   if it is not discoverable in the path""",
                   default=None)
    p.add_argument("--samtools", dest="samtools",
                   required=False,
                   help="""Optionally, specify the samtools executable file
                   if it is not discoverable in the path""",
                   default=None)
    p.add_argument("--syri", dest="syri",
                   required=False,
                   help="""Optionally, specify the syri executable file
                   if it is not discoverable in the path""",
                   default=None)
    p.add_argument("--plotsr", dest="plotsr",
                   required=False,
                   help="""Optionally, specify the plotsr executable file
                   if it is not discoverable in the path""",
                   default=None)
    
    args = p.parse_args()
    validate_args(args)
    
    # Determine the most parsimonious ordering of FASTA files
    ## TBD if this is needed
    
    # Let user know the order of FASTA files
    print("Input FASTA files in the following order:")
    for index, fastaFile in enumerate(args.fastaFiles):
        print(f"{index + 1}: {fastaFile}")
    
    # Validate that all files have the same sequence IDs
    firstSeqIDs = set()
    for index, fastaFile in enumerate(args.fastaFiles):
        # Extract the sequence IDs from the input FASTA file
        inputSeqIDs = set()
        with open(fastaFile, "r") as fileIn:
            for line in fileIn:
                if line.startswith(">"):
                    seqID = line.split()[0][1:]
                    inputSeqIDs.add(seqID)
        
        # Check if the input sequence IDs match the first file's sequence IDs
        if index == 0:
            firstSeqIDs = inputSeqIDs
        else:
            if firstSeqIDs != inputSeqIDs:
                raise ValueError(f"Input FASTA file '{fastaFile}' has sequence IDs not found in " + 
                                f"the first input file '{args.fastaFiles[0]}'")
    
    # Set up the working directory structure
    alignmentsDir = os.path.join(args.outputDirectory, "alignments")
    os.makedirs(alignmentsDir, exist_ok=True)
    
    syriDir = os.path.join(args.outputDirectory, "syri")
    os.makedirs(syriDir, exist_ok=True)
    
    # Run minimap2 of input files against one another
    alignmentFiles = []
    for i in range(0, len(args.fastaFiles)-1):
        # Get details of the two FASTA files to align
        fastaFile1 = args.fastaFiles[i]
        fastaFile2 = args.fastaFiles[i+1]
        
        fastaPrefix1 = os.path.basename(fastaFile1).rsplit(".", maxsplit=1)[0]
        fastaPrefix2 = os.path.basename(fastaFile2).rsplit(".", maxsplit=1)[0]
        
        # Format output file and flag names
        prefix = f"{fastaPrefix1}_{fastaPrefix2}"
        outputFileName = os.path.join(alignmentsDir, f"{prefix}.sorted.bam")
        flagName = os.path.join(alignmentsDir, f"{prefix}.is.ok.flag")
        
        # Run if the output file or flag file doesn't exist
        if not os.path.exists(flagName) or not os.path.exists(outputFileName):
            run_sorted_minimap2(fastaFile1, fastaFile2, outputFileName, args.preset,
                                args.minimap2, args.samtools, args.threads)
            open(flagName, "w").close()
        
        # Store the file for syri processing
        alignmentFiles.append([prefix, fastaFile1, fastaFile2, outputFileName])
    
    # Run syri on the BAM files generated by minimap2
    syriFiles = []
    for prefix, fastaFile1, fastaFile2, bamFile in alignmentFiles:
        # Format output file and flag names
        syriOutputName = os.path.join(syriDir, f"{prefix}syri.out")
        flagName = os.path.join(syriDir, f"{prefix}.is.ok.flag")
        
        # Run if the output file or flag file doesn't exist
        if not os.path.exists(syriOutputName) or not os.path.exists(flagName):
            run_syri(bamFileName, fastaFile1, fastaFile2, prefix, syriDir, args.syri)
            open(flagName, "w").close()
        
        # Store the file for plotsr processing
        syriFiles.append(syriOutputName)
    
    # Generate a genomes.txt file for plotsr
    genomeTxtFile = os.path.join(args.outputDirectory, "genomes.txt")
    with open(genomeTxtFile, "w") as fileOut:
        # Write the header line
        fileOut.write("#file\tname\ttags\n")
        # Write details for each FASTA file
        for fastaFile in args.fastaFiles:
            fastaPrefix = os.path.basename(fastaFile).rsplit(".", maxsplit=1)[0]
            fileOut.write(f"{fastaFile}\t{fastaPrefix}\tlw:1.5\n")
    
    # Run plotsr to generate the plots
    pngPlotFile = os.path.join(args.outputDirectory, "plotsr.png")
    pngFlagFile = os.path.join(args.outputDirectory, "plotsr.png.is.ok.flag")
    if not os.path.exists(pngPlotFile) or not os.path.exists(pngFlagFile):
        run_plotsr(syriFiles, genomeTxtFile, pngPlotFile, plotsrPath)
        open(pngFlagFile, "w").close()
    else:
        print(f"png output file already exists; skipping.")
    
    pdfPlotFile = os.path.join(args.outputDirectory, "plotsr.pdf")
    pdfFlagFile = os.path.join(args.outputDirectory, "plotsr.pdf.is.ok.flag")
    if not os.path.exists(pdfPlotFile) or not os.path.exists(pdfFlagFile):
        run_plotsr(syriFiles, genomeTxtFile, pdfPlotFile, plotsrPath)
        open(pdfFlagFile, "w").close()
    else:
        print(f"pdf output file already exists; skipping.")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

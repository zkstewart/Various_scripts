#! python3
# binge_to_msa.py
# Takes a BINge clustering output file, obtains each cluster as a FASTA file, then
# aligns them using MAFFT, outputting the aligned sequences.

import os, argparse, shutil, subprocess, sys
from pyfaidx import Fasta

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) # 2 dirs up is where we find GFF3IO
from Function_packages import ZS_SeqIO

MAFFT_ALIGN_PY = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "Fasta_related", "Alignment", "Auto_aligners", "mafftAlign.py")

class FastaCollection:
    '''
    Wrapper for pyfaidx Fasta objects which allows multiple to be combined
    and queried as one logical entity.
    
    Parameters:
        fastaFiles -- a list of strings pointing to the locations of FASTA files
                      which are to be loaded in using pyfaidx.Fasta
    '''
    def __init__(self, fastaFiles):
        self.fastaFiles = fastaFiles
        self.records = []
        self.pipePrefixes = []
        
        self._parse_fastas()
    
    def _parse_fastas(self):
        for fastaFile in self.fastaFiles:
            self.records.append(Fasta(fastaFile))
            # Extract pipe prefixes from FASTA titles
            pipePrefixes = set()
            with open(fastaFile, "r") as fileIn:
                for line in fileIn:
                    if line.startswith(">"):
                        seqPrefix = line[1:].split(None, 1)[0]
                        if "|" in seqPrefix:
                            pipePrefix = seqPrefix.split("|")[0] + "|"
                            pipePrefixes.add(pipePrefix)
            self.pipePrefixes.append(pipePrefixes)
    
    def __getitem__(self, key):
        for i, records in enumerate(self.records):
            try:
                return records[key]
            except:
                for pipePrefix in self.pipePrefixes[i]:
                    try:
                        return records[pipePrefix + key]
                    except:
                        pass
        raise KeyError(f"'{key}' not found in collection")
    
    def __contains__(self, key):
        try:
            self[key] # __getitem__ raises exception if the key isn't found
            return True
        except:
            return False
    
    def __iter__(self):
        for records in self.records:
            yield from records
    
    def __repr__(self):
        return (f"<FastaCollection object;num_records='{len(self.records)}';" +
                f"fastaFiles={self.fastaFiles}"
        )

def validate_args(args):
    # Validate input directory location
    if not os.path.isdir(args.workingDirectory):
        raise FileNotFoundError(f"-i '{args.workingDirectory} does not exist!")
    
    # Validate BINge clustering output file
    args.bingeFile = os.path.join(args.workingDirectory, "analysis", args.analysisFolder, "BINge_clustering_result.tsv")
    if not os.path.isfile(args.bingeFile):
        raise FileNotFoundError(f"BINge clustering output file '{args.bingeFile}' does not exist! " +
                                "Please ensure that the BINge clustering has been performed in the specified working directory.")
    
    # Validate MAFFT arguments
    if args.mafft is None:
        args.mafft = shutil.which("mafft")
        if args.mafft is None:
            raise FileNotFoundError(f"'mafft' not discoverable in your system PATH and was not specified as an argument.")
    else:
        if not os.path.isfile(args.mafft):
            raise FileNotFoundError(f"'mafft' was not found at the location indicated (--mafft {args.mafft})")
    if args.maxiterate < 0:
        raise ValueError(f"--maxiterate must be a non-negative integer.")
    if args.threads < 1:
        raise ValueError(f"--threads must be an integer greater than or equal to 1.")
    
    # Validate output file locations
    if os.path.exists(args.outputDirectory):
        print(f"Output directory '{args.outputDirectory}' already exists; will attempt to resume program operations ...")
    else:
        os.makedirs(args.outputDirectory)
        print(f"Created output directory '{args.outputDirectory}' as " + 
              "part of argument validation")

def parse_binge_clusters(bingeFile):
    '''
    Reads in the output file of BINge as a dictionary assocating clusters to their
    sequence members.
    
    Parameters:
        bingeFile -- a string pointing to the location of a BINge cluster output file.
    Returns:
        clusterDict -- a dictionary with structure like:
                       {
                             0: [seqid1, seqid2, ...],
                             1: [ ... ],
                             ...
                         }
    '''    
    clusterDict = {}
    lineNum = 0
    with open(bingeFile, "r") as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n ").split("\t")
            
            # Handle header lines
            if lineNum == 0:
                assert line.startswith("#BINge clustering information file"), \
                    ("BINge file is expected to start with a specific comment line! " +
                    "Your file is hence not recognised as a valid BINge cluster file.")
                lineNum += 1
            elif lineNum == 1:
                assert sl == ["cluster_num", "sequence_id", "cluster_type"], \
                    ("BINge file is expected to have a specific header line on the second line! " +
                     "Your file is hence not recognised as a valid BINge cluster file.")
                lineNum += 1
            
            # Handle content lines
            else:
                clustNum, seqID, clusterType = int(sl[0]), sl[1], sl[2]
                clusterDict.setdefault(clustNum, [])
                clusterDict[clustNum].append(seqID)
    return clusterDict

def run_mafftAlignPy(inputLocations, outputDirectory, mafftExe, fastaExtensions,
                     algorithm=None, maxiterate=0, threads=1, codons=False):
    '''
    Parameters:
        inputLocations -- a string indicating the location of the files
                          and/or directories containing files to be aligned
        outputDirectory -- a string indicating the location where the aligned files
                           should be written
        mafftExe -- a string indicating the location of the MAFFT executable
        fastaExtensions -- a list of file extensions to consider as FASTA files
        algorithm -- a string indicating the MAFFT algorithm to use; if None, will
                     use the default algorithm
        maxiterate -- an integer indicating the number of iterations to use for
                      alignment; default == 0 (no iterations)
    '''
    # Validate input values
    assert isinstance(inputLocations, list)
    assert isinstance(outputDirectory, str)
    assert isinstance(mafftExe, str)
    assert isinstance(fastaExtensions, list)
    assert maxiterate >= 0
    assert threads >= 1
    assert isinstance(codons, bool)
    
    # Construct the cmd for subprocess
    cmd = [
        "python", MAFFT_ALIGN_PY, "-i", *inputLocations,
        "-o", outputDirectory, "--mafft", mafftExe,
        "--fastaExtensions", *fastaExtensions,
        "--alg", algorithm if algorithm else "auto",
        "--maxiterate", str(maxiterate), "--threads", str(threads),
        "--codons" if codons else ""
    ]
    if cmd[-1] == "":
        cmd.pop()
    
    # Run the command
    run_mafftAlignPy = subprocess.Popen(" ".join(cmd), shell = True,
                                        stdout = subprocess.PIPE,
                                        stderr = subprocess.PIPE)
    
    maffpyout, mafftpyerr = run_mafftAlignPy.communicate()
    if not "Program completed successfully!" in maffpyout.decode("utf-8").strip():
        raise Exception(("run_mafftAlignPy encountered an error; have a look " +
                        f'at the stdout ({maffpyout.decode("utf-8")}) and stderr ' + 
                        f'({mafftpyerr.decode("utf-8")}) to make sense of this.'))

def main():
    usage = """%(prog)s will receive a BINge analysis directory and produce MSA files
    for each cluster in the analysis. Depending on the size of your results, you may
    want to filter the clusters to only those you are interested in.
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="workingDirectory",
                   required=True,
                   help="Specify the location where BINge clustering has been performed")
    p.add_argument("-o", dest="outputDirectory",
                   required=True,
                   help="Output directory where results are to be written")
    # Optional
    p.add_argument("--analysis", dest="analysisFolder",
                   required=False,
                   help="""Specify the analysis folder to filter; if not provided,
                   the most recent analysis folder will be viewed""",
                   default="most_recent")
    p.add_argument("--threads", dest="threads",
                   required=False,
                   type=int,
                   help="""Optionally, specify the number of threads to use for alignment;
                   default == 1""",
                   default=1)
    p.add_argument("--mafft", dest="mafft",
                   required=False,
                   help="""Optionally, specify the mafft executable file
                   if it is not discoverable in the path""",
                   default=None)
    p.add_argument("--alg", dest="algorithm",
                   required=False,
                   choices=["auto", "einsi", "linsi", "ginsi", "fftns1", "fftns2", "fftnsi"],
                   help="""Optionally, specify the MAFFT algorithm to use for alignment;
                   default == 'auto'""",
                   default=None)
    p.add_argument("--maxiterate", dest="maxiterate",
                   required=False,
                   type=int,
                   help="""Optionally, specify the --maxiterate value to use for alignment;
                   default == 0""",
                   default=0)
    
    args = p.parse_args()
    validate_args(args)
    
    # Prepare output directory
    clusterFilesDir = os.path.join(args.outputDirectory, "raw")
    os.makedirs(clusterFilesDir, exist_ok=True)
    
    aaAlignedDir = os.path.join(args.outputDirectory, "aligned_aa")
    os.makedirs(aaAlignedDir, exist_ok=True)
    
    cdsAlignedDir = os.path.join(args.outputDirectory, "aligned_cds")
    os.makedirs(cdsAlignedDir, exist_ok=True)
    
    # Parse the BINge clusters
    bingeDict = parse_binge_clusters(args.bingeFile)
    
    # Locate the amino acid and CDS FASTA files
    sequencesDir = os.path.join(args.workingDirectory, "sequences")
    aaFiles = [ os.path.join(sequencesDir, f) for f in os.listdir(sequencesDir) if f.endswith(".aa") ]
    if len(aaFiles) == 0:
        raise FileNotFoundError(f"No amino acid FASTA files found in '{sequencesDir}'! " +
                                "Please ensure that your BINge working directory is correct and intact.")
    
    cdsFiles = [ os.path.join(sequencesDir, f) for f in os.listdir(sequencesDir) if f.endswith(".cds") ]
    if len(cdsFiles) == 0:
        raise FileNotFoundError(f"No CDS FASTA files found in '{sequencesDir}'! " +
                                "Please ensure that your BINge working directory is correct and intact.")
    
    # Parse the FASTA files
    aaFasta = FastaCollection(aaFiles)
    cdsFasta = FastaCollection(cdsFiles)
    
    # Produce each raw cluster FASTA file
    warnedOnce = False
    for clustNum, seqIDs in bingeDict.items():
        aaClusterFile = os.path.join(clusterFilesDir, f"cluster_{clustNum}.aa")
        cdsClusterFile = os.path.join(clusterFilesDir, f"cluster_{clustNum}.cds")
        
        if os.path.isfile(aaClusterFile) or os.path.isfile(cdsClusterFile):
            if not warnedOnce:
                print(f"WARNING: Raw FASTA file for cluster #{clustNum} already exists and I " +
                      "will not overwrite it. If you want to regenerate any raw FASTA files, " +
                      "please remove the existing files.")
                warnedOnce = True
            continue
        
        with open(aaClusterFile, "w") as aaFileOut, open(cdsClusterFile, "w") as cdsFileOut:
            for seqID in seqIDs:
                aaSeq = aaFasta[seqID]
                cdsSeq = cdsFasta[seqID]
                
                aaFileOut.write(f">{seqID}\n{str(aaSeq)}\n")
                cdsFileOut.write(f">{seqID}\n{str(cdsSeq)}\n")
    
    # Align the CDS clusters with mafft
    run_mafftAlignPy([clusterFilesDir], cdsAlignedDir, args.mafft, [".cds"],
                     algorithm=args.algorithm, maxiterate=args.maxiterate,
                     threads=args.threads, codons=True)
    
    # Translate the CDS MSAs to amino acids
    warnedOnce = False
    for msaFile in os.listdir(cdsAlignedDir):
        if msaFile.endswith(".cds"):
            # Get the file names
            cdsFileName = os.path.join(cdsAlignedDir, msaFile)
            aaFileName = os.path.join(aaAlignedDir, msaFile.replace(".cds", ".aa"))
            
            # Skip if the MSA already exists and is not empty
            if os.path.exists(aaFileName) and os.path.getsize(aaFileName) > 0:
                if not warnedOnce:
                    print(f"WARNING: Aligned FASTA file for cluster #{clustNum} already exists and I " +
                        "will not overwrite it. If you want to re-run this or any other alignments, " +
                        "please remove the existing files.")
                    warnedOnce = True
                continue
            
            # Load the CDS MSA as a FASTA object
            FASTA_obj = ZS_SeqIO.FASTA(cdsFileName, isAligned=True)
            
            # Translate the MSA
            translatedFASTA_obj = FASTA_obj.translate()
            
            # Write output file
            translatedFASTA_obj.write(aaFileName, asAligned=True)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

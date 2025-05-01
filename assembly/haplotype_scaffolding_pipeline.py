#! python3
# haplotype_scaffolding_pipeline.py
# Script to take in a hifiasm haplotype assembly and a reference haplotype assembly
# to scaffold the hifiasm assembly using the reference assembly as a guide.

import os, argparse, sys, shutil, re, subprocess, platform
from Bio import SeqIO

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) # 3 dirs up is where we find GFF3IO
from Function_packages import ZS_SeqIO, ZS_AlignIO

# Define functions
def validate_args(args):
    def _not_found_error(program):
        raise FileNotFoundError(f"{program} not discoverable in your system PATH and was not specified as an argument.")
    def _specified_wrong_error(program, path):
        raise FileNotFoundError(f"{program} was not found at the indicated location '{path}'")
    
    # Validate input file locations
    for inputFile in args.inputGenome:
        if not os.path.isfile(inputFile):
            raise FileNotFoundError(f"Input genome file '{inputFile}' not found.")
    for referenceFile in args.referenceGenome:
        if not os.path.isfile(referenceFile):
            raise FileNotFoundError(f"Reference genome file '{referenceFile}' not found.")
    
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
    
    if args.ragtag == None:
        args.ragtag = shutil.which("ragtag.py")
        if args.ragtag == None:
            _not_found_error("ragtag.py")
    else:
        if not os.path.isfile(args.ragtag):
            _specified_wrong_error("ragtag.py", args.ragtag)
    
    # Validate output file location
    args.outputDirectory = os.path.abspath(args.outputDirectory)
    if os.path.isdir(args.outputDirectory) and os.listdir(args.outputDirectory) != []:
        print(f"Output directory '{args.outputDirectory}' already exists; I'll write output files here.")
        print("But, I won't overwrite any existing files, so beware that if a previous run had issues, " +
              "you may need to delete/move files first.")
    if not os.path.isdir(args.outputDirectory):
        os.makedirs(args.outputDirectory)
        print(f"Output directory '{args.outputDirectory}' has been created as part of argument validation.")

def run_ragtag(inputFile, referenceFile, outputDir, ragtagPath, threads=1):
    '''
    Runs RagTag to scaffold the input file using the reference file as a guide.
    
    Parameters:
        inputFile -- a string indicating the location of the file to be scaffolded
        referenceFile -- a string indicating the location of the reference file
        outputDir -- a string indicating the location to write ragtag outputs
        ragtagPath -- a string indicating the location of the ragtag.py file
        threads -- (OPTIONAL) an integer indicating the number of threads when running ragtag
    '''
    # Format ragtag command
    cmd = [
        ragtagPath, "scaffold", "-t", str(threads), "-o", outputDir,
        referenceFile, inputFile
    ]
    
    # Run ragtag
    if platform.system() != "Windows":
        run_ragtag = subprocess.Popen(" ".join(cmd), shell = True,
                                      stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
    else:
        run_ragtag = subprocess.Popen(cmd, shell = True,
                                      stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
    ragtagout, ragtagerr = run_ragtag.communicate()
    if not "INFO: Finished running" in ragtagerr.decode("utf-8"):
        raise Exception('ragtag error text below\n' + ragtagerr.decode("utf-8"))

def find_longest_seq(fastaFile):
    '''
    Parse a FASTA file (without loading it into memory) to find the longest sequence.
    
    Parameters:
        fastaFile -- a string indicating the location of the FASTA file
    Returns:
        longestSeqID -- a string indicating the ID of the longest sequence
        numSeqs -- an integer indicating the number of sequences in the FASTA file
    '''
    numSeqs = 0
    with open(fastaFile, "r") as fastaFile:
        records = SeqIO.parse(fastaFile, "fasta")
        longestSeq = [None, 0] # [seq, length]
        for record in records:
            numSeqs += 1
            length = len(record)
            if length > longestSeq[1]:
                longestSeq[0] = record.id
                longestSeq[1] = length
    return longestSeq[0], numSeqs

## Main
def main():
    MULTILINE_LENGTH = 70
    
    # User input
    usage = """%(prog)s automates the scaffolding of a haplotype assembly using a reference assembly.
    It first ensures that haplotypes are equivalently assigned (e.g., hap1 is, to the best of the program's ability,
    consistently hap1 rather than hap2) and then uses ragtag to scaffold the haplotype assembly where
    necessary. Note 1: Contigs for hap1 and hap2 of the reference genome must be equivalently named.
    Note 2: Defaults for --minQueryAlign and --minAlignLen are based on paf2dotplot.
    """
    p = argparse.ArgumentParser(description=usage)
    # Reqs
    p.add_argument("-i", dest="inputGenome",
                   required=True,
                   nargs=2,
                   help="Input hap1 and hap2 genome FASTA files requiring scaffolding")
    p.add_argument("-r", dest="referenceGenome",
                   required=True,
                   nargs=2,
                   help="Input hap1 and hap2 genome FASTA files for use as a reference")
    p.add_argument("-o", dest="outputDirectory",
                   required=True,
                   help="Specify location to write output files to")
    # Opts (behavioural)
    p.add_argument("--minQueryAlign", dest="minQueryAlign",
                   required=False,
                   type=int,
                   help="""Optionally, specify the minimum cumulative alignment length for
                   a query to be considered a match (i.e., the sum of all alignments
                   that meet at least --minAlignLen); default == 400000""",
                   default=400000)
    p.add_argument("--minAlignLen", dest="minAlignLen",
                   required=False,
                   type=int,
                   help="""Optionally, specify the minimum alignment length to
                   consider a match; default == 10000""",
                   default=10000)
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
    p.add_argument("--ragtag", dest="ragtag",
                   required=False,
                   help="""Optionally, specify the ragtag.py file
                   if it is not discoverable in the path""",
                   default=None)
    
    args = p.parse_args()
    validate_args(args)
    
    # Set up the working directory structure
    referenceDir = os.path.join(args.outputDirectory, "reference")
    os.makedirs(referenceDir, exist_ok=True)
    
    inputDir = os.path.join(args.outputDirectory, "input")
    os.makedirs(inputDir, exist_ok=True)
    
    hap1Dir = os.path.join(args.outputDirectory, "hap1")
    os.makedirs(hap1Dir, exist_ok=True)
    
    hap2Dir = os.path.join(args.outputDirectory, "hap2")
    os.makedirs(hap2Dir, exist_ok=True)
    
    # Symlink files to the working directory locations
    ref1File = os.path.join(referenceDir, "hap1.fasta")
    if not os.path.exists(ref1File):
        os.symlink(args.referenceGenome[0], ref1File)
    
    ref2File = os.path.join(referenceDir, "hap2.fasta")
    if not os.path.exists(ref2File):
        os.symlink(args.referenceGenome[1], ref2File)
    
    input1File = os.path.join(inputDir, "hap1.fasta")
    if not os.path.exists(input1File):
        os.symlink(args.inputGenome[0], input1File)
    
    input2File = os.path.join(inputDir, "hap2.fasta")
    if not os.path.exists(input2File):
        os.symlink(args.inputGenome[1], input2File)
    
    # Index the reference files
    if not os.path.exists(f"{ref1File}.fai"):
        ZS_SeqIO.StandardProgramRunners.samtools_faidx(ref1File, args.samtools)
    
    if not os.path.exists(f"{ref2File}.fai"):
        ZS_SeqIO.StandardProgramRunners.samtools_faidx(ref2File, args.samtools)
    
    # Run minimap2 of input files against the reference files
    minimap2FlagName = os.path.join(inputDir, "minimap2_was_successful.flag")
    if not os.path.exists(minimap2FlagName):
        for inputIndex, queryFile in enumerate([input1File, input2File]):
            for refIndex, refFile in enumerate([ref1File, ref2File]):
                outputFileName = os.path.join(inputDir, f"i{inputIndex+1}_vs_r{refIndex+1}.paf")
                runner = ZS_AlignIO.Minimap2(queryFile, refFile, args.preset, args.minimap2, args.threads)
                runner.minimap2(outputFileName, force=True) # allow overwriting since the flag was not created
        open(minimap2FlagName, "w").close()
    else:
        print(f"Minimap2 alignment has already been performed; skipping.")
    
    # Parse minimap2 PAF files
    i1r1Dict = {}
    i1r2Dict = {}
    i2r1Dict = {}
    i2r2Dict = {}
    irDictList = [[i1r1Dict, i1r2Dict], [i2r1Dict, i2r2Dict]]
    refContigs = set()
    
    for inputIndex in range(2):
        for refIndex in range(2):
            pafDict = irDictList[inputIndex][refIndex]
            pafFile = os.path.join(inputDir, f"i{inputIndex+1}_vs_r{refIndex+1}.paf")
            
            with open(pafFile, "r") as fileIn:
                for line in fileIn:
                    # Extract relevant data
                    sl = line.rstrip("\r\n ").split("\t")
                    qid, qlen, qstart, qend, strand, tid, tlen, tstart, tend, \
                        numresidues, lenalign, mapq = sl[0:12]
                    
                    # Convert to integers
                    qlen, qstart, qend, tlen, tstart, tend, numresidues, lenalign, mapq = \
                        map(int, [qlen, qstart, qend, tlen, tstart, tend, numresidues, lenalign, mapq])
                    refContigs.add(tid)
                    
                    # Skip if the alignment doesn't meet length minimum
                    if lenalign < args.minAlignLen:
                        continue
                    
                    # Store the alignment length
                    pafDict.setdefault(tid, {})
                    pafDict[tid].setdefault(qid, {"numresidues": 0, "lenalign": 0})
                    pafDict[tid][qid]["numresidues"] += numresidues
                    pafDict[tid][qid]["lenalign"] += lenalign
    refContigs = sorted(list(refContigs), key=lambda x: int(" ".join(re.findall(r'\d+', x))))
    
    # Adjust the PAF dictionaries to discount multimapping & filter queries that fail the minimum criteria
    for inputIndex in range(2):
        for refIndex in range(2):
            pafDict = irDictList[inputIndex][refIndex]
            
            # Iterate over each reference contig
            for tid, qidDict in pafDict.items():
                toDelete = []
                # Iterate over each query ID
                for qid in list(qidDict.keys()):
                    # Calculate the identity before adjusting anything
                    identity = qidDict[qid]["numresidues"] / qidDict[qid]["lenalign"]
                    pafDict[tid][qid]["identity"] = identity
                    
                    # Sum up alignment length the query made to all other contigs
                    otherAlignLen = 0
                    for tid2, qidDict2 in pafDict.items():
                        if tid == tid2 or (qid not in qidDict2):
                            continue
                        otherAlignLen += qidDict2[qid]["lenalign"]
                    
                    # Discount the alignment length made to other contigs
                    pafDict[tid][qid]["lenalign"] -= otherAlignLen
                    
                    # If the cumulative query alignment length is less than the minimum, remove it
                    if pafDict[tid][qid]["lenalign"] < args.minQueryAlign:
                        toDelete.append(qid)
                # Remove any query IDs that were discounted
                for qid in toDelete:
                    del pafDict[tid][qid]
    
    # Ensure that query contigs do not map to multiple reference contigs
    for inputIndex in range(2):
        for refIndex in range(2):
            pafDict = irDictList[inputIndex][refIndex]
            foundContigs = set()
            
            # Iterate over each reference contig
            for tid, qidDict in pafDict.items():
                for qid in qidDict.keys():
                    # Behaviour one: error out
                    if qid in foundContigs:
                        raise ValueError(f"Query ID '{qid}' maps to multiple reference contigs; can't resolve this!")
                    # Behaviour two: remove the query ID from all but the best reference contig
                    ## TBD if the behaviour one condition ever occurs
    
    # Associate input haplotype sequences to a reference haplotype
    hap1Dict = {}
    hap2Dict = {}
    for contig in refContigs:
        # Get the alignment lengths for each input/reference pair
        i1r1 = irDictList[0][0][contig]
        i1r2 = irDictList[0][1][contig]
        i2r1 = irDictList[1][0][contig]
        i2r2 = irDictList[1][1][contig]
        
        # Get the keys which will allow comparison
        i1 = set(i1r1.keys()).intersection(set(i1r2.keys()))
        i2 = set(i2r1.keys()).intersection(set(i2r2.keys()))
        
        # Optimise the pairwise assignment of the haplotypes to the references
        "Variable name refers to what the assignment would be e.g., i1 assigned to r1, i2 assigned to r2"
        i1r1_i2r2_identity = sum([ i1r1[k]["identity"] - i1r2[k]["identity"] for k in i1 ]) + sum([ i2r2[k]["identity"] - i2r1[k]["identity"] for k in i2 ])
        #i1r1_i2r2_len = sum([ i1r1[k]["lenalign"] - i1r2[k]["lenalign"] for k in i1 ]) + sum([ i2r2[k]["lenalign"] - i2r1[k]["lenalign"] for k in i2 ])
        #i1r1_i2r2 = sum([ (i1r1[k]["lenalign"] * i1r1[k]["identity"]) - (i1r2[k]["lenalign"] * i1r2[k]["identity"]) for k in i1 ]) + \
        #            sum([ (i2r2[k]["lenalign"] * i2r2[k]["identity"]) - (i2r1[k]["lenalign"] * i2r1[k]["identity"]) for k in i2 ])
        
        if i1r1_i2r2_identity > 0: # optimise for identity as this will give fewer variants
            # Assign haplotype 1 to reference 1 and haplotype 2 to reference 2
            hap1Dict[contig] = list(i1r1.keys())
            hap2Dict[contig] = list(i2r2.keys())
        else:
            # Assign haplotype 1 to reference 2 and haplotype 2 to reference 1
            hap1Dict[contig] = list(i2r1.keys())
            hap2Dict[contig] = list(i1r2.keys())
    
    # Print estimated haplotype assignments
    print("# Based on minimap2 alignments, the following haplotype assignments were made:")
    print("## Haplotype 1:")
    for contig in refContigs:
        print(f"### {contig}: {', '.join(hap1Dict[contig])}")
    print("# Haplotype 2:")
    for contig in refContigs:
        print(f"### {contig}: {', '.join(hap2Dict[contig])}")
    
    # Load in sequences as a FastaCollection object
    inputSequences = ZS_SeqIO.FastaCollection([input1File, input2File])
    refSequences = ZS_SeqIO.FastaCollection([ref1File, ref2File])
    
    # Generate scaffolds for each haplotype chromosome
    for hapDir, hapDict, refFile in zip([hap1Dir, hap2Dir], [hap1Dict, hap2Dict], [ref1File, ref2File]):
        for refContig, inputContigs in hapDict.items():
            outputFileName = os.path.join(hapDir, f"{refContig}.fasta")
            outputFlagName = os.path.join(hapDir, f"{refContig}.is.ok.flag")
            if not os.path.exists(outputFlagName):
                # Output 1:1 relationships
                if len(inputContigs) == 1:
                    inputContig = inputContigs[0]
                    inputSeq = inputSequences[inputContig]
                    with open(outputFileName, "w") as fileOut:
                        fileOut.write(f">{refContig}\n{str(inputSeq)}\n")
                # Ragtag scaffold multiple input contigs
                else:
                    # Create a working directory for the files
                    chrDir = os.path.join(hapDir, refContig)
                    os.makedirs(chrDir, exist_ok=True)
                    
                    # Write the input sequences to a file
                    rawSequencesFile = os.path.join(chrDir, f"{refContig}.input.fasta")
                    with open(rawSequencesFile, "w") as fileOut:
                        for inputContig in inputContigs:
                            inputSeq = inputSequences[inputContig]
                            fileOut.write(f">{inputContig}\n{str(inputSeq)}\n")
                    
                    # Write the reference sequence to a file
                    refSequenceFile = os.path.join(chrDir, f"{refContig}.ref.fasta")
                    with open(refSequenceFile, "w") as fileOut:
                        refSeq = refSequences[refContig]
                        fileOut.write(f">{refContig}\n{str(refSeq)}\n")
                    
                    # Run ragtag to scaffold the sequences
                    run_ragtag(rawSequencesFile, refSequenceFile, os.path.join(chrDir, "ragtag_output"),
                               args.ragtag, threads=1)
                    
                    # Rewrite the output file to the haplotype directory
                    ragtagOutputFile = os.path.join(chrDir, "ragtag_output", f"ragtag.scaffold.fasta")
                    seqID, numIDs = find_longest_seq(ragtagOutputFile) # detect if ragtag couldn't scaffold them all together
                    if numIDs > 1:
                        print(f"{numIDs} IDs were found in the ragtag output '{ragtagOutputFile}'; " + 
                              f"we will just use the longest sequence '{seqID}'")
                    
                    with open(ragtagOutputFile, "r") as fileIn, open(outputFileName, "w") as fileOut:
                        records = SeqIO.parse(fileIn, "fasta")
                        for record in records:
                            if record.id == seqID:
                                fileOut.write(f">{refContig}\n")
                                fileOut.write(str(record.seq) + "\n")
                
                # Create a flag to indicate that the file was created
                open(outputFlagName, "w").close()
            else:
                print(f"contig '{refContig}' for '{os.path.basename(hapDir)}' already exists; skipping.")
    
    # Generate final output files
    for hapIndex, hapDir in enumerate([hap1Dir, hap2Dir]):
        combinedOutputFile = os.path.join(args.outputDirectory, f"hap{hapIndex+1}.fasta")
        combinedOutputFlag = os.path.join(args.outputDirectory, f"hap{hapIndex+1}.is.ok.flag")
        
        if not os.path.exists(combinedOutputFlag):
            # Generate a combined output file
            with open(combinedOutputFile, "w") as fileOut:
                for contig in refContigs:
                    contigFile = os.path.join(hapDir, f"{contig}.fasta") # file is assumed to still exist
                    with open(contigFile, "r") as fileIn:
                        records = SeqIO.parse(fileIn, "fasta")
                        for record in records:
                            # Write multiline output
                            fileOut.write(f">{contig}\n")
                            sequence = str(record.seq)
                            fileOut.write("\n".join(
                                [
                                    sequence[i:i+MULTILINE_LENGTH]
                                    for i in range(0, len(sequence), MULTILINE_LENGTH)
                                ]
                            ) + "\n")
            # Create a flag to indicate that the file was created
            open(combinedOutputFlag, "w").close()
        else:
            print(f"final output '{combinedOutputFile}' already exists; skipping.")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

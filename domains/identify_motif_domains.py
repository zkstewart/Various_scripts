#! python3

import os, sys, subprocess

sys.path.append(r"C:\git\Various_scripts\Function_packages")
import ZS_Utility, ZS_BlastIO, ZS_SeqIO, ZS_AlignIO, ZS_HmmIO

def execute_through_wsl(cmd):
    '''
    Parameters:
        cmd -- a string of the command to be executed through WSL.
    '''
    # Explode the command into a list
    cmd = cmd.split(" ")
    
    # Make sure the first argument will be handled right
    cmd[0] = ZS_Utility.wsl_which(cmd[0])
    
    # Change any paths to WSL paths
    cmd = [
        c if "\\" not in c
        else ZS_Utility.convert_to_wsl_if_not_unix(c)
        for c in cmd
    ]
    
    # Format for execution through WSL
    cmd = ["wsl", "~", "-e", *cmd]
    
    # Run the command
    run_cmd = subprocess.Popen(cmd, shell = True,
                               stdout = subprocess.PIPE,
                               stderr = subprocess.PIPE)
    out, err = run_cmd.communicate()

# Specify input files and locations
motifFile = os.path.abspath("motif1.fasta")
dbDir = os.path.abspath("databases")

# Specify the program locations
makeblastdb = "/home/stewarz2/anaconda3/bin/makeblastdb"
blastp = "/home/stewarz2/anaconda3/bin/blastp"
msaCurate = r"C:\git\Various_scripts\Fasta_related\Alignment\msa_curate.py"
hmmerDir = "/mnt/c/bio/hmmer-3.4/bin"

################################################################################
# STEP 1: BLAST motif sequence to database sequences

# Set up output directory
blastResultsDir = os.path.abspath("blast_results")
os.makedirs(blastResultsDir, exist_ok=True)

dbFiles = [ os.path.join(dbDir, f) for f in os.listdir(dbDir) if f.endswith(".aa") ]
for dbFile in dbFiles:
    prefix = os.path.basename(dbFile).split(".")[0]
    
    # Create the blast database if it doesn't exist
    if not os.path.isfile(f"{dbFile}.phr"):
        execute_through_wsl(f"{makeblastdb} -in {dbFile} -dbtype prot -out {dbFile}")
    
    # Run the blastp command
    outputFileName = os.path.join(blastResultsDir, f"{prefix}.outfmt6")
    if not os.path.isfile(outputFileName):
        execute_through_wsl(f"{blastp} -query {motifFile} -db {dbFile} -outfmt 6 -out {outputFileName}")

################################################################################
# STEP 2: Parse BLAST files for ALL hits
blastDict = {}
resultsFiles = [ os.path.join(blastResultsDir, f) for f in os.listdir(blastResultsDir) if f.endswith(".outfmt6") ]
for resultFile in resultsFiles:
    prefix = os.path.basename(resultFile).split(".")[0]
    
    # Parse the file
    resultsParser = ZS_BlastIO.BLAST_Results(resultFile)
    resultsParser.evalue = 100
    resultsParser.parse()
    
    # Store data
    if resultsParser.results != {}:
        blastKey = list(resultsParser.results.keys())[0]
        blastDict[prefix] = resultsParser.results[blastKey]

################################################################################
# STEP 3: Build a single MSA from all the hits

# Set up output directory
fullDir = os.path.abspath("full_msa")
os.makedirs(fullDir, exist_ok=True)

# Figure out what our output file will be called
outputPrefix = os.path.basename(motifFile).split(".")[0]
fullMsaFile = os.path.join(fullDir, f"{outputPrefix}.full_msa.fasta")

# Do step 3 if the output file does not exist
if not os.path.exists(fullMsaFile):
    # Set up a FASTA object to store the sequences
    hitsFasta = ZS_SeqIO.FASTA(motifFile)
    
    # Iterate through the database files and grab the FULL sequences
    for dbFile in dbFiles:
        prefix = os.path.basename(dbFile).split(".")[0]
        
        # Get the hit data
        if prefix in blastDict:
            hitData = blastDict[prefix]
            
            # Parse the FASTA file
            records = ZS_SeqIO.FASTA(dbFile)
            
            # Extract the aligned regions and add to the output FASTA object
            for hitID, identity, qstart, qend, hitstart, hitend, evalue in hitData:
                seq = records[hitID].seq[hitstart-1:hitend]
                hitsFasta.add(ZS_SeqIO.FastASeq(hitID, seq))
    
    # Align the sequences with MAFFT
    tmpFile = os.path.join(fullDir, "tmp.fasta")
    hitsFasta.write(tmpFile)
    
    aligner = ZS_AlignIO.MAFFT(ZS_Utility.wsl_which("mafft"), algorithm="einsi", maxiterate=5)
    alignedFasta = aligner.align(tmpFile)
    
    os.unlink(tmpFile)
    
    # Write the aligned FASTA object to file
    alignedFasta.write(fullMsaFile, asAligned=True)
# Otherwise just load the file in
else:
    alignedFasta = ZS_SeqIO.FASTA(fullMsaFile, isAligned=True)

################################################################################
# STEP 4: Produce a HMM profile from the aligned sequences

# Figure out what our output file will be called
hmmFile = os.path.join(fullDir, f"{outputPrefix}.hmm")

# Specify the edited MSA file
"Done manually"
editedFile = os.path.join(fullDir, f"{outputPrefix}.full_msa.edit.fasta")

# Do step 4 if the output file does not exist
if not os.path.exists(hmmFile):
    # Create the HMM profile
    hmm = ZS_HmmIO.HMM(hmmerDir)
    hmm.create(editedFile, hmmFile, hmmName=outputPrefix)

################################################################################
# STEP 5: Query the HMM profile against all database sequences
EVALUE = 1e-10

# Set up output directory
hmmsearchDir = os.path.abspath("full_hmmsearch_results")
os.makedirs(hmmsearchDir, exist_ok=True)

# Iterate through the database files and query the HMM profile
hmmDict = {}
for dbFile in dbFiles:
    prefix = os.path.basename(dbFile).split(".")[0]
    
    # Figure out what our output file will be called
    outputFileName = os.path.join(hmmsearchDir, f"{prefix}.hmmsearch")
    
    # Run the query if the output file does not exist
    if not os.path.isfile(outputFileName):
        hmmer = ZS_HmmIO.HMMER(hmmerDir, hmmFile, threads=16)
        domDict = hmmer.run(dbFile, outputFileName)
    
    # Parse the file in with E-value cutoff
    domDict = ZS_HmmIO.hmmer_parse(outputFileName, EVALUE)
    
    hmmDict[prefix] = domDict

################################################################################
# STEP 6: Extract the aligned regions from the HMM hits

EXPECTED_LEN = 32 # Expected length of the motif, hardcoded in

# Figure out what our output file will be called
hmmHitMsaFile = os.path.join(fullDir, f"{outputPrefix}.final.fasta")

# Do step 6 if the output file does not exist
if not os.path.exists(hmmHitMsaFile):
    # Set up a FASTA object to store the aligned regions
    hitsFasta = ZS_SeqIO.FASTA(None)
    
    # Iterate through the database files and grab the aligned regions
    for dbFile in dbFiles:
        prefix = os.path.basename(dbFile).split(".")[0]
        
        # Get the hit data
        if prefix in hmmDict:
            hitData = hmmDict[prefix]
            
            # Parse the FASTA file
            records = ZS_SeqIO.FASTA(dbFile)
            
            # Extract the aligned regions and add to the output FASTA object
            for hitID, hitList in hitData.items():
                for domainID, start, end, evalue in hitList:
                    hitLen = end - start + 1
                    seq = records[hitID].seq[start-1:end + (EXPECTED_LEN - hitLen)] # Add the expected length to the end
                    hitsFasta.add(ZS_SeqIO.FastASeq(hitID, seq))
    
    # Align the sequences with MAFFT
    tmpFile = os.path.join(fullDir, "tmp.fasta")
    hitsFasta.write(tmpFile)
    
    aligner = ZS_AlignIO.MAFFT(ZS_Utility.wsl_which("mafft"), algorithm="einsi", maxiterate=5)
    alignedFasta = aligner.align(tmpFile)
    
    os.unlink(tmpFile)
    
    # Drop anything that ended up being too gappy
    alignedFasta = ZS_AlignIO.MSA.drop_gappy_seqs(alignedFasta, allowedGappiness=0.30, inPlace=True)
    
    # Write the aligned FASTA object to file
    alignedFasta.write(hmmHitMsaFile, asAligned=True)
# Otherwise just load the file in
else:
    alignedFasta = ZS_SeqIO.FASTA(hmmHitMsaFile, isAligned=True)

################################################################################
# STEP 7: Produce a final HMM profile from all of the aligned sequences

# Figure out what our output file will be called
finalHmmFile = os.path.join(fullDir, f"{outputPrefix}.final.hmm")

# Do step 7 if the output file does not exist
if not os.path.exists(finalHmmFile):
    # Create the HMM profile
    hmm = ZS_HmmIO.HMM(hmmerDir)
    hmm.create(hmmHitMsaFile, finalHmmFile, hmmName=outputPrefix)

################################################################################
# STEP 8: Produce a final HMM search result from all of the aligned sequences
EVALUE = 1e-10 # derived through manual inspection of HMMER results

# Set up output directory
finalHmmsearchDir = os.path.abspath("final_hmmsearch_results")
os.makedirs(finalHmmsearchDir, exist_ok=True)

# Iterate through the database files and query the HMM profile
finalHmmDict = {}
for dbFile in dbFiles:
    prefix = os.path.basename(dbFile).split(".")[0]
    
    # Figure out what our output file will be called
    outputFileName = os.path.join(finalHmmsearchDir, f"{prefix}.hmmsearch")
    
    # Run the query if the output file does not exist
    if not os.path.isfile(outputFileName):
        hmmer = ZS_HmmIO.HMMER(hmmerDir, finalHmmFile, threads=16)
        domDict = hmmer.run(dbFile, outputFileName)
    
    # Parse the file in with E-value cutoff
    domDict = ZS_HmmIO.hmmer_parse(outputFileName, EVALUE)
    
    finalHmmDict[prefix] = domDict

################################################################################
# STEP 9: Format a report of the final HMM search results

# Figure out what our output file will be called
finalReportFile = os.path.join(finalHmmsearchDir, f"{outputPrefix}.final.tsv")

# Do step 9 if the output file does not exist
if not os.path.exists(finalReportFile):
    # Open the file for writing
    with open(finalReportFile, "w") as reportFileOut:
        # Write the header
        reportFileOut.write("species\tsequence_id\tdomain_id\tstart\tend\tevalue\n")
        
        # Iterate through the database files and output any significant hits
        for dbFile in dbFiles:
            prefix = os.path.basename(dbFile).split(".")[0]
            
            # Get the hit data
            if prefix in finalHmmDict:
                hitData = finalHmmDict[prefix]
                
                # Parse the FASTA file
                records = ZS_SeqIO.FASTA(dbFile)
                
                # Write output FASTA
                with open(os.path.join(finalHmmsearchDir, f"{prefix}.hits.fasta"), "w") as fastaFileOut:
                    for hitID, hitList in hitData.items():
                        for domainID, start, end, evalue in hitList:
                            reportFileOut.write(f"{prefix}\t{hitID}\t{domainID}\t{start}\t{end}\t{evalue}\n")
                            fastaFileOut.write(str(records[hitID]) + "\n")

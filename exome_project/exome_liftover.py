#! python3
# exome_liftover.py
# Program to enable discovery of exon sequences
# from genome sequences on the basis of exome
# sequencing alignments

import sys, argparse, os
sys.path.append(os.path.dirname(os.path.dirname(__file__))) # 2 dirs up is where we find dependencies
from Function_packages import ZS_SeqIO, ZS_HmmIO

def validate_args(args):
    # Validate input data location
    if not os.path.isdir(args.alignmentsDir):
        print('I am unable to locate the directory where the alignments files are (' + args.alignmentsDir + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    for exe in ["hmmpress", "hmmbuild", "hmmsearch"]:
        if not os.path.isfile(os.path.join(args.hmmerDir, exe)) and not os.path.isfile(os.path.join(args.hmmerDir, exe + ".exe")):
            print("{0} does not exist at {1}".format(exe, args.hmmerDir))
            print('Make sure you\'ve typed the location correctly and try again.')
            quit()
    if not os.path.isfile(args.genomeFile):
        print('I am unable to locate the genome FASTA file (' + args.genomeFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate numeric arguments
    if args.Evalue < 0:
        print('Evalue should be greater than 0')
        quit()
    # Handle file output
    if os.path.isdir(args.outputDir):
        print("""Output directory already exists; note that this program will NOT
              overwrite existing files and instead will opt to resume any steps not
              completed. If files are corrupted or truncated expect unhandled errors""")
    else:
        try:
            os.mkdir(args.outputDir)
            print("Created '{0}' directory as part of argument validation".format(args.outputDir))
        except:
            print("Wasn't able to create '{0}' directory; does '{1}' actually exist?".format(args.outputDir, os.path.dirname(args.outputDir)))

def predict_exon_from_domDict():
    pass

if __name__ == "__main__":
    usage = """%(prog)s receives a directory full of aligned FASTA files as part of the
    Oz Mammals genome project. Its goal is to transform these alignments into HMMs that can
    then be queried against a genome of interest to locate the relevant exon sequence from
    said genome.
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-a", dest="alignmentsDir", required=True,
                help="Specify the directory where aligned FASTA files are located")
    p.add_argument("-hmm", dest="hmmerDir", required=True,
                help="Specify the directory where HMMER executables are located")
    p.add_argument("-g", dest="genomeFile", required=True,
                help="Specify the location of the genome FASTA file to find exons in")
    p.add_argument("-o", dest="outputDir", required=True,
                help="Output directory location (working and final files go here)")
    # Opts
    p.add_argument("-e", dest="Evalue", required=False, type=float,
                help="Optionally, specify the E-value cut-off for HMMER results (default==1e-20)",
                default=1e-20)
    
    args = p.parse_args()
    validate_args(args)
    
    # Locate all files
    files = [os.path.join(args.alignmentsDir, file) for file in os.listdir(args.alignmentsDir)]

    # Create HMMs from aligned files
    hmmsDir = os.path.join(args.outputDir, "hmms")
    hmmsList = []
    os.makedirs(hmmsDir, exist_ok=True)
    for f in files:
        hmmName = os.path.join(hmmsDir, os.path.basename(f).rsplit(".", maxsplit=1)[0] + ".hmm")
        hmm = ZS_HmmIO.HMM(args.hmmerDir)
        
        # If HMM exists, load it in
        if os.path.isfile(hmmName):
            hmm.load_HMM_file(hmmName)
        # Create HMM if it doesn't exist
        else:
            hmm.load_FASTA_from_file(f)
            hmm.create_HMM(hmmName, hmmBuildExtraArgs="--dna")
        hmmsList.append(hmm)
    
    # Use our HMMs to query the genome for possible exon hits
    domtbloutsDir = os.path.join(args.outputDir, "domtblouts")
    os.makedirs(domtbloutsDir, exist_ok=True)
    for hmm in hmmsList:
        # Derive domtblout name
        domtbloutName = os.path.join(
            domtbloutsDir,
            "{0}.domtblout".format(os.path.basename(hmm.hmmFile).rsplit(".", maxsplit=1)[0])
        )
        
        # If domtblout doesn't exist, run HMMER
        if not os.path.isfile(domtbloutName):
            hmmer = ZS_HmmIO.HMMER(hmm)
            hmmer.load_FASTA_from_file(args.genomeFile)
            hmmer.set_output_name(domtbloutName)
            hmmer.set_Evalue(args.Evalue)
            hmmer.run_search()
            domDict = hmmer.domDict # This is what we want out of HMMER
        # If it does exist, simply load it in
        else:
            domDict = ZS_HmmIO.hmmer_parse(domtbloutName, args.Evalue)

        # Perform liftover operation
        ## TBD, need to define a function above and use it with hmmer.domDict
        ## predict_exon_from_domDict(...)
    print("Program completed successfully!")

#! python3
# dge_pipe.py
# Pipeline script for performing differential gene expression analysis

import os, argparse, subprocess
from pathlib import Path

# Define functions
def validate_args(args):
    def fail_validation():
        print("Make sure you\'ve typed the file name or location correctly and try again.")
        quit()
    
    # Validate that all arguments are provided
    for key, value in vars(args).items():
        if value == None:
            print(key + " argument was not specified. Fix this and try again.")
            quit()
    
    # Validate input file locations
    if not os.path.isfile(args.genomeFile):
        print("I am unable to locate the genome FASTA file ({0})".format(args.genomeFile))
        fail_validation()
    if not os.path.isfile(args.annotationFile):
        print("I am unable to locate the GFF3 annotation file ({0})".format(args.annotationFile))
        fail_validation()
    if not os.path.isdir(args.rnaseqDirectory):
        print("I am unable to locate the directory containing RNAseq reads ({0})".format(
            args.rnaseqDirectory))
        fail_validation()
    for prefix in args.rnaseqPrefixes:
        f1 = Path(args.rnaseqDirectory, "{0}_1.fastq".format(prefix)).as_posix()
        f2 = Path(args.rnaseqDirectory, "{0}_2.fastq".format(prefix)).as_posix()
        if not os.path.isfile(f1):
            print("I am unable to locate one of the specified RNAseq files ({0})".format(f1))
            fail_validation()
        if not os.path.isfile(f2):
            print("I am unable to locate one of the specified RNAseq files ({0})".format(f2))
            fail_validation()
    if not os.path.isdir(args.trimmomaticDir):
        print("I am unable to locate the directory containing the Trimmomatic .jar file ({0})".format(
            args.trimmomaticDir))
        fail_validation()
    if not os.path.isfile(Path(args.trimmomaticDir, args.trimmomaticJar).as_posix()):
        print("I am unable to locate the Trimmomatic .jar file in the specified directory ({0})".format(
            Path(args.trimmomaticDir, args.trimmomaticJar).as_posix()))
        fail_validation()
    if not os.path.isdir(args.starDir):
        print("I am unable to locate the directory containing the STAR executable file ({0})".format(
            args.starDir))
        fail_validation()
    if not os.path.isfile(Path(args.starDir, args.starExe).as_posix()):
        print("I am unable to locate the STAR executable file in the specified directory ({0})".format(
            Path(args.starDir, args.starExe).as_posix()))
        fail_validation()
    if not os.path.isdir(args.python2Dir):
        print("I am unable to locate the directory containing the Python 2 executable file ({0})".format(
            args.python2Dir))
        fail_validation()
    if os.path.isfile(Path(args.python2Dir, "python").as_posix()):
        args.python2Exe = "python"
    elif os.path.isfile(Path(args.python2Dir, "python2.7").as_posix()):
        args.python2Exe = "python2.7"
    elif os.path.isfile(Path(args.python2Dir, "python2").as_posix()):
        args.python2Exe = "python2"
    else:
        print("I am unable to locate the Python 2 executable file in the specified directory ({0})".format(
            Path(args.starDir).as_posix()))
        print("I expect there to be a file called \"python\", \"python2.7\", or \"python2\"")
        fail_validation()
    
    

def setup_working_directory(baseDir, species, genomeFile, rnaseqDirectory, rnaPrefixes):
    # Symbolic link dirs
    genomeLocation = Path(baseDir, "genome").as_posix()
    os.mkdir(genomeLocation)
    os.symlink(genomeFile, Path(genomeLocation, "{0}.fasta".format(species)))

    rnaseqLocation = Path(baseDir, "rnaseq_reads").as_posix()
    os.mkdir(rnaseqLocation)
    for prefix in rnaPrefixes:
        os.symlink(Path(rnaseqDirectory, "{0}_1.fastq".format(prefix)).as_posix(),
                    Path(rnaseqLocation, "{0}_1.fastq".format(prefix)).as_posix())

    # Working dirs
    trimWorkDir = Path(baseDir, "trimmomatic")
    starWorkDir = Path(baseDir, "star_map")
    countWorkDir = Path(baseDir, "htseq_count")

    os.mkdir(trimWorkDir)
    os.mkdir(starWorkDir)
    os.mkdir(countWorkDir)

    return genomeLocation, rnaseqLocation, trimWorkDir, starWorkDir, countWorkDir

def qsub(scriptName):
    p = subprocess.Popen(["qsub", scriptName])
    stdout, stderr = p.communicate()
    return stdout

def generate_trim_script(scriptName, workDir, species, trimDir, trimJar, rnaseqLocation, rnaPrefixes, threads, walltime=80, mem="50G"):
    trimScript = r"""#!/bin/bash -l
    #PBS -N trim_{species}
    #PBS -l walltime={walltime}:00:00
    #PBS -l mem={mem}
    #PBS -l ncpus={threads}

    cd {workDir}

    ### MANUAL SETUP BELOW
    ## SETUP: Load modules
    module load java/1.8.0_92

    ## SETUP: Specify trimmomatic location
    TRIMDIR={trimDir}
    TRIMJAR={trimJar}

    ## SETUP: Specify file prefixes
    SPECIES={species}
    RNADIR={rnaseqLocation}

    FILEPREFIXES="{rnaPrefixes}"
    # Note: FILEPREFIXES assumes a file name like ${PREFIX}_1.fastq with paired ${PREFIX}_2.fastq
    ### MANUAL SETUP END

    ### AUTOMATIC SETUP START
    ## SETUP: Paired-end trimmomatic settings
    COMMAND="ILLUMINACLIP:{trimDir}/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25"
    ### AUTOMATIC SETUP END

    ### RUN PROGRAM
    ## STEP 1: Run Trimmomatic
    for file in $FILEPREFIXES; do java -jar $TRIMDIR/$TRIMJAR PE -threads {threads} -trimlog $SPECIES.logfile $RNADIR/$/{file/}_1.fastq $RNADIR/$/{file/}_2.fastq -baseout $/{file/}.trimmed.fq.gz $/{COMMAND/}; done

    ## STEP 2: Unzip files
    for file in $FILEPREFIXES; do gunzip $/{file/}.trimmed_1P.fq.gz $/{file/}.trimmed_2P.fq; done

    ## STEP 3: Combine files
    for file in $FILEPREFIXES; do cat $/{file/}.trimmed_1P.fq.gz >> $/{SPECIES/}_1P.fq; cat $/{file/}.trimmed_2P.fq.gz >> $/{SPECIES/}_2P.fq; done
    touch 
    """.format(species=species, walltime=walltime, mem=mem, threads=threads, workDir=workDir, trimDir=trimDir,
            trimJar=trimJar, rnaseqLocation=rnaseqLocation, rnaPrefixes=" ".join(rnaPrefixes))
    
    with open(scriptName, "w") as fileOut:
        fileOut.write(trimScript)

def generate_star_script(scriptName, species, starDir, genomeLocation, rnaseqLocation, threads, previousJob, walltime=80, mem="90G"):
    starScript = r"""#!/bin/bash -l
    #PBS -N star_{species}
    #PBS -l walltime={walltime}:00:00
    #PBS -l mem={mem}
    #PBS -l ncpus={threads}
    #PBS -W depend=afterok:{previousJob}

    cd $PBS_O_WORKDIR

    ## MANUAL SETUP BELOW
    # SETUP: Specify STAR location
    STARDIR={starDir}

    # SETUP: Specify prefix
    SPECIES={species}

    # SETUP: Specify computational resources
    CPUS={threads}
    ## MANUAL SETUP END

    ## AUTO SETUP BELOW
    # SETUP: Genome location
    GENDIR={genomeLocation}
    GENFILE={species}.fasta

    # SETUP: Trimmed reads location
    RNADIR={rnaseqLocation}
    RNAFILE1={species}_1P.fq
    RNAFILE2={species}_2P.fq
    ## AUTO SETUP END

    ## RUN PROGRAM
    # STEP 1: Copy genome here. Need to do this since STAR can only tolerate 1 index per directory...
    cp $GENDIR/$GENFILE .

    # STEP 2: Generate index
    $STARDIR/source/STAR --runThreadN $CPUS --runMode genomeGenerate --genomeDir $PBS_O_WORKDIR --genomeFastaFiles $GENFILE

    # STEP 3: Run 2-pass procedure
    $STARDIR/source/STAR --runThreadN $CPUS --genomeDir $PBS_O_WORKDIR --readFilesIn $RNAFILE1 $RNAFILE2 --twopassMode Basic
    """.format(species=species, walltime=walltime, mem=mem, threads=threads, starDir=starDir,
            genomeLocation=genomeLocation, rnaseqLocation=rnaseqLocation, previousJob=previousJob)

    with open(scriptName, "w") as fileOut:
        fileOut.write(starScript)

def generate_htseq_script(scriptName, species, py2Dir, annotationFile, samFile, previousJob, walltime=10, mem="50G"):
    htseqScript = r"""#!/bin/bash -l
    #PBS -N htseq_{species}
    #PBS -l walltime={walltime}:00:00
    #PBS -l mem={mem}
    #PBS -l ncpus=1
    #PBS -W depend=afterok:{previousJob}

    cd $PBS_O_WORKDIR

    ## MANUAL SETUP BELOW
    # SETUP: Specify mapped reads SAM file
    SAMFILE={samFile}

    # SETUP: Specify Python2 location
    PY2DIR={py2Dir}

    # SETUP: GFF annotation file
    GFF={annotationFile}
    ## MANUAL SETUP END

    ## AUTO SETUP BELOW
    # SETUP: Output file name
    OUT={species}.htseq.counts
    ## AUTO SETUP END

    ## RUN PROGRAM
    $PY2DIR/python2.7 -m HTSeq.scripts.count -r name -a 0 -t gene -i ID $SAMFILE $GFF > $OUT
    """.format(species=species, walltime=walltime, mem=mem, py2Dir=py2Dir, samFile=samFile,
            annotationFile=annotationFile, previousJob=previousJob)

    with open(scriptName, "w") as fileOut:
        fileOut.write(htseqScript)

# Argument input
usage = """%(prog)s does things.
It will establish a working environment within the directory where this script is located. As such, you should
perform a single DGE analysis within the present directory.

Also note that skipping steps should only be done if they were done previously by this script, or if you have set
up the working directory to emulate how it would be laid out were this script run previously.
"""

p = argparse.ArgumentParser(description=usage)
p.add_argument("-g", dest="genomeFile",
                help="Specify the location and file name of the genome file")
p.add_argument("-a", dest="annotationFile",
                help="Specify the location and file name of the gene annotation GFF3 file")
p.add_argument("-td", dest="rnaseqDirectory",
                help="Specify the location of the RNAseq reads file(s)")
p.add_argument("-r", dest="rnaseqPrefixes", nargs="+",
                help="""Specify one or more prefixes of the RNAseq FASTQ file(s).
                File(s) should be named like \{PREFIX\}_1.fastq and \{PREFIX\}_2.fastq""")
p.add_argument("-sp", dest="species",
                help="Provide a short, descriptive name for the species being analysed")
p.add_argument("-td", dest="trimmomaticDir",
                help="Specify the location of the Trimmomatic jar file")
p.add_argument("-te", dest="trimmomaticJar",
                help="Optionally specify the name of the Trimmomatic jar file",
                default = "trimmomatic-0.36.jar")
p.add_argument("-sd", dest="starDir",
                help="Specify the location of the STAR executable file")
p.add_argument("-se", dest="starExe",
                help="Optionally specify the name of the STAR executable file",
                default = "STAR")
p.add_argument("-p2", dest="python2Dir",
                help="Specify the location of the Python 2 directory which has HTSeq installed",
                default = "STAR")
p.add_argument("-skip_trim", dest="skip_trim", action="store_true",
                help="Optionally skip trimming if already performed")
p.add_argument("-skip_map", dest="skip_map", action="store_true",
                help="Optionally skip mapping if already performed")
p.add_argument("-skip_count", dest="skip_count", action="store_true",
                help="Optionally skip counting if already performed")
args = p.parse_args()
#validate_args(args)

## Hard-coded testing
class Arg:
    def __init__(self):
        self.genomeFile = r""
        self.annotationFile = r""
        self.rnaseqDirectory = r""
        self.rnaseqPrefixes = []
        self.species = ""
        self.trimmomaticDir = r""
        self.trimmomaticJar = r""
        self.starDir = r""
        self.python2Dir = r""
        self.skip_trim = False
        self.skip_map = False
        self.skip_count = False

args = Arg()

# Setup directory
baseDir = os.getcwd()
genomeLocation, rnaseqLocation, trimWorkDir, starWorkDir, countWorkDir = setup_working_directory(baseDir,
        args.species, args.genomeFile, args.rnaseqDirectory, args.rnaPrefixes)

# Trim reads
if not args.skip_trim:
    trimScriptName = Path(trimWorkDir, "run_trim.sh").as_posix()
    generate_trim_script(trimScriptName, trimWorkDir, args.species, args.trimmomaticDir,
            args.trimmomaticJar, args.rnaseqDirectory, args.rnaPrefixes, args.cpus)
    trimJob = qsub(trimScriptName)
else:
    trimJob = ""

# Map reads
if not args.skip_map:
    starScriptName = Path(starWorkDir, "run_star.sh").as_posix()
    generate_star_script(starScriptName, args.species, args.starDir, genomeLocation,
            rnaseqLocation, args.cpus, trimJob)
    starJob = qsub(starScriptName)
else:
    starJob = ""

# Count mapped reads
if not args.skip_count:
    htseqScriptName = Path(countWorkDir, "run_htseq.sh").as_posix()
    generate_htseq_script(htseqScriptName, args.species, args.python2Dir, args.annotationFile,
            Path(starWorkDir, "Aligned.out.sam").as_posix(), starJob)
    countJob = qsub(htseqScriptName)
else:
    countJob = ""

# DGE analysis

# DGE visualisation

# GOseq analysis

# GOseq visualisation

# Visualisation


#! python3
# Script to receive a SAM or BAM file and reduce
# its header SQ lines to only list contigs
# that actually have alignments in the file.

import os, argparse, shutil, subprocess

# Define functions
def validate_args(args):
    # Validate input file location
    if not os.path.isfile(args.bamFile):
        print(f'I am unable to locate the BAM file ({args.bamFile})')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate output file location
    if args.newOutputFile != None:
        if os.path.isfile(args.newOutputFile):
            print(f'A file already exists at "{args.newOutputFile}"')
            print('This program will not overwrite existing files; specify a different name and try again.')
            quit()

def validate_programs_locatable(programList):
    '''
    Raises an error if any program is not locatable.
    
    Parameters:
        programList -- a list containing one or more strings corresponding to
                       program executable names which are expected to be locatable
                       via the system PATH
    '''
    for program in programList:
        which = shutil.which(program)
        if which == None:
            raise FileNotFoundError(
                f"{program} not found in system PATH!"
            )

def tmp_file_name_gen(prefix, suffix):
    '''
    Hidden function for use by Class methods.
    Params:
        prefix -- a string for a file prefix e.g., "tmp"
        suffix -- a string for a file suffix e.g., "fasta". Note that we don't
                    use a "." in this, since it's inserted between prefix and suffix
                    automatically.
    Returns:
        tmpName -- a string for a file name which does not exist in the current dir.
    '''
    ongoingCount = 1
    while True:
        if not os.path.isfile("{0}.{1}".format(prefix, suffix)):
            return "{0}.{1}".format(prefix, suffix)
        elif os.path.isfile("{0}.{1}.{2}".format(prefix, ongoingCount, suffix)):
            ongoingCount += 1
        else:
            return "{0}.{1}.{2}".format(prefix, ongoingCount, suffix)

def get_aligned_contigs(amFile):
    '''
    Parameters:
        amFile -- a string indicating the location of a (S/B)AM file
    Returns:
        contigSet -- a set containing all contigs mentioned in the AM file
    '''
    # Run samtools pipe for getting contig IDs
    samtoolsProcess = subprocess.Popen(f"samtools view {amFile} | cut -f 3 | sort | uniq", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = samtoolsProcess.communicate()
    stdout, stderr = stdout.decode(), stderr.decode()
    if stderr != "":
        raise Exception(f"samtools view contigs died with stderr == {stderr}")
    
    # Get and return contig set
    contigSet = set(stdout.rstrip("\n").split("\n"))
    return contigSet

def write_reheader_file(amFile, contigSet, reheaderFile):
    '''
    Parameters:
        amFile -- a string indicating the location of a (S/B)AM file
        contigSet -- a set containing all contigs mentioned in the AM file
        reheaderFile -- a string indicating the file name to write the
                        file that can be used for reheader-ing
    '''
    # Run samtools to get header data
    samtoolsProcess = subprocess.Popen(f"samtools view -H {amFile}", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = samtoolsProcess.communicate()
    stdout, stderr = stdout.decode(), stderr.decode()
    if stderr != "":
        raise Exception(f"samtools view header died with stderr == {stderr}")
    
    # Extract header lines, filtered by contigSet
    foundContigs = set()
    headerLines = []
    for line in stdout.rstrip("\n").split("\n"):
        # Handle @SQ lines
        if line.startswith("@SQ"):
            # Extract contig
            sqTag, snTag, lnTag = line.split("\t")
            contig = snTag.split(":", maxsplit=1)[1] # removes SN: prefix
            
            # Hold onto line if it's in our contigSet
            if contig in contigSet:
                headerLines.append(line)
                foundContigs.add(contig)
        
        # Pass every other line
        else:
            headerLines.append(line)
    
    # End program if something is desperately wrong
    if len(foundContigs) == 0:
        print("ERROR: Somehow, we found NO @SQ lines?!?")
        
        if len(contigSet) == 0:
            print("This is likely because contigSet is empty.")
            print("In other words, your BAM file itself is probably empty!")
        else:
            print("Something must be very wrong with your file, and I don't know how.")
        
        print("Program will exit now to prevent erroneous behaviour.")
        quit()
    
    # Raise a warning if something looks amiss
    if foundContigs != contigSet:
        print("WARNING: Reheader file might be missing @SQ lines!")
        print("Lines that are missing are:")
        print("\n".join([ contigSet.difference(foundContigs) ]))
        print("\n".join([ foundContigs.difference(contigSet) ]))
    
    # Write output file
    with open(reheaderFile, "w") as fileOut:
        fileOut.write("{0}\n".format("\n".join(
            headerLines
        )))

def reheader_bam_file(amFile, reheaderFile, outputFile):
    '''
    Parameters:
        amFile -- a string indicating the location of a BAM file
        reheaderFile -- a string indicating the location of the SAM file
                        containing a header for reheadering
        outputFile -- a string indicating the name for the output
    '''
    # Get a temporary file name for reheadering
    tmpFile = tmp_file_name_gen("tmp_samtools", "sam")
    
    # Run samtools reheader
    samtoolsProcess = subprocess.Popen(f"samtools reheader -P {reheaderFile} {amFile} > {tmpFile}", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = samtoolsProcess.communicate()
    stdout, stderr = stdout.decode(), stderr.decode()
    if stderr != "":
        raise Exception(f"samtools reheader died with stderr == {stderr}")
    
    # Move file to output location
    shutil.move(tmpFile, outputFile)

def main():
    # User input
    usage = """%(prog)s will reduce the header of a BAM file to just
    list @SQ lines that are actually present inside the file. It will do
    this in-place by default, but it can optionally write a new output file
    instead if desired.
    """
    p = argparse.ArgumentParser(description=usage)
    ## Required
    p.add_argument("-i", dest="bamFile",
                   required=True,
                   help="Input SAM/BAM file")
    ## Optional
    p.add_argument("--writeNewFile", dest="newOutputFile",
                   required=False,
                   help="""Optionally, specify a file name here if your intent is
                   not to reheader the file in-place, but instead to write a new
                   output file""")
    args = p.parse_args()
    validate_args(args)
    
    # Validate that PATH programs are locatable
    validate_programs_locatable(["samtools"])
    
    # Parse the file to get the IDs of contigs mentioned in alignments
    contigSet = get_aligned_contigs(args.bamFile)
    
    # Create temp file for reheader
    reheaderFile = tmp_file_name_gen("tmp_reheader", "sam")
    write_reheader_file(args.bamFile, contigSet, reheaderFile)
    
    # Reheader in-place if desired
    if args.newOutputFile == None:
        reheader_bam_file(args.bamFile, reheaderFile, args.bamFile)
    # Write new file otherwise
    else:
        reheader_bam_file(args.bamFile, reheaderFile, args.newOutputFile)
    
    # Clean up temp file
    os.unlink(reheaderFile)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

#! python3
# merge_invariant_sites.py
# A script to receive a variant-sites only VCF
# and a VCF containing invariant sites, and merges
# the invariant site predictions into the first VCF.
# This will allow for nucleotide diversity as from pixy
# to be calculated effectively.

import os, argparse, gzip
from contextlib import contextmanager, ExitStack

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.variantsVCF):
        print('I am unable to locate the variant calls VCF file (' + args.variantsVCF + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.invariantsVCF):
        print('I am unable to locate the invariant calls VCF file (' + args.invariantsVCF + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate numeric arguments
    if args.qualCutoff < 1:
        print('QUAL value cutoff must be a positive integer or float.')
        quit()
    # Validate output file locations
    if os.path.exists(args.tmpDir):
        print('Temporary dir location already exists (' + args.outputFileName + ')')
        print("I don't want to have to deal with checking for existing files and all that.")
        print("So, I'm going to ask you to specify a different temporary directory. Sorry and thanks.")
        quit()
    else:
        try:
            os.mkdir(args.tmpDir)
            print("Temporary directory created successfully at " + args.tmpDir)
        except:
            print('I am unable to create the temporary directory (' + args.tmpDir + ')')
            print('Make sure you specify its location as being in an existing location and try again.')
    if os.path.exists(args.outputFileName):
        print('File already exists at output location (' + args.outputFileName + ')')
        print('Make sure you specify a unique file name and try again.')
        quit()

@contextmanager
def open_vcf_file(filename):
    if filename.endswith(".gz"):
        with gzip.open(filename, "rt") as f:
            yield f
    else:
        with open(filename) as f:
            yield f

def parse_vcf_line(line):
    sl = line.rstrip("\r\n\t ").split("\t")
    return *sl[0:9], sl[9:]

## Main
def main():
    # User input
    usage = """%(prog)s reads in the primary VCF containing only variant calls,
    and a second VCF containing invariant sites. The invariant sites are then
    merged into the primary VCF, and the output VCF is written to the specified
    output file location.
    """
    p = argparse.ArgumentParser(description=usage)
    # Req
    p.add_argument("-i1", dest="variantsVCF",
                   required=True,
                   help="""Input the first VCF file containing variant calls that
                   you want variant sites to be merged into""")
    p.add_argument("-i2", dest="invariantsVCF",
                   required=True,
                   help="""Input the second VCF file containing invariant calls that
                   you want to merge into -i1""")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file name for the filtered SNPs")
    # Opts
    p.add_argument("--tmp", dest="tmpDir",
                   required=False,
                   help="""Optionally, specify a temporary directory to store
                   intermediate files; default == invariants_tmp""",
                   default="invariants_tmp")
    p.add_argument("--qual", dest="qualCutoff",
                   required=False,
                   type=float,
                   help="""Optionally, specify the QUAL value to filter SNPs by;
                   default == 30.0""",
                   default=30.0)
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse through invariant VCF file and write to temporary files
    with open_vcf_file(args.invariantsVCF) as fileIn, ExitStack() as stack:
        handles = {}
        for line in fileIn:
            if not line.startswith("#"):
                # Extract line details
                contig, pos, varid, ref, alt, qual, filt, info, fmt, samples = parse_vcf_line(line)
                
                # Skip if not invariant or below quality cutoff
                if alt != "." or qual == "." or (float(qual) < args.qualCutoff):
                    continue
                
                # Set up file handle for each contig
                if contig not in handles:
                    handles[contig] = stack.enter_context(open(os.path.join(args.tmpDir, f"{contig}.vcf"), "w"))
                
                # Write to file
                handles[contig].write(line)
    
    # Merge files together
    with open(args.variantsVCF, "r") as fileIn, open(args.outputFileName, "w") as fileOut, ExitStack() as stack:
        handles = {}
        lastContig = None
        for line in fileIn:
            # Write out primary file's header lines
            if line.startswith("#"):
                fileOut.write(line)
                continue
            
            # Extract details from current line
            contig, pos, varid, ref, alt, qual, filt, info, fmt, samples = parse_vcf_line(line)
            
            # Finish writing invariant file lines if we've moved onto a new contig
            if lastContig != None and contig != lastContig:
                firstIter = True
                while True:
                    # Get the current line from the invariant file
                    if firstIter:
                        invariantLine = handles[f"{lastContig}_last"]
                        firstIter = False
                    else:
                        invariantLine = handles[lastContig].readline()
                    
                    # Break if we've reached the end of the file
                    if invariantLine == "":
                        break
                    
                    # Otherwise, write the line
                    fileOut.write(invariantLine)
            lastContig = contig
            
            # Iterate through invariant file until current position is skipped
            firstIter = True
            while True:
                # Set up file handle for each contig
                if contig not in handles:
                    handles[contig] = stack.enter_context(open(os.path.join(args.tmpDir, f"{contig}.vcf"), "r"))
                
                # Get the current line from the invariant file
                if firstIter:
                    invariantLine = handles[f"{contig}_last"] if f"{contig}_last" in handles else handles[contig].readline()
                    firstIter = False
                else:
                    invariantLine = handles[contig].readline()
                handles[f"{contig}_last"] = invariantLine # set for 'finish writing invariant ...' block
                
                # Break if we've reached the end of the file
                if invariantLine == "":
                    break
                
                # Extract details from the invariant line
                invariantValues = parse_vcf_line(invariantLine)
                
                # Write the line if the position precedes the current variant position
                if int(invariantValues[1]) < int(pos):
                    fileOut.write(invariantLine)
                
                # Skip if the position is the same as the current variant position
                elif int(invariantValues[1]) == int(pos):
                    continue
                
                # If the position is greater than the current variant position, break
                else:
                    break
            
            # Write current position
            fileOut.write(line)
        
        # Finish writing invariant file lines if we've ended the primary VCF
        firstIter = True
        while True:
            # Get the current line from the invariant file
            if firstIter:
                invariantLine = handles[f"{lastContig}_last"]
                firstIter = False
            else:
                invariantLine = handles[lastContig].readline()
            
            # Break if we've reached the end of the file
            if invariantLine == "":
                break
            
            # Otherwise, write the line
            fileOut.write(invariantLine)
    
    # Clean up temporary files
    try:
        os.unlink(args.tmpDir)
    except PermissionError:
        print(f"I'm unable to delete the temporary directory at '{args.tmpDir}'. " +
              "Please delete it manually.")
    
    print("Program finished successfully!")

if __name__ == "__main__":
    main()

#! python3
# msa_extract_seqs.py
# Program to extract all sequences from a directory of MSA FASTA files into
# a single output FASTA. In the future, criteria might be provided to determine
# which sequences are extracted.

# Import external packages
import os, argparse
from Bio import SeqIO

# Define functions for later use
## Validate arguments
def validate_args(args):
        # Validate that necessary arguments have been provided
        if args.inputLocation == None:
                print('-i argument must be provided; fix your input and try again.')
                quit()
        # Validate the FASTA file input location depending on type of input
        if len(args.inputLocation) == 1:
                # Check that specified value is a path
                if not os.path.isdir(args.inputLocation[0]):
                        print('One value was provided for -i, which means you should have provided a directory containing FASTA files.')
                        print('The provided value "' + args.inputLocation[0] + '" is not a directory; either it does not exist or it is a file; fix your input and try again.')
                        quit()
        else:
                # Check that the specified values are files
                for file in args.inputLocation:
                        if not os.path.isfile(file):
                                print('Multiple values were provided for -i, which means you should have provided individual FASTA files.')
                                print('The provided value "' + file + '" is not a file; either it does not exist or it is a directory; fix your input and try again.')
                                quit()
        # Reformat prefixes for use
        if args.prefixes != None:
                args.prefixes = tuple(args.prefixes)
        # Ensure that the output location is sensible
        outputDir = os.path.dirname(args.outputFileName)
        if outputDir != '':
                outputDir = os.path.abspath(outputDir)
                if not os.path.isdir(outputDir):
                        print('The output directory "' + outputDir + '" does not exist; create this location before specifying this directory as the output file location.')
                        quit()
        # Prevent file overwrites
        if os.path.isfile(args.outputFileName):
                print('There is already a file named "' + args.outputFileName + '". Move/delete/rename this file and try again.')
                quit()
        return args

def main():
    #### USER INPUT SECTION
    usage = """%(prog)s will identify FASTA files provided by the user individually
    or within a specified directory and extract all sequences from these files, remove
    '-' characters as found in MSA FASTAs, and produce a single concatenated FASTA file.
    While intended for MSAs, this can also be used for normal FASTA files. Sequences with
    prefix codes allow this program to have different behaviour if they are detected;
    retrieve_all means if we find a single sequence with one of the specified prefixes we
    grab all sequences, remove_all is the inverse; retrieve_only means we'll only grab sequences
    with the prefix, and remove_only is the inverse.
    """

    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", "-input", dest="inputLocation", nargs="+",
                help="Input a single directory or multiple file locations")
    p.add_argument("-o", "-output", dest="outputFileName",
                help="Specify the file to write output sequences in FASTA to")
    # Opts
    p.add_argument("-p", "-prefix", dest="prefixes", nargs="+",
                help="Optionally specify one or more prefixes present in your sequence IDs to look for; -b will control what happens when we find one")
    p.add_argument("-b", "-behaviour", dest="behaviour", choices=['retrieve_all', 'remove_all', 'retrieve_only', 'remove_only'], default = 'retrieve_all',
                help="Optionally control what happens when we find a sequence with prefix within a FASTA file; this isn't necessary if you don't specify any prefixes")

    args = p.parse_args()
    args = validate_args(args)

    # Find FASTA files depending on how inputLocation was specified
    msaFileNameList = []
    if len(args.inputLocation) == 1:
            # Scan through files and detect our files of interest
            for file in os.listdir(args.inputLocation[0]):
                    file = os.path.join(args.inputLocation[0], file)
                    if not os.path.isfile(file):
                            continue
                    with open(file, 'r') as fileIn:
                            for line in fileIn:
                                    if line.startswith('>'):
                                            msaFileNameList.append(os.path.abspath(file))
                                    break
            # Ensure that we found some input files
            if msaFileNameList == []:
                    print('I did not find any FASTA files in the provided input directory "' + args.inputLocation[0] + '"')
                    print('Make sure you specified the correct location, or make sure all files are located at this directory. Program will exit now.')
                    quit()
    else:
            # Make sure that the provided files all exist
            for file in args.inputLocation:
                    if not os.path.isfile(file):
                            print('Input file "' + file + '" either does not exist or is not a file.')
                            print('Make sure you spelled the file name/location correctly and try again.')
                            quit()
                    with open(file, 'r') as fileIn:
                            for line in fileIn:
                                    if not line.startswith('>'):
                                            print('Input file "' + file + '" does not appear to be FASTA formatted i.e., it lacks the ">" character at its start.')
                                            print('Make sure you spelled the correct file or fix this file and try again.')
                                            quit()
                                    break
                    msaFileNameList.append(os.path.abspath(file))

    # Read through files, extract sequences, and format an output file
    with open(args.outputFileName, 'w') as fileOut:
            for fasta in msaFileNameList:
                    # Scan for prefixes if relevant
                    if args.prefixes != None:
                            foundIDs = []
                            records = SeqIO.parse(open(fasta, 'r'), 'fasta')
                            for record in records:
                                    if record.description.startswith(args.prefixes):
                                            foundIDs.append(record.description)
                    # Extract output depending on behaviour, if relevant
                    records = SeqIO.parse(open(fasta, 'r'), 'fasta')
                    for record in records:
                            if args.prefixes == None:
                                    fileOut.write('>' + record.description + '\n' + str(record.seq).replace('-', '') + '\n')
                            else:
                                    if args.behaviour == 'retrieve_all' and foundIDs != []:
                                            fileOut.write('>' + record.description + '\n' + str(record.seq).replace('-', '') + '\n')
                                    elif args.behaviour == 'remove_all' and foundIDs == []:
                                            fileOut.write('>' + record.description + '\n' + str(record.seq).replace('-', '') + '\n')
                                    elif args.behaviour == 'retrieve_only' and record.description in foundIDs:
                                            fileOut.write('>' + record.description + '\n' + str(record.seq).replace('-', '') + '\n')
                                    elif args.behaviour == 'remove_only' and record.description not in foundIDs:
                                            fileOut.write('>' + record.description + '\n' + str(record.seq).replace('-', '') + '\n')

    # All done!
    print('Program completed successfully!')

if __name__ == "__main__":
    main()

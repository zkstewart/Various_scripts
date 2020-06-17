#! python3
# run_mmseqs2.py
# Wrapper script to make it easier to run a search with MMseqs2

import argparse, os
from itertools import groupby
from Bio import SeqIO
from collections import OrderedDict

# Various functions for program operations
def validate_args(args):
        # Validate that relevant arguments were provided
        for key, value in vars(args).items():
                if key != "recomputeFasta":
                        if value == None:
                                print('"' + key + '" arg was not specified.')
                                quit()
                elif args.sortMethod == "recompute":
                        if value == None:
                                print('"' + key + '" arg was not specified (and sorting method is "recompute").')
                                quit()
        # Validate input file locations
        if not os.path.isfile(args.inputFile):
                print('I am unable to locate the input file (' + args.inputFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if args.sortMethod == "recompute":
                if not os.path.isfile(args.recomputeFasta):
                        print('I am unable to locate the recompute FASTA file (' + args.recomputeFasta + ')')
                        print('Make sure you\'ve typed the file name or location correctly and try again.')
                        quit()
        # Validate output file location
        if os.path.exists(args.outputFile):
                print('The specified output file name already exists (' + args.outputFile + ')')
                print('Specify a different name or move/delete/rename the existing file and try again.')
                quit()

def mms2sort_evalue(inputFile, outputFile, delInputFile=False):
        grouper = lambda x: x.split('\t')[0]
        with open(inputFile, 'r') as fileIn, open(outputFile, 'w') as fileOut:
                for key, group in groupby(fileIn, grouper):
                        group = list(group)
                        # Skip empty groups if input file is poorly formatted
                        """Note that this is a self insult. I don't know why mms2sort_recompute
                        adds in blank lines. It makes no sense to my current mind (tired) so...
                        """
                        if key.rstrip("\r\n") == "":
                                continue
                        # Sort group if relevant
                        if len(group) > 1:
                                for i in range(len(group)):
                                        group[i] = group[i].rstrip('\n').split('\t')
                                try:
                                        group.sort(key = lambda x: (float(x[10]),-float(x[11])))
                                except:
                                        print('"' + key + '"')
                                        print(group)
                                        quit()
                                for i in range(len(group)):
                                        group[i] = '\t'.join(group[i])
                        else:
                                group[0] = group[0].rstrip('\n')
                        # Put in output
                        for entry in group:
                                fileOut.write(entry + '\n')
        if delInputFile == True:
                os.unlink(inputFile)

def mms2sort_recompute(inputFile, outputFile, recomputeFasta):
        """
        Recompute requires two separate steps occur.
        1) The file needs to be sorted by query ID (i.e., column 1)
        2) After (1), the file needs to be sorted by E-value (i.e., column 11)
        This function performs (1), and feeds in the intermediate result to
        mms2sort_evalue for (2).
        """
        print("Starting step 1: initial sort by sequence IDs")
        # Parse recompute FASTA for sequence ID order
        records = SeqIO.parse(open(recomputeFasta, 'r'), "fasta")
        recomputeDict = OrderedDict()
        for record in records:
                recomputeDict[record.id] = []

        # Generate temporary file name
        tmpFileName = temp_file_name_gen("tmp_mms2sort_recompute", "out")

        # Parse convertalis file and store in memory
        """
        Yes, this is memory-intensive, but my use case is most often on a HPC.
        I just want it to run quickly.
        """
        with open(inputFile, 'r') as fileIn:
                for line in fileIn:
                        queryID = line.split()[0]
                        recomputeDict[queryID].append(line.rstrip("\r\n"))
        
        # Write output file with sorting by query ID
        with open(tmpFileName, 'w') as fileOut:
                for seqID, linesList in recomputeDict.items():
                        fileOut.write('\n'.join(linesList) + '\n')
        
        # Perform step 2, and call for clean up at the end
        print("Starting step 2: Sort by E-value")
        mms2sort_evalue(tmpFileName, outputFile, delInputFile=True)

def temp_file_name_gen(prefix, suffix):
        import os
        ongoingCount = 1
        while True:
                if not os.path.isfile(prefix + "." + suffix):
                        return prefix + "." + suffix
                elif os.path.isfile(prefix + "." + str(ongoingCount) + "." + suffix):
                        ongoingCount += 1
                else:
                        return prefix + "." + str(ongoingCount) + "." + suffix

# Main call
def main():
        #### USER INPUT SECTION
        usage = """Wrapper script to sort output MMseqs2 results. This functionality
        should eventually be part of run_mmseqs2.py, but it will require significant refactoring
        I don't presently have time for."""
        # Required
        p = argparse.ArgumentParser(description=usage)
        p.add_argument("-i", dest="inputFile", type = str,
                help="Specify the input convertalis file")
        p.add_argument("-o", dest="outputFile", type = str,
                help="Specify the output sorted file")
        p.add_argument("--sort_method", dest="sortMethod", type = str, choices = ['evalue', 'recompute'], default = 'recompute',
                help="Specify the method of sorting - use recompute if the entire file is unsorted, or evalue if sequence blocks are unsorted")
        p.add_argument("--recompute_fasta", dest="recomputeFasta", type = str,
                help="If sorting by \"recompute\", specify the query FASTA file to use as a sorting guide")
        args = p.parse_args()
        validate_args(args)

        if args.sortMethod == "recompute":
                mms2sort_recompute(args.inputFile, args.outputFile, args.recomputeFasta)
        elif args.sortMethod == "evalue":
                print("Starting sort by E-value")
                mms2sort_evalue(args.inputFile, args.outputFile)

        # Done!
        print('Program completed successfully!')

if __name__ == '__main__':
        main()

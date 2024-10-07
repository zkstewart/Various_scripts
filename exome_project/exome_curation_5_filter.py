#! python3
# exome_curation_5_filter.py
# Program to enable manual curation of exome sequencing
# to occur for the Oz Mammals Genomics initiative as part
# of Matthew Phillips and Andrew Baker (et. al.'s) group.

import sys, argparse, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) # 2 dirs up is where we find dependencies
from Function_packages import ZS_SeqIO
from Fasta_related import Alignment

def validate_args(args):
    # Validate input data location
    if not os.path.isdir(args.alignmentsDir):
        print('I am unable to locate the directory where the alignments files are (' + args.alignmentsDir + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Ensure paired args are provided
    if args.filterID != None and len(args.filterValues) == 0:
        print("If filterID is specified, you must also specify one or more filterValues")
        quit()
    # Handle file output
    if os.path.isfile(args.outputFileName):
        print(args.outputFileName + ' already exists. Specify a new file name or move/rename the existing file.')
        quit()

def main():
    usage = """%(prog)s receives a directory of MSAs that have been chunked using
    exome_curation_4_chunk.py and filters + concats these into a single MSA that can
    be used for phylogenetics.
    
    Note: This should be step 5 in the Oz Mammals project!
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-a", dest="alignmentsDir",
                   required=True,
                   help="Specify the directory where aligned FASTA files are located")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Specify the file name for the output MSA FASTA")
    # Opts
    p.add_argument("--rows_to_drop", dest="rowsToDrop",
                   required=False,
                   nargs="+",
                   default=["Codons", "GeneName"],
                   help="""Optionally, specify one or more space-separated IDs to drop these
                   sequences from the MSA (default == 'Codons GeneName')""")
    p.add_argument("--filter_id", dest="filterID",
                   required=False,
                   help="""Optionally specify a single row ID which contains values that will
                   mark columns to be dropped (default == 'Codons')""",
                   default="Codons")
    p.add_argument("--filter_values", dest="filterValues",
                   required=False,
                   nargs="+",
                   help="""Optionally, specify one or more space-separated values that mark
                   columns to be dropped (default == '4 5')""",
                   default=["4", "5"])
    
    args = p.parse_args()
    validate_args(args)
    
    # Locate all files
    files = [os.path.join(args.alignmentsDir, file) for file in os.listdir(args.alignmentsDir)]
    
    # Load FASTA files
    fastaObjs = []
    for file in files:
        f = ZS_SeqIO.FASTA(file, isAligned=True)
        fastaObjs.append(f)
    
    # Filter columns if relevant
    if args.filterID != None:
        for x in range(len(fastaObjs)):
            FASTA_obj = Alignment.drop_aligned_FASTA_columns_by_value(fastaObjs[x], args.filterID, args.filterValues)
            
            # Store the updated FASTA objects back in their list
            fastaObjs[x] = FASTA_obj
    
    # Drop any rows if relevant
    if len(args.rowsToDrop) != 0:
        for FASTA_obj in fastaObjs:
            Alignment.drop_FASTA_rows_by_id(FASTA_obj, args.rowsToDrop)
    
    # Concatenate FASTA objects
    for x in range(1, len(fastaObjs)):
        fastaObjs[0].concat(fastaObjs[x]) # it all gets concatenated into the first FASTA object
    
    # Write output
    fastaObjs[0].write(args.outputFileName, asAligned=True)
    
    # Done!
    print('Program completed successfully!')
    
if __name__ == "__main__":
    main()

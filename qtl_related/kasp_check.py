#! python3
# kasp_check.py
# Checks whether a (potential) KASP marker flags as positive
# in a WGS FASTQ file.

import os, argparse, gzip, codecs, screed
from contextlib import contextmanager

def get_codec(fileName):
    try:
        f = codecs.open(fileName, encoding='utf-8', errors='strict')
        for line in f:
            break
        f.close()
        return "utf-8"
    except:
        try:
            f = codecs.open(fileName, encoding='utf-16', errors='strict')
            for line in f:
                break
            f.close()
            return "utf-16"
        except UnicodeDecodeError:
            print(f"'{fileName}' is neither utf-8 nor utf-16 encoded; please convert to one of these formats.")

@contextmanager
def read_gz_file(filename):
    if filename.endswith(".gz"):
        with gzip.open(filename, "rt") as f:
            yield f
    else:
        with open(filename, "r", encoding=get_codec(filename)) as f:
            yield f

def parse_markers(markersFasta):
    # Store marker sequences in dict
    markerSeqs = {}
    for record in screed.open(markersFasta):
        prefix, suffix = record.name.split("_", maxsplit=1)
        markerSeqs.setdefault(prefix, {})
        markerSeqs[prefix][suffix] = record.sequence
    
    # Validate marker sequences
    for markerID, seqsDict in markerSeqs.items():
        if not "common" in seqsDict:
            raise ValueError(f"Marker '{markerID}' should have a 'common' sequence suffix")
        if not len(seqsDict) == 3:
            raise ValueError(f"Marker '{markerID}' should have 3 values in it, not {len(seqsDict)}")
    
    # Format the full markers to find
    markers = {}
    for markerID, seqsDict in markerSeqs.items():
        for seqID, seq in seqsDict.items():
            if seqID != "common":
                fullSeq = seq + seqsDict["common"]
                rcSeq = screed.rc(fullSeq)
                markers[f"{markerID}_{seqID}"] = [fullSeq, rcSeq]
    
    return markers

def find_marker(markers, fastqFile):
    # Convert markers for quick lookup
    markerSet = { x:markerID for markerID, seqs in markers.items() for x in seqs }
    sizes = set([ len(x) for x in markerSet.keys() ])
    minSize = min(sizes)
    
    # Look for markers in the FASTQ file
    foundMarkers = set()
    for record in screed.open(fastqFile):
        for i in range(0, len(record.sequence)-minSize+1):
            for size in sizes:
                s = record.sequence[i:i+size]
                if s in markerSet:
                    foundMarkers.add(markerSet[s])
    
    return foundMarkers

def validate_args(args):
    # Validate input files
    args.fastqFile = os.path.abspath(args.fastqFile)
    if not os.path.isfile(args.fastqFile):
        raise FileNotFoundError(f"-i FASTQ file '{args.fastqFile}' is not a file!")
    
    args.markersFasta = os.path.abspath(args.markersFasta)
    if not os.path.isfile(args.markersFasta):
        raise FileNotFoundError(f"-m markers FASTA file '{args.markersFasta}' is not a file!")
    
    # Validate output file location
    args.outputFileName = os.path.abspath(args.outputFileName)
    if os.path.isfile(args.outputFileName):
        raise FileExistsError(f"-o file '{args.outputFileName}' already exists!")

def main():
    usage = """%(prog)s checks whether one or more KASP markers flag as positive given an input
    FASTQ file. If your marker length is too short there's a good chance this will flag as positive
    on every file. Onus is on you to use this correctly.
    """
    # Establish main parser
    p = argparse.ArgumentParser(description=usage)
    
    # Set arguments shared by subparsers
    p.add_argument("-i", dest="fastqFile",
                   required=True,
                   help="Specify the FASTQ (optionally .gz'd) file to assess")
    p.add_argument("-m", dest="markersFasta",
                   required=True,
                   help="Specify the location of the marker(s) FASTA file")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="""Specify the location to write results to""")
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse markers
    markers = parse_markers(args.markersFasta)
    
    # Search through FASTQ file for markers
    foundMarkers = find_marker(markers, args.fastqFile)
    
    # Print and output results
    with open(args.outputFileName, "w") as fileOut:
        for marker in foundMarkers:
            fileOut.write(f"{marker}\n")
            print(f"# Found '{marker}'")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

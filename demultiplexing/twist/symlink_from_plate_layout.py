#!/usr/bin/env python3

import argparse
import os
import re

ROW_LABELS = ["A", "B", "C", "D", "E", "F", "G", "H"]
DIR_REGEX = re.compile(r"Plate\d{1,3}-Pool-\d{1}")
WELL_REGEX = re.compile(r"Well_\d{1,2}")

def validate_args(args):
    # Validate input file
    args.inputTSV = os.path.abspath(args.inputTSV)
    if not os.path.isfile(args.inputTSV):
        raise FileNotFoundError(f"-t '{inputFile}' is not a file or does not exist.")
    
    # Validate input directory
    args.demultiplexDirectory = os.path.abspath(args.demultiplexDirectory)
    if not os.path.isdir(args.demultiplexDirectory):
        raise NotADirectoryError(f"-d '{args.demultiplexDirectory}' is not a directory or does not exist.")

def parse_tsv(fileName):
    EXPECTED_FIRST = ["plate", "well", "sample"]
    
    layoutDict = {}
    with open(fileName, "r") as fileIn:
        firstLine = True
        for line in fileIn:
            sl = line.rstrip().split("\t")
            if firstLine:
                if sl != EXPECTED_FIRST:
                    raise ValueError(f"First line of '{fileName}' should be {EXPECTED_FIRST} not {sl}")
                firstLine = False
            else:
                plateNum, well, sample = sl
                if " " in sample:
                    oldSample = sample
                    sample = sample.replace(" ", "_")
                    print(f"# Sample with whitespace '{oldSample}' changed to be '{sample}'")
                
                poolNum = ROW_LABELS.index(well[0]) # first character should be the row letter
                wellNum = well[1:]
                
                plateKey = f"Plate{plateNum}"
                poolKey = f"Pool-{poolNum+1}"
                wellKey = f"Well_{wellNum}"
                
                layoutDict.setdefault(plateKey, {})
                layoutDict[plateKey].setdefault(poolKey, {})
                layoutDict[plateKey][poolKey][wellKey] = sample
    return layoutDict

def main():
    ##### USER INPUT SECTION
    usage = """%(prog)s parses the sequencing plate layout as derived from
    parse_sequencing_plate_layout.py or similar, and creates symlinks
    from the demultiplexed output of demultiplex_twist_pools.sh
    whereby each plate row (A to H) is a pool of 12 samples (for the plate
    columns). Symlinks will be created in the PWD.
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-t", dest="inputTSV",
                   required=True,
                   help="Specify the input plate layout TSV")
    p.add_argument("-d", dest="demultiplexDirectory",
                   required=True,
                   help="Specify the directory containing demultiplexed outputs")
    p.add_argument("--fwd", dest="fwdSuffix",
                   required=False,
                   help="""Optionally specify the file ending for the forward read;
                   default='_R1.fastq.gz'""",
                   default="_R1.fastq.gz")
    p.add_argument("--rvs", dest="rvsSuffix",
                   required=False,
                   help="""Optionally specify the file ending for the reverse read;
                   default='_R2.fastq.gz'""",
                   default="_R2.fastq.gz")
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse the layout TSV
    layoutDict = parse_tsv(args.inputTSV)
    
    # Search through the demultiplexed outputs
    located = {}
    for location in os.listdir(args.demultiplexDirectory):
        platePoolLocation = os.path.join(args.demultiplexDirectory, location)
        if "Plate" in location and "Pool-" in location and os.path.isdir(platePoolLocation):
            # Interpret the directory
            result = DIR_REGEX.search(location)
            if result is None:
                raise ValueError(f"Location '{location}' failed to match the directory regex for some reason")
            result = result.group()
            
            plate, pool = result.split("-", maxsplit=1)
            located.setdefault(plate, {})
            located[plate].setdefault(pool, {})
            
            # Search the directory for well outputs
            for sublocation in os.listdir(platePoolLocation):
                if "Well_" in sublocation:
                    result = WELL_REGEX.search(sublocation)
                    if result is None:
                        raise ValueError(f"Location '{sublocation}' failed to match the welll regex for some reason")
                    well = result.group()
                    
                    # Store forward or reverse read
                    located[plate][pool].setdefault(well, {"fwd": None, "rvs": None})
                    if sublocation.endswith(args.fwdSuffix):
                        located[plate][pool][well]["fwd"] = os.path.join(platePoolLocation, sublocation)
                    elif sublocation.endswith(args.rvsSuffix):
                        located[plate][pool][well]["rvs"] = os.path.join(platePoolLocation, sublocation)
                    else:
                        raise ValueError(f"Location '{sublocation}' should end with '{args.fwdSuffix}' or '{args.rvsSuffix}")
    
    # Check that the layout and demultiplexed outputs match
    layoutKeys = set([
        (plate, pool, well)
        for plate, poolDict in layoutDict.items()
        for pool, wellDict in poolDict.items()
        for well in wellDict.keys()
    ])
    locatedKeys = set([
        (plate, pool, well)
        for plate, poolDict in located.items()
        for pool, wellDict in poolDict.items()
        for well in wellDict.keys()
    ])
    layoutOnly = layoutKeys.difference(locatedKeys)
    locatedOnly = locatedKeys.difference(layoutKeys)
    
    if len(layoutOnly) > 0:
        raise ValueError(f"Layout TSV indicated samples not found in demultiplex directory, including: {layoutOnly}")
    if len(locatedOnly) > 0:
        raise ValueError(f"Demultiplex directory contains files not indicated by the layout TSV, including: {locatedOnly}")
    
    # Symlink everything
    for plate, pool, well in layoutKeys:
        sample = layoutDict[plate][pool][well]
        reads = located[plate][pool][well]
        
        fwdSrc, rvsSrc = reads["fwd"], reads["rvs"]
        fwdDest = os.path.join(os.getcwd(), f"{sample}{args.fwdSuffix}")
        rvsDest = os.path.join(os.getcwd(), f"{sample}{args.rvsSuffix}")
        
        os.symlink(fwdSrc, fwdDest)
        os.symlink(rvsSrc, rvsDest)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

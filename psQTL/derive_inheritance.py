#! python3
# derive_inheritance.py
# Script to take in a psQTL ED.tsv.gz file alongside a
# VCF (from which the ED file was derived) and produce
# a new VCF file with the derived inheritance information
# of highly segregating variants/CNVs.

import os, argparse, re, gzip, codecs, sys
from contextlib import contextmanager

#sys.path.append("/mnt/c/git/Various_scripts")
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Function_packages import ZS_VCFIO

# Various functions for program operations
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.vcfFile):
        raise FileNotFoundError(f'I am unable to locate the VCF file ({args.vcfFile})')
    if not os.path.isfile(args.edFile):
        raise FileNotFoundError(f'I am unable to locate the Euclidean Distance file ({args.edFile})')
    if not os.path.isfile(args.metadataFile):
        raise FileNotFoundError(f'I am unable to locate the metadata file ({args.metadataFile})')
    
    # Validate numeric arguments
    if args.cutoff < 0:
        raise ValueError(f"-c value '{args.cutoff}' must be >= 0")
    if args.power < 1:
        raise ValueError(f"--power value '{args.power}' must be >= 1!")
    if args.missingFilter < 0 or args.missingFilter > 1:
        raise ValueError(f"--missing value '{args.missingFilter}' must be between 0 and 1!")
    
    # Validate output file location
    if os.path.exists(args.outputFileName):
        raise FileExistsError(f'The specified output file name ({args.outputFileName}) already exists')

def get_codec(fileName):
    try:
        f = codecs.open(fileName, encoding='utf-8', errors='strict')
        for line in f:
            break
        return "utf-8"
    except:
        try:
            f = codecs.open(fileName, encoding='utf-16', errors='strict')
            for line in f:
                break
            return "utf-16"
        except UnicodeDecodeError:
            print(f"Can't tell what codec '{fileName}' is!!")

@contextmanager
def read_gz_file(filename):
    if filename.endswith(".gz"):
        with gzip.open(filename, "rt") as f:
            yield f
    else:
        with open(filename, "r", encoding=get_codec(filename)) as f:
            yield f

def parse_metadata(metadataFile):
    '''
    Copied from psQTL pipeline.
    
    Parameters:
        metadataFile -- a string indicating the path to a metadata file
    Returns:
        metadataDict -- a dictionary with structure like:
                        {
                            "bulk1": [ "sample1", "sample2", ... ],
                            "bulk2": [ "sample3", "sample4", ... ]
                        }
    '''
    ACCEPTED_BULK1 = ['bulk1', '1', 'bulk 1', 'b1']
    ACCEPTED_BULK2 = ['bulk2', '2', 'bulk 2', 'b2']
    
    metadataDict = {}
    foundSamples = set()
    with read_gz_file(metadataFile) as fileIn:
        for line in fileIn:
            l = line.rstrip("\r\n ")
            # Skip blank lines
            if l == "":
                continue
            
            # Skip comment lines
            if line.startswith("#"):
                continue
            
            # Split line based on delimiter
            if "\t" in line:
                sl = l.split("\t")
            elif "," in line:
                sl = l.split(",")
            else:
                raise ValueError(f"Metadata file is not tab or comma-delimited; offending line is '{l}'")
            
            # Parse out relevant information
            try:
                sample, pop = sl
                sample, pop = sample.strip(), pop.strip().lower() # make lowercase for easier comparison
            except ValueError:
                raise ValueError(f"Metadata file does not have two columns; offending line is '{l}'")
            
            # Validate that the population is one of the expected values
            if not pop in ACCEPTED_BULK1 + ACCEPTED_BULK2:
                raise ValueError(f"'{pop}' is not in the expected format for denoting bulks; offending line is '{l}'")
            
            # Unify the population names
            if pop in ACCEPTED_BULK1:
                pop = 'bulk1'
            elif pop in ACCEPTED_BULK2:
                pop = 'bulk2'
            
            # Validate that the sample is non-redundant
            if sample in foundSamples:
                raise ValueError(f"Sample '{sample}' is listed more than once in the metadata file")
            foundSamples.add(sample)
            
            # Store the sample in the appropriate population
            metadataDict.setdefault(pop, set())
            metadataDict[pop].add(sample)
    
    # Make sure that the metadata file has the expected number of populations
    if len(metadataDict) == 0:
        raise ValueError("Metadata file is empty; please provide at least two samples from different bulks")
    if len(metadataDict) == 1:
        foundBulk = list(metadataDict.keys())[0]
        raise ValueError(f"Metadata file should have 2 populations ('bulk1' and 'bulk2'); I only found '{foundBulk}'")
    
    # Reformat sets and lists then return
    for pop in metadataDict:
        metadataDict[pop] = list(metadataDict[pop])
    return metadataDict

def parse_regions(rawRegions):
    # Validate regions
    regions = []
    regionsRegex = re.compile(r"^([^:]+):(\d+)-(\d+)$")
    for region in rawRegions:
        reMatch = regionsRegex.match(region)
        
        # Handle chr:start-end format
        if reMatch != None:
            contigID, start, end = reMatch.groups()
            start = int(start)
            end = int(end)
            
            # Validate start and end positions
            if start < 0:
                raise ValueError(f"--region start position '{start}' is < 0!")
            if start >= end:
                raise ValueError(f"--region start position '{start}' is >= end position '{end}'!")
            # Store region
            regions.append([contigID, start, end])
        # Handle invalid format
        elif ":" in region:
            raise ValueError(f"Invalid region input '{region}'; you included a ':' but did " + 
                             "not format the region as 'chr:start-end'!")
        # Handle chr format
        else:
            raise ValueError(f"Invalid region input '{region}'; you did not include a ':' " +
                             "to specify a range within the chromosome!")
    return regions

def parse_ed_file(edFile, regions, metadataDict, power, missingFilter, cutoff):
    edDict = {}
    BULK1_ALLELES = len(metadataDict["bulk1"]) * 2
    BULK2_ALLELES = len(metadataDict["bulk2"]) * 2
    
    with read_gz_file(edFile) as fileIn:
        header = fileIn.readline()
        for line in fileIn:
            chrom, posi, variant, b1alleles, b2allales, edist = line.rstrip("\r\n ").split("\t")
            posi = int(posi)
            b1alleles = int(b1alleles)
            b2allales = int(b2allales)
            edist = float(edist) ** power
            
            # Check if the variant is significantly segregating
            if edist < cutoff:
                continue
            
            # Check if missing filter would pass
            if (b1alleles / BULK1_ALLELES) < missingFilter or (b2allales / BULK2_ALLELES) < missingFilter:
                continue
            
            # Check if the variant is in the specified regions
            skipThis = True
            for region in regions:
                if chrom == region[0] and posi >= region[1] and posi <= region[2]:
                    skipThis = False
                    break
            else:
                continue
            if skipThis:
                continue
            
            # Store the variant
            edDict.setdefault(chrom, {})
            edDict[chrom][posi] = edist
    return edDict

def gt_code_to_allele(gtCode, ref_alt):
    return "/".join([ref_alt[gt] for gt in gtCode])

def gt_code_format(gtCode):
    return "/".join(map(str, gtCode))

# Main call
def main():
    #### USER INPUT SECTION
    usage = """%(prog)s ..."""
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-v", dest="vcfFile",
                   required=True,
                   help="Input VCF file")
    p.add_argument("-e", dest="edFile",
                   required=True,
                   help="Input Euclidean Distance file")
    p.add_argument("-m", dest="metadataFile",
                   required=True,
                   help="Input metadata file")
    p.add_argument("-r", dest="regions",
                   required=True,
                   nargs="+",
                   help="""Specify which regions to check. Give the
                   chromosome(s) to plot by their individual genome contig identifiers
                   (e.g., 'chr1') with the option to specify a range within the chromosome
                   with chr:start-end format (e.g., 'chr1:1000000-2000000').
                   """,
                   default=[])
    p.add_argument("-c", dest="cutoff",
                   required=True,
                   type=float,
                   help="""Specify the cutoff for the Euclidean Distance value
                   to determine if a variant is significantly segregating""")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output FASTA file name for representative sequences")
    # Optional arguments
    p.add_argument("--showNumGT", dest="showNumGT",
                   required=False,
                   action="store_true",
                   help="""Provide this flag if you would like to show the genotype numbers,
                   rather than the nucleotides.""",
                   default=False)
    p.add_argument("--power", dest="power",
                   required=False,
                   type=int,
                   help="""Optionally, specify the power to raise Euclidean distances to
                   reduce noise (default: 4)""",
                   default=4)
    p.add_argument("--missing", dest="missingFilter",
                   type=float,
                   required=False,
                   help="""Optionally, specify the proportion of missing data that is
                   tolerated in both bulk populations before a variant is filtered out
                   (recommended: 0.5)""",
                   default=0.5)
    p.add_argument("--parents1", dest="bulk1Parents",
                   required=False,
                   nargs="+",
                   help="Specify one or more parent IDs for bulk 1",
                   default=[])
    p.add_argument("--parents2", dest="bulk2Parents",
                   required=False,
                   nargs="+",
                   help="Specify one or more parent IDs for bulk 2",
                   default=[])
    args = p.parse_args()
    validate_args(args)
    
    # Parse metadata file
    metadataDict = parse_metadata(args.metadataFile)
    metadataSamples = set([ # used for error checking later
        sample
        for samples in metadataDict.values()
        for sample in samples
    ])
    
    # Parse regions
    args.regions = parse_regions(args.regions)
    
    # Parse ED file
    edDict = parse_ed_file(args.edFile, args.regions, metadataDict,
                           args.power, args.missingFilter, args.cutoff)
    
    # Write results to output file
    with open(args.outputFileName, "w") as fileOut:
        # Get consistently ordered sample IDs
        parentIDs = set(args.bulk1Parents + args.bulk2Parents)
        bulk1Samples = [ x for x in metadataDict["bulk1"] if x not in parentIDs ]
        bulk2Samples = [ x for x in metadataDict["bulk2"] if x not in parentIDs ]
        bulk1Samples.sort()
        bulk2Samples.sort()
        
        # Write header
        fileOut.write("contig\tposition\tED" + "\t")
        if len(args.bulk1Parents) > 0:
            fileOut.write("\t".join(args.bulk1Parents) + "\t")
        fileOut.write("\t".join(bulk1Samples) + "\t")
        if len(args.bulk2Parents) > 0:
            fileOut.write("\t".join(args.bulk2Parents) + "\t")
        fileOut.write("\t".join(bulk2Samples) + "\n")
        
        # Iterate through VCF file
        vcfIterator = ZS_VCFIO.SimpleGenotypeIterator(args.vcfFile)
        firstLine = True
        for values in vcfIterator:
            # Grab header line containing sample IDs
            if firstLine:
                vcfSamples = set(values)
                firstLine = False
                
                # Check that the metadata file matches the VCF file
                if vcfSamples != metadataSamples:
                    vcfDiff = vcfSamples.difference(metadataSamples)
                    metadataDiff = metadataSamples.difference(vcfSamples)
                    
                    print("ERROR: Metadata file does not match the VCF file!")
                    
                    if len(vcfDiff) > 0:
                        print("In your VCF, the following samples exist which are " + 
                            "absent from the metadata: ", ", ".join(vcfDiff))
                    if len(metadataDiff) > 0:
                        print("In your metadata, the following samples exist which are " + 
                            "absent from the VCF: ", ", ".join(metadataDiff))
                    quit()
            # Handle content lines
            else:
                contig, pos, ref, alt, snpDict = values
                ref_alt = [ref, *alt]
                
                # Skip if the variant is not in the ED dict
                if not contig in edDict or not pos in edDict[contig]:
                    continue
                
                # Format values for output
                b1ParentGTs = [
                    "./." if sample not in snpDict else
                    gt_code_to_allele(snpDict[sample], ref_alt) if not args.showNumGT else
                    gt_code_format(snpDict[sample])
                    for sample in args.bulk1Parents
                ]
                b2ParentGTs = [
                    "./." if sample not in snpDict else
                    gt_code_to_allele(snpDict[sample], ref_alt) if not args.showNumGT else
                    gt_code_format(snpDict[sample])
                    for sample in args.bulk2Parents
                ]
                
                b1GTs = [
                    "./." if sample not in snpDict else 
                    gt_code_to_allele(snpDict[sample], ref_alt) if not args.showNumGT else
                    gt_code_format(snpDict[sample])
                    for sample in bulk1Samples
                ]
                b2GTs = [
                    "./." if sample not in snpDict else 
                    gt_code_to_allele(snpDict[sample], ref_alt) if not args.showNumGT else
                    gt_code_format(snpDict[sample])
                    for sample in bulk2Samples
                ]
                
                # Write formatted output
                row = [contig, str(pos), str(edDict[contig][pos])]
                for gtList in [b1ParentGTs, b1GTs, b2ParentGTs, b2GTs]:
                    row.extend(gtList)
                fileOut.write("\t".join(row) + "\n")
    
    # Done!
    print('Program completed successfully!')

if __name__ == '__main__':
    main()

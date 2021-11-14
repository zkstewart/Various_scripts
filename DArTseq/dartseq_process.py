#! python3
# Script to perform preprocessing operations on DArTseq reads,
# including barcode removal, demultiplexing, and contamination
# cleaning via stacks process_radtags

import os, argparse, subprocess

# Define functions
def validate_args(args):
    # Validate input file location
    if not os.path.isfile(args.multiplexedFastq):
        print('I am unable to locate the FASTQ file (' + args.multiplexedFastq + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()

def parse_metadata_csv(metadataCsv, speciesIdCol, barcodeCol):
    header = None
    speciesIds = []
    barcodes = []
    
    with open(metadataCsv, "r") as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n ").split(",")
            # Handle header line
            if header == None:
                header = sl
                # Raise errors if columns don't exist
                if speciesIdCol not in header:
                    raise Exception("Species ID column does not exist")
                elif barcodeCol not in header:
                    raise Exception("Barcode column does not exist")
                # Get column indices if they exist
                else:
                    speciesIdIndex = header.index(speciesIdCol)
                    barcodeIndex = header.index(barcodeCol)
            # Handle content lines
            else:
                speciesIds.append(sl[speciesIdIndex])
                barcodes.append(sl[barcodeIndex])
    return speciesIds, barcodes

def create_barcode_file(speciesIds, barcodes, outputFileName="process_radtags_barcodes.tsv"):
    # Validations
    if os.path.isfile(outputFileName):
        raise FileExistsError("File name <{0}> already exists".format(outputFileName))
    assert len(speciesIds) == len(barcodes)
    
    # Write to file
    with open(outputFileName, "w") as fileOut:
        for i in range(len(speciesIds)):
            speciesId = speciesIds[i]
            barcode = barcodes[i]
            fileOut.write("{0}\t{1}\n".format(speciesId, barcode))

def run_process_radtags(multiplexedFastq, enzyme, barcodesFileName="process_radtags_barcodes.tsv", outputDirectory="process_radtags_out"):
    # Format cmd
    cmd = "process_radtags -f {0} -b {1} -o {2} -e {3} --inline_null -c -q -r --disable-rad-check".format(multiplexedFastq, barcodesFileName, outputDirectory, enzyme)
    run = subprocess.Popen(cmd, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    out, err = run.communicate()
    out = out.decode("utf-8")
    print("stdout\n{0}\n".format(out))
    print("-----\nstderr\n{0}\n".format(err))

def main():
    # User input
    usage = """%(prog)s does things...
    
    It assumes process_radtags can be found in environment variables.
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-mfq", dest="multiplexedFastq",
                   required=True,
                   help="Input multiplexed fastq file")
    p.add_argument("-csv", dest="metadataCsv",
                   required=True,
                   help="Input Metadata CSV file")
    p.add_argument("-e", dest="enzyme",
                   required=True,
                   choices=[
                        'aciI', 'ageI', 'aluI', 'apaLI', 'apeKI', 'apoI', 'aseI', 'bamHI',
                        'bbvCI', 'bfaI', 'bfuCI', 'bgIII', 'bsaHI', 'bspDI', 'bstYI', 'btgI',
                        'cac8I', 'claI', 'csp6I', 'ddeI', 'dpnII', 'eaeI', 'ecoRI', 'ecoRV',
                        'ecoT22I', 'haeIII', 'hinP1I', 'hindIII', 'hpaII', 'hpyCH4IV', 'kpnI', 'mluCI',
                        'mseI', 'mslI', 'mspI', 'ncoI', 'ndeI', 'ngoMIV', 'nheI', 'nlaIII',
                        'notI', 'nsiI', 'nspI', 'pacI', 'pspXI', 'pstI', 'rsaI', 'sacI',
                        'sau3AI', 'sbfI', 'sexAI', 'sgrAI', 'speI', 'sphI', 'taqI', 'xbaI', 
                        'xhoI'
                    ],
                   help="Restriction enzyme used for cutting (only 1 supported for this script)")
    p.add_argument("-id", dest="speciesIdCol",
                   required=True,
                   help="Column name where species ID is located")
    p.add_argument("-b", dest="barcodeCol",
                   required=True,
                   help="Column name where column string is located")
    args = p.parse_args()
    validate_args(args)

    # Parse CSV file columns
    speciesIds, barcodes = parse_metadata_csv(args.metadataCsv, args.speciesIdCol, args.barcodeCol)

    # Create barcodes file
    create_barcode_file(speciesIds, barcodes)
    
    # Run process_radtags
    run_process_radtags(args.multiplexedFastq, args.enzyme)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

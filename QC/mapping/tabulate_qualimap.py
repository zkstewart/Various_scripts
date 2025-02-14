#! python3
# tabulate_qualimap.py
# Simple script to read in the output directories from
# Qualimap and tabulate the results for easy understanding

import os, argparse, re
import pandas as pd
from xlsxwriter.utility import xl_col_to_name
from statistics import median

basicStatsRegex = re.compile(r"Total Sequences<\/td><td>(\d+)<\/td><\/tr><tr><td>Sequences flagged as poor quality<\/td><td>(\d+)<\/td>")
passFailRegex = re.compile(r"alt=\"\[(\w+?)\]\"\/><a href=\"#M\d+\">([\w\s]+?)<\/a>")

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isdir(args.inputDir):
        raise FileNotFoundError(f"I am unable to locate the parent directory where Qualimap subdirs are ({args.inputDir}). " + 
                                "Make sure you've typed the file name or location correctly and try again.")
    
    # Validate output file location
    args.outputFileName = os.path.abspath(args.outputFileName)
    if not args.outputFileName.endswith(".xlsx"):
        args.outputFileName += ".xlsx"
    
    if os.path.exists(args.outputFileName):
        raise FileExistsError(f"File already exists at output location ({args.outputFileName}). " +
                              "Make sure you specify a unique file name and try again.")
    elif not os.path.isdir(os.path.dirname(args.outputFileName)):
        raise FileNotFoundError(f"Output directory does not exist ({os.path.dirname(args.outputFileName)}). " +
                                "Make sure you specify a parent directory for your output file which already exists.")

def parse_genome_results(genomeResultsFile):
    '''
    Parameters:
        genomeResultsFileir -- a string with the path to a Qualimap 'genome_results.txt' file
    Returns:
        genomeResultsDict -- a dictionary with structure like:
                             {
                                 "median_insert_size": int,
                                 "mean_mapping_quality": float,
                                 "general_error_rate": float,
                                 "mean_coverage": float,
                                 "std_coverage": float,
                                 "1x": float,
                                 "5x": float,
                                 "10x": float,
                                 "15x": float,
                                 "20x": float,
                                 "25x": float
                             }
    '''
    genomeResultsDict = {}
    with open(genomeResultsFile, "r") as fileIn:
        for line in fileIn:
            l = line.rstrip("\r\n ")
            if "median insert size" in l:
                genomeResultsDict["median_insert_size"] = int(l.split(" ")[-1])
            if "mean mapping quality" in l:
                genomeResultsDict["mean_mapping_quality"] = float(l.split(" ")[-1])
            if "general error rate" in l:
                genomeResultsDict["general_error_rate"] = float(l.split(" ")[-1])
            if "mean coverageData" in l:
                genomeResultsDict["mean_coverage"] = float(l.split(" ")[-1][:-1].replace(",", "")) # drop the X at the end
            if "std coverageData" in l:
                genomeResultsDict["std_coverage"] = float(l.split(" ")[-1][:-1].replace(",", ""))
            if "reference with a coverageData >= 1X" in l:
                genomeResultsDict["1x"] = float(l.split("There is a ")[1].split("%")[0])
            if "reference with a coverageData >= 5X" in l:
                genomeResultsDict["5x"] = float(l.split("There is a ")[1].split("%")[0])
            if "reference with a coverageData >= 10X" in l:
                genomeResultsDict["10x"] = float(l.split("There is a ")[1].split("%")[0])
            if "reference with a coverageData >= 15X" in l:
                genomeResultsDict["15x"] = float(l.split("There is a ")[1].split("%")[0])
            if "reference with a coverageData >= 20X" in l:
                genomeResultsDict["20x"] = float(l.split("There is a ")[1].split("%")[0])
            if "reference with a coverageData >= 25X" in l:
                genomeResultsDict["25x"] = float(l.split("There is a ")[1].split("%")[0])
    return genomeResultsDict

def parse_cov_ref_file(covRefFile):
    '''
    Parameters:
        covRefFile -- a string with the path to a Qualimap 'coverage_across_reference.txt' file
    Returns:
        covRefDict -- a dictionary with structure like:
                      {
                          TBD
                      }
    '''
    covRefDict = {}
    
    coverages = []
    with open(covRefFile, "r") as fileIn:
        for line in fileIn:
            if line.startswith("#"):
                continue
            
            pos, cov, std = line.rstrip("\r\n ").split("\t")
            coverages.append(float(cov))
    
    covRefDict["minimum_cov"] = min(coverages)
    covRefDict["median_cov"] = median(coverages)
    covRefDict["maximum_cov"] = max(coverages)
    
    return covRefDict

def get_results_from_qm_files(qmFiles):
    '''
    Parameters:
        qmFiles -- a dictionary with structure like:
                    {
                        "sample1": ["path/to/sample1_qualimap/genome_results.txt", "path/to/sample1_qualimap/raw_data_qualimapReport/coverage_across_reference.txt"],
                        "sample2": ["path/to/sample2_qualimap/genome_results.txt", "path/to/sample2_qualimap/raw_data_qualimapReport/coverage_across_reference.txt"],
                        ...
                    }
    Returns:
        TBD
    '''
    qmDict = {}
    
    # Parse all files
    for sample, (genomeResultsFile, covRefFile) in qmFiles.items():
        genomeResultsDict = parse_genome_results(genomeResultsFile)
        covRefDict = parse_cov_ref_file(covRefFile)
        
        # Combine dicts
        qmDict[sample] = {**covRefDict, **genomeResultsDict}
    
    # Convert to table
    qmTable = pd.DataFrame(qmDict).T
    
    # Order columns
    qmTable = qmTable[["median_insert_size", "mean_mapping_quality", "general_error_rate",
                       "minimum_cov", "median_cov", "mean_coverage", "maximum_cov", "std_coverage",
                       "1x", "5x", "10x", "15x", "20x", "25x"]]
    return qmTable

def write_table_to_excel(table, sheetName, outputFileName):
    writer = pd.ExcelWriter(outputFileName, engine = "xlsxwriter")
    
    table.to_excel(writer, sheet_name = sheetName)
    workbook = writer.book
    worksheet = writer.sheets[sheetName]
    
    # Set column width based on header length
    columnLength = max([len(str(x)) for x in table.index])
    worksheet.set_column(0, 0, columnLength+1)
    
    for i, header in enumerate(table.columns):
        columnLength = max(table[header].astype(str).map(len).max(), len(header))
        worksheet.set_column(i+1, i+1, columnLength+1)
    
    # Figure out the shape/column indices where data is stored
    # columnNumbers = table.columns.get_indexer(table.columns).tolist()
    # columnLetters = list(map(xl_col_to_name, [i+1 for i in columnNumbers]))
    # numRows, numCols = table.shape
    
    writer.close()

## Main
def main():
    # User input
    usage = """%(prog)s accepts a directory containing subdirectories with Qualimap outputs
    and parses their contents to create a table containing an overview of their statistics.
    The output file will be written in Excel format; if your -o value does not end with 
    '.xlsx', the program will append it for you.
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="inputDir",
                   required=True,
                   help="Input directory containing Qualimap result subdirs")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file name for the reads count table")
    
    args = p.parse_args()
    validate_args(args)
    
    # Locate all files
    qmFiles = {}
    for directory in [os.path.join(args.inputDir, d) for d in os.listdir(args.inputDir) ]:
        if os.path.isdir(directory):
            # Get the genome_results.txt file
            genomeResultFile = os.path.join(directory, "genome_results.txt")
            if not os.path.exists(genomeResultFile):
                continue
            
            # Get the coverage_across_reference.txt file
            coverageFile = os.path.join(directory, "raw_data_qualimapReport", "coverage_across_reference.txt")
            if not os.path.exists(coverageFile):
                continue
            
            # Store in dictionary
            qmFiles[os.path.basename(directory)] = [genomeResultFile, coverageFile]
    
    # Parse all files
    qmTable = get_results_from_qm_files(qmFiles)
    
    # Write content
    try:
        write_table_to_excel(qmTable, "qualimap_qc", args.outputFileName)
    except Exception as e:
        os.unlink(args.outputFileName)
        raise e
    
    # Alert user to program success
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

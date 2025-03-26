#! python3
# tabulate_fastqc.py
# Simple script to read in the output directories from
# FASTQC and tabulate the results for easy understanding

import os, argparse, re
import pandas as pd
from xlsxwriter.utility import xl_col_to_name

passFailRegex = re.compile(r"alt=\"\[(\w+?)\]\"\/><a href=\"#M\d+\">([\w\s]+?)<\/a>")

# Stats regex for older FASTQC versions
basicStatsRegex = re.compile(r"Total Sequences<\/td><td>(\d+)<\/td><\/tr><tr><td>Sequences flagged as poor quality<\/td><td>(\d+)<\/td>")

# Stats regex for newer FASTQC versions
basicStatsRegex2 = re.compile(r"Total Sequences<\/td><td>(\d+)<\/td><\/tr><tr><td>Total Bases<\/td><td>[\.\s\w]+?<\/td><\/tr><tr><td>Sequences flagged as poor quality<\/td><td>(\d+)<\/td>")

# Define functions
def validate_args(args):
    # Validate input file locations
    if not os.path.isdir(args.fqcDir):
        raise FileNotFoundError(f"I am unable to locate the parent directory where FASTQC subdirs are ({args.fqcDir}). " + 
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

def parse_fcq_file(fqcFile):
    '''
    Parameters:
        fqcFile -- a string with the path to a FASTQC output file
    Returns:
        statsDict -- a dictionary with the following keys:
                     total_reads, low_quality_reads, per_base_sequence_quality, per_tile_sequence_quality,
                     per_sequence_quality_scores, per_base_sequence_content, per_sequence_gc_content,
                     per_base_n_content, sequence_length_distribution, sequence_duplication_levels,
                     overrepresented_sequences, adapter_content
    '''
    statsDict = {}
    with open(fqcFile, "r") as fileIn:
        contents = fileIn.read()
        
        # Parse results from file
        try:
            numReads, poorQuality = basicStatsRegex.search(contents).groups()
        except:
            numReads, poorQuality = basicStatsRegex2.search(contents).groups()
        passFail = passFailRegex.findall(contents)
        
        # Add in basic stats
        statsDict = {}
        statsDict["total_reads"] = int(numReads)
        statsDict["low_quality_reads"] = int(poorQuality)
        
        # Convert pass/fail results to dictionary
        statsDict.update({
            category.lower().replace(" ", "_"): flag
            for flag, category in passFail
        })
    
    return statsDict

def get_results_from_fqc_files(fqcFiles):
    '''
    Parameters:
        fqcFiles -- a dictionary with structure like:
                    {
                        "sample1": ["path/to/sample1_1_fastqc.html", "path/to/sample1_2_fastqc.html"],
                        "sample2": ["path/to/sample2_1_fastqc.html", "path/to/sample2_2_fastqc.html"],
                        ...
                    }
    Returns:
        forwardTable -- a pandas DataFrame with the following columns:
                        sample, total_reads, low_quality_reads, per_base_sequence_quality, per_tile_sequence_quality,
                        per_sequence_quality_scores, per_base_sequence_content, per_sequence_gc_content,
                        per_base_n_content, sequence_length_distribution, sequence_duplication_levels,
                        overrepresented_sequences, adapter_content
        reverseTable -- a pandas DataFrame with the same columns as forwardTable
    '''
    forwardDict = {}
    reverseDict = {}
    
    # Parse all files
    for sample, (fHtml, rHtml) in fqcFiles.items():
        # Parse forward/reverse reports
        fDict = parse_fcq_file(fHtml)
        rDict = parse_fcq_file(rHtml)
        
        # Store in larger dictionary
        forwardDict[sample] = fDict
        reverseDict[sample] = rDict
    
    # Convert to table
    forwardTable = pd.DataFrame(forwardDict).T
    reverseTable = pd.DataFrame(reverseDict).T
    
    return forwardTable, reverseTable

def write_tables_to_excel(forwardTable, reverseTable, outputFileName):
    writer = pd.ExcelWriter(outputFileName, engine = "xlsxwriter")
    for table, sheetName in zip([forwardTable, reverseTable], ["Forward Reads", "Reverse Reads"]):
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
        columnNumbers = table.columns.get_indexer(table.columns).tolist()
        columnLetters = list(map(xl_col_to_name, [i+1 for i in columnNumbers]))
        numRows, numCols = table.shape
        
        # Colour PASS/WARNING/FAIL values
        passFormat = workbook.add_format({"bg_color": "#00CC66"})
        warningFormat = workbook.add_format({"bg_color": "#FFCC00"})
        failFormat = workbook.add_format({"bg_color": "#FF0000"})
        
        worksheet.conditional_format(f"{columnLetters[0]}2:{columnLetters[-1]}{str(numRows+1)}",
                                     {"type": "cell", "criteria": "==",
                                      "value": '"PASS"', "format": passFormat})
        worksheet.conditional_format(f"{columnLetters[0]}2:{columnLetters[-1]}{str(numRows+1)}",
                                     {"type": "cell", "criteria": "==",
                                      "value": '"WARNING"', "format": warningFormat})
        worksheet.conditional_format(f"{columnLetters[0]}2:{columnLetters[-1]}{str(numRows+1)}",
                                     {"type": "cell", "criteria": "==",
                                      "value": '"FAIL"', "format": failFormat})
    writer.close()

## Main
def main():
    # User input
    usage = """%(prog)s accepts a directory containing subdirectories with FASTQC output files
    and parses their contents to create a table containing an overview of their statistics.
    The output file will be written in Excel format; if your -o value does not end with 
    '.xlsx', the program will append it for you.
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="fqcDir",
                   required=True,
                   help="Input directory containing FASTQC result subdirs")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output file name for the reads count table")
    
    args = p.parse_args()
    validate_args(args)
    
    # Locate all files
    fqcFiles = {}
    for directory in [os.path.join(args.fqcDir, f) for f in os.listdir(args.fqcDir) ]:
        if os.path.isdir(directory):
            # Get the html files
            htmlFiles = [ f for f in os.listdir(directory) if f.endswith(".html") ]
            if len(htmlFiles) == 0:
                raise FileNotFoundError(f"No FASTQC result files found in '{directory}'!")
            
            # Sort forward and reverse files
            if len(htmlFiles) > 1:
                commonPrefix = os.path.commonprefix(htmlFiles)
                htmlFiles.sort(key=lambda x: int(x.replace(commonPrefix, "")[0]))
            
            # Store in dictionary
            fqcFiles[os.path.basename(directory)] = [ os.path.join(directory, f) for f in htmlFiles ]
    
    # Parse all files
    forwardTable, reverseTable = get_results_from_fqc_files(fqcFiles)
    
    # Write content
    try:
        write_tables_to_excel(forwardTable, reverseTable, args.outputFileName)
    except Exception as e:
        os.unlink(args.outputFileName)
        raise e
    
    # Alert user to program success
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

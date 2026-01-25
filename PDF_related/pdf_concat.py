#!/usr/bin/env python3
# pdf_concat.py
# A script to provide simple capabilities for combining PDFs with page selection
# and rotation

import os, argparse, PyPDF2
from collections import Counter

def validate_args(args):
    args.instructions = parse_pdf(args.pdfs)
    if os.path.exists(args.outputName):
        raise FileExistsError(f"-o file '{args.outputName}' already exists and will not be overwritten")

def parse_pdf(pdfs):
    '''
    Parses out the argument input format which specifies page numbers and rotations.
    This format is expected to occur in one of four variants:
      1) 'Just the file name' [filename.pdf]
         - Take a file as-is and concatenate it in its entirety.
      2) 'Subset of the file' [filename.pdf:1-5,7,9,10-end]
         - Take a file, extract specific pages, and concatenate those.
      3) 'Rotation of the file' [filename.pdf:+90]
         - Take a file, rotate all pages by the specified degrees, and concatenate.
      4) 'Subset and rotate' [filename.pdf:1-3,5-end:-45]
         - Take a file, rotate a subset of its pages, and concatenate those.
    '''
    instructions = {}
    for x, inputStr in enumerate(pdfs):
        sections = inputStr.split(":")
        if len(sections) > 3:
            raise ValueError(f"-i value '{inputStr}' can only have up to 3 sections (delimited by ':') " +
                             f"but it instead has {len(sections)}")
        
        for i, section in enumerate(sections):
            # Parse first mandatory section: file name
            if i == 0:
                fileName = section
                if not os.path.isfile(fileName):
                    raise FileNotFoundError(f"-i file '{fileName}' does not exist; check the name or location for typos")
                instructions[x] = {"file": fileName, "pages": None, "rotation": None}
                documentLength = len(PyPDF2.PdfReader(fileName).pages) # to be used in following sections of this input
            
            # Parse remaining optional sections
            else:
                # Rotation
                if section.startswith("+") or section.startswith("-"):
                    direction = section[0]
                    degrees = section[1:]
                    if not degrees.isdigit():
                        raise ValueError(f"The section of -i '{inputStr}' starting with a '{direction}' " +
                                         f"does not have an integer degrees value ('{degrees}')")
                    
                    degrees = int(degrees)
                    if degrees == 0:
                        print(f"Note: rotation of '{section}' has 0 degrees which means no change is applied")
                        continue
                    elif degrees > 360:
                        raise ValueError(f"Degrees value must not exceed 360; you gave a value of '{degrees}'")
                    
                    if degrees % 90 != 0:
                        raise ValueError(f"Degrees value should be a multiple of 90; you gave a value of '{degrees}'")
                    
                    # Store rotation instruction
                    if instructions[x]["rotation"] != None:
                        raise ValueError(f"The section of -i '{inputStr}' specifies rotation direction more than once")
                    instructions[x]["rotation"] = eval(section)
                
                # Subsetting
                else:
                    # Find which pages have been selected
                    pages = section.replace(" ", "").split(",")
                    selectedPages = []
                    for p in pages:
                        originalp = p # for error reporting
                        
                        # Impute the document length
                        if "end" in p:
                            p = p.replace("end", str(documentLength))
                        
                        # Handle ranges
                        if "-" in p:
                            # Ensure range is formatted correctly
                            startEnd = p.split("-")
                            if len(startEnd) != 2:
                                raise ValueError(f"The section of -i '{inputStr}' specifies a page range '{originalp}' " + 
                                                 "with more than one '-' delimiter; expected format is as 'start-end' values")
                            start, end = startEnd
                            
                            # Ensure values are integers
                            try:
                                start = int(start)
                            except ValueError:
                                raise ValueError(f"The section of -i '{inputStr}' specifies a non-integer page value '{start}'")
                            try:
                                end = int(end)
                            except ValueError:
                                raise ValueError(f"The section of -i '{inputStr}' specifies a non-integer page value '{end}'")
                            
                            # Ensure start-end range is valid
                            reverse = False
                            if start > end:
                                start, end = end, start # do this so we can apply a single check for start/end validity
                                reverse = True # note how to handle start/end when storing page ranges
                                
                            if start < 1:
                                raise ValueError(f"The section of -i '{inputStr}' specifies a start ('{start}') which is less than 1 (i.e., the first page)")
                            if end > documentLength:
                                raise ValueError(f"The section of -i '{inputStr}' specifies an end ('{end}') which exceeds " +
                                                 f"the number of pages in the '{sections[0]}' file (i.e., '{documentLength}' pages)")
                            
                            # Store page ranges with handling for reversing
                            if reverse:
                                selectedPages.extend(range(end-1, start-2, -1)) # reversed conversion of 1-based inclusive to 0-based inclusive range
                            else:
                                selectedPages.extend(range(start-1, end)) # convert 1-based inclusive to 0-based inclusive range
                        
                        # Handle individual pages
                        else:
                            try:
                                pageNum = int(p)
                            except ValueError:
                                raise ValueError(f"The section of -i '{inputStr}' specifies a non-integer page value '{pageNum}'")

                            selectedPages.append(pageNum-1) # convert 1-based to 0-based
                    
                    # Store subsetting instruction
                    if instructions[x]["pages"] != None:
                        raise ValueError(f"The section of -i '{inputStr}' specifies page selections more than once")
                    instructions[x]["pages"] = selectedPages
        
        # Impute pages if needed and warn about duplicate pages
        if instructions[x]["pages"] == None:
            instructions[x]["pages"] = list(range(0, documentLength))
        
        if len(set(instructions[x]["pages"])) != len(instructions[x]["pages"]):
            duplicates = Counter(instructions[x]["pages"])
            duplicates = [ pagenum for pagenum, count in duplicates.items() if count > 1 ]
            
            print(f"WARNING: the section of -i '{inputStr}' has page number duplications for: " +
                  f"{duplicates}. If this is intended, ignore this warning. Otherwise, you may need " +
                  "to fix your input instructions, as these pages will be written out as duplicates.")
    
    return instructions

def main():
    usage = """%(prog)s enables simple editing capabilities for PDFs while
    concatenating them into a joined document. Input files can be provided in
    four ways to control the behaviour of what to do with that file.
    1) Provide the file name alone ('filename.pdf'), to concatenate the entire file.
    2) Provide the file name, followed by a colon, then instructions for which pages
    to obtain ('filename.pdf:1-5,7,9,14-end'). As indicated, ranges can be given and
    'end' can be specified in lieu of the full length of the document.
    3) Provide the file name, followed by a colon, then instructions for how to rotate
    the document ('filename.pdf:+90'). The rotation angle should be given as a '+' or '-'
    value to indicate clockwise or counter-clockwise rotation, respectively.
    4) A combination of options 2 and 3 ('filename.pdf:1-5,9-end:-180') to subset
    and rotate.
    
    Note that you can provide the same file more than once, if you want to apply
    rotation to a subset of pages but not another subset of pages. For example:
    'filename1.pdf:1-5 filename1.pdf:6-end:+90' would take the first 5 pages
    of the document as-is, then rotate all pages from 6 onwards 90 degrees clockwise.
    You can also change the ordering of pages like 'filename1.pdf:5-1' to reverse the ordering
    of the first five pages.
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="pdfs",
                   required=True,
                   nargs="+",
                   help="Input PDF file names to join")
    p.add_argument("-o", dest="outputName",
                   required=True,
                   help="Output file name of concatenated document")
    args = p.parse_args()
    validate_args(args)
    
    # Iterate through each input PDF and concatenate as specified
    writer = PyPDF2.PdfWriter()
    for instruction in args.instructions.values():
        with open(instruction["file"], "rb") as fileIn:
            # Interpret rotation angle as number of rotations
            
            reader = PyPDF2.PdfReader(fileIn)
            for pageNum in instruction["pages"]:
                page = reader.pages[pageNum]
                writer.add_page(page)
                
                if instruction["rotation"] != None:
                    writer.pages[-1].rotate(instruction["rotation"])
    
    # Write to a file then save the compiled PDF
    with open(args.outputName, "wb") as fileOut:
        writer.write(fileOut)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()

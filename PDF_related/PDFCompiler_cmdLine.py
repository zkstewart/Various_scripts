#! python3
# pdfCompiler - little program to combine PDF files wholesale. More refined page number
#               selection will be a thing to do with another program

# Import required packages
import PyPDF2, os, argparse

def validate_args(args):
    # Validate input files exist
    for file in args.inputNames:
        if not os.path.isfile(file):
            print("File doesn't exist ({0})".format(file))
            quit()

##### USER INPUT SECTION
usage = """%(prog)s produces a single concatenated PDF from multiple input PDFs
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-i", dest="inputNames", nargs="+",
    help="PDF file names")
p.add_argument("-o", dest="outputName",
    help="output file name")
args = p.parse_args()

# Create the new PDF to add new PDFs into
compiledPDF = PyPDF2.PdfFileWriter()

# Loop to add each PDFs contents into one main PDF
for pdf in args.inputNames:
    pdfFile = open(pdf, 'rb')
    pdfReader = PyPDF2.PdfFileReader(pdfFile)
    # Loop through each page and add to file without rotating
    for pageNum in range(pdfReader.numPages):
        page = pdfReader.getPage(pageNum)
        compiledPDF.addPage(page)

# Write to a file then save the compiled PDF
with open(args.outputName, 'wb') as pdfOutput:
    compiledPDF.write(pdfOutput)
    pdfOutput.close()

print("Program completed successfully!")

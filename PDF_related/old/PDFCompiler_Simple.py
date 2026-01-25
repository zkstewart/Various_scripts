#! python3
# pdfCompiler - little program to combine PDF files wholesale. More refined page number
#               selection will be a thing to do with another program

# Import required packages
import PyPDF2, os

# Loop to give user chance to enter PDF names to combine
print('Enter the prefix of the PDFs you wish to combine. These files should be sortable numerically (since that\'s what Python is gonna do)')

prefix = input()

# Specify the output file name
print('Enter the name you want the compiled PDF to be called. Do not include the extension, and make sure not to use illegal characters (i.e. \\/:?"<>|)')
outputName = input()

###### Core functionality: PDF compilation

# Create the new PDF to add new PDFs into
compiledPDF = PyPDF2.PdfFileWriter()

# Get PDFs
dirFiles = os.listdir('.')
pdfFiles = []
for file in dirFiles:
        if file.startswith(prefix):
                pdfFiles.append(file)

# Loop to add each PDFs contents into one main PDF
for pdf in pdfFiles:
        pdfFile = open(pdf, 'rb')
        pdfReader = PyPDF2.PdfFileReader(pdfFile)
        # Loop through each page and add to file
        for pageNum in range(pdfReader.numPages):
                page = pdfReader.getPage(pageNum)
                compiledPDF.addPage(page)

# Write to a file then save the compiled PDF
pdfOutput = open(outputName + '.pdf', 'wb')
compiledPDF.write(pdfOutput)
pdfOutput.close()

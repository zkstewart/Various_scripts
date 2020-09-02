import PyPDF2, os, re

# Global value for page numbers later
numCount = 1
# Create regex for later command
numOrNot = re.compile(r'\d{1,6}-\d{1,6}')

print('Enter the name of the PDF file you want to split into individual pages. Don\' include the extension')
while True:
        try:
                pdfName = input()
                if os.path.isfile(pdfName + '.pdf') == False:
                        raise Exception
                print('PDF located successfully')
                print('')
                break
        except:
                print('PDF failed to load. If you misspelt the name, try again. If the name is correct but still won\'t work, make sure the file is in the same folder as this script')
                continue

pdfFile = open(pdfName + '.pdf', 'rb')
pdfFileReader = PyPDF2.PdfFileReader(pdfFile)
pdfWriter = PyPDF2.PdfFileWriter()

# Specify the output file name
print('Enter the prefix of the name which you want the output PDF to be called. When saving the new pages, they will be saved in the format "prefix_page1" onwards')
outputFileName = input()
print('')

# Specify which page numbers to grab
print('Enter the page number range to split into separate PDFs in the format "num1-num2" e.g. "4-24", or type "all" to split each page of the file')
print('As a reminder, this might create a lot of files. You might want to run this within a new folder to contain it all')
pageRange = input()
while True:
        if len(numOrNot.findall(pageRange)) == 1:
                pageRange = numOrNot.search(pageRange).group()
                splitPageRange = pageRange.split(sep='-')
                firstPage = int(splitPageRange[0])
                secondPage = int(splitPageRange[1])
                break
        elif pageRange.lower() == 'all':
                pageRange = pageRange.lower()
                break
        elif pageRange != 'all' and len(numOrNot.findall(pageRange)) != 1:
                print('Something went wrong. Make sure to type a range in the format "num1-num2" e.g. "4-24", or type "all" to split each page of the file')
                continue
       
# Loop through the PDF to extract the noted pages
if len(numOrNot.findall(pageRange)) == 1:
    numCount = firstPage
    for pageNum in range(firstPage-1, secondPage):                      # firstPage needs to be minused by one since for loops start at 0, whereas the PDF starts at page 1
        pageObj = pdfFileReader.getPage(pageNum)
        pdfWriter = PyPDF2.PdfFileWriter()
        pdfWriter.addPage(pageObj)
        pdfOutputFile = open(outputFileName + '_page' + str(numCount) + '.pdf', 'wb')
        pdfWriter.write(pdfOutputFile)
        pdfOutputFile.close()
        numCount += 1
elif pageRange == 'all':
    for pageNum in range(pdfFileReader.getNumPages()):
            pageObj = pdfFileReader.getPage(pageNum)
            pdfWriter = PyPDF2.PdfFileWriter()
            pdfWriter.addPage(pageObj)
            pdfOutputFile = open(outputFileName + '_page' + str(numCount) + '.pdf', 'wb')
            pdfWriter.write(pdfOutputFile)
            pdfOutputFile.close()
            numCount += 1

# Close the initial PDF
pdfFile.close()

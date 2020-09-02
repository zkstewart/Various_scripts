import PyPDF2, os, re

# Create regex for later command
numOrNot = re.compile(r'(\d{1,6})-(\d{1,6})')

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
print('Enter the name which you want the output PDF to be called. Do not include the file extension, and make sure not to use illegal characters (i.e. \\/:?"<>|)')
outputFileName = input()
print('')

# Specify which page numbers to grab
print('Enter the page number range to extract in the format "num1-num2" e.g. "4-24". You can specify multiple ranges separated by a comma, as well as specify single pages e.g. "1, 4-6".')
while True:
        firstPage = []
        secondPage = []
        try:
                pageRange = input()
                for entry in pageRange.split(','):
                        entryRange = entry.replace(' ', '')
                        numbers = numOrNot.findall(entryRange)
                        if len(numbers)==1:
                                firstPage.append(int(numbers[0][0]))
                                secondPage.append(int(numbers[0][1]))
                        else:
                                firstPage.append(int(entryRange))
                                secondPage.append(int(entryRange))
                break
        except KeyboardInterrupt:
                quit()
        except:
                print('Something went wrong. Make sure to type numbers in the format "num1-num2" e.g. "4-24", without use of quotation marks')
                
       
# Loop through the PDF to extract the noted pages
for i in range(len(firstPage)):
        for x in range(firstPage[i]-1, secondPage[i]):            # page numbers needs to be -1'd since for loops start at 0, whereas the PDF starts at page 1
                pageObj = pdfFileReader.getPage(x)        
                pdfWriter.addPage(pageObj)


# Make the new PDF
pdfOutputFile = open(outputFileName + '.pdf', 'wb')
pdfWriter.write(pdfOutputFile)
pdfOutputFile.close()
pdfFile.close()

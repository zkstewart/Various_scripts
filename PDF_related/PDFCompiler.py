#! python3
# pdfCompiler - little program to combine PDF files wholesale. More refined page number
#               selection will be a thing to do with another program

# Import required packages
import PyPDF2, os

# Global value
fileNames = []

# Loop to give user chance to enter PDF names to combine
print('Enter the name of the PDFs you wish to combine. When you have selected all the PDFs you want to combine, type \'go\' and press ENTER')
print('Note that the order in which you enter them here will determine the order in which they appear in the combined PDF')

while True:
        try:
                name = input()
                if name.lower() == 'go':
                        break
                if os.path.isfile(name + '.pdf') == False:
                        raise Exception
                fileNames += [name + '.pdf']
                print(name + '.pdf' + ' located successfully')
                print('')
                print('Enter the next PDFs name, or type \'go\' to compile the files you have specified')
                continue
        except:
                print('PDF failed to load. If you misspelt the name, try again. If the name is correct but still won\'t work, make sure the file is in the same folder as this script')
                continue

# Specify the output file name
print('Enter the name you want the compiled PDF to be called. Do not include the extension, and make sure not to use illegal characters (i.e. \\/:?"<>|)')
outputName = input()

# Offer chance to do some ognTSM rotations on the files
print('Do any of the files need to be rotated? Type \'yes\' or \'no\' then press ENTER')
print('Note that if you say yes, you will need to specify the degree of rotation for each file. You can leave 0 for no rotation which will be explained later if you choose this option')

while True:
        yesno = input()
        if yesno != 'yes' and yesno != 'no':
                print('Was that a yes, or a no? Type it in then press ENTER')
        else:
                break

###### Core functionality: PDF compilation

# Create the new PDF to add new PDFs into
compiledPDF = PyPDF2.PdfFileWriter()

# Loop to add each PDFs contents into one main PDF
if yesno.lower() == 'yes':
        for pdf in fileNames:
                print('How many times do you want to rotate ' + pdf + ' to the right/clockwise 90 degrees? If you do not wish to rotate this file, enter 0')
                while True:
                        try:
                                pdfRotateTimes = int(input())
                        except:
                                print('You appear to have entered something that isn\'t a number. Try again')
                                continue
                        if pdfRotateTimes <= 3:
                                break
                        if pdfRotateTimes > 3 and pdfRotateTimes < 100:
                                print('I dont think you actually want to rotate the file ' + str(pdfRotateTimes) + ' times. Enter a number from 0 to 3')
                                continue
                        elif pdfRotateTimes >= 100:
                                print('Rotating it this many times would be extremely painful ... for you')
                                print('(p.s. Enter a number from 0 to 3)')
                                continue

                pdfFile = open(pdf, 'rb')
                pdfReader = PyPDF2.PdfFileReader(pdfFile)
                # Loop through each page and rotate before adding to file
                for pageNum in range(pdfReader.numPages):
                        page = pdfReader.getPage(pageNum)
                        if pdfRotateTimes == 1:
                                page.rotateClockwise(90)
                        elif pdfRotateTimes == 2:
                                page.rotateClockwise(180)
                        elif pdfRotateTimes == 3:
                                page.rotateClockwise(270)
                        compiledPDF.addPage(page)
                
if yesno.lower() == 'no':
        for pdf in fileNames:
                pdfFile = open(pdf, 'rb')
                pdfReader = PyPDF2.PdfFileReader(pdfFile)
                # Loop through each page and add to file without rotating
                for pageNum in range(pdfReader.numPages):
                        page = pdfReader.getPage(pageNum)
                        compiledPDF.addPage(page)

# Write to a file then save the compiled PDF
pdfOutput = open(outputName + '.pdf', 'wb')
compiledPDF.write(pdfOutput)
pdfOutput.close()

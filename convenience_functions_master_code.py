#! python3

# convenience_functions_master_code.py

# Combination of multiple functions as a go-to script for convenience
# functions (i.e., scripts that run as a double-click program rather than
# on the command-line).

# Load packages for main
import os

# Define functions
def lineremover(functionList):
        # Set up
        import pyperclip
        # Core loop
        while True:
                print('Press any button to continue the loop (\\n remover). Enter \'help\' or \'exit\' if wanted.')
                textInput = input()
                # Help and exit
                if textInput.lower() == 'help':
                        detailed_help(functionList)
                elif textInput.lower() == 'exit':
                        return 'success'
                # Main function
                else:
                        plainText = pyperclip.paste()
                        noReturn = plainText.replace('\r', '')  # Handle Windows-specific formatting
                        strippedText = noReturn.replace('\n', '')
                        pyperclip.copy(strippedText)
                        print('Copied to clipboard')

def hyphenremover(functionList):
        # Set up
        import pyperclip
        # Core loop
        while True:
                print('Press any button to continue the loop (- remover). Enter \'help\' or \'exit\' if wanted.')
                textInput = input()
                # Help and exit
                if textInput.lower() == 'help':
                        detailed_help(functionList)
                elif textInput.lower() == 'exit':
                        return 'success'
                # Main function
                else:
                        plainText = pyperclip.paste()
                        strippedText = plainText.replace('-', '')
                        pyperclip.copy(strippedText)
                        print('Copied to clipboard')

def expasyjoiner(functionList):
        # Set up
        import pyperclip
        # Core loop
        while True:
                print('Press any button to continue the loop (expasy sequence joiner). Enter \'help\' or \'exit\' if wanted.')
                textInput = input()
                # Help and exit
                if textInput.lower() == 'help':
                        detailed_help(functionList)
                elif textInput.lower() == 'exit':
                        return 'success'
                # Main function
                else:
                        plainText = pyperclip.paste()
                        expasySeq = plainText.split(sep=' ')
                        outSeq = ''
                        for amino in expasySeq:
                                if amino == 'Met':
                                        outSeq += 'M'
                                elif amino == 'Stop':
                                        outSeq += '*'
                                else:
                                        outSeq += amino
                        pyperclip.copy(outSeq)
                        print('Copied to clipboard')

# Define text for functions
def detailed_help(functionList):
        import textwrap
        lineremover = r'''
        The _lineremover_ function will remove '\r' and '\n' characters from
        text in the clipboard.
        '''
        hyphenremover = r'''
        The _hyphenremover_ function will remove '-' characters from
        text in the clipboard.
        '''
        expasyjoiner = r'''
        The _expasyjoiner_ function will parse Expasy amino acid characters
        in the clipboard and provide the raw protein sequence to the clipboard.
        '''
        printList = str(functionList).replace("'", "")
        printList = eval(printList)
        for entry in printList:
                entry = textwrap.dedent(entry)
                entry = entry.strip('\n').replace('\n', ' ')
                for line in textwrap.wrap(entry, width=50):
                        print(line)
                print('')

# Define main loop
def core_loop(functionList):
        while True:
                print('Specify function to run. Type \'help\' at any time to see list of functions. Type \'exit\' or ctrl+c to quit.')
                # Input and ctrl+c check
                try:
                        cmdInput = input()
                except KeyboardInterrupt:
                        print('Received Keyboard Interrupt. Exiting now.')
                        quit()
                # Help and exit
                if cmdInput.lower() == 'help':
                        detailed_help(functionList)
                elif cmdInput.lower() == 'exit':
                        quit()
                # Enact function
                elif cmdInput.lower() in functionList:
                        while True:
                                try:
                                        eval(cmdInput + '(functionList)')
                                        break
                                except KeyboardInterrupt:
                                        quit()
                                except:
                                        print('Something went wrong. Try again.')
                else:
                        print('I didn\'t recognise that function. Make sure you typed it correctly, or type \'help\' to see a list of available functions.\n')

# Start program
functionList = ['lineremover', 'hyphenremover', 'expasyjoiner']          
core_loop(functionList)

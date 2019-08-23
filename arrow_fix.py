#! python3
# fasta_arrow_fix
# Remove '|arrow' seqid suffix and make all sequence characters uppercase

import os, argparse

# Define functions
def arrow_fix(fasta_file, output_file_name):
        with open(fasta_file, 'r') as file_in, open(output_file_name, 'w') as file_out:
                for line in file_in:
                        if line.startswith('>'):
                                file_out.write(line.split('|')[0] + '\n')
                        else:
                                file_out.write(line.upper())

# Main call
def main():
        def validate_args(args):
                # Validate input file locations
                if args.input == None:
                        print('No fasta argument was provided. Fix this and try again.')
                        quit()
                if not os.path.isfile(args.input):
                        print('I am unable to locate the fasta file (' + args.input + ')')
                        print('Make sure you\'ve typed the file name or location correctly and try again.')
                        quit()
                # Handle file overwrites
                if args.output == None:
                        print('No output file name argument was provided. Fix this and try again.')
                        quit()
                if os.path.isfile(args.output):
                        print('There is already a file named "' + args.output + '".')
                        print('Either specify a new file name or move/delete this older file again try again.')
                        quit()
        ##### USER INPUT SECTION
        usage = """%(prog)s reads in a provided fasta file and undoes the changes arrow induces
        """
        p = argparse.ArgumentParser(description=usage)
        p.add_argument("-i", dest="input",
                          help="fasta file name")
        p.add_argument("-o", dest="output",
                          help="output file name")
        args = p.parse_args()
        validate_args(args)
        # Run arrow fix function
        arrow_fix(args.input, args.output)

if __name__ == '__main__':
        main()

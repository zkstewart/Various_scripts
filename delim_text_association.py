#! python3
# delim_text_association.py
# Extensible code for the association of column values in multiple text files
# to one another.

# Imports
import os, argparse

# Definition of class types
class TabDelimTable:
        def __init__(self, file_location, file_type):
                '''key_only should be provided for the first file which is
                a single-column file; key_value should be provided for all subsequent
                files which should be two-column files of key:value pairs'''
                self.file_location = file_location
                if file_type == 'key_value':
                        self.read_key_value()
                elif file_type == 'key_only':
                        self.read_key_only()
        
        def read_key_value(self, delimiter='\t'):
                self.kv_dict = {}
                first_line = True
                with open(self.file_location, 'r') as file_in:
                        for l in file_in:
                                sl = l.rstrip('\r\n').split(delimiter)
                                if first_line == True:
                                        self.header = sl[1]
                                        first_line = False
                                        continue
                                self.kv_dict[sl[0]] = sl[1]
        
        def read_key_only(self):
                self.k_list = []
                first_line = True
                with open(self.file_location, 'r') as file_in:
                        for l in file_in:
                                if first_line == True:
                                        self.header = l.rstrip('\r\n')
                                        first_line = False
                                        continue
                                self.k_list.append(l.rstrip('\r\n'))

# Definition of script functions
def associate_klist_to_kvdict(k_list, kv_dict, delimiter, second_file_name, blank_character):
        k_skip = {}
        v_list = []
        if blank_character == True:
                blank = ''
        else:
                blank = '.'
        for k in k_list:
                # Split keys into individual entries
                if delimiter in k:
                        k = k.split(delimiter)
                else:
                        k = [k]
                # Generate a corresponding string of values associated with each key
                tmp_v_list = []
                for indiv_k in k:
                        if indiv_k in kv_dict:
                                tmp_v_list.append(kv_dict[indiv_k])
                        else:
                                if indiv_k not in k_skip:
                                        print('"' + indiv_k + '" was not found in "' + second_file_name +
                                              '"\nThis might be a problem.')
                                        k_skip.setdefault(indiv_k)
                                tmp_v_list.append(blank)
                v_string = delimiter.join(tmp_v_list)
                v_list.append(v_string)
        return v_list

def main():
        def validate_args(args):
                # Disallow blank arguments
                if args.file_locations == None:
                        print('file_locations argument was not specified')
                        quit()
                elif args.output_file_name == None:
                        print('output_file_name argument was not specified')
                        quit()
                elif args.delimiter == None:
                        print('delimiter argument was not specified')
                        quit()
                elif args.output_delimiter == None:
                        print('output_delimiter argument was not specified')
                        quit()
                # Validate inputs
                for file_location in args.file_locations:
                        if not os.path.isfile(file_location):
                                print('"' + file_location + '" is not recognised',
                                      'as a file which exists. Make sure the file name',
                                      'was entered correctly; if the file is not in the same'
                                      'directory as this Python script file, make sure to specify'
                                      'the full path to the file in question.')
                                quit()
                # Validate output
                if os.path.isfile(args.output_file_name):
                        print('"' + args.output_file_name + '" refers to a file',
                              'that already exists. To prevent unfortunate file overwrites,',
                              'please specify a new file name or location to save output',
                              'to.')
                        quit()
                # Fix potential user errors with delimiters
                '''Because of how the command-line help is written, it's probable
                that a user might find themself trying "-od \t" to produce a tab character
                but, because it is not escaped by default by argparse, it will produce
                rubbish output. We can fix it here quickly and prevent issues
                '''
                if args.delimiter == '\\t':
                        args.delimiter = '\t'
                if args.output_delimiter == '\\t':
                        args.output_delimiter = '\t'
                return args
        
        ## Argument parsing
        usage = """%(prog)s is a flexible program for associating data in columns
across multiple delimited text files. "Key"s obtained from a single file will
be associated with other files containing "Key:Value" pairs. Note the below
prerequisites for successful program execution:
        -All files are assumed to contain headers
        -The first file should have only one column; all other files should have
        two columns
        """
        
        p = argparse.ArgumentParser(description=usage, formatter_class=argparse.RawDescriptionHelpFormatter)
        p.add_argument("-i", dest="file_locations", nargs="+",
                       help="""Input text files for column association; the first
                       file specified should contain a single column of "key" values
                       to be associated to "key:value" pairs in the other file(s)""")
        p.add_argument("-o", dest="output_file_name",
                       help="""Output file will be produced at the specified
                       location""")
        p.add_argument("-d", dest="delimiter",
                       help="""Specify the delimiter used in text files for parsing (e.g.,
                       ';' or ',' sans apostrophe""")
        p.add_argument("-od", dest="output_delimiter",
                       help="""Specify the delimiter to use for output file writing; this
                       should be different to that used for parsing. (e.g., '\\t')""")
        p.add_argument("-b", dest="blank_character", action='store_true',
                       help="By default missing key:value pairs will be marked as '.'; specify -b argument to have no character inserted for marked values", default=False)
        args = p.parse_args()
        args = validate_args(args)
        
        # Read in primary file for keys
        keys = TabDelimTable(args.file_locations[0], 'key_only')
        # Associate keys to key:value files
        v_list_storage = [[keys.header] + keys.k_list]
        for i in range(1, len(args.file_locations)):
                key_value = TabDelimTable(args.file_locations[i], 'key_value')
                v_list = associate_klist_to_kvdict(keys.k_list, key_value.kv_dict, args.delimiter, args.file_locations[i], args.blank_character)
                v_list_storage.append([key_value.header] + v_list)
        # Write output file
        with open(args.output_file_name, 'w') as file_out:
                for i in range(len(v_list_storage[0])):
                        for x in range(len(v_list_storage)):
                                if x != len(v_list_storage) - 1:
                                        file_out.write(v_list_storage[x][i] + args.output_delimiter)
                                else:
                                        file_out.write(v_list_storage[x][i] + '\n')

if __name__ == '__main__':
        main()
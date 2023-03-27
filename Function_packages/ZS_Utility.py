#! python3
# ZS_Utility.py
# Various statics that need to be shared amongst the other
# function packages.

import os, codecs

def tmp_file_name_gen(prefix, suffix):
    '''
    Function for use when creating temporary files.
    Params:
        prefix -- a string for a file prefix e.g., "tmp"
        suffix -- a string for a file suffix e.g., "fasta". Note that we don't
                    use a "." in this, since it's inserted between prefix and suffix
                    automatically.
    Returns:
        tmpName -- a string for a file name which does not exist in the current dir.
    '''
    ongoingCount = 1
    while True:
        if not os.path.isfile("{0}.{1}".format(prefix, suffix)):
            return "{0}.{1}".format(prefix, suffix)
        elif os.path.isfile("{0}.{1}.{2}".format(prefix, ongoingCount, suffix)):
            ongoingCount += 1
        else:
            return "{0}.{1}.{2}".format(prefix, ongoingCount, suffix)

def get_codec(fileName):
    try:
        f = codecs.open(fileName, encoding='utf-8', errors='strict')
        for line in f:
            pass
        return "utf-8"
    except:
        try:
            f = codecs.open(fileName, encoding='utf-16', errors='strict')
            for line in f:
                pass
            return "utf-16"
        except UnicodeDecodeError:
            print(f"BlastIO class can't tell what codec '{fileName}' is!!")

if __name__ == "__main__":
    pass

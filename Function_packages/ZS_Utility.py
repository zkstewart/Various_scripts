#! python3
# ZS_Utility.py
# Various statics that need to be shared amongst the other
# function packages.

import os, codecs, platform, re

def tmp_file_name_gen(prefix, suffix):
    '''
    Function for use when creating temporary files.
    Parameters:
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
            print(f"Can't tell what codec '{fileName}' is!!")

def base_subprocess_cmd(exe):
    '''
    Helper function to begin formatting a cmd list for use with subprocess
    that saves a bit of code repetition and handles platform dependence of
    behaviour i.e., with Windows, will setup a cmd that instructs WSL to
    execute. Saves a bit of code repetion is all.
    '''
    if platform.system() == "Windows":
        cmd = ["wsl", "~", "-e", convert_windows_to_wsl_path(exe)]
    else:
        cmd = [exe]
    return cmd

def convert_windows_to_wsl_path(windowsPath):
    '''    
    Provides simple functionality to infer the WSL path from
    a windows path, provided as a string.
    
    Parameters:
        windowsPath -- a string indicating the full path to a file
                        or directory of interest; this MUST include
                        the root character e.g., 'D:\\' or 'C:\\'
    Returns:
        wslPath -- a string indicating the inferred full path to the
                    given file or directory using WSL formatting
    '''
    # Get abspath if not already done
    windowsPath = os.path.abspath(windowsPath)
    
    # Check that the path is something we can work with
    driveRegex = re.compile(r"^([A-Za-z]{1}):\\")
    assert driveRegex.match(windowsPath) != None, \
        f"'{windowsPath}' is not recognised as a full, root drive inclusive path"
    
    # assert os.path.exists(windowsPath), \
    #     f"'{windowsPath}' is not recognised as an existing path"
    "We don't actually need this in this function; sometimes the file shouldn't exist yet"
    
    # If it is, convert it
    driveLetter = driveRegex.match(windowsPath).group(1)
    wslPath = "/{0}".format("/".join(
        [
            "mnt",
            driveLetter.lower(),
            *windowsPath.split("\\")[1:]
        ]
    ))
    
    return wslPath

if __name__ == "__main__":
    pass

#! python3
# ZS_Utility.py
# Various statics that need to be shared amongst the other
# function packages.

import os, codecs, platform, re, subprocess, shutil

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

def wsl_which(program):
    '''
    A function to expand upon shutil.which to work with WSL. Emulates its behaviour
    by directly reaching into the WSL shell and calling the Unix 'which' command.
    If this program is running in Unix, then shutil.which will be called.
    '''
    if platform.system() != 'Windows':
        return shutil.which(program)
    else:
        cmd = ["wsl", "~", "-e", "which", program]
        run_wsl_which = subprocess.Popen(cmd, shell = True,
                                         stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        whichout, whicherr = run_wsl_which.communicate()
        if whicherr.decode("utf-8") != "":
            return None
        elif whichout.decode("utf-8").rstrip("\r\n ") != "":
            return whichout.decode("utf-8").rstrip("\r\n ")
        else:
            raise Exception(("wsl_which encountered an unhandled situation; have a look " +
                             "at the stdout '" + whichout.decode("utf-8") + "' and stderr '" + 
                             whicherr.decode("utf-8") + "' to make sense of this."))

def wsl_exists(program):
    '''
    A function to expand upon os.path.exists to work with WSL. Emulates its behaviour
    by directly reaching into the WSL shell and calling the Unix 'ls' command.
    If this program is running in Unix, then os.path.exists will be called.
    
    Parameters:
        program -- a string indicating the full path to the program of interest.
        isFolder -- OPTIONAL; a boolean indicating whether the program is a folder
                    or not. Default is False (i.e., a file or program).
    '''
    if platform.system() != 'Windows':
        return os.path.exists(program)
    else:
        cmd = ["wsl", "~", "-e", "ls", program]
        run_wsl_exists = subprocess.Popen(cmd, shell = True,
                                         stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        existsout, existserr = run_wsl_exists.communicate()
        if "No such file or directory" in existserr.decode("utf-8"):
            return False
        elif existsout.decode("utf-8").rstrip("\r\n") != "":
            return True
        else:
            raise Exception(("wsl_exists encountered an unhandled situation; have a look " +
                             "at the stdout '" + existsout.decode("utf-8") + "' and stderr '" + 
                             existserr.decode("utf-8") + "' to make sense of this."))

def wsl_isfile(program):
    '''
    A function to expand upon os.path.isfile to work with WSL. Emulates its behaviour
    by directly reaching into the WSL shell and calling the Unix 'ls' command.
    If this program is running in Unix, then os.path.isfile will be called.
    
    Parameters:
        program -- a string indicating the full path to the program of interest.
    '''
    if platform.system() != 'Windows':
        return os.path.isfile(program)
    else:
        cmd = ["wsl", "~", "-e", "ls", program]
        run_wsl_exists = subprocess.Popen(cmd, shell = True,
                                         stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        existsout, existserr = run_wsl_exists.communicate()
        if "No such file or directory" in existserr.decode("utf-8"):
            return False
        elif existsout.decode("utf-8").rstrip("\r\n ").endswith(program):
            return True
        else:
            return False

def wsl_isdir(program):
    '''
    A function to expand upon os.path.isdir to work with WSL. Emulates its behaviour
    by directly reaching into the WSL shell and calling the Unix 'ls' command.
    If this program is running in Unix, then os.path.isdir will be called.
    
    Parameters:
        program -- a string indicating the full path to the program of interest.
    '''
    if platform.system() != 'Windows':
        return os.path.isdir(program)
    else:
        cmd = ["wsl", "~", "-e", "ls", program]
        run_wsl_exists = subprocess.Popen(cmd, shell = True,
                                         stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        existsout, existserr = run_wsl_exists.communicate()
        if "No such file or directory" in existserr.decode("utf-8"):
            return False
        elif existsout.decode("utf-8").rstrip("\r\n ").endswith(program):
            return False
        else:
            return True

def convert_to_wsl_if_not_unix(fileLocation):
    if platform.system() != 'Windows':
        return fileLocation
    else:
        if os.path.basename(fileLocation) == fileLocation: # WSL won't interpret relative paths correctly
            return convert_windows_to_wsl_path(os.path.abspath(fileLocation))
        elif is_unix_path(fileLocation):
            return fileLocation
        else:
            return convert_windows_to_wsl_path(fileLocation)

def base_subprocess_cmd(exe):
    '''
    Helper function to begin formatting a cmd list for use with subprocess
    that saves a bit of code repetition and handles platform dependence of
    behaviour i.e., with Windows, will setup a cmd that instructs WSL to
    execute. Saves a bit of code repetion is all.
    '''
    if platform.system() == "Windows":
        cmd = ["wsl", "~", "-e", convert_to_wsl_if_not_unix(exe)]
    else:
        cmd = [exe]
    return cmd

def is_unix_path(path):
    '''
    Helper function to check if a path is a unix path. Relies on a simple
    heuristic of checking if the path does not contain backslashes.
    
    Parameters:
        path -- a string indicating the path to check.
    Returns:
        isUnix -- a boolean indicating whether the path is a unix path.
    '''
    return path.count("\\") == 0

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

def convert_wsl_to_windows_path(wslPath):
    '''    
    Provides simple functionality to infer the WSL path from
    a windows path, provided as a string.
    
    Parameters:
        wslPath -- a string indicating the WSL path to a file
                   or directory of interest; this MUST include
                   the mnt directory e.g., "/mnt/d"
    Returns:
        windowsPath -- a string indicating the inferred full path to the
                       given file or directory using windows formatting
    '''
    # Check that the path is something we can work with
    driveRegex = re.compile(r"/mnt/(\w)/")
    assert driveRegex.match(wslPath) != None, \
        f"'{wslPath}' is not recognised as a properly formatted WSL path"
    
    # If it is, convert it
    driveLetter = driveRegex.match(wslPath).group(1)
    windowsPath = "{0}:\\{1}".format(driveLetter.upper(),
                               "\\".join(wslPath.split("/")[3:])
                               )
    
    return windowsPath

if __name__ == "__main__":
    pass

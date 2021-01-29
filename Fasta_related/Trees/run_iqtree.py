#! python3
# run_iqtree.py
# Wrapper script to make it easier to run batch IQTree analysis

import argparse, os, subprocess, shutil

# Validate arguments
def validate_args(args):
    # iqt2suffixes is used to obtain just the files in the directory that (should) correspond to fasta files i.e., those not made by iqtree2 or this script
    iqt2Suffixes = (".bionj", ".ckp.gz", ".contree", ".iqtree", ".log", ".mldist", ".model.gz", ".nwk", ".splits.nex", ".treefile", ".uniqueseq.phy")


    # Validate input argument
    if not os.path.isfile(args.input) and not os.path.isdir(args.input):
            print('I am unable to locate a file or directory as provided (' + args.input + ')')
            print('Make sure you\'ve typed the file name or location correctly and try again.')
            quit()
    # Format files list from user-provided input
    if os.path.isfile(args.input):
        args.files = [os.path.abspath(args.input)]
    else:
        args.files = []
        for file in os.listdir(args.input):
            if file.endswith(iqt2Suffixes): # Necessary for job skipping
                continue
            args.files.append(os.path.join(os.path.abspath(args.input), file)) # Make sure the full path is available just in case
    # Validate IQTree location
    if args.iqtreeDir != '':
            if not os.path.isfile(os.path.join(args.iqtreeDir, 'iqtree2')) and not os.path.isfile(os.path.join(args.iqtreeDir, 'iqtree2.exe')):
                    print('I cannot find "iqtree2" or "iqtree2.exe" at the location provided (' + args.iqtreeDir + ')')
                    quit()
    else:
            print('You haven\'t specified a location for the IQTree executable. If this is in your PATH that\'s OK. Otherwise, I\'ll likely crash soon.')
    # Validate numeric arguments
    if args.cpus < 1:
            print('CPUs cannot be less than 1. Specify any number >= 1 and try again.')
            quit()
    return args

# IQTree call
def execute_iqtree(iqtreeDir, inputFile, cpus, bootstraps):
        # Format command
        if bootstraps > 0:
            cmd = "{0} -s {1} -T {2} -B {3}".format(os.path.join(iqtreeDir, 'iqtree2'), inputFile, cpus, bootstraps)
        else:
            cmd = "{0} -s {1} -T {2}".format(os.path.join(iqtreeDir, 'iqtree2'), inputFile, cpus)
        print("# " + cmd)
        # Run command
        exe_iqtree = subprocess.Popen(cmd, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        iqtree_out, iqtree_err = exe_iqtree.communicate()
        if iqtree_err.decode("utf-8") != '':
                raise Exception('IQTree error text below\n' + iqtree_err.decode("utf-8"))

def iqtree_to_newark(inputFile):
    shutil.copyfile("{0}.treefile".format(inputFile), "{0}.nwk".format(inputFile))

# Main call
def main():
        #### USER INPUT SECTION
        usage = """Wrapper script to perform MMseqs2 search. Provide the arguments below."""
        # Required
        p = argparse.ArgumentParser(description=usage)
        p.add_argument("-i", dest="input", type = str,
                          help="Specify a file or directory containing ONLY files to be used as input")
        p.add_argument("-c", dest="cpus", type = int, default = 1,
                          help="Specify the number of CPUs/threads to provide as an argument (default == 1)")
        p.add_argument("-b", dest="bootstraps", type = int, default = 1000,
                          help="Specify the number of bootstraps to provide as an argument (default == 1000; enter 0 for no bootstrapping)")
        p.add_argument("-e", dest="iqtreeDir", type = str, default = "",
                          help="Specify the directory where the IQTree executable is located; if it is accessible from your PATH, you can leave this blank")
        p.add_argument("-s", dest="skipCompleted", action = "store_true", default = False,
                          help="Optionally skip running IQTree on FASTAs when .nwk files already exist for that file (enables job resumption)")
        args = p.parse_args()
        args = validate_args(args)

        # Call IQTree
        for file in args.files:
            alreadyCompleted = False
            if os.path.isfile("{0}.nwk".format(file)):
                alreadyCompleted = True
            
            if args.skipCompleted == True and alreadyCompleted == True:
                print("Skipping file \"{0}\"...".format(file))
            else:
                execute_iqtree(args.iqtreeDir, file, args.cpus, args.bootstraps)
                iqtree_to_newark(file)

        # Done!
        print('Program completed successfully!')

if __name__ == '__main__':
        main()

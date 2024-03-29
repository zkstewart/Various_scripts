#! python3
# run_mmseqs2.py
# Wrapper script to make it easier to run a search with MMseqs2

import argparse, os

# Various functions for program operations

## Validate arguments
def validate_args(args):
        # Validate query file location
        if not os.path.isfile(args.query):
                print('I am unable to locate the query file (' + args.query + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        args.querydir = os.path.dirname(os.path.abspath(args.query))    # Provides an easy way to handle directories without needing to specify these as separate arguments
        args.query = os.path.basename(args.query)                       # Finalises the separation
        # Validate target file location (if relevant)
        if args.profile_terminate == False:
                if not os.path.isfile(args.target):
                        print('I am unable to locate the target file (' + args.target + ')')
                        print('Make sure you\'ve typed the file name or location correctly and try again.')
                        quit()
                args.targetdir = os.path.dirname(os.path.abspath(args.target))
                args.target = os.path.basename(args.target)
        # Validate output file location
        args.outputdir = os.path.dirname(os.path.abspath(args.output))
        args.output = os.path.basename(args.output)
        if not os.path.isdir(args.outputdir):
                print('The specified output directory does not exist (i.e., "' + args.outputdir + ')".')
                print('Create this directory first and then try again.')
                quit()
        # Validate temporary directory
        if args.tmpDir == "":
            print("-tmp needs to be given a value with at least one character")
            quit()
        elif os.path.basename(args.tmpDir) != args.tmpDir:
            print("-tmp needs to be a basename only; it must be written to the current working directory")
            quit()
        # Validate mmseqs location
        if args.mmseqs2dir != '':
                if not os.path.isfile(os.path.join(args.mmseqs2dir, 'mmseqs')) and not os.path.isfile(os.path.join(args.mmseqs2dir, 'mmseqs.exe')):
                        print('I cannot find "mmseqs" or "mmseqs.exe" at the location provided (' + args.mmseqs2dir + ')')
                        quit()
        else:
                print('You haven\'t specified a location for the MMseqs executable. If this is in your PATH that\'s OK. Otherwise, I\'ll likely crash soon.')
        # Validate numeric arguments
        if args.evalue < 0:
                print('E-value cannot be a negative value. Specify any number >= 0 and try again.')
                quit()
        elif args.threads < 1:
                print('CPUs cannot be less than 1. Specify any number >= 1 and try again.')
                quit()
        elif args.num_iterations < 1:
                print('Num_iterations cannot be less than 1. Specify any number >= 1 and try again.')
                quit()
        elif args.alt_ali < 0:
                print('Alt_ali cannot be a negative value. Specify any number >= 0 and try again.')
                quit()
        # Handle conflicting arguments
        if args.alt_ali > 0 and args.num_iterations > 1:
                print('MMseqs doesn\'t support alternative alignments for profile searches. You\'re going to need to choose one or the other - there\'s no way around this.')
                quit()
        return args

## Log file related functions
def temp_file_name_gen(prefix):
        import os
        ongoingCount = 1
        while True:
                if not os.path.isfile(prefix):
                        return prefix
                elif os.path.isfile(prefix + str(ongoingCount)):
                        ongoingCount += 1
                else:
                        return prefix + str(ongoingCount)

def log_update(logName, text):
        with open(logName, 'a') as logFile:
                logFile.write(text + '\n')

## MMseqs2 related functions
def makemms2db(mmseqs2dir, query):
        import os, subprocess
        # Format command
        dbname = query + '_seqDB'
        cmd = os.path.join(mmseqs2dir, 'mmseqs') + ' createdb "' + query + '" "' + dbname + '"'
        # Query DB generation
        print("# " + cmd)
        run_makedb = subprocess.Popen(cmd, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        makedbout, makedberr = run_makedb.communicate()
        if makedberr.decode("utf-8") != '':
                raise Exception('Make MMseqs2 query db error text below\n' + makedberr.decode("utf-8"))

def makemms2profile(mmseqs2dir, query, tmpdir, iterations=4, sensitivity=7.5):
        import os, subprocess
        # Format command pipeline
        cmd1 = '{0} createdb "{1}" "{2}"'.format(os.path.join(mmseqs2dir, 'mmseqs'), query, query + '_seqDB')
        cmd2 = '{0} search "{1}" "{1}" "{2}" "{3}" --num-iterations {4} -s {5}'.format(os.path.join(mmseqs2dir, 'mmseqs'), query + '_seqDB', query + '_selfsearchDB', tmpdir, iterations, sensitivity)
        cmd3 = '{0} result2profile "{1}" "{1}" "{2}" "{3}"'.format(os.path.join(mmseqs2dir, 'mmseqs'), query + '_seqDB', query + '_selfsearchDB', query + '_profileDB')
        # Run profile creation pipeline
        print("# " + cmd1)
        run_makedb = subprocess.Popen(cmd1, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        makedbout, makedberr = run_makedb.communicate()
        if makedberr.decode("utf-8") != '':
                raise Exception('MMseqs2 createdb error text below\n' + makedberr.decode("utf-8"))

        print("# " + cmd2)
        run_search = subprocess.Popen(cmd2, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        searchout, searcherr = run_search.communicate()
        if searcherr.decode("utf-8") != '':
                errList = searcherr.decode("utf-8").rstrip("\n").split("\n")
                for err in errList:
                        if not err.startswith("posix_madvise"): # mmseqs2 version e1a1c raises an error I don't understand, but it doesn't seem to be an exception
                                raise Exception('MMseqs2 search error text below\n' + searcherr.decode("utf-8"))

        print("# " + cmd3)
        run_r2s = subprocess.Popen(cmd3, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        r2sout, r2serr = run_r2s.communicate()
        if r2serr.decode("utf-8") != '':
                raise Exception('MMseqs2 result2profile error text below\n' + r2serr.decode("utf-8"))

def indexmms2(mmseqs2dir, query, tmpdir, threads):
        import os, subprocess
        # Format command
        dbname = query + '_seqDB'
        cmd = os.path.join(mmseqs2dir, 'mmseqs') + ' createindex "' + dbname + '" "' + tmpdir + '" --threads ' + str(threads)
        # Run query index
        print("# " + cmd)
        run_index = subprocess.Popen(cmd, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        indexout, indexerr = run_index.communicate()
        if indexerr.decode("utf-8") != '':
                raise Exception('Indexing MMseqs2 query db error text below\n' + indexerr.decode("utf-8"))
        # Create output flag
        with open(dbname + ".idxComplete", "w") as fileOut:
            fileOut.write("Index completed")

def runmms2(mmseqs2dir, queryDB, targetDB, tmpdir, searchName, searchType, params):
        '''
        params = [evalue, threads, num_iterations, sensitivity, alt_ali]
        '''
        assert searchType in ["auto", "blastn", "tblastx"]
        searchType = "0" if searchType == "auto" else "2" if searchType == "blastn" else "3"
        
        import os, subprocess
        # Format command
        cmd = '{} search "{}" "{}" "{}" "{}" --search-type {} -e {} --threads {} --num-iterations {} -s {} --alt-ali {}'.format(os.path.join(mmseqs2dir, 'mmseqs'), queryDB, targetDB, searchName, tmpdir, searchType, *params)
        print("# " + cmd)
        # Run query
        run_mms2 = subprocess.Popen(cmd, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        mms2out, mms2err = run_mms2.communicate()
        if mms2err.decode("utf-8") != '':
                errList = mms2err.decode("utf-8").rstrip("\n").split("\n")
                for err in errList:
                        if not err.startswith("posix_madvise"): # mmseqs2 version e1a1c raises an error I don't understand, but it doesn't seem to be an exception
                                raise Exception('MMseqs2 search error text below\n' + mms2err.decode("utf-8"))

def mms2tab(mmseqs2dir, queryDB, targetDB, searchName, threads):
        import os, subprocess
        # Create tab-delim BLAST-like output
        cmd = os.path.join(mmseqs2dir, 'mmseqs') + ' convertalis "' + queryDB + '" "' + targetDB + '" "' + searchName + '" "' + searchName + '.m8" --threads ' + str(threads)
        print("# " + cmd)
        run_mms2 = subprocess.Popen(cmd, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        mms2out, mms2err = run_mms2.communicate()
        if mms2err.decode("utf-8") != '':
                raise Exception('MMseqs2 tabular output generation error text below\n' + mms2err.decode("utf-8"))

def mms2sort(searchName):
        from itertools import groupby
        # Get file names
        outName = searchName.rsplit('.', maxsplit=1)[0] + '_sorted.' + searchName.rsplit('.', maxsplit=1)[1]
        # Parse file
        grouper = lambda x: x.split('\t')[0]
        with open(searchName, 'r') as fileIn, open(outName, 'w') as fileOut:
                for key, group in groupby(fileIn, grouper):
                        group = list(group)
                        # Sort group if relevant
                        if len(group) > 1:
                                for i in range(len(group)):
                                        group[i] = group[i].rstrip('\n').split('\t')
                                group.sort(key = lambda x: (float(x[10]),-float(x[11])))
                                for i in range(len(group)):
                                        group[i] = '\t'.join(group[i])
                        else:
                                group[0] = group[0].rstrip('\n')
                        # Put in output
                        for entry in group:
                                fileOut.write(entry + '\n')

def mms2profiletab(queryFasta, searchName, outName):
        # Determine profile name
        profileName = queryFasta.rsplit(".", maxsplit=1)[0]
        # Write modified file
        with open(searchName, 'r') as fileIn, open(outName, 'w') as fileOut:
                for line in fileIn:
                        sl = line.split("\t")
                        sl[0] = profileName
                        fileOut.write("\t".join(sl))

def mms2_index_exists(fileNamePrefix, directory):
        indexName = fileNamePrefix + '.idxComplete'
        if indexName in os.listdir(directory):
                return True
        else:
                return False

# Main call
def main():
        #### USER INPUT SECTION
        usage = """Wrapper script to perform MMseqs2 search. Provide the arguments below."""
        # Required
        p = argparse.ArgumentParser(description=usage)
        p.add_argument("-q", "--query", dest="query", type = str, required = True,
                          help="Specify the file to be used as a query")
        p.add_argument("-t", "--target", dest="target", type = str, required = True,
                          help="Specify the file to be used as a target; this can be the same as the query.")
        p.add_argument("-o", "--output", dest="output", type = str, required = True,
                          help="Specify the prefix of the output files (can include the path to direct intermediate and output files to a specific location e.g., /home/example_dir/mmseqsout)")
        p.add_argument("-tmp", "--tmpDir", dest="tmpDir", type = str, required = True,
                          help="Specify location to write temporary files to; this will be used when resuming a run, so keep it UNIQUE to a single run")
        p.add_argument("-st", "--searchtype", dest="searchtype", default="auto", choices=["auto", "blastn", "tblastx"],
                          help="Specify what type of search to run; (default==auto)")
        # Optional
        p.add_argument("-m", "--mmseqs2dir", dest="mmseqs2dir", type = str, default = "",
                          help="Specify the directory where the MMseqs2 executable is located; if it is accessible from your PATH, you can leave this blank")
        p.add_argument("-e", "--evalue", dest="evalue", type = float, default = 10,
                          help="Specify the E-value cut-off to provide as an argument (default == 10)")
        p.add_argument("-c", "--cpus", dest="threads", type = int, default = 1,
                          help="Specify the number of CPUs/threads to provide as an argument (default == 1)")
        p.add_argument("-n", "--num_iterations", dest="num_iterations", type = int, default = 4,
                          help="Specify the number of iterations to provide as an argument (default == 4)")
        p.add_argument("-s", "--sensitivity", dest="sensitivity", type = float, choices = [1,2,3,4,5,5.7,6,7,7.5], default = 7,
                          help="Specify the sensitivity number to be provided as an argument (default == 7)")
        p.add_argument("-a", "--alt-ali", dest="alt_ali", type = int, default = 0,
                          help="Specify the number of alternative alignments (similar to BLAST's HSPs) to be provided as an argument (default == 0)")
        p.add_argument("-i", "--index_skip", dest="skip_index", action="store_true", default=False,
                          help="Optionally specify whether you want to skip the indexing step")
        p.add_argument("-bs", "--blast_sort", dest="blast_sort", action="store_true", default=False,
                          help="Optionally specify whether you want the output file to be sorted by E-value similar to BLAST")
        p.add_argument("-r", "--resume", dest="resume", action="store_true", default=False,
                          help="Optionally specify whether you want the program to check for files and skip processing steps")
        p.add_argument("-p", "--profile_query", dest="profile_query", action="store_true", default=False,
                          help="Optionally specify whether you want to automatically convert the query into a profile")
        p.add_argument("-pt", "--profile_terminate", dest="profile_terminate", action="store_true", default=False,
                          help="Optionally specify whether you want to exit after query profile generation")
        p.add_argument("-iiq", "--id_is_query", dest="id_is_query", action="store_true", default=False,
                          help="Optionally specify whether you want to convert the output sequence ID to be the query file sans suffix")
        p.add_argument("-tai", "--target_as_is", dest="target_as_is", action="store_true", default=False,
                          help="""Optionally specify whether the target name should be treated as-is i.e., I won't assume a DB of it has _seqDB suffix.
                          This is mostly useful for using a database made by mmseqs databases command, not by this wrapper script.""")
        args = p.parse_args()
        args = validate_args(args)
        print('Program arguments appear to be OK. If you have any errors, try deleting the mms2tmp folder and try again.')
        
        # Make temporary folder
        tmpdir = os.path.join(args.outputdir, args.tmpDir)
        if not os.path.isdir(tmpdir):
                os.mkdir(tmpdir)

        # Make a log file for peace of mind when running BIG jobs to know that things are still happening
        logName = temp_file_name_gen(os.path.join(args.outputdir, args.output + '.mms2log'))
        open(logName, 'w').close()

        # Alert user to resume status
        if args.resume:
                print('You\'ve specified that you want to resume the run. I will attempt to do that.')
                log_update(logName, 'You\'ve specified that you want to resume the run. I will attempt to do that.')
        
        # Query, target, and output dir details
        querydir = os.listdir(args.querydir)
        if args.profile_terminate == False:
                targetdir = os.listdir(args.targetdir)
        outputdir = os.listdir(args.outputdir)

        # Format parameters for later MMseqs2 search function
        params = [args.evalue, args.threads, args.num_iterations, args.sensitivity, args.alt_ali]
        
        # Figure out our target name based on target_as_is parameter
        ## This name is only used for our if checks in this main function for resuming
        if args.target_as_is:
            targetName = args.target
        else:
            targetName = args.target + '_seqDB'

        # Make query db
        if args.profile_query == False:
                if args.resume == False or (args.query + '_seqDB' not in querydir):
                        print('Running query DB generation...')
                        log_update(logName, 'Running query DB generation...')
                        makemms2db(args.mmseqs2dir, os.path.join(args.querydir, args.query))
                else:
                        print('Skipping query DB generation...')
                        log_update(logName, 'Skipping query DB generation...')
        else:
                if args.resume == False or (args.query + '_profileDB' not in querydir):
                        print('Running query profile DB generation...')
                        log_update(logName, 'Running query profile DB generation...')
                        makemms2profile(args.mmseqs2dir, os.path.join(args.querydir, args.query), tmpdir)
                else:
                        print('Skipping query profile DB generation...')
                        log_update(logName, 'Skipping query profile DB generation...')
        if args.profile_terminate == True:
                print('Terminating program after query profile DB generation...')
                log_update(logName, 'Terminating program after query profile DB generation...')
                quit()
        
        # Make target db
        if args.resume == False or (targetName not in targetdir):
                print('Running target DB generation...')
                log_update(logName, 'Running target DB generation...')
                makemms2db(args.mmseqs2dir, os.path.join(args.targetdir, args.target)) # Note that we don't use targetName, since mms2db adds _seqDB automatically
        else:
                print('Skipping target DB generation...')
                log_update(logName, 'Skipping target DB generation...')
        
        # Indexing
        if not args.skip_index:
                # Index query db
                if args.profile_query == False:
                        index = mms2_index_exists(args.query + '_seqDB', args.querydir)
                        if args.resume == False or (index == False):
                                print('Indexing query DB...')
                                log_update(logName, 'Indexing query DB...')
                                indexmms2(args.mmseqs2dir, os.path.join(args.querydir, args.query), tmpdir, args.threads)
                        else:
                                print('Skipping query DB indexing...')
                                log_update(logName, 'Skipping query DB indexing...')
                else:
                        print('Skipping query DB indexing (not needed for profile)...')
                        log_update(logName, 'Skipping query DB indexing (not needed for profile)...')
                
                # Index target db
                index = mms2_index_exists(targetName, args.targetdir)
                if args.resume == False or (index == False):
                        print('Indexing target DB...')
                        log_update(logName, 'Indexing target DB...')
                        indexmms2(args.mmseqs2dir, os.path.join(args.targetdir, args.target), tmpdir, args.threads)
                else:
                        print('Skipping target DB indexing...')
                        log_update(logName, 'Skipping target DB indexing...')
        else:
                print('Indexing step is being skipped due to argument flag...')
                log_update(logName, 'Indexing step is being skipped due to argument flag...')
        
        # Run MMseqs2 search
        if args.resume == False or (args.output + '_mms2SEARCH.m8' not in outputdir):
                '''Because MMseqs2 creates output for each thread (mms2SEARCH.0 -> mms2SEARCH.n)
                we can't reliably check if a single _mms2SEARCH file exists like in older versions
                of this program. We should be able to simple check for the tabular file since it
                should always exist almost immediately after a search operation completes
                '''
                print('Running MMseqs2 search...')
                log_update(logName, 'Running MMseqs2 search...')
                if args.profile_query == False:
                        query = os.path.join(args.querydir, args.query + "_seqDB")
                else:
                        query = os.path.join(args.querydir, args.query + "_profileDB")
                target = os.path.join(args.targetdir, targetName)
                runmms2(args.mmseqs2dir, query, target, tmpdir, os.path.join(args.outputdir, args.output + '_mms2SEARCH'), args.searchtype, params)
        else:  
                print('Skipping MMseqs2 search...[If you want to re-run the search, delete the previous file (' + os.path.join(args.outputdir, args.output) + '_mms2SEARCH) and the mms2tmp directory]')
                log_update(logName, 'Skipping MMseqs2 search...[If you want to re-run the search, delete the previous file (' + os.path.join(args.outputdir, args.output) + '_mms2SEARCH) and the mms2tmp directory]')
        
        # Generate tabular output
        if args.resume == False or (args.output + '_mms2SEARCH.m8' not in outputdir):
                print('Generating MMseqs2 tabular output...')
                log_update(logName, 'Generating MMseqs2 tabular output...')
                if args.profile_query == False:
                        query = os.path.join(args.querydir, args.query + "_seqDB")
                else:
                        query = os.path.join(args.querydir, args.query + "_profileDB")
                target = os.path.join(args.targetdir, targetName)
                mms2tab(args.mmseqs2dir, query, target, os.path.join(args.outputdir, args.output + '_mms2SEARCH'), args.threads)
        else:
                print('Skipping MMseqs2 table generation...')
                log_update(logName, 'Skipping MMseqs2 table generation...')

        # Sort if necessary
        if args.resume == False or (args.output + '_mms2SEARCH_sorted.m8' not in outputdir):
                if args.blast_sort:
                        print('Sorting MMseqs2 output file...')
                        log_update(logName, 'Sorting MMseqs2 output file...')
                        mms2sort(os.path.join(args.outputdir, args.output + '_mms2SEARCH.m8'))
        else:
                print('Skipping MMseqs2 sorting...')
                log_update(logName, 'Skipping MMseqs2 sorting...')

        # Extra tabular modification for profile searches
        if args.id_is_query == True:
                # Determine input and outname names depending upon whether sorting occurred
                if args.profile_query == True:
                        if args.blast_sort:
                                inputName = os.path.join(args.outputdir, args.output + '_mms2SEARCH_sorted.m8')
                                finalName = os.path.join(args.outputdir, args.output + '_mms2PROFILESEARCH_sorted.m8')
                        else:
                                inputName = os.path.join(args.outputdir, args.output + '_mms2SEARCH.m8')
                                finalName = os.path.join(args.outputdir, args.output + '_mms2PROFILESEARCH.m8')
                else:
                        if args.blast_sort:
                                inputName = os.path.join(args.outputdir, args.output + '_mms2SEARCH_sorted.m8')
                                finalName = os.path.join(args.outputdir, args.output + '_mms2FAMILYSEARCH_sorted.m8')
                        else:
                                inputName = os.path.join(args.outputdir, args.output + '_mms2SEARCH.m8')
                                finalName = os.path.join(args.outputdir, args.output + '_mms2FAMILYSEARCH.m8')
                # Perform modification if necessary
                if args.resume == False or (finalName not in outputdir):
                        print('Converting MMseqs2 tabular output to ID-is-query output...')
                        log_update(logName, 'Converting MMseqs2 tabular output to ID-is-query output...')
                        mms2profiletab(args.query, inputName, os.path.join(args.outputdir, finalName))
                else:
                        print('Skipping MMseqs2 ID-is-query output modification...')
                        log_update(logName, 'Skipping MMseqs2 profile-based output modification...')

        # Done!
        print('Program completed successfully!')

if __name__ == '__main__':
        main()

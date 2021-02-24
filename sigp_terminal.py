import os, platform, subprocess, time, hashlib, random, shutil, pyperclip

def run_signalp_sequence(signalpdir, cygwindir, organism, tmpDir, seqID, protString):
        # Determine whether seqId and protString values are the proper type
        if not type(seqID) == str and not type(protString) == str:
                if not type(seqID) == list and not type(protString) == list:
                        print('run_signalp_sequence: seqID and protString inputs should both be str or both be list; this isn\'t true here, so I cannot procede.')
                        print('Fix the code leading up to this function call.')
                        quit()
                # If they are lists, ensure they have the same length
                else:
                        if len(seqID) != len(protString):
                                print('run_signalp_sequence: seqID and protString inputs are lists of nonequivalent length; I cannot procede unless this is true.')
                                print('Fix the code leading up to this function call.')
                                quit()
        # Generate temporary file for sequence
        if type(seqID) == list:
                while True:
                        tmpFileName = tmp_file_name_gen(os.path.join(tmpDir, 'tmp_sigpInput_' + ''.join([sid[0:5] for sid in seqID])[0:25] + '_'), '.fasta', ''.join([prot[0:10] for prot in protString]) + str(time.time()))
                        if os.path.isfile(tmpFileName):
                                continue
                        else:
                                break

        else:
                while True:
                        tmpFileName = tmp_file_name_gen(os.path.join(tmpDir, 'tmp_sigpInput_' + seqID + '_'), '.fasta', protString + str(time.time()))
                        if os.path.isfile(tmpFileName):
                                continue
                        else:
                                break
        with open(tmpFileName, 'w') as fileOut:
                if type(seqID) == list:
                        for i in range(len(seqID)):
                                fileOut.write('>' + seqID[i].lstrip('>') + '\n' + protString[i] + '\n')      # lstrip any > characters just in case they're already present
                else:
                        fileOut.write('>' + seqID.lstrip('>') + '\n' + protString + '\n')
        # Run signalP
        if type(seqID) == list:
                sigpResultFile = tmp_file_name_gen(os.path.join(tmpDir, 'tmp_sigpResults_' + ''.join([sid[0:5] for sid in seqID])[0:25] + '_'), '.txt', ''.join([prot[0:10] for prot in protString]) + str(time.time()))
        else:
                sigpResultFile = tmp_file_name_gen(os.path.join(tmpDir, 'tmp_sigpResults_' + seqID + '_'), '.txt', protString + str(time.time()))
        signalp_unthreaded(signalpdir, cygwindir, organism, tmpDir, tmpFileName, sigpResultFile)
        # Join and parse signalP results files
        sigPredictions = {}
        with open(sigpResultFile, 'r') as fileIn:
                for line in fileIn:
                        if line.startswith('#'):
                                continue
                        sl = line.split('\t')
                        sigPredictions[sl[0]] = [int(sl[3]), int(sl[4]), float(sl[5])] # [start, stop, score]
        # Clean up temporary files
        os.remove(tmpFileName)
        os.remove(sigpResultFile)
        # Return signalP prediction dictionary
        return sigPredictions

def signalp_unthreaded(signalpdir, cygwindir, organism, tmpDir, fastaFile, sigpResultFile):
        # Get the full fasta file location
        fastaFile = os.path.abspath(fastaFile)
        # Format signalP script text
        sigpTmpDir = tmp_file_name_gen(os.path.join(signalpdir, 'tmp_sigp_run_'), '', str(time.time()) + sigpResultFile)
        scriptText = '"{0}" -t {1} -f short -n "{2}" -T "{3}" "{4}"'.format(os.path.join(signalpdir, 'signalp'), organism, sigpResultFile, sigpTmpDir, fastaFile)
        # Generate a script for use with cygwin (if on Windows)
        if platform.system() == 'Windows':
                sigpScriptFile = os.path.join(tmpDir, tmp_file_name_gen('tmp_sigpScript_', '.sh', scriptText + str(time.time())))
                with open(sigpScriptFile, 'w') as fileOut:
                        fileOut.write(scriptText.replace('\\', '/'))
        # Run signalP depending on operating system
        if platform.system() == 'Windows':
                cmd = os.path.join(cygwindir, 'bash') + ' -l -c "' + sigpScriptFile.replace('\\', '/') + '"'
                runsigP = subprocess.Popen(cmd, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
                sigpout, sigperr = runsigP.communicate()
                os.remove(sigpScriptFile)       # Clean up temporary file
        else:
                os.putenv("PYTHONPATH",os.pathsep.join([os.getenv("PYTHONPATH",""),signalpdir]))
                runsigP = subprocess.Popen(scriptText, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
                sigpout, sigperr = runsigP.communicate()
        # Process output
        okayLines = ['is an unknown amino amino acid', 'perl: warning:', 'LC_ALL =', 'LANG =', 'are supported and installed on your system', '# temporary directory will not be removed']
        for line in sigperr.decode("utf-8").split('\n'):
                # If sigperr indicates null result, create an output file we can skip later
                if line.rstrip('\n') == '# No sequences predicted with a signal peptide':
                        with open(sigpResultFile, 'w') as fileOut:
                                fileOut.write(line)
                        break
                # Check if this line has something present within okayLines
                okay = 'n'
                for entry in okayLines:
                        if entry in line or line == '':
                                okay = 'y'
                                break
                if okay == 'y':
                        continue
                # If nothing matches the okayLines list, we have a potentially true error
                else:
                        raise Exception('SignalP error occurred when processing file name ' + fastaFile + '. Error text below\n' + sigperr.decode("utf-8"))
        # Clean up tmp dir
        if os.path.isdir(sigpTmpDir):
                shutil.rmtree(sigpTmpDir)

def tmp_file_name_gen(prefix, suffix, hashString):
        # Main function
        tmpHash = hashlib.sha256(bytes(str(hashString) + str(time.time()) + str(random.randint(0, 100000)), 'utf-8') ).hexdigest()       # This should always give us something unique even if the string for hashString is the same across different runs
        while True:
                if os.path.isfile(prefix + tmpHash + suffix):
                        tmpHash += 'X'
                else:
                        return prefix + tmpHash + suffix

if __name__ == "__main__":
        signalpdir = r"D:\Bioinformatics\Protein_analysis\signalp-4.1f.CYGWIN\signalp-4.1"
        tmpDir = r"D:\Bioinformatics\Protein_analysis\signalp-4.1f.CYGWIN\signalp-4.1\tmp"
        cygwindir = ""
        organism = "euk"
        while True:
                try:
                        while True:
                            print("Copy or enter the sequence(s) to check for sigP")
                            text = input()
                            if text == "":
                                    text = pyperclip.paste()
                            # Continue if blank
                            if text == "":
                                    print("Nothing was entered; try again.")
                                    continue
                            # Clean up text
                            text = text.replace("\r", "")
                            lines = text.rstrip("\n").split("\n")
                            # Get fasta-format
                            seqs = []
                            seqIDs = []
                            if ">" in text:
                                for i in range(0, len(lines), 2):
                                    seqIDs.append(lines[i][1:])
                                    seqs.append(lines[i+1])
                            else:
                                for i in range(0, len(lines)):
                                    seqIDs.append("seq{0}".format(i+1))
                                    seqs.append(lines[i])
                            # Run sigP
                            predictions = run_signalp_sequence(signalpdir, cygwindir, organism, tmpDir, seqIDs, seqs)
                            # Format results
                            out = []
                            for seqID in seqIDs:
                                    if seqID in predictions:
                                            prediction = predictions[seqID]
                                            out.append("{0}\t{1}-{2}, {3}".format(seqID, prediction[0], prediction[1], prediction[2]))
                                    else:
                                            out.append("{0}\t.".format(seqID))
                            # Present results
                            print("Done! Results printed to terminal below and saved to your clipboard")
                            for line in out:
                                    print(line)
                            pyperclip.copy("\n".join(out))
                except KeyboardInterrupt:
                        print("Keyboard interrupt received; program exiting")
                        quit()
                except:
                        print("Something went wrong, try it again?")

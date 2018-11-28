#! python3
# msa_curate.py
# Program to assist in the curation of multiple sequence alignments. This is 
# in very early development, but I expect this program to eventually have
# a few different "modes" which will enable different curation behaviour
# to take place.

# Import external packages
import os, argparse, shutil, sys

# Define functions for later use
## Validate arguments
def validate_args(args):
        # Validate that necessary arguments have been provided
        if args.inputLocation == None:
                print('-i argument must be provided; fix your input and try again.')
                quit()
        if args.outputLocation == None:
                print('-o argument must be provided; fix your input and try again.')
                quit()
        if args.mode == None:
                print('-m argument must be provided; fix your input and try again.')
                quit()
        # Validate the FASTA file input location depending on type of input
        if len(args.inputLocation) == 1:
                # Check that specified value is a path
                if not os.path.isdir(args.inputLocation[0]):
                        print('One value was provided for -i, which means you should have provided a directory containing FASTA files.')
                        print('The provided value "' + args.inputLocation[0] + '" is not a directory; either it does not exist or it is a file; fix your input and try again.')
                        quit()
        else:
                # Check that the specified values are files
                for file in args.inputLocation:
                        if not os.path.isfile(file):
                                print('Multiple values were provided for -i, which means you should have provided individual FASTA files.')
                                print('The provided value "' + file + '" is not a file; either it does not exist or it is a directory; fix your input and try again.')
                                quit()
        # Ensure that the output location is sensible
        if os.path.isfile(args.outputLocation):
                print('The specified output location "' + args.outputLocation + '" is a file. You should be specifying a directory (that may or may not exist); fix your input and try again.')
                quit()
        elif not os.path.isdir(args.outputLocation):
                pathSplit = list(os.path.split(args.outputLocation))
                if pathSplit[0] == '':
                        pathSplit[0] = os.getcwd()
                if not os.path.isdir(pathSplit[0]):
                        print('The specified output location "' + args.outputLocation + '" is not a directory.')
                        print('The location one directory up i.e., "' + pathSplit[0] + '" is also not a directory.')
                        print('This program won\'t create new directories that aren\'t contained within an existing one.')
                        print('Either specify a new location or create one of the above two mentioned locations and try again.')
                        quit()
        # Handle arguments specific to modes
        if args.mode == 'outliers':
                # Ensure that necessary prerequisites are installed
                try:
                        from Bio import AlignIO
                except:
                        print('biopython is a prerequisite of the "outliers" function; install this then try again.')
                        quit()
                try:
                        import hdbscan
                except:
                        print('hdbscan is a prerequisite of the "outliers" function; install this then try again.')
                        quit()
                try:
                        import numpy
                except:
                        print('numpy is a prerequisite of the "outliers" function; install this then try again.')
                        quit()
                # Replace value if necessary
                if args.rscriptdir == None:
                        args.rscriptdir = ''
                # Validate program execution is successful
                program_execution_check(os.path.join(args.rscriptdir, 'Rscript'))
        if args.mode == 'start_trim':
                placeholder = True      # Currently there are no special handling for these functions, but this might change if the functions evolve
        if args.mode == 'full_trim':
                # Validate that propTrim is sensible
                if not 0 <= args.propTrim <= 1:
                        print('-p value must be greater than or equal to 0, and less than or equal to 1; it is a proportion value.')
                        print('Fix your input and try again.')
                        quit()
                # Validate that dropProp is sensible
                if not 0 <= args.dropProp <= 1:
                        print('-d value must be greater than or equal to 0, and less than or equal to 1; it is a proportion value.')
                        print('Fix your input and try again.')
                        quit()
        if args.mode == 'full_trim' or args.mode == 'start_trim':
                # Validate that lengthMin is sensible
                if 0 > args.lengthMin:
                        print('-l value must be greater than or equal to 0; it is either a proportion value like "0 <= lengthMin <= 1" or an absolute integer value > 1.')
                        print('Fix your input and try again.')
                        quit()
        return args

def program_execution_check(cmd):
        import subprocess
        run_cmd = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
        cmdout, cmderr = run_cmd.communicate()
        if cmderr.decode("utf-8") != '' and not cmderr.decode("utf-8").startswith('Usage'):
                print('Failed to execute program "' + cmd + '". Is this executable in the location specified/discoverable in your PATH, or does the executable even exist? I won\'t be able to run properly if I can\'t execute this program.')
                print('---')
                print('stderr is below for debugging purposes.')
                print(cmderr.decode("utf-8"))
                print('Program closing now.')
                quit()

## MSA curation-related
### full_trim
def msa_trim(msaFastaIn, pctTrim, minLength, outType, msaFastaOut, indivSeqDrop, skipOrDrop):
        '''
        msaFastaIn is the path to the aligned MSA FASTA file to be trimmed
        pctTrim refers to the minimum proportion of sequences present in a single column to demarcate the start and end of an alignment
        minLength refers to the minimum length of a MSA after trimming before we decide to not trim at all; if this value is less than 1,
        we assume it's a ratio, otherwise it is an absolute length. 
        outType influences whether this function returns a Biopython MSA object ("obj") or an output file ("file")
        msaFastaOut is only relevant when outType == file; otherwise it will be ignored
        '''
        # Set up
        import os, copy
        from Bio import AlignIO
        from Bio.Seq import Seq
        from Bio.Alphabet import SingleLetterAlphabet
        from Bio.Align import MultipleSeqAlignment
        # Ensure outType is sensible
        if outType.lower() not in ['obj', 'file', 'both']:
                print('msa_trim: This function requires an outType to be specified with a specific format.')
                print('Your provided value "' + outType + '" should instead be "obj", to return the Biopython MSA object, "file" to produce an output file which uses the string provided as msaFastaOut, or "both" to do both aforementioned things.')
                print('Format this correctly and try again.')
                quit()
        if (outType.lower() == 'file' or outType.lower() == 'both') and not type(msaFastaOut) == str:
                print('msa_trim: You specified a file output but didn\'t provide a string for msaFasta out - the file output name.')
                print('Format this correctly and try again.')
                quit()
        # Process minLength and ensure it is sensible
        try:
                int(minLength)
        except:
                print('msa_trim: minLength must be an integer or capable of conversion to integer.')
                print('Format this correctly and try again.')
                quit()
        if minLength < 0:
                print('msa_trim: minLength must be greater than 0.')
                print('Format this correctly and try again.')
                quit()
        # Process indivSeqDrop and ensure it is sensible
        if indivSeqDrop != None:
                try:
                        float(indivSeqDrop)
                        if not 0 <= float(indivSeqDrop) <= 1:
                                print('msa_trim: indivSeqDrop appears to be a float, but it is not a value from 0 to 1.')
                                print('Format this correctly and try again.')
                                quit()
                        indivSeqDrop = float(indivSeqDrop)
                except:
                        print('msa_trim: indivSeqDrop was not specified as None, but is also not capable of conversion to float.')
                        print('Format this correctly and try again.')
                        quit()
        # Process skipOrDrop and ensure it is sensible
        if skipOrDrop.lower() not in ['skip', 'drop']:
                print('msa_trim: skipOrDrop must equal "skip" or "drop"; I don\'t recognise ' + skipOrDrop + '.')
                print('Format this correctly and try again.')
                quit()
        # Load in fasta file as MSA object
        msa = AlignIO.read(msaFastaIn, 'fasta')
        # Loop through aligned columns and find the first position from the 5' end that meets our pctTrim value
        for i in range(len(msa[0].seq)):
                col = msa[:,i]
                pctBases = 1 - (col.count('-') / len(col))
                if pctBases <= pctTrim:
                        continue
                break
        # Same but for 3' end
        for x in range(len(msa[0].seq), 0, -1):
                col = msa[:,x-1]
                pctBases = 1 - (col.count('-') / len(col))
                if pctBases <= pctTrim:
                        continue
                break
        # Check our values to ensure they're sensible
        if i >= x:      # If i >= x, that means we'd be trimming the sequence to 1bp or a negative value; in other words, we can't trim it at this pctTrim as printed below
                if skipOrDrop.lower() == 'skip':
                        print('#"' + os.path.basename(msaFastaIn) + '" can\'t be trimmed at this pctTrim value since no columns contain this proportion of sequences; no trimming will be performed.')
                        return msa              # If the user isn't expecting a returned object this should just disappear; if they want a file out, we won't modify it
                elif skipOrDrop.lower() == 'drop':
                        print('#"' + os.path.basename(msaFastaIn) + '" can\'t be trimmed at this pctTrim value since no columns contain this proportion of sequences; msa will be dropped.')
                        return None
        # Compare our MSA length post-trimming to our specified cut-offs to determine whether we're doing anything to this sequence or not
        seqLen = x - i          # This works out fine in 1-based notation
        if minLength < 1:
                ratio = seqLen / len(msa[0])
                if ratio < minLength:
                        if skipOrDrop.lower() == 'skip':
                                print('#"' + os.path.basename(msaFastaIn) + '" trimming reduces length more than minLength proportion cut-off; no trimming will be performed.')
                                return msa      # We're not going to make any changes if trimming shortens it too much
                        elif skipOrDrop.lower() == 'drop':
                                print('#"' + os.path.basename(msaFastaIn) + '" trimming reduces length more than minLength proportion cut-off; msa will be dropped.')
                                return None
        else:
                if seqLen < minLength:
                        if skipOrDrop.lower() == 'skip':
                                print('#"' + os.path.basename(msaFastaIn) + '" trimming reduces length more than absolute minLength cut-off; no trimming will be performed.')
                                return msa      # As above
                        elif skipOrDrop.lower() == 'drop':
                                print('#"' + os.path.basename(msaFastaIn) + '" trimming reduces length more than absolute minLength cut-off; msa will be dropped.')
                                return None
        # Check MSA columns to see what the most common character is for each position i.e., whether they are gap or sequence
        charType = ''
        for z in range(i, x):
                col = msa[:,z]
                pctBases = 1 - (col.count('-') / len(col))
                if pctBases >= 0.50:
                        charType += 'N'
                else:
                        charType += '-'
        # Trim our MSA object
        origMsa = copy.deepcopy(msa)    # Since we're going to be making changes to the msa from here on but still might want to return the unedited msa, we need to create a backup
        newMsa = MultipleSeqAlignment([])
        for y in range(len(msa)):       # Don't overwrite i from above! I made this mistake...
                msa[y].seq = Seq(str(msa[y].seq)[i:x], SingleLetterAlphabet())
                # Optionally remove sequences that don't appear to "fit" the alignment [i.e., contain a lot of gaps] according to user-specified gap proportion cut-off
                if indivSeqDrop != None:
                        indivGapRepresentation = ''
                        for z in range(len(msa[y].seq)):
                                if msa[y].seq[z] != '-' and charType[z] != '-':         # i.e., if it agrees with consensus, we don't count it as a gap
                                        indivGapRepresentation += 'N'
                                elif msa[y].seq[z] == '-' and charType[z] == '-':       # i.e., if it agrees with consensus, we don't count it as a gap
                                        indivGapRepresentation += 'N'
                                else:
                                        indivGapRepresentation += '-'                   # i.e., if it differs from consensus, it's treated as a gap/gap-inducer
                        gapCount = indivGapRepresentation.count('-')
                        gapProp = gapCount / len(indivGapRepresentation)
                        if gapProp > indivSeqDrop:
                                continue
                        newMsa.append(msa[y])
        # If we dropped sequences, make sure we don't have any blank columns now
        if indivSeqDrop != None and len(newMsa) > 1:
                for a in range(len(newMsa[0].seq), 0, -1):
                        col = newMsa[:,a-1]
                        if set(col) == {'-'}:
                                for b in range(len(newMsa)):
                                        newMsa[b].seq = Seq(str(newMsa[b].seq)[0:a-1] + str(newMsa[b].seq)[a:], SingleLetterAlphabet())
        # If we dropped sequences, ensure that our newMsa still has more than one entry in it
        if indivSeqDrop != None:
                if len(newMsa) < 2:
                        if skipOrDrop.lower() == 'skip':
                                print('#"' + os.path.basename(msaFastaIn) + '" removing gappy sequences according to indivSeqDrop cut-off means we do not have >= 2 sequences in this msa; no trimming will be performed.')
                                return origMsa
                        elif skipOrDrop.lower() == 'drop':
                                print('#"' + os.path.basename(msaFastaIn) + '" removing gappy sequences according to indivSeqDrop cut-off means we do not have >= 2 sequences in this msa; msa will be dropped.')
                                return None
                msa = newMsa
        # Return results either as the MSA object, as an output file, or as both
        if outType.lower() == 'file' or outType.lower() == 'both':
                with open(msaFastaOut, 'w') as fileOut:
                        fileOut.write(msa.format('fasta'))
        if outType.lower() == 'obj' or outType.lower() == 'both':
                return msa

def msa_score(msa, scoringMethod, cutoff):
        # Set up
        import os
        from Bio import AlignIO
        import Bio.Align
        # Define functions integral to this one
        def blosum_score(pair):
                # Set up
                from Bio.SubsMat import MatrixInfo
                blosum = MatrixInfo.blosum62    # blosum62 seems to be the most commonly used matrix, I'll lean on consensus w/r/t which one is probably best here
                # Make sure pair is formatted correctly
                pair = (pair[0].upper(), pair[1].upper())
                if pair not in blosum:
                        pair = tuple(reversed(pair))
                # Return the score
                return blosum[pair]
        def sumpairs_score(pair):
                if pair[0].upper() == pair[1].upper():
                        return 1
                else:
                        return 0
        # Determine what input we are receiving
        if type(msa) == Bio.Align.MultipleSeqAlignment:
                fileType = 'obj'
        elif type(msa) == str:
                if not os.path.isfile(msa):
                        print('msa_score: You\'ve provided a string input but it isn\'t detected as a file.')
                        print('Did you type the file name correctly, or do you need to provide the full path to it? Fix your input and try again.')
                        quit()
                fileType = 'file'
        else:
                print('msa_score: You\'ve provided an input type I am not compatible with.')
                print('The input should be a Bio.Align.MultipleSeqAlignment object or a string indicating the file location. Fix your input and try again.')
                quit()
        # If we're working with a file input, load it as an obj
        if fileType == 'file':
                msa = AlignIO.read(msa, 'fasta')
        # Loop through msa columns and perform pairwise scoring
        overallScores = []
        for i in range(len(msa[0].seq)):
                col = msa[:,i]
                colScores = []
                # Column scores based on possible pairs exclusive of '-'
                for x in range(len(col)-1):
                        for y in range(x+1, len(col)):
                                pair = (col[x].upper(), col[y].upper())
                                if pair[0] != '-' and pair[1] != '-':
                                        if scoringMethod == 'blosum':
                                                colScores.append(blosum_score(pair))
                                        elif scoringMethod == 'sp':
                                                colScores.append(sumpairs_score(pair))
                                else:
                                        continue
                # Average column score with penalty for gaps
                if colScores != []:
                        colScore = (sum(colScores) / len(colScores)) * (1 - (col.count('-') / len(col)))
                else:
                        colScore = 0    # If there is only a single letter, it receives no score
                overallScores.append(colScore)
        # Compare our overall score against cutoff and determine if this alignment is "good enough"
        finalScore = sum(overallScores) / len(overallScores)
        if finalScore >= cutoff:
                return True
        else:
                return False

def odseq_outlier_detect(msaFileNameList, rScriptDir, tmpDir, threshold, distMetric, bootStraps):
        # Set up
        import os, subprocess, pathlib
        ## Ensure input parameters are sensible
        # rScriptDir
        if rScriptDir == None:
                rScriptDir = ''         # We'll assume Rscript is locatable in PATH if unspecified
        elif rScriptDir != '' and not (os.path.isfile(os.path.join(rScriptDir, 'Rscript.exe')) or os.path.isfile(os.path.join(rScriptDir, 'Rscript'))):
                print('odseq_outlier_detect: rScriptDir does not appear to contain the Rscript file.')
                print('Fix your input and try again.')
                quit()
        # tmpDir
        if tmpDir == None:
                tmpDir = '.'            # We'll just put the file in the current directory if it isn't specified
        elif tmpDir != '' and not os.path.isdir(tmpDir):
                print('odseq_outlier_detect: tmpDir is not an existing directory or not able to be located.')
                print('Fix your input and try again.')
                quit()
        # threshold
        try:
                threshold = float(threshold)
        except:
                print('odseq_outlier_detect: threshold needs to be a float or capable of conversion to float.')
                print('Fix your input and try again.')
                quit()
        # distMetric
        if distMetric.lower() not in ['linear', 'affine']:
                print('odseq_outlier_detect: distMetric needs to be an option available in the list below. Fix your input and try again.')
                print(['linear', 'affine'])
                quit()
        # bootStraps
        try:
                bootStraps = int(bootStraps)
        except:
                print('odseq_outlier_detect: bootStraps needs to be an integer or capable of conversion to integer.')
                print('Fix your input and try again.')
                quit()
        # msaFileNameList
        if type(msaFileNameList) == str:
                msaFileNameList = [msaFileNameList]
        elif type(msaFileNameList) != list:
                print('odseq_outlier_detect: msaFileNameList type is not recognisable. It should be a list, but instead it is ' + str(type(msaFileNameList)) + '.')
                print('Fix your input and try again.')
                quit()
        if msaFileNameList == []:
                print('odseq_outlier_detect: msaFileNameList is empty. I don\'t know what to do in this situation since it shouldn\'t happen.')
                print('Code your call to this function properly to skip it.')
                quit()
        # Ensure that the msaFileNames are locatable
        ongoingCount = 0
        for fileName in msaFileNameList:
                if not os.path.isfile(fileName):
                        print('odseq_outlier_detect: index ' + str(ongoingCount) + ' in msaFileNameList is not able to be located. You might need to specify the full path to the file.')
                        print('Fix your input and try again.')
                        quit()
                msaFileNameList[ongoingCount] = os.path.abspath(fileName)
                ongoingCount += 1
        # Create script file
        scriptText = ['library("msa")', 'library("odseq")']
        for fileName in msaFileNameList:
                fileName = pathlib.Path(fileName).as_posix()
                scriptText.append('filename = "' + fileName + '"')
                scriptText.append('alig <- readAAMultipleAlignment(filename)')
                scriptText.append('y <- odseq(alig, threshold = {}, distance_metric = "{}", B = {})'.format(threshold, distMetric.lower(), bootStraps))
                scriptText.append('print(filename)')
                scriptText.append('print(y)')
        scriptFile = file_name_gen(os.path.join(tmpDir, 'tmp_odseq_script'), '.R')
        with open(scriptFile, 'w') as fileOut:
                fileOut.write('\n'.join(scriptText))
        # Format cmd
        cmd = '"' + os.path.join(rScriptDir, 'Rscript') + '" ' + scriptFile
        run_odseq = subprocess.Popen(cmd, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        odseqout, odseqerr = run_odseq.communicate()
        odseqout = odseqout.decode("utf-8")
        if 'FALSE' not in odseqout and 'TRUE' not in odseqout:  # Rscript writes module loading details to stderr so we need to check that it worked properly in other ways
                raise Exception('Rscript ODseq error text below' + str(odseqerr.decode("utf-8")))
        # Parse ODseq results
        odseqTable = odseqout.split('[1]')
        if odseqTable[0] == '':
                del odseqTable[0]
        odseqDict = {}
        for i in range(len(odseqTable)):
                # Extract the chunk of information from this ODseq result
                chunk = odseqTable[i].split('\n')
                fileName = os.path.basename(chunk[0].strip(' "'))
                # Clean up information-less lines
                del chunk[0]
                if chunk[-1] == '':
                        del chunk[-1]
                # Assemble the boolean results in a list [Note that this is an ordered list which corresponds to the input MSA, hence why there should be no need to parse file names for reassociation - we already have msaFileNameList]
                odseqResults = []
                for x in range(len(chunk)):
                        if x % 2 == 1:                  # i.e., if we're looking at an odd number in the list, it should be a boolean result line rather than sequence ID line
                                sl = chunk[x].split()   # It's probably overly cautious, but it does let us use sequences named 'TRUE' and 'FALSE' (why would you do this...) which a simple if == statement would not
                                odseqResults += sl
                # Convert 'TRUE' to True, 'FALSE' to False
                for x in range(len(odseqResults)):
                        if odseqResults[x] == 'TRUE':
                                odseqResults[x] = True
                        elif odseqResults[x] == 'FALSE':
                                odseqResults[x] = False
                        else:
                                print('odseq_outlier_detect: unrecognised output in odseqResults (' + str(odseqResults[x]) + '). What\'s going on?')
                                quit()
                # Add to our dictionary using filename as key
                odseqDict[fileName] = odseqResults
        # Ensure everything worked fine
        if len(odseqDict) != len(msaFileNameList):
                print('odseq_outlier_detect: length of odseqDict != length of msaFileNameList. Inspect the below stderr report to see what went wrong and try to fix it.')
                print(odseqerr.decode("utf-8"))
                quit()
        # Clean up tmp file
        os.unlink(scriptFile)
        # Return results
        return odseqDict

def protein_alphabet_reduce(proteinList, reduceNum):
        # Set up
        elevenLetter = {'E':'D','L':'I','M':'I','Q':'K','R':'K','T':'S','V':'I','W':'F','Y':'F'}
        fifteenLetter = {"L":"L","V":"L","I":"L","M":"L","C":"C","A":"A","G":"G","S":"S","T":"T","P":"P","F":"F","Y":"F","W":"W","E":"E","D":"D","N":"N","Q":"Q","K":"K","R":"K","H":"H"}
        # Ensure reduceNum is sensible and hold onto the correct reduce dict OR return our unmodified protein list
        if reduceNum == None or reduceNum == 'n' or type(reduceNum) == bool:
                return proteinList      # No modification is necessary
        elif int(reduceNum) == 11:
                reduceDict = elevenLetter
        elif int(reduceNum) == 15:
                reduceDict = fifteenLetter
        else:
                print('I didn\'t recognise the reduceNum value provided to the protein_alphabet_reduce function. It should be None, \'n\', 11, or 15.')
                print('I\'m just going to treat this as None... if you don\'t want this behaviour, fix your input.')
                return proteinList      # No modification is necessary
        # Ensure our proteinList is a list; if a single str is provided, make it a list (then return the str back from the function later)
        listAtEnd = True
        if type(proteinList) == str:
                proteinList = [proteinList]
                listAtEnd = False
        # Main function
        for i in range(len(proteinList)):
                newseq = ''
                for letter in proteinList[i]:
                        if letter in reduceDict:
                                newseq += reduceDict[letter]
                        else:
                                newseq += letter
                proteinList[i] = newseq
        # Return our modified list
        if listAtEnd == False:
                proteinList = proteinList[0]
        return proteinList

def alfree_matrix(fastaFile, wordSize, reduceNum, alfAlgorithm):
        # Set up
        from alfpy import word_pattern, word_vector, word_distance
        from alfpy.utils import seqrecords, distmatrix
        from alfpy.utils.data import seqcontent
        # Read in unclustered domains file
        unclustDoms = open(fastaFile)
        records = seqrecords.read_fasta(unclustDoms)
        unclustDoms.close()
        # Extract details from records using alfpy-provided functions
        seqList = records.seq_list
        lengthList = records.length_list
        idList = records.id_list
        # Optional reduction of protein alphabet
        seqList = protein_alphabet_reduce(seqList, reduceNum)
        # Compute distance matrix for word sizes
        matrices = []           # Currently I'm not returning multiple matrices, but this exists to enable multiple word sizes to be used and combined
        wordSizes = [wordSize]  # As above, not used now, but maybe in the future this will be relevant
        for num in wordSizes:
                p = word_pattern.create(seqList, word_size=wordSize)
                if alfAlgorithm == 'canberra':
                        weightmodel = word_vector.WeightModel(seqcontent.get_weights('protein'))
                        counts = word_vector.CountsWeight(lengthList, p, weightmodel)
                else:
                        counts = word_vector.Counts(lengthList, p)
                dist = word_distance.Distance(counts, alfAlgorithm)
                matrices.append(distmatrix.create(idList, dist))
        # Return value
        return matrices, idList

def msa_outlier_detect(msaFileNameList, statsSave, removeIdentical, scoreType):
        # Set up
        import hdbscan, statistics
        from Bio import AlignIO
        from Bio.Align import MultipleSeqAlignment
        import numpy as np
        # Set default HDBSCAN parameters for detecting outliers in MSA
        allowSingle = True      # Making this user tunable is probably incorrect.
        clustSelect = 'eom'     # This method is pretty imprecise as it stands,
        minSize = 2             # but it seems like these settings are the only
        minSample = 2           # ones that really "work"
        # Define functions integral to this one
        def sumpairs_score(pair):
                if pair[0].upper() == pair[1].upper():
                        return 1
                else:
                        return 0
        def pairwise_sumpairs_matrix(msa):
                spScore = []
                for a in range(len(msa)):
                        spScore.append([])
                        for b in range(len(msa)):
                                colScores = []
                                gapCount = 0
                                if a == b:
                                        spScore[-1].append(0.00)
                                        continue
                                for i in range(len(msa[a].seq)):
                                        pair = (msa[a][i], msa[b][i])
                                        if pair[0] != '-' and pair[1] != '-':
                                                colScores.append(sumpairs_score(pair))
                                        elif pair[0] == '-' and pair[1] == '-':         # Don't penalise an alignment that has a gap induced by another sequence
                                                continue
                                        elif set(msa[a][i:]) == {'-'}:                  # Don't penalise an alignment that has ended
                                                continue
                                        else:                                           # If this induces a gap in another alignment or is gapped itself, penalise it
                                                gapCount += 1
                                # Average column score with penalty for gaps
                                if colScores != []:
                                        colScore = 1 - (sum(colScores) / len(colScores)) * (1 - (gapCount / len(msa[a])))       # 1 - (..) since we want this to be a measure of dissimilarity in line with what alfpy does
                                        #colScore = 1 - (sum(colScores) / len(colScores))
                                else:
                                        colScore = 1
                                spScore[-1].append(colScore)
                return spScore
        
        # Ensure msaFileNameList is correctly formatted
        if type(msaFileNameList) == str:
                msaFileNameList = [msaFileNameList]
        elif type(msaFileNameList) != list:
                print('msa_outlier_detect: msaFileNameList type is not recognisable. It should be a list, but instead it is ' + str(type(msaFileNameList)) + '.')
                print('Fix your input and try again.')
                quit()
        if msaFileNameList == []:
                print('msa_outlier_detect: msaFileNameList is empty. I don\'t know what to do in this situation since it shouldn\'t happen.')
                print('Code your call to this function properly to skip it.')
                quit()
        # Ensure statsSave is correctly formatted
        if statsSave != True and statsSave != False:
                print('msa_outlier_detect: statsSave should equal True or False, not ' + str(statsSave) + '.')
                print('Fix your input and try again.')
                quit()
        # Ensure removeIdentical is correctly formatted
        if removeIdentical != True and removeIdentical != False:
                print('msa_outlier_detect: removeIdentical should equal True or False, not ' + str(removeIdentical) + '.')
                print('Fix your input and try again.')
                quit()
        # Ensure scoreType is correctly formatted
        if scoreType.lower() not in ['sp', 'alfree']:
                print('msa_outlier_detect: scoreType should equal "sp" or "alfree", not ' + str(scoreType) + '.')
                print('Fix your input and try again.')
                quit()
        # Loop through MSA files and try to identify outliers
        outlierDict = {}
        ongoingCount = 0
        for fileName in msaFileNameList:
                msa = AlignIO.read(fileName, 'fasta')
                origMsaLen = len(msa)
                # Remove identical sequences if relevant [these can skew our pairwise scoring matrix if highly similar sequence bits slip through]
                if removeIdentical == True:
                        madeChanges = False
                        newMsa = MultipleSeqAlignment([])
                        identicalPairs = []
                        for y in range(len(msa)):
                                for z in range(len(msa)):
                                        if y == z:
                                                continue
                                        elif msa[y].seq == msa[z].seq:
                                                identicalPairs.append([y, z])
                        if identicalPairs != []:
                                # Populate identicalPairs list with full groups
                                for pair in identicalPairs:
                                        for n in range(len(identicalPairs)):
                                                if pair[0] in identicalPairs[n] or pair[1] in identicalPairs[n]:
                                                        if pair[0] not in identicalPairs[n]:
                                                                identicalPairs[n].append(pair[0])
                                                        elif pair[1] not in identicalPairs[n]:
                                                                identicalPairs[n].append(pair[1])
                                # Sort then remove redundancy in identicalPairs to render our groups
                                for i in range(len(identicalPairs)):
                                        identicalPairs[i].sort()
                                identicalPairs.sort()
                                removeGroups = []
                                for pair in identicalPairs:
                                        if pair not in removeGroups:
                                                removeGroups.append(pair)
                                # Remove all but one entry from each identical group from our newMsa
                                for y in range(len(msa)):
                                        found = False
                                        for group in removeGroups:
                                                if y in group[1:]:      # This means we can allow the first number in each removeGroup
                                                        found = True
                                                        break
                                        if found == False:
                                                newMsa.append(msa[y])
                                if len(newMsa) > 2:     # We don't want to work with a structure with less than 3 sequences since our pairwise scoring becomes less impactful. It's a bit arbitrary in some respects,
                                        madeChanges = True
                                        msa = newMsa    # but this should help a domain model to not be overwhelmed by identical sequences while still letting us meaningfully use means and stdev for outlier detection.
                # Perform pairwise scoring
                if scoreType.lower() == 'sp':
                        spScore = pairwise_sumpairs_matrix(msa)
                        spdm = np.array(spScore)
                else:
                        # Make a temp file with our altered MSA
                        tmpName = file_name_gen('tmp_alfree', '.fasta')
                        with open(tmpName, 'w') as fileOut:
                                fileOut.write(msa.format('fasta'))
                        spScore, idList = alfree_matrix(tmpName, 1, None, 'google')
                        spdm = spScore[0].data
                        spScore = spdm.tolist()
                        # Clean up tmp file
                        os.unlink(tmpName)
                # Cluster with HDBSCAN
                clusterer = hdbscan.HDBSCAN(metric='precomputed', cluster_selection_method = clustSelect, min_cluster_size = int(minSize), min_samples = int(minSample), allow_single_cluster = allowSingle)
                clusterer.fit(spdm)
                # Calculate basic statistics from pairwise scoring excluding HDBSCAN detected outliers
                if statsSave == True:
                        spMeanList = []
                        spAllMeansList = []
                        for i in range(len(spScore)):
                                spMean = statistics.mean(spScore[i][:i] + spScore[i][i+1:])     # Exclude self-match
                                spAllMeansList.append(spMean)
                                if clusterer.labels_[i] == -1:
                                        continue
                                spMeanList.append(spMean)
                        spMeansMean = statistics.mean(spMeanList)
                        spMeansPsdt = statistics.stdev(spMeanList)     # If len(spMeanList) == 1, stdev == 0. This makes it more likely we remove a legitimate sequence. TESTING: See if I should make this 10% of mean or something like that?
                # Convert HDBSCAN groups into boolean list of True == outlier, False == not outlier
                outlierList = []
                for i in range(len(clusterer.labels_)):
                        label = clusterer.labels_[i]
                        if label == -1 or (-1 not in clusterer.labels_ and 1 in clusterer.labels_):     # If the second condition is True, HDBSCAN didn't find any outliers but it did find at least 2 separate clusters. In this case, it seems
                                if statsSave == True:                                                   # appropriate to treat every sequence as a potential outlier, and use our statistical distribution to pick out huge outliers
                                        '''Note: This serves as a _really_ rough heuristic check to justify HDBSCAN's clustering decision. It involves 
                                        a basic comparison using pstdev to see if this row's mean SP score is anomalous compared to others. At the start
                                        is an additional hard cut-off check to make sure we don't remove something "a little bit" different to a group of 
                                        sequences that are otherwise very similar. I added this cut-off check as a result of manual inspection of results
                                        to prevent a mistake from happening to a specific cluster I was testing. The test scenario was this [0.26951219512195124,
                                        0.33048780487804874, 0.2, 0.22073170731707314, 0.2182926829268293, 0.21707317073170732]. 0.33 was detected as an outlier,
                                        but 0.33 is still really similar for a SP score. 0.2 * 2 == 0.4, and this helps to rescue our example. 0.5 seems to
                                        be a point where the cluster is no longer highly homogenous.
                                        The second hard cut-off check was derived from [0.5294117647058824, 0.47500000000000003, 0.46029411764705885, 0.5042016806722689,
                                        0.4613445378151261]. It exceeds our 0.5 cut-off so we want to be less lenient with it, but the first sequence is still
                                        quite similar to the others from manual inspection. 0.1 seems to be a good point where, even if the minimum mean is 0.5,
                                        a sequence with distance 0.6 to the others still looks "normal" in a MSA. 0.7 as the max for this cut-off is a bit
                                        arbitrary and it shouldn't really happen, but it's just to prevent any weirdness from happening (e.g., a cluster
                                        of 0.9 distances should not group with a 1.0 distance [even though a 0.9 cluster shouldn't exist])
                                        '''
                                        if spAllMeansList[i] < 0.5 and spAllMeansList[i] < min(spAllMeansList) * 2:
                                                outlierList.append(False)
                                        elif spAllMeansList[i] < 0.6 and spAllMeansList[i] < min(spAllMeansList) + 0.15:        # This is another hard cut-off case derived from real data
                                                outlierList.append(False)                                                       # tldr; 0.45 and 0.6 are compatible in a cluster
                                        elif spAllMeansList[i] < 0.7 and spAllMeansList[i] < min(spAllMeansList) + 0.1:
                                                outlierList.append(False)
                                        elif spAllMeansList[i] > spMeansMean + (1.5*spMeansPsdt):
                                                outlierList.append(True)
                                        else:
                                                outlierList.append(False)
                                else:
                                        outlierList.append(True)
                        else:
                                outlierList.append(False)
                # Add in previously deleted identical values if removeIdentical is specified and we made changes
                if removeIdentical == True:
                        if madeChanges == True:
                                for x in range(origMsaLen):
                                        found = False
                                        for group in removeGroups:
                                                if x in group[1:]:
                                                        found = group
                                        if found != False:
                                                outlierList.insert(x, outlierList[found[0]])    # This took me a bit of mental effort to devise, then a bit more to understand why it worked, but it's quite simple
                assert origMsaLen == len(outlierList)                                                   # When we find an index that was removed, we just insert an identical copy of its remaining sequence result, the index of which is given by found[0]
                outlierDict[os.path.basename(fileName)] = outlierList                                   # ... but we still check to make sure it worked correctly, since errors here will cause errors later anyway
                ongoingCount += 1
        return outlierDict

def outlier_dict_merge(dict1, dict2, strictness):
        # Ensure that strictness value is sensible
        if strictness.lower() not in ['relaxed', 'strict']:
                print('outlier_dict_merge: strictness value is not recognised; it should be "relaxed" or "strict", not "' + strictness.lower() + '".')
                print('Fix the code for this section.')
                quit()
        # Ensure dict inputs are compatible
        if dict1.keys() != dict2.keys():
                print('outlier_dict_merge: dict1 and dict2 don\'t have identical keys. This shouldn\'t be true if they were produced using the same MSA file name list.')
                print('Fix your input and try again.')
                quit()
        # Merge into a new output dict
        mergedDict = {}
        for key in dict1.keys():
                value1 = dict1[key]
                value2 = dict2[key]
                if len(value1) != len(value2):
                        print('outlier_dict_merge: values within dict1 and dict2 don\'t have identical length. This shouldn\'t be true if they were produced using the same MSA file name list.')
                        print('Fix your input and try again. A bit of debug info is below.')
                        print('Key = ' + str(key) + '. Value1 = ' + str(value1) + '. Value2 = ' + str(value2) + '.')
                        quit()
                mergedList = []
                for i in range(len(value1)):
                        if strictness.lower() == 'relaxed':
                                if value1[i] == True and value2[i] == True:     # If both outlier detection methods agree that it is an outlier, we mark it as outlier when strictness == relaxed
                                        mergedList.append(True)
                                else:                                           # If they disagree whether it's an outlier or both agree it is not an outlier, we marked it as not outlier
                                        mergedList.append(False)
                        elif strictness.lower() == 'strict':
                                if value1[i] == True or value2[i] == True:     # If either outlier detection methods identify it as an outlier, we mark it as outlier when strictness == strict
                                        mergedList.append(True)
                                else:
                                        mergedList.append(False)
                mergedDict[key] = mergedList
        return mergedDict

def curate_msa_from_outlier_dict(outlierDict, msaFileNameList):
        # Set up
        from Bio import AlignIO
        from Bio.Align import MultipleSeqAlignment
        from Bio.Seq import Seq
        from Bio.Alphabet import SingleLetterAlphabet
        # Ensure msaFileNameList is correctly formatted
        if type(msaFileNameList) == str:
                msaFileNameList = [msaFileNameList]
        elif type(msaFileNameList) != list:
                print('curate_msa_from_outlier_dict: msaFileNameList type is not recognisable. It should be a list, but instead it is ' + str(type(msaFileNameList)) + '.')
                print('Fix your input and try again.')
                quit()
        if msaFileNameList == []:
                print('curate_msa_from_outlier_dict: msaFileNameList is empty. I don\'t know what to do in this situation since it shouldn\'t happen.')
                print('Code your call to this function properly to skip it.')
                quit()
        # Ensure outlierDict and msaFileNameList are compatible
        if len(outlierDict) != len(msaFileNameList):
                print('curate_msa_from_outlier_dict: msaFileNameList length is not identical to outlierDict length. This shouldn\'t be true if they were produced using the same MSA file name list.')
                print('Fix your input and try again. A bit of debug info is below.')
                print('Len outlierDict = ' + str(len(outlierDict)) + '. Len msaFileNameList = ' + str(len(msaFileNameList)) + '.')
                quit()
        # Loop through msaFileNameList and make modifications to MSAs in place
        for i in range(len(msaFileNameList)):
                # Parse MSA and make sure it corresponds to outlierDict correctly
                msa = AlignIO.read(msaFileNameList[i], 'fasta')
                outlierList = outlierDict[os.path.basename(msaFileNameList[i])]
                if len(msa) != len(outlierList):
                        print('curate_msa_from_outlier_dict: MSA file "' + msaFileNameList[i] + '" length is not the same as outlierDict entry. This shouldn\'t be true if they were produced using the same MSA file name list.')
                        print('Fix your input and try again.')
                        quit()
                newMsa = MultipleSeqAlignment([])
                for y in range(len(msa)):
                        if outlierList[y] == False:
                                newMsa.append(msa[y])
                # If we dropped sequences, ensure that our newMsa still has more than one entry in it [Note: this should technically never happen]
                if len(newMsa) < 2:
                        print('curate_msa_from_outlier_dict: MSA file "' + msaFileNameList[i] + '" after outlier removal has less than two entries. This shouldn\'t be possible, so something is wrong with the code.')
                        print('A bit of debug info is below.')
                        print('Len msa = ' + str(len(msa)) + '. Len newMsa = ' + str(len(newMsa)) + '. outlistList = ' + str(outlierList) + '.')
                        quit()
                # If we dropped sequences, make sure we don't have any blank columns now
                if len(newMsa) < len(msa):
                        for a in range(len(newMsa[0].seq), 0, -1):
                                col = newMsa[:,a-1]
                                if set(col) == {'-'}:
                                        for b in range(len(newMsa)):
                                                newMsa[b].seq = Seq(str(newMsa[b].seq)[0:a-1] + str(newMsa[b].seq)[a:], SingleLetterAlphabet())
                # Produce our updated MSA fasta file
                with open(msaFileNameList[i], 'w') as fileOut:
                        fileOut.write(newMsa.format('fasta'))

## start_trim
def msa_start_find(msaFastaIn, minLength, outType, msaFastaOut, skipOrDrop):
        # Set up
        import copy
        from collections import Counter
        from Bio import AlignIO
        from Bio.Seq import Seq
        from Bio.Alphabet import SingleLetterAlphabet
        # Define functions integral to this one
        def common_start_find(msa):
                # Loop through each sequence and find its start site
                startSiteCount = [[i, 0] for i in range(len(msa[0]))]
                for i in range(len(msa)):
                        for x in range(len(msa[i])):
                                if msa[i][x] == '-':
                                        continue
                                else:
                                        startSiteCount[x][1] += 1
                                        break
                # Identify the most common start site(s)
                startSiteCount.sort(key = lambda x: -x[1])
                commonStarts = []
                for startSite in startSiteCount:
                        if commonStarts == []:
                                commonStarts.append([startSite[0]])
                                prevStartNum = startSite[1]
                        else:
                                if startSite[1] == prevStartNum:
                                        commonStarts.append([startSite[0]])
                                else:
                                        break
                startSiteCount.sort()   # We want this sorted by index again
                # Associate the amino acid for these starts
                for i in range(len(commonStarts)):
                        col = msa[:,commonStarts[i][0]].replace('-', '')
                        aaCount = Counter(col)
                        maxCount = list(aaCount.most_common(1)[0])[1]
                        chars = ''
                        for aa, count in aaCount.items():
                                if count == maxCount:
                                        chars += aa
                        commonStarts[i].append(chars)
                # Extract information from commonStarts in more useable format
                commonStartAA = ''
                for commonPair in commonStarts:
                        if commonPair[1] not in commonStartAA:
                                commonStartAA += commonPair[1]
                return startSiteCount, commonStarts, commonStartAA
        # Ensure outType is sensible
        if outType.lower() not in ['obj', 'file', 'both']:
                print('msa_trim: This function requires an outType to be specified with a specific format.')
                print('Your provided value "' + outType + '" should instead be "obj", to return the Biopython MSA object, "file" to produce an output file which uses the string provided as msaFastaOut, or "both" to do both aforementioned things.')
                print('Format this correctly and try again.')
                quit()
        if (outType.lower() == 'file' or outType.lower() == 'both') and not type(msaFastaOut) == str:
                print('msa_trim: You specified a file output but didn\'t provide a string for msaFasta out - the file output name.')
                print('Format this correctly and try again.')
                quit()
        # Process minLength and ensure it is sensible
        try:
                int(minLength)
        except:
                print('msa_start_find: minLength must be an integer or capable of conversion to integer.')
                print('Format this correctly and try again.')
                quit()
        if minLength < 0:
                print('msa_start_find: minLength must be greater than 0.')
                print('Format this correctly and try again.')
                quit()
        # Load in fasta file as MSA object
        msa = AlignIO.read(msaFastaIn, 'fasta')
        # Count start sites and identify common starts
        startSiteCount, commonStarts, commonStartAA = common_start_find(msa)
        # Loop through each sequence and identify its start position with ranked checking system
        prevMsa = None
        while True:
                # Loop exit condition: no changes were made
                if prevMsa != None:
                        if prevMsa.format('fasta') == prevMsa.format('fasta'):
                                break
                prevMsa = copy.deepcopy(msa)
                # Loop through MSA, identifying start sites for changing and doing this
                for i in range(len(msa)):
                        msaSeq = str(msa[i].seq)
                        # Identify the sequence start site
                        for x in range(len(msaSeq)):
                                if msaSeq[x] == '-':
                                        continue
                                else:
                                        startSite = x
                                        break
                        # Identify start site candidates with associated evidence
                        startCandidates = []
                        for x in range(len(msaSeq)):
                                if (startSiteCount[x][1] - 1 > 0 and msaSeq[x] in commonStartAA) or msaSeq[x] == 'M' or x == startSite: # i.e., if this position is the same AA as others at a common start site OR it's an M OR it's the original start site
                                        col = prevMsa[:,x]      # We want to work with the MSA as a snapshot of what it was at the start of the iteration rather than have changes effect results internally
                                        # Calculate the proportion of start sites that support this one
                                        if x == startSite:
                                                startProp = (startSiteCount[x][1] - 1) / len(msa)       # -1 since we don't want to include this sequence in its own proportion calculation
                                        else:
                                                startProp = (startSiteCount[x][1]) / len(msa)
                                        # Calculate the proportion of common start AAs at this position
                                        tmpCount = 0
                                        for letter in list(set(commonStartAA + 'M')):                   # We always want M to be considered as a common start at this point
                                                tmpCount += col.count(letter)
                                        aaProp = (tmpCount - 1) / len(msa)                              # -1 as before
                                        # Calculate the distance of this from a common start site
                                        dist = None
                                        for commonPair in commonStarts:
                                                if dist == None:
                                                        dist = abs(x-commonPair[0])
                                                else:
                                                        tmpDist = abs(x-commonPair[0])
                                                        if tmpDist < dist:
                                                                dist = tmpDist
                                        # Store results
                                        #startCandidates.append([x, startProp, aaProp, dist])
                                        startCandidates.append([x, startProp, dist, aaProp])
                        ### Sort candidates based on weighting of ranks in order of priority startProp > aaProp > dist > length
                        # Sort candidates based on weighting of ranks in order of priority startProp > dist > aaProp > length
                        #startCandidates.sort(key = lambda x: (-x[1], -x[2], x[3], x[0]))
                        startCandidates.sort(key = lambda x: (-x[1], x[2], -x[3], x[0]))
                        # Generate new sequence & store as MSA
                        msaSeq = Seq('-' * startCandidates[0][0] + msaSeq[startCandidates[0][0]:], SingleLetterAlphabet())
                        #newSeq = SeqRecord(Seq(msaSeq, SingleLetterAlphabet()), id=msa[i].id, name=msa[i].name, description=msa[i].description)
                        msa[i].seq = msaSeq
                # New common start details for iteration
                startSiteCount, commonStarts, commonStartAA = common_start_find(msa)
        # Update the MSA to remove any empty columns
        for a in range(len(msa[0].seq), 0, -1):
                col = msa[:,a-1]
                if set(col) == {'-'}:
                        for b in range(len(msa)):
                                msa[b].seq = Seq(str(msa[b].seq)[0:a-1] + str(msa[b].seq)[a:], SingleLetterAlphabet())
        # Compare our MSA length post-trimming to our specified cut-offs to determine whether we're doing anything to this sequence or not
        seqLen = len(msa[0])
        if minLength < 1:
                ratio = seqLen / len(msa[0])
                if ratio < minLength:
                        if skipOrDrop.lower() == 'skip':
                                print('#"' + os.path.basename(msaFastaIn) + '" trimming reduces length more than minLength proportion cut-off; no trimming will be performed.')
                                return msa      # We're not going to make any changes if trimming shortens it too much
                        elif skipOrDrop.lower() == 'drop':
                                print('#"' + os.path.basename(msaFastaIn) + '" trimming reduces length more than minLength proportion cut-off; msa will be dropped.')
                                return None
        else:
                if seqLen < minLength:
                        if skipOrDrop.lower() == 'skip':
                                print('#"' + os.path.basename(msaFastaIn) + '" trimming reduces length more than absolute minLength cut-off; no trimming will be performed.')
                                return msa      # As above
                        elif skipOrDrop.lower() == 'drop':
                                print('#"' + os.path.basename(msaFastaIn) + '" trimming reduces length more than absolute minLength cut-off; msa will be dropped.')
                                return None
        # Return results either as the MSA object, as an output file, or as both
        if outType.lower() == 'file' or outType.lower() == 'both':
                with open(msaFastaOut, 'w') as fileOut:
                        fileOut.write(msa.format('fasta'))
        if outType.lower() == 'obj' or outType.lower() == 'both':
                return msa

## General purpose funtions
def file_name_gen(prefix, suffix):
        import os
        ongoingCount = 2
        while True:
                if not os.path.isfile(prefix + '1' + suffix):
                        return prefix + '1' + suffix
                elif os.path.isfile(prefix + str(ongoingCount) + suffix):
                        ongoingCount += 1
                else:
                        return prefix + str(ongoingCount) + suffix

# Customise help messages based on current user input
modeHandling = None
if '-m outliers' in ' '.join(sys.argv) or '-mode outliers' in ' '.join(sys.argv):
        modeHandling = 'outliers'
elif '-m start_trim' in ' '.join(sys.argv) or '-mode start_trim' in ' '.join(sys.argv):
        modeHandling = 'start_trim'
elif '-m full_trim' in ' '.join(sys.argv) or '-mode full_trim' in ' '.join(sys.argv):
        modeHandling = 'full_trim'

#### USER INPUT SECTION
usage = """%(prog)s will curate multiple sequence alignments (MSAs)...

If a single -i value is provided, this program assumes that any file in this
directory with contents starting with the '>' character is a FASTA MSA file; if
this proves untrue, errors will occur. If multiple values are provided to -i,
it is assumed that you have specified FASTA MSA files individually.

Mode information: 'outliers' will automatically detect and remove outlier sequences
from provided FASTA MSA files; a detailed log will be produced in the output
directory indicating what changes were made. 'start_trim' will attempt to trim protein
sequences to equalise their start amino acid relative to each other; this is useful
when trying to identify the correct CDS start positions in a MSA. 'full_trim' will
attempt to trim all sequences in a MSA to capture the most "central" or "conserved"
region within the MSA; this function is likely performed better by some other
available programs, have a look online before trying this.

Note that 'start_trim' is designed specifically to handle proteins, whereas the
other systems will handle the underlying sequence agnostically to whether it is
protein or DNA.

Tips and tricks: "full_trim" with stricter dropDrop (-d) value can be used as an
alternative to "outliers"
"""

# Reqs
p = argparse.ArgumentParser(description=usage, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("-i", "-input", dest="inputLocation", nargs="+",
               help="Input a single directory or multiple file locations")
p.add_argument("-o", "-output", dest="outputLocation",
               help="Specify the location to write output files to")
p.add_argument("-m", "-mode", dest="mode", choices=['outliers', 'start_trim', 'full_trim'],
               help="Specify the program mode to run in; details are provided above")
## Opts
p.add_argument("-t", "-threads", dest="threads", type=int, default=1,
               help="Optionally specify the number of threads to run functions with (default == 1).")
## Opts with visibility settings
### outliers
p.add_argument("-r", "-rscriptdir", dest="rscriptdir", type=str,
               help="""outliers: Specify the R directory that contains Rscript.exe. If this is already in your PATH,
               you can leave this blank."""
               if (modeHandling == 'outliers' or modeHandling == None) else argparse.SUPPRESS)
p.add_argument("-s", "-strictness", dest="strictness", choices=['relaxed', 'strict'], default='relaxed',
               help="""outliers: Specify the strictness with which outliers will be removed; 'relaxed' settings
               require two different lines of evidence to suggest a sequence is an outlier, whereas 'strict'
               settings only require one, thus removing more sequences"""
               if (modeHandling == 'outliers' or modeHandling == None) else argparse.SUPPRESS)
### full_trim
p.add_argument("-p", "-propTrim", dest="propTrim", type=float, default=0.7,
               help="""full_trim: Specify the proportion of a MSA column that must contain
               sequences (rather than gaps) to mark the start and end positions for trimming
               (default == 0.7)"""
               if (modeHandling == 'full_trim' or modeHandling == None) else argparse.SUPPRESS)
p.add_argument("-d", "-dropProp", dest="dropProp", type=float, default=0.5,
               help="""full_trim: Specify the proportion of a sequence allowed to disagree with
               consensus (i.e., whether a column is a gap or a sequence) before removal of the sequence
               (default == 0.5)"""
               if (modeHandling == 'full_trim' or modeHandling == None) else argparse.SUPPRESS)
### start_trim & full_trim
p.add_argument("-l", "-lengthMin", dest="lengthMin", type=float, default=0.25,
               help="""full_trim: Specify the minimum length a MSA must be for trimming to occur;
               providing an integer > 1 will enforce a minimum length as a value, whereas providing
               a 0 < float >= 1 value enforce a minimum proportion relative to the original MSA
               (default == 0.25)"""
               if (modeHandling == 'full_trim' or modeHandling == 'start_trim' or modeHandling == None) else argparse.SUPPRESS)
## signalP functionality on/off

args = p.parse_args()
## HARDCODED TESTING
args.inputLocation = ['E:\genome\small_peptide_annot\orthofinder\manually_curated_aligns']
args.outputLocation = r'E:\genome\small_peptide_annot\orthofinder\manually_curated_aligns\outliers'
args.mode = 'start_trim'
#args.rscriptdir
args = validate_args(args)

# Find FASTA files depending on how inputLocation was specified
msaFileNameList = []
if len(args.inputLocation) == 1:
        # Scan through files and detect our files of interest
        for file in os.listdir(args.inputLocation[0]):
                file = os.path.join(args.inputLocation[0], file)
                if not os.path.isfile(file):
                        continue
                with open(file, 'r') as fileIn:
                        for line in fileIn:
                                if line.startswith('>'):
                                        msaFileNameList.append(os.path.abspath(file))
                                break
        # Ensure that we found some input files
        if msaFileNameList == []:
                print('I did not find any FASTA files in the provided input directory "' + args.inputLocation[0] + '"')
                print('Make sure you specified the correct location, or make sure all files are located at this directory. Program will exit now.')
                quit()
else:
        # Make sure that the provided files all exist
        for file in args.inputLocation:
                if not os.path.isfile(file):
                        print('Input file "' + file + '" either does not exist or is not a file.')
                        print('Make sure you spelled the file name/location correctly and try again.')
                        quit()
                with open(file, 'r') as fileIn:
                        for line in fileIn:
                                if not line.startswith('>'):
                                        print('Input file "' + file + '" does not appear to be FASTA formatted i.e., it lacks the ">" character at its start.')
                                        print('Make sure you spelled the correct file or fix this file and try again.')
                                        quit()
                                break
                msaFileNameList.append(os.path.abspath(file))

# Make the output directory if it doesn't exist
if not os.path.isdir(args.outputLocation):
        os.mkdir(args.outputLocation)
else:
        # Ensure that file overwrites won't happen
        for file in msaFileNameList:
                if os.path.isfile(os.path.join(args.outputLocation, os.path.basename(file))):
                        print('"' + os.path.basename(file) + '" file already exists in output directory "' + os.path.abspath(args.outputLocation) + '"')
                        print('This program will not overwrite existing files; delete/rename/move this existing file or specify a new output directory location and try again.')
                        quit()

# Make output files for later modification in place
finalMsaFileList = []
for file in msaFileNameList:
        shutil.copy(file, args.outputLocation)
        finalMsaFileList.append(os.path.join(args.outputLocation, os.path.basename(file)))

# Enact relevant function
## Outliers
if args.mode == 'outliers':
        odSeqResults = odseq_outlier_detect(msaFileNameList, args.rscriptdir, args.outputLocation, 0.01, 'affine', 1000)        # Values are somewhat arbitrary; 0.01 refers to ODseq threshold and seems to be appropriate, 'affine' means we will penalise gaps, and 1000 is the number of bootstrap replicates - it seems to be fast enough so the large number isn't a concern
        spOutResults = msa_outlier_detect(msaFileNameList, True, True, 'alfree')                                                # The first True here means we use basic distribution statistics to justify HDBSCAN's outlier prediction; it helps to temper HDBSCAN and should thus be turned on
        mergedOutlierResults = outlier_dict_merge(odSeqResults, spOutResults, args.strictness)                                  # The second True above removes identical sequences for the purpose of calculating distribution statistics; this helps to prevent a domain being overwhelmed by identicality (it's a word, look it up)
        curate_msa_from_outlier_dict(mergedOutlierResults, finalMsaFileList)
elif args.mode == 'full_trim':
        for msaFileName in finalMsaFileList:
                msa_trim(msaFileName, args.propTrim, args.lengthMin, 'file', msaFileName, args.dropProp, 'skip')
elif args.mode == 'start_trim':
        for msaFileName in finalMsaFileList:
                msa_start_find(msaFileName, args.lengthMin, 'file', msaFileName, 'skip')

        
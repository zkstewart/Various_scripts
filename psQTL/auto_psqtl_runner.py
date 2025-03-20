#! python3
# auto_psqtl_runner.py
# Script for running psQTL evaluation on the simulated data.

import os, sys, subprocess
from multiprocessing import Process, Pipe, Queue

class BasicProcess(Process):
    def __init__(self, *args, **kwargs):
        Process.__init__(self)
        self.args = args
        self.kwargs = kwargs
        self._exception_receiver, self._exception_sender = Pipe()
        self.exception = None
    
    def run(self):
        try:
            self.task(*self.args, **self.kwargs)
        except Exception as e:
            self._exception_sender.send(e)
    
    def task(self, *args, **kwargs):
        # Override this method in a subclass
        raise NotImplementedError
    
    def check_errors(self):
        if self._exception_receiver.poll():
            self.exception = self._exception_receiver.recv()
        if self.exception:
            raise self.exception

def psqtl_init(location, psqtlDir):
    # Format and run the command
    cmd = ["python", os.path.join(psqtlDir, "psQTL_prep.py"), "initialise", "-d", location,
           "--meta", os.path.join(location, "metadata.txt"), "--fvcf", os.path.join(location, "variants.vcf")]
    run = subprocess.Popen(" ".join(cmd), shell = True,
                           stdout = subprocess.PIPE,
                           stderr = subprocess.PIPE)
    out, err = run.communicate()
    if err.decode("utf-8") != "":
        raise Exception(("ERROR: prep-init encountered an error; have a look " +
                        f'at the stdout ({out.decode("utf-8")}) and stderr ' + 
                        f'({err.decode("utf-8")}) to make sense of this.'))

def psqtl_proc(location, psqtlDir):
    # Format and run the command    
    cmd = ["python", os.path.join(psqtlDir, "psQTL_proc.py"), "call", "-d", location,
           "--ignoreIdentical"]
    run = subprocess.Popen(" ".join(cmd), shell = True,
                           stdout = subprocess.PIPE,
                           stderr = subprocess.PIPE)
    out, err = run.communicate()
    if err.decode("utf-8") != "":
        raise Exception(("ERROR: proc-call encountered an error; have a look " +
                        f'at the stdout ({out.decode("utf-8")}) and stderr ' + 
                        f'({err.decode("utf-8")}) to make sense of this.'))

def psqtl_plot(location, genomeFile, psqtlDir):
    params = location.split("/", maxsplit=7)[-1]
    
    # Format and run the command    
    cmd = ["python", os.path.join(psqtlDir, "psQTL_post.py"), "plot", "-d", location,
           "-i", "call", "-f", genomeFile, "-o", os.path.join(location, f"{params.replace('/', '_')}.png"),
           "--power", "4", "-p", "scatter", "line"]
    run = subprocess.Popen(" ".join(cmd), shell = True,
                           stdout = subprocess.PIPE,
                           stderr = subprocess.PIPE)
    out, err = run.communicate()
    if err.decode("utf-8") != "":
        raise Exception(("ERROR: post-plot encountered an error; have a look " +
                        f'at the stdout ({out.decode("utf-8")}) and stderr ' + 
                        f'({err.decode("utf-8")}) to make sense of this.'))

class psQTLProcess(BasicProcess):
    '''
    Parameters:
        location -- a string indicating the location of the VCF and metadata files.
        genomeFile -- a string indicating the location of the genome file.
        psqtlDir -- a string indicating the location of the psQTL scripts.
    '''
    def task(self, location, genomeFile, psqtlDir):
        psqtl_init(location, psqtlDir)
        psqtl_proc(location, psqtlDir)
        psqtl_plot(location, genomeFile, psqtlDir)

# Other functions
def psqtl_pipe(locations, genomeFile, psqtlDir, threads):
    '''
    Parameters:
        locations -- a list of directories containing the VCF and metadata files.
        genomeFile -- a string indicating the location of the genome file.
        psqtlDir -- a string indicating the location of the psQTL scripts.
        threads -- an integer indicating how many threads to run GMAP with.
    '''
    for i in range(0, len(locations), threads): # only process n (threads) directories at a time
        processing = []
        for x in range(threads): # begin processing n files
            if i+x < len(locations): # parent loop may excess if n > the number of directories
                location = locations[i+x]
                
                psqtlWorkerThread = psQTLProcess(location, genomeFile, psqtlDir)
                psqtlWorkerThread.start()
                processing.append(psqtlWorkerThread)
        
        # Gather results
        for psqtlWorkerThread in processing:
            psqtlWorkerThread.join()
            psqtlWorkerThread.check_errors()

def main():
    genomeFile = "genome.fasta"
    psqtlDir = "/home/stewarz2/scripts/psQTL"
    
    # Locate the analysis directories
    locations = []
    for root, dirs, files in os.walk(os.getcwd()):
        if "metadata.txt" in files and "variants.vcf" in files:
            locations.append(root)
    
    # Run the psQTL pipeline
    psqtl_pipe(locations, genomeFile, psqtlDir, 12)

if __name__ == "__main__":
    main()

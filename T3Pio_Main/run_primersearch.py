import argparse
import subprocess
import shutil
import sys
import os
import glob
import concurrent.futures
from pathlib import Path


mismatch_percent = 6


def runPrimerSearch(seq_file, primer_file, output_file):
    '''
    check if primersearch is on path first, then run primersearch on the given sequence file and primer list file
    return the valid primersearch output file, otherwise return None if the command was not run successfully
    '''
    
    if shutil.which('primersearch') is None:
        print (f"primersearch is not found on PATH")
        sys.exit()
        
    command = ['primersearch','-seqall',seq_file,'-infile',primer_file,
               '-mismatchpercent',f'{mismatch_percent}','-outfile',output_file]
    process = subprocess.run(command, capture_output=True, text=True)
    if process.returncode == 0:
        return(output_file)
    

def parse_argument():
    # note
    # usage: python3 run_primersearch.py 
    #       -i input_folder -f input_file_type -o output_folder -p primer_list -t num_threads
    parser = argparse.ArgumentParser(prog = 'run_primersearch.py')
    
    parser.add_argument('-i', '--input', metavar = '', required = True, help = 'Specify input folder')
    parser.add_argument('-o', '--output', metavar = '', required = True, help = 'Specify output folder')
    parser.add_argument('-f', '--format', metavar = '', required = False, help = 'Specify input file format')
    parser.add_argument('-p', '--primers', metavar = '', required = True, help = 'Specify primers list file')
    parser.add_argument('-t', '--thread', metavar = '', required = False, help = 'Specify number of threads')
    
    return parser.parse_args()
    
if __name__ == "__main__":
    
    args = parse_argument()
    
    input_folder = args.input  
    output_folder = args.output
    input_format = args.format
    primerlist_file = args.primers
    num_thread = args.thread
    
    fasta_list = [] #list of all the sequence files
    if os.path.isdir(input_folder):
        if input_format:
            fasta_list.extend(list(glob.glob(f"{input_folder}/*.{input_format}")))
        else:
            fasta_list.extend(list(glob.glob(f"{input_folder}/*.fasta"))) #default to fasta format
    else:
        fasta_list.append(input_folder) #then input_folder is just a file
        
    def concurrent_work(fasta_file):
        
        file_base = Path(fasta_file).stem # basename of the fasta file without .fasta extension
        print (f"running primersearch on: {file_base}")
        runPrimerSearch(fasta_file, primerlist_file, f"{output_folder}/{file_base}.ps")
    
    if num_thread:
        num_workers = int(num_thread)
    else:
        num_workers = 20
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
        results = executor.map(concurrent_work, fasta_list)
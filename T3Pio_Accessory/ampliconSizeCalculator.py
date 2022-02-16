def primerSearchRunner(target,primerFile,mismatch):

    primerSearchFile = primerFile.split('/')[-1]

    primerSearchFile = primerSearchFile+'.ps'

    return_code = subprocess.check_output(['primersearch','-seqall',target,'-infile',primerFile,'-mismatchpercent',mismatch,'-outfile',primerSearchFile])

    return primerSearchFile






def parallelRuner(primerFiles):

    pool = multiprocessing.Pool(processes=args.numCores)

    parallelStatic = partial(parallelFunctions,target=args.target,mismatch=args.mismatch)

    result_list = pool.map(parallelStatic,primerFiles)

    pool.close()

    return result_list


import Bio.Emboss.PrimerSearch
import argparse
import glob
import subprocess
import multiprocess
from functools import partial

parser = argparse.ArgumentParser(description='This script is designed to intake a target sequence and primer sequences and run primersearch. It will return a file with amplicon sizes and associated primer pair.')

parser.add_argument('target',type=str,help='Enter the path to the target sequence file')

parser.add_argument('primerDir',type=str,help='Enter the path to the directory containing the primer sequence files')

parser.add_argument('mismatch',type=str,help='Enter the mismatch percent for primersearch to use')

parser.add_argument('numCores',type=int,help='Enter the number of cores for the process to use')

parser.add_argument('outFile',type=str,help='Enter the name of the outFile')

args = parser.parse_args()


if __name__ == '__main__':

    
    primerFiles = sorted(glob.glob(args.primerDir+'*'))

    

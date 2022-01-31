def primersearchRunner(genome,primerFile,outDir):

    outName = primerFile.split('/')[-1]
    
    outName = outName[:9]

    outFile = outDir+'/'+outName

    return_code = subprocess.check_output(['primersearch','-seqall',genome,'-infile',primerFile,'-mismatchpercent','6','-outfile',outFile])


import argparse
import subprocess
import os
import glob

parser = argparse.ArgumentParser(description='Runs primersearch')

parser.add_argument('genome',type=str,help='Path to genome')

parser.add_argument('primerDir',type=str,help='Path to primer dir')

parser.add_argument('outDir',type=str,help='Directory name for output files')

args = parser.parse_args()



if __name__ == '__main__':
    
    if not os.path.exists(args.outDir):
        os.makedirs(args.outDir,exist_ok=True)

    primerFiles = sorted(glob.glob(args.primerDir+'*'))

    for primers in primerFiles:

        primersearchRunner(args.genome,primers,args.outDir)

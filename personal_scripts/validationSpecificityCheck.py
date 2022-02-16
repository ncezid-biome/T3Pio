class OutputRecord:

    def __init__(self):

        self.amplifiers={}

class Amplifier:

    def __init__(self):

        self.hit_info = ''
        self.length = 0

    
def readPrimersearch(handle):

    record = OutputRecord()

    for line in handle:
        if not line.strip():
            continue
        elif line.startswith('Primer name'):
            name = line.split()[-1]
            record.amplifiers[name] = []
        elif line.startswith('Amplimer'):
            amplifier = Amplifier()
            record.amplifiers[name].append(amplifier)
        elif line.startswith('\tSequence: '):
            amplifier.hit_info = line.replace('\tSequence: ','')
        elif line.startswith('\tAmplimer length: '):
            length = line.split()[-2]
            amplifier.length = int(length)
        else:
            amplifier.hit_info += line
        
    for name in record.amplifiers:
        for amplifier in record.amplifiers[name]:
            amplifier.hit_info = amplifier.hit_info.rstrip()

    return record


def primersearchRunner(assembly,primerFile,outDir):

    outfileName = primerFile.split('/')[-1].split('.')[0]

    outfileName = outfileName[:9]

    outfileName = outDir+outfileName

    return_code = subprocess.check_output(['primersearch','-seqall',assembly,'-infile',primerFile,'-mismatchpercent','6','-outfile',outfileName])


    return outfileName


def fileReader(File):

    with open(File,'r') as f:

        fileInfo = f.readlines()

    f.close

    return fileInfo


def parallelFunctions(primerFile,assembly,outDir):

    primersearchFile = primersearchRunner(assembly,primerFile,outDir)

    primersearchInfo = fileReader(primersearchFile)

    subprocess.check_output(['rm',primersearchFile])

    primersearchRecord = readPrimersearch(primersearchInfo)

    return primersearchRecord




def parallelRunner(primerFiles):

    pool = multiprocessing.Pool(processes=args.numCores)

    parallelStatic = partial(parallelFunctions,assembly=assembly,outDir=outDir)

    result_list = pool.map(parallelStatic,primerFiles)

    pool.close()

    return result_list


import argparse
import multiprocessing
from functools import partial
import glob
import subprocess
import os


parser = argparse.ArgumentParser(description='Runs primersearch on assembly files and returns collalated results as well as list of primers hitting non-targeted genomes')

parser.add_argument('primerFilePath',type=str,help='Path to directory contatining primer files')

parser.add_argument('assemblyPath',type=str,help='Path to directory containing assembly files')

parser.add_argument('outDir',type=str,help='Name of the outDir to store results')

parser.add_argument('numCores',type=int,help='Number of cores to use with multiprocess')

#parser.add_argument('assemblyID',type=str,help='Tab delim file identifing the assemblies and their genome ID')

#parser.add_argument('target',type=str,help='Primer targeted species')

args = parser.parse_args()


if __name__ == '__main__':


    primerFiles = sorted(glob.glob(args.primerFilePath+'*'))

    assemblies = sorted(glob.glob(args.assemblyPath+'*'))

    
    for assembly in assemblies:

        assemblyName = assembly.split('/')[-1].split('.')[0]

        outDir = args.outDir+'/'+assemblyName+'/'

        
        if not os.path.exists(outDir):

            os.makedirs(outDir,exist_ok=True)


        results = parallelRunner(primerFiles)

        singleHitOut = open(outDir+'singleHitInfo','w')

        multiHitOut = open(outDir+'multiHitInfo','w')

        for result in results:

            for k,v in result.amplifiers.items():

                if len(v) == 1: 

                    print(k,file=singleHitOut)

                    print(v[0].hit_info,file=singleHitOut)
                    print('\t'+str(v[0].length),file=singleHitOut)

                elif len(v) > 1:

                    print(k,file=multiHitOut)

                    for i in v:

                        print(i.hit_info,file=multiHitOut)
                        print('\t'+str(i.length),file=multiHitOut)

        singleHitOut.close()
        multiHitOut.close()

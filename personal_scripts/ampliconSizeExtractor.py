def fastaToDict(fastaFile):

    fastaName = fastaFile.split('/')[-1].split('.')[0]

    fastaDict = SeqIO.to_dict(SeqIO.parse(fastaFile,'fasta'))

    return fastaDict


def primersearchRunner(assemblyFile,primerFile,mismatch,outDir,primersearch):

    psFile = assemblyFile.split('/')[-1].split('.')[0]

    psFile = psFile+'.ps'
    
    return_code= subprocess.run([primersearch,'-seqall',assemblyFile,'-infile',primerFile,'-mismatchpercent',mismatch,'-outfile',outDir+'/'+psFile])

    return(outDir+'/'+psFile)



def fileOpener(psFile):

    with open(psFile,'r') as f:

        psContents = f.readlines()

    f.close()

    return psContents

def ampSizeFileWriter(psClass,outFileBaseName):

    outFile = open(outFileBaseName+'--AmpliconSizes','w')

    for k,v in psClass.amplifiers.items():
        if len(v) > 0:
            lengths = []
        
            outFile.write(k+'\t')
        
            for i in v:
                if v.index(i) != len(v)-1:
                    outFile.write(str(i.length)+'\t')
                elif v.index(i) == len(v)-1:
                    outFile.write(str(i.length)+'\n')

    outFile.close()


def ampliconSequenceWriter(psClass,assemDict,outFileBaseName):

    outFile = open(outFileBaseName+'--Amplicons','w')

    for k,v in psClass.amplifiers.items():
        if len(v) > 0:
            for i in v:
            
                i = i.hit_info.split()
                leftHit = int(i[i.index('forward')+3])
                rightHit = int(i[i.index('reverse')+3].strip('[').strip(']'))
                leftPrime = (i[i.index('forward')-2])
                leftLen = len(re.sub(r'\[[^\]]+\]','N',leftPrime))
                rightPrime = (i[i.index('reverse')-2])
                rightLen = len(re.sub(r'\[[^\]]+\]','N',rightPrime))
                targetSequence = str(assemDict[i[0]].seq)
                
                amplicon = targetSequence[leftHit+leftLen-1:len(targetSequence)-rightHit-rightLen+1]
            
                outFile.write('>'+i[0]+'--'+k+'\n')
                outFile.write(amplicon+'\n')


    outFile.close()
            

def parallelFunctions(assembly,primerFile,outDir,keep,mismatch,primersearch):

    outFileBaseName = assembly.split('/')[-1].split('.')[0]

    outFileBaseName = outDir+'/'+outFileBaseName

    psFile = primersearchRunner(assembly,primerFile,mismatch,outDir,primersearch)

    psContents = fileOpener(psFile)

    psClass = Bio.Emboss.PrimerSearch.read(psContents)

    assemDict = fastaToDict(assembly)

    ampSizeFileWriter(psClass,outFileBaseName)
    
    ampliconSequenceWriter(psClass,assemDict,outFileBaseName)

    if keep == 'n':
        
        subprocess.run(['rm',psFile])


def parallelRunner(assemblies):

    pool = multiprocessing.Pool(processes=args.cores)

    parallelStatic = partial(parallelFunctions,primerFile=args.primerFile,outDir=args.outDir,keep=args.keepInterumFiles,mismatch=args.mismatchPercent,primersearch=primersearch)

    result_list = pool.map(parallelStatic,assemblies)

    pool.close()




import argparse
import Bio.Emboss.PrimerSearch
import subprocess
from Bio import SeqIO
import multiprocessing
from functools import partial
import glob
import os
import re

parser = argparse.ArgumentParser(description='This script is designed to take in a directory of assemblies and and a primersearch formated file of primers and run EMBOSS-6.4.0/primersearch on the assemblies and return the amplicons and associated sizes for each of the assemblies')

parser.add_argument('assemblyDir',type=str,help='Enter the path to the directory containing the assemblies')

parser.add_argument('primerFile',type=str,help='Enter the path to the file containing the primers')

parser.add_argument('outDir',type=str,help='Enter the path and name of the directory where the results will be stored')

parser.add_argument('keepInterumFiles',type=str,help='Enter y if you want to keep the output primersearch files or n if you want them removed from the output directory')

parser.add_argument('mismatchPercent',nargs='?',type=str,default='6',help='Enter a new mismatch percent for primersearch to use. The default is set at 6 percent')

parser.add_argument('cores',nargs='?',type=int,default=1,help='Change the default number of cores to use when processing inputs')

parser.add_argument('embossPath',nargs='?',type=str,default='/apps/x86_64/emboss/EMBOSS-6.4.0_bin/bin/primersearch',help='If EBMOSS gets moved to a different location input that new location here')


args = parser.parse_args()


if __name__ == '__main__':

    assemblyFiles = sorted(glob.glob(args.assemblyDir+'*.fasta'))

    primersearch = args.embossPath

    if not os.path.exists(args.outDir+'/'):
        os.makedirs(args.outDir+'/',exist_ok=True)

    parallelRunner(assemblyFiles)

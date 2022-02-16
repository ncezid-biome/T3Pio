def outFileWrite(fastqList,outDir):

    r1 = fastqList[0].split('/')[-1]
    r2 = fastqList[1].split('/')[-1]
    
    outFile1 = outDir+'/'+r1
    outFile2 = outDir+'/'+r2

    system("cat "+fastqList[0]+' '+fastqList[2]+" >> "+outFile1)
    system("cat "+fastqList[1]+' '+fastqList[3]+" >> "+outFile2)


import glob
import argparse
import subprocess
from os import system

parser = argparse.ArgumentParser(description='This script takes a list of directories containing lane 1 and 2 for R1 and R2 and outputs single combined R1s and R2s')

parser.add_argument('InputDirList',type=str,help='txt file containing paths to reads directories')

parser.add_argument('OutputDir',type=str,help='Directory path where the combined reads will go')

args = parser.parse_args()

if __name__ =='__main__':

    with open(args.InputDirList,'r') as f:
        dirList = f.readlines()
    f.close()

    for i in dirList:

        reads = sorted(glob.glob(i.strip('\n')+'/*'))

        outFileWrite(reads,args.OutputDir)

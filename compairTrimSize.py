import argparse
import fnmatch
import glob
import os

parser = argparse.ArgumentParser(description = 'Creates a list of isolates that have appropriately sized trimmed files for assembly')

parser.add_argument('trueReadPairs',type=str,help='File containing a list of reads that have passed previous qc')

parser.add_argument('trimSize',type=str,help='File containing the size of each trimmed reads file')

args = parser.parse_args()

def pairTrimSizeComp(pair, size):
    checkedTrimList = []
    for i in pair:
        iso = i.split('/')
        iso = iso[len(iso)-1]
        iso = iso.split('_L0')[0]
        if 'APHI' in iso:
            try:
                aphiReads = sorted(glob.glob('trimmed/APHI/*'+iso+'*.gz'))
                r1Size = os.path.getsize(aphiReads[0])
                r2Size = os.path.getsize(aphiReads[1])
                readTotalSize = r1Size + r2Size
                readTSizeMB = readTotalSize >> 20
                if readTSizeMB >= 100:
                    holder = str(readTSizeMB) + '\t' + aphiReads[0] + '\t' + aphiReads[1] + '\t' + iso
                    checkedTrimList.append(holder)
                else:
                    continue
            except IndexError:
                continue
        else:
            try:
                name = iso.split('_')[0]
                readInfo = fnmatch.filter(size, '*'+name+'*')
                readPath = readInfo[0].split('\t')[1]
                readPath = readPath.strip('\n')
                reads = sorted(glob.glob(readPath+'/*'+name+'*.gz'))
                r1Size = os.path.getsize(reads[0])
                r2Size = os.path.getsize(reads[1])
                readTSize = r1Size + r2Size
                readTSizeMB = readTSize >> 20
                if readTSizeMB >= 100:
                    holder = str(readTSizeMB)+'\t' + reads[0] + '\t' + reads[1] + '\t' + iso
                    checkedTrimList.append(holder)
                else:
                    continue
            except IndexError:
                continue
    return checkedTrimList

                

with open(args.trueReadPairs,'r') as f:
    readPairs = f.readlines()

with open(args.trimSize,'r') as g:
    sizeFile = g.readlines()

trimmedForAssem = pairTrimSizeComp(readPairs, sizeFile)

h = open('InfoFileForAssembly.txt','a')

for j in trimmedForAssem:
    print(j,file=h)

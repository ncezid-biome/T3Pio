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

def justContigs(contigList):

    """ Creates list of contigs from Cleared contig files generated from blast and coverage process

    Params
    ------
    contigList: List
         List of contig info consisting of contig name/length of contig/percent ID/coverage of reference

    Returns
    -------
    contigs: List
         List of contig names

    """

    contigs = []

    for contigInfo in contigList:

        contig=contigInfo[0].split('-')[0]

        contigs.append(contig)

    return contigs



def multiHitPrimers(contigs,psObj):

    """Takes list of contig names and psObj where primers have hit multiple targets in a genome and assess for target specificity
       and amplicon length

    Params
    ------
    contigs: List
         List of contig names
    psObj: Amplifier objects
         Amplifier objects from OutputRecord object dictionary

    Returns
    -------
    quality: Boolean
         True/False value

    """

    quality = True

    for obj in psObj:

        if obj.hit_info.split()[0] not in contigs:
            
            if obj.length < 1000:
                quality = False

    return quality
        
def multiHitWriter(primer,psObj):

    fileName = open(primer+'--MultiHitInfo','a')

    for obj in psObj:

        print(obj.hit_info,file=fileName)
        print(obj.length,file=fileName)

    fileName.close()

def contigCompare(contigs,psInfo):

    badPrimers = []

    offTargetLarge = []

    for k,v in psInfo.amplifiers.items():

        if len(v) == 1:
            contig = v[0].hit_info.split()[0]
            
            if contig not in contigs:
                if v[0].length < 1000:
                    badPrimers.append(k)
                elif v[0].length >= 1000:
                    offTargetLarge.append(K)

        elif len(v) > 1:
            
            quality = multiHitPrimers(contigs,v)

            multiHitWriter(k,v)

            if quality == False:

                badPrimers.append(k)

    return(badPrimers,offTargetLarge)
            

def fileOpener(inFile):

    with open(inFile,'r') as f:

        fileContents = f.readlines()

    f.close()


    return fileContents

def ParallelFunctions(psFile,contigs):

    psFile = fileOpener(psFile)

    psInfo = readPrimersearch(psFile)

    badPrimers = contigCompare(contigs,psInfo)

    return badPrimers
    

def ParallelRunner(psFiles):
        
    pool = multiprocessing.Pool(processes=args.numCores)

    parallelStatic = partial(ParallelFunctions,contigs=contigs)

    result_list = pool.map(parallelStatic,psFiles)

    pool.close()

    return result_list

import argparse
from csv import reader
import multiprocessing 
from functools import partial
import glob


parser = argparse.ArgumentParser(description='Hopefully a primersearch parser')

parser.add_argument('psFilePath',type=str,help='Path to a primersearch files')

parser.add_argument('psName',type=str,help='Metagenome name associated with psFiles')

parser.add_argument('contigFile',type=str,help='Path to file containing cleared contigs')

parser.add_argument('numCores',type=int,help='Number of cores for multiprocess')

parser.add_argument('outfile',type=str,help='Name of the outfile')

args = parser.parse_args()


if __name__ == '__main__':

    psFiles = sorted(glob.glob(args.psFilePath+'*'+args.psName+'*'))

    with open(args.contigFile,'r') as f:

        csv_reader = reader(f)

        contigList = list(csv_reader)

    f.close()

    contigs = justContigs(contigList)

    results = ParallelRunner(psFiles)

    f = open(args.outfile,'w')

    for primers in results:

        for primer in primers[0]:

            print(primer,file=f)

    f.close()


    f = open(args.outfile+'OffTargetLargeHits','w')

    for primers in results:

        for primer in primers[1]:

            print(primer,file=f)

    f.close()

#'Length of sal amplicons 
#flag long/short amps'

#'Flag off target primers that pass size filter' 

#'align multi hits/ harvest data'

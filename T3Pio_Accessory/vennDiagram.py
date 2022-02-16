def fileReader(primerFile):

    with open(primerFile,'r') as f:

        primerInfo = f.readlines()

    print(primerFile)

    return primerInfo

def listCreator(primerFiles):
    
    primerDict = {}

    primerList = []

    primerOrthList = []

    for primers in primerFiles:

        primerInfo = fileReader(primers)

        for primer in primerInfo:
            
            primer = primer.split('\t')

            primerDict[primer[0]] = [primer[1],primer[2].strip('\n')]

            primerList.append(primer[0])

            primerOrthList.append(primer[0][:9])


    return(primerDict,primerList,primerOrthList)

def listCompare(primerlist1,primerlist2):

    list1Unique = list(set(primerlist1) - set(primerlist2))

    list2Unique = list(set(primerlist2) - set(primerlist1))

    sharedList = list(set(primerlist1) & set(primerlist2))


    return(list1Unique,list2Unique,sharedList)


def fileWriter(primerList,outFile):

    out = open(outFile,'w')

    for info in primerList:

        print(info,file=out)

def primerFileWriter(sharedList,primerDict,outfile):

    out = open(outfile,'w')

    for info in sharedList:

        primers = primerDict[info]

        print(info+'\t'+primers[0]+'\t'+primers[1],file=out)

import argparse
import glob
import os

#list(set(list1).intersection(list2))
#list(set(list1).differences(list2))

parser = argparse.ArgumentParser(description='This script is designed to find the unique and shared primers and orthogroups between two directories containing primer files')

parser.add_argument('primerDir1',type=str, help='Enter path to the first directory')

parser.add_argument('primerDir2',type=str,help='Enter path to the second directory')

parser.add_argument('sharedPrimers',type=str,help='Enter name of the out directory')

parser.add_argument('uniqueList',type=str,help='Enter the name of the file for the unique primer/orth numbers to go')

args = parser.parse_args()


if __name__ == '__main__':

    primers1 = sorted(glob.glob(args.primerDir1+'*'))

    primers2 = sorted(glob.glob(args.primerDir2+'*'))

#    if not os.path.exists(args.sharedPrimerDir+'/'):
#        os.makedirs(args.sharedPrimerDir+'/',exist_ok=True) 


    primerlist1 = listCreator(primers1)

    primerlist2 = listCreator(primers2)

    primerInfo = listCompare(primerlist1[1],primerlist2[1])

    orthInfo = listCompare(primerlist1[2],primerlist2[2])

    uniqueList1File = fileWriter(primerInfo[0],args.uniqueList+'primers1')

    uniqueList1File = fileWriter(primerInfo[1],args.uniqueList+'primers2')

    uniqueList3File = fileWriter(orthInfo[0],args.uniqueList+'orth1')

    uniqueList4File = fileWriter(orthInfo[1],args.uniqueList+'orth2')

    sharedOrths = fileWriter(orthInfo[2],'sharedOrths')

    if primerlist1[2] >= primerlist2[2]:
    
        primerFileWriter(primerInfo[2],primerlist1[0],args.sharedPrimers)

    else:
        
        primerFileWriter(primerInfo[2],primerlist2[0],args.sharedPrimers)

        
    print(len(primerlist1[0]))
    print(len(primerlist2[0]))

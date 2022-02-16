import argparse
import csv

parser = argparse.ArgumentParser(description='This script is for finding the variable bases in isolate binding sites')

parser.add_argument('Primers',type=str,help='Provide a file with primers')

parser.add_argument('Isolates',type=str,help='Provide a file with isolate sequences')

args = parser.parse_args()

class Primers:

    def __init__(self,name,forward,reverse,forwardDegen,reverseDegen,forwardDegenLoc,reverseDegenLoc):
        
        self.name = name
        self.forward = forward
        self.reverse = reverse
        self.forwardDegen = forwardDegen
        self.reverseDegen = reverseDegen
        self.forwardDegenLoc = forwardDegenLoc
        self.reverseDegenLoc = reverseDegenLoc

class Isolate:

    def __init__(self,primerInfo,name,fowardBases,reverseBases):
        
        self.primerInfo = primerInfo
        self.name = name
        self.forwardBases = forwardBases
        self.reverseBases = reverseBases

def DegeneracyLocation(primer):
    
    reverseSequence = primer[::-1]
    
    degeneracyLocationList = []

    degenerateBaseList = []
    
    counter = 0
    
    for base in reverseSequence:
        
        if base not in ('ATGC'):

            degeneracyLocationList.append(counter)
            
            degenerateBaseList.append(base)

        counter +=1
    
    return(degeneracyLocationList,degenerateBaseList)



def PrimerFileParser(primerFile):
    
    primerObjectsList = []

    with open(primerFile,'r') as p:
        
        primerInfo = p.readlines()
        
    p.close()
        
    for primers in primerInfo:
        
        splitPrimers = primers.split('\t')
        
        forwardDegenInfo = DegeneracyLocation(splitPrimers[1])
        
        reverseDegenInfo = DegeneracyLocation(splitPrimers[2])

        primerClass = Primers(splitPrimers[0],splitPrimers[1],splitPrimers[2],forwardDegenInfo[1],reverseDegenInfo[1],forwardDegenInfo[0],reverseDegenInfo[0])

        primerObjectsList.append(primerClass)

    return(primerObjectsList)



def IsolateVariableBases(sequence,primerPositions):
    
    variableBaseList = []
    
    for position in primerPositions:
        
        base = sequence[position]
        
        variableBaseList.append(base)
    
    return(variableBaseList)


primerObjects = PrimerFileParser(args.Primers)


isolateList = []

with open(args.Isolates) as tsvFile:
    tsvreader = csv.reader(tsvFile,delimiter='\t')
    for line in tsvreader:
        isolateList.append(line)

isolateList = isolateList[1:]

isolateDegenInfoList = []

for sequenceData in isolateList:

    for primerObj in primerObjects:
        
        try:
            sequence = sequenceData[primerObjects.index(primerObj)+1]
            
            if sequence != '0' or sequence != 0:
        
                forwardSeq = sequence[:len(primerObj.forward)]

                forwardSeq = forwardSeq[::-1]
        
                reverseSeq = sequence[::-1]
         
                reverseSeq = reverseSeq[:len(primerObj.reverse)]
            
                reverseSeq = reverseSeq[::-1]
            
                forwardBases = IsolateVariableBases(forwardSeq,primerObj.forwardDegenLoc)

                reverseBases = IsolateVariableBases(reverseSeq,primerObj.reverseDegenLoc)
                
                isolateClass = Isolate(primerObj,sequenceData[0],forwardBases,reverseBases)

                isolateDegenInfoList.append(isolateClass)

        except IndexError:
            
            continue
        

#f = open('NucleotideCodex.tsv','w')

masterList = []

for i in isolateDegenInfoList:

    printingList = []

    printingList.append(i.name)
    printingList.append(i.primerInfo.name)
    
    if not i.primerInfo.forwardDegen:
        printingList.append('-')
    else:
        printingList.append(i.primerInfo.forwardDegen)

    if not i.primerInfo.reverseDegen:
        printingList.append('-')
    else:
        printingList.append(i.primerInfo.reverseDegen)

    if not i.primerInfo.forwardDegenLoc:
        printingList.append('-')
    else:
        printingList.append(i.primerInfo.forwardDegenLoc)

    if not i.primerInfo.reverseDegenLoc:
        printingList.append('-')
    else:
        printingList.append(i.primerInfo.reverseDegenLoc)

    if not i.forwardBases:
        printingList.append('-')
    else:
        printingList.append(i.forwardBases)

    if not i.reverseBases:
        printingList.append('-')
    else:
        printingList.append(i.reverseBases)

    masterList.append(printingList)

with open('NucleotideCodex.csv','w') as f:
    writer = csv.writer(f)
    writer.writerows(masterList)

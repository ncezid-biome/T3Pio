import argparse
import glob
import re
import numpy
import csv

parser = argparse.ArgumentParser(description='Find differences in primer binding sites')

parser.add_argument('primerSearch',type=str,help='Input primerSearch Files')

parser.add_argument('fnaFiles',type=str,help='Input files used with primerSearch')

parser.add_argument('primerFiles',type=str,help='Input primer files')


args = parser.parse_args()


fnaDict = {}
primerDict ={}


fnaFiles = sorted(glob.glob(args.fnaFiles+'*'))
primerFiles = sorted(glob.glob(args.primerFiles+'*'))
searchFiles = sorted(glob.glob(args.primerSearch+'*'))

for i in fnaFiles:
    with open(i,'r') as f:
        fnaInfo = f.readlines()
    f.close()
    fnaInfo = [s.strip('\n') for s in fnaInfo]
    fnaInfo = ''.join(fnaInfo)
    delimiters = 'bp','>'
    regex = '|'.join(map(re.escape,delimiters))
    fnaInfo = re.split(regex,fnaInfo)
    i = 1
    while i < len(fnaInfo):
        name = fnaInfo[i].split(' ')[0]
        fnaDict[name] = fnaInfo[i+1]
        i += 2

for i in primerFiles:
    with open(i,'r') as f:
        primerInfo = f.readlines()
    f.close()
    for j in primerInfo:
        j =j.split('\t')
        pName = j[0]
        lSize = len(j[1])
        rSize = j[2].strip('\n')
        rSize = len(rSize)
        pSizes =[lSize,rSize]
        primerDict[pName]=pSizes

outsideDict={}
seqNames = []

for i in searchFiles:
    with open(i,'r') as f:
        searchInfo = f.readlines()
    f.close()
    for j in searchInfo:
        if j.startswith('Primer name'):
            pName = j.split(' ')[2]
            pName = pName.strip('\n')
            outsideDict.setdefault(pName,[])
            del searchInfo[searchInfo.index(j)]
        if j.startswith('\tSequence:'):
            insideDict={}
            seq = j.split(' ')[1]
            seq2 = seq.split('.')[0]
            seqNames.append(seq2)
            forwardHit = searchInfo.index(j)+2
            forwardHit = int(searchInfo[forwardHit].split(' ')[5])
            reverseHit = (searchInfo.index(j)+3)
            reverseHit = searchInfo[reverseHit].split(' ')[5]
            reverseHit = int(reverseHit.replace('[','').replace(']',''))
            hitList = [forwardHit,reverseHit]
            insideDict[seq] = hitList
            outsideDict[pName].append(insideDict)
            del searchInfo[searchInfo.index(j)]


seqSet = set(seqNames)
seqList = list(seqSet)

bindArray = numpy.zeros((len(seqList)+1,len(primerDict)+1),dtype=object)

rowNum = 1
for i in seqList:
    bindArray[rowNum,0] = i
    rowNum += 1
colList = []
colNum=1
for k in primerDict.keys():
    bindArray[0,colNum] = k
    colNum+=1
    colList.append(k)

#
with open('isolatesInStock.txt','r') as f:
    a = f.readlines()
f.close()
aSet = set(a)
aList = list(aSet)
with open('primerZero.txt','r') as f:
    b = f.readlines()
f.close()
joCheckerDict = {}

ampArray = numpy.zeros((len(aList)+2,len(primerDict)+1),dtype=object)

rowNum = 1
for i in aList:
    i = i.strip('\n')
    ampArray[rowNum,0] = i
    rowNum += 1
colNum=1
for k in primerDict.keys():
    ampArray[0,colNum] = k
    colNum+=1


for k1,v1 in primerDict.items():
    print(k1)
    count = 1
    ampList = []
    isoSeqList = outsideDict[k1]
    for i in isoSeqList:
        for k2,v2 in i.items():
            seq = fnaDict[k2]
            seq = seq.replace('-','')
            lBind = seq[v2[0]-1:v2[0]+v1[0]-1]
            seqLen = len(seq)
            rSpot = seqLen - v2[1]
            rBind = seq[rSpot-v1[1]+1:rSpot+1]
            print(lBind+'\t'+rBind)
            sName = k2.split('.')[0]
            rowIndexNum = seqList.index(sName)
            colIndexNum = colList.index(k1)
            bindArray[rowIndexNum+1,colIndexNum+1] = lBind+'\t'+rBind
            for k in aList:
                kn = k.strip('\n')
                if kn in sName:
                    ampList.append(seq[v2[0]-1:rSpot+1])
                    rowIndexNumAmp = aList.index(k)
                    colIndexNumAmp = colList.index(k1)
                    ampArray[rowIndexNumAmp+1,colIndexNumAmp+1] = (seq[v2[0]+v1[0]:rSpot-v1[1]]).strip('\n')
                    iso = ampArray[rowIndexNumAmp+1,0]
                    og = ampArray[0,colIndexNumAmp+1]
                    fastaOut = open('joFastaFiles/'+iso+'.fasta','a')
                    fastaOut.write('>'+iso+'-'+og+'\n'+seq[v2[0]+v1[0]:rSpot-v1[1]]+'\n')
    ampSet = set(ampList)
    ampSetList = list(ampSet)
    joCheckerDict[k1]=ampSetList
    count += 1






'''
numCount =1 
for k,v in joCheckerDict.items():
    numAmps = len(v)
    ampArray[len(aList)+1,numCount]=numAmps
    numCount += 1

with open('AmpliconsForSelectedIsolates','w') as f:
    writer = csv.writer(f)
    writer.writerows(ampArray)

with open('landingSites','w') as f:
    writer = csv.writer(f)
    writer.writerows(bindArray)

diffBindList = []

oneAndZerosList = []

upper = 1

while upper < len(bindArray):
    first = bindArray[upper]
    lower = upper + 1
    while lower < len(bindArray):
        second = bindArray[lower]
        leftDiff = 0
        rightDiff = 0
        diffSites = 0
        spot = 1
        while spot < len(first):
            if first[spot] != 0 and second[spot]:
                fLeft = first[spot].split('\t')[0]
                fRight =  first[spot].split('\t')[1]
                sLeft = second[spot].split('\t')[0]
                sRight = second[spot].split('\t')[1]
                lCount = sum(1 for a,b in zip(fLeft,sLeft) if a != b)
                rCount = sum(1 for c,d in zip(fRight, sRight) if c != d)
                leftDiff += lCount
                rightDiff += rCount
                if lCount > 0 or rCount > 0:
                    diffSites += 1
            spot += 1
        pairDiff = [first[0],second[0],str(leftDiff),str(rightDiff),str(diffSites)]
        diffBindList.append(pairDiff)
        oneAndTwo = [first[0],second[0]]
        oneAndZerosList.append(oneAndTwo)
        print(pairDiff)
        lower+=1
    upper+=1





bitArray = numpy.zeros((len(oneAndZerosList)+1,len(primerDict)+2),dtype=object)

rowNum = 1
for i in oneAndZerosList:
    bitArray[rowNum,0] = i[0]
    bitArray[rowNum,1] = i[1]
    rowNum += 1

colList = []
colNum=2
for k in primerDict.keys():
    bitArray[0,colNum] = k
    colNum+=1
    colList.append(k)

upper = 1
count = 1

while upper < len(bindArray):
    first = bindArray[upper]
    lower = upper + 1
    while lower < len(bindArray):
        second = bindArray[lower]
        spot = 1
        pgSpot = 2
        while spot < len(first):
            if first[spot] != 0 and second[spot]:
                fLeft = first[spot].split('\t')[0]
                fRight =  first[spot].split('\t')[1]
                sLeft = second[spot].split('\t')[0]
                sRight = second[spot].split('\t')[1]
                lCount = sum(1 for a,b in zip(fLeft,sLeft) if a != b)
                rCount = sum(1 for c,d in zip(fRight, sRight) if c != d)
                if lCount > 0 or rCount > 0:
                    bitArray[count,pgSpot] = 1
            spot += 1
            pgSpot += 1
        print(bitArray[count,0])
        lower+=1
        count += 1
    upper+=1

with open('RevisedPrimerBindings2.0','w') as f:
    writer = csv.writer(f)
    writer.writerows(bitArray)



#######################################################
with open('isolatesInStock.txt','r') as f:
    a = f.readlines()

aSet = set(a)
aList = list(aSet)

for i in oneAndZerosList:
    if i[0] not in aList or i[1] not in aList:
        del oneAndZerosList[oneAndZerosList.index(i)]

with open('BindingSiteDifferences','w') as f:
    writer = csv.writer(f)
    writer.writerows(diffBindList)
'''
                    

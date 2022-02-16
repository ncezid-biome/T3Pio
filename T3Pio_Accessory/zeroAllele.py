import numpy
import glob
import os
import csv

with open('faaFiles/Results_Jan10/Orthogroups.csv') as f:
    reader = csv.reader(f)
    data = [r for r in reader]

isolates = (data[0][0].split('\t')[1:])

orthogroups = sorted(glob.glob('correctFilteredPrimerSearch/*'))

OGComparitorDict = {}

zeroList = []


for i in isolates:
    OGComparitorDict[i] = [0]*len(orthogroups)

amp = 0
while amp < len(orthogroups):
    ampDict = {}
    OG = orthogroups[amp].split('/')[1]
    OG = OG.split('.')[0]
    print(OG)
    with open(orthogroups[amp],'r') as a:
        ampInfo = a.readlines()
    a.close()
    ampList = []
    cnt = 0
    while cnt < len(ampInfo):
        ampName=ampInfo[cnt].split('.')[0]
        ampName = ampName.strip('>')
        ampDict[ampName] = ampInfo[cnt+1]
        allele = ampInfo[cnt+1].strip('\n')
        ampList.append(ampInfo[cnt+1])
        cnt += 2
    ampSet = set(ampList)
    ampSet = list(ampSet)
    for k in ampDict.keys():
        alleleNum = ampSet.index(ampDict[k]) + 1
        isoAlleleList = OGComparitorDict[k]
        del isoAlleleList[amp]
        isoAlleleList = isoAlleleList.insert(amp, alleleNum)
    amp += 1


sameAlleles = []
sampleToSampleList = []

for k1 in list(OGComparitorDict.keys()):
    for k2 in list(OGComparitorDict.keys()):
        if k1 != k2:
            if OGComparitorDict.get(k1) == OGComparitorDict.get(k2): #logic to go through and find identical isolates and tag them
                if k2 not in sameAlleles:
                    sameAlleles.append(k2)
                    sampleToSampleList.append(k1 + ' is the same as ' + k2)

print(len(sameAlleles))
f=open('primerPairSameNineFive.txt','a')
for line in sampleToSampleList:
    print(line,file=f)

#numpy.savetxt('pipelineAlleles.csv',OGComparitorDict,delimiter=',')

with open('pipelineAlleles.csv','w') as f:
    for k in OGComparitorDict.keys():
        f.write('%s.%s\n'%(k,OGComparitorDict[k]))


#(pd.DataFrame.from_dict(data=OGComparitorDict,orient='index').to_csv('pipelineAleles.csv',header=False))

#with open('pipelineAlleles.csv','w',newline='') as csv_file:
#    writer = csv.writer(csv_file)
#    for k,v in OGComparitorDict.items():
#        writer.writerow([k,v])

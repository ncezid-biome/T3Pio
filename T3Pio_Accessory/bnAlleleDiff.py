import numpy as np
import csv
import subprocess

allele_Dict = np.load('PrimerPairNintyFiveAllele.npy').item()

#fixAlleles = allele_Dict['Sal_JRB_S06004']
#del allele_Dict['Sal_JRB_S06004']
#allele_Dict['Sal_JRA_S06004'] = fixAlleles

totalAlleleMatrix = np.zeros(((len(allele_Dict)+1),(len(allele_Dict)+1)),dtype=object)

p = open('Sean_Allele_3_Col_Comp','a')
slIsoList = []
orderIsoList = []
print('ISO1'+'\t'+'ISO2'+'\t'+'SL',file=p)
x=1
for k1 in allele_Dict.keys():
    y=1
    totalAlleleMatrix[x][0] = k1
    totalAlleleMatrix[0][x] = k1
    for k2 in allele_Dict.keys():
        list1 = allele_Dict[k1]
        list2 = allele_Dict[k2]
        size = 0
        count = 0
        while size < len(list1):
            if list1[size] != 0 and list2[size] != 0:
                if list1[size] != list2[size]:
                    count += 1
                size += 1
            else:
                size += 1
        totalAlleleMatrix[x][y] = count
        iso1 = k1.strip('>')
        iso2 = k2.strip('>')
        if k2 not in slIsoList:
            print(iso1 + '\t' + iso2 + '\t' + str(count),file=p)
        y += 1
    slIsoList.append(k1)
    orderIsoList.append(iso1)
    x += 1

np.savetxt('Total_Allele_Pairwise_Comp',totalAlleleMatrix,delimiter = '\t',fmt='%s')



bn_allele_Dict = {}

f = open('BNCoreCalls.csv')

bnAlleles = csv.reader(f)

for row in bnAlleles:
    isolate = row[0]
    isolate = isolate.replace('.','_')
    alleles = row[1:]
    bn_allele_Dict[isolate] = alleles

np.save('bn_correct_alleles.npy',allele_Dict)


totalAlleleMatrix = np.zeros(((len(allele_Dict)+1),(len(allele_Dict)+1)),dtype=object)

usedIsoList = []
f = open('BioNum_3_Correct_Comp','a')

print('ISO1'+'\t'+'ISO2'+'\t'+'BN',file=f)

x=1
for k1 in orderIsoList:
    y=1
    totalAlleleMatrix[x][0] = k1
    totalAlleleMatrix[0][x] = k1
    for k2 in orderIsoList:
        list1 = bn_allele_Dict[k1]
        list2 = bn_allele_Dict[k2]
        size = 0
        count = 0
        while size < len(list1):
            if list1[size] != 0 and list2[size] != 0:
                if list1[size] != list2[size]:
                    count += 1
                size += 1
            else:
                size += 1
        totalAlleleMatrix[x][y] = count
        if k2 not in usedIsoList:
            print(k1 + '\t' + k2 + '\t' + str(count),file=f)
        y += 1
    usedIsoList.append(k1)
    x += 1


np.savetxt('BioNum_Allele_Correct_Comp',totalAlleleMatrix,delimiter = '\t',fmt='%s')


from __future__ import print_function
import argparse
import csv

parser = argparse.ArgumentParser(description='Creates a simple pairwise compairison')

parser.add_argument('csvFile',type=str,help='Input tsv file containing isolates and associated alleles')

parser.add_argument('outFileName',type=str,help='Name you want the outfile to be')

args = parser.parse_args()

alleleList = []

with open(args.csvFile) as f:
    reader = csv.reader(f)
    for row in reader:
        alleleList.append(row)

compList = []

primerInfo = alleleList[0]
diffCompList = []

#del alleleList[0]

for a in alleleList:
    for b in alleleList:
        print(a[0] + ' for ' + b[0])
        diffCount = 0
        step = 1
        if a[0] != b[0]:
            if a[0] not in compList and b[0] not in compList:
                while step < len(a):
                    if int(a[step].strip(' ')) != 0 and int(b[step].strip(' ')) != 0:
                        if int(a[step].strip(' ')) != int(b[step].strip(' ')):
                            diffCount += 1
                    step += 1
                diffInfo = [a[0],b[0],diffCount]    
                diffCompList.append(diffInfo)
    compList.append(a[1])

#f = open(args.outFileName,'a')

with open(args.outFileName,'w') as f:
    writer = csv.writer(f)
    writer.writerows(diffCompList)




#for i in diffCompList:
#    print(i[0]+'\t'+i[1]+'\t'+str(i[2]),file=f)

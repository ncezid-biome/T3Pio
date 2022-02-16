import argparse
import pandas as pd
import statistics

parser = argparse.ArgumentParser(description='This is to take in a csv file and give stats about it')

parser.add_argument('csvFile',type=str,help='csv file')

parser.add_argument('columnStart',type=int,help='where to start stat finding')

args = parser.parse_args()

myDf = pd.read_csv(args.csvFile,dtype=None)
myDf1 = myDf.values.tolist()

statList = []
isoStat = {}
outbreaks = []

for row in myDf1:
    start = args.columnStart
    numZero = 0
    iso = row[start-1]
    isoStat[iso]=0
    outbreaks.append(row[start-3])
    while start < len(row):
        allele = int(row[start])
        if allele == 0:
            numZero += 1
        start += 1
    statList.append(numZero)
    isoStat[iso]=numZero

average = sum(statList)/len(statList)

f = open(args.csvFile+'BelowAverageIsolates.txt','a')


sdev = statistics.stdev(statList)
sdev2 = sdev+sdev

for k,v in isoStat.items():
    if v < (average-sdev):
        print(k+'\t'+'NumberZeros: '+v,file=f)

print(average)

numOutBs = set(outbreaks)
print(len(numOutBs))

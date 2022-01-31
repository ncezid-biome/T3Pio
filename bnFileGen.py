import pandas as pd
import argparse
import csv

parser = argparse.ArgumentParser(description='This is to take a bn core MLST sheet and pull out specific isolates')

parser.add_argument('WantedIsolates',type=str,help='Text list of desired isolates')

parser.add_argument('bnCsv',type=str,help='CSV sheet from BN')

args = parser.parse_args()

with open (args.WantedIsolates,'r') as i:
    isoNames = i.readlines()
i.close()

bnDf = pd.read_csv(args.bnCsv,dtype=None)
bnDf1 = bnDf.values.tolist()

trueBNList = []

for bn in bnDf1:
    for iso in isoNames:
        bn[4] = bn[4].replace('-','_')
        if bn[4] in iso:
            bn[4] = iso
            trueBNList.append(bn)
            print(bn[4])

with open('BNMasterFileOutput.csv','w') as c:
    writer = csv.writer(c)
    writer.writerows(trueBNList)

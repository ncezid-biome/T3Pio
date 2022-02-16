def sequenceDict(ref):

    seqDict = SeqIO.to_dict(SeqIO.parse(ref,'fasta'))

    return(seqDict)



def coverageCheck(seqDict,blastFile,cov):

    clearedContigs = []

    with open(blastFile) as b:
        blastInfo = b.readlines()
        b.close()

    for blast in blastInfo:
        
        try:

            blast = blast.split('\t')
        
            if int(blast[3])/len(seqDict[blast[1]]) > cov:

                contigMetrics = [blast[1],blast[3],blast[2],int(blast[3])/len(seqDict[blast[1]])]
                clearedContigs.append(contigMetrics)

        except IndexError:
            continue

    return(clearedContigs)


import argparse
from Bio import SeqIO
import csv

parser = argparse.ArgumentParser(description='Takes blast6 output files and calculates coverage percent across subject sequence')

parser.add_argument('blastFile',type=str,help='Path to blast output file')

parser.add_argument('ref',type=str,help='Path to reference genome')

parser.add_argument('cov',type=float,help='Coverage threshold')

parser.add_argument('outFile',type=str,help='Outfile name')

args = parser.parse_args()


if __name__ == '__main__':

    seqDict = sequenceDict(args.ref)

    clearedContigs = coverageCheck(seqDict,args.blastFile,args.cov)

    with open(args.outFile,'w',newline="") as f:

        writer = csv.writer(f)

        writer.writerows(clearedContigs)
    

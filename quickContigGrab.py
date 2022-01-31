from Bio import SeqIO
import glob
import argparse

parser = argparse.ArgumentParser(description='Takes cleared contigs and makes a multifasta file')

parser.add_argument('clearedFileDir',type=str,help='Input the directory containing the cleared contig files')

parser.add_argument('ref',type=str,help='Input path to the reference file')

parser.add_argument('outFile',type=str,help='Name of the output multifasta file')

args = parser.parse_args()

if __name__ == '__main__':

    clearedFiles = sorted (glob.glob(args.clearedFileDir+'*'))

    contigs = []

    for cFile in clearedFiles:

        with open(cFile,'r') as f:
            contigInfo = f.readlines()

        f.close()

        for info in contigInfo:
            
            contigs.append(info.split(',')[0])

    contigSet = set(contigs)

    contigList = list(contigSet)

    seqDict = SeqIO.to_dict(SeqIO.parse(args.ref,'fasta'))

    f = open(args.outFile,'w')
    
    for k in seqDict.keys():
        
        if k not in contigList:

            seq = str(seqDict[k].seq)
    
            print('>'+k,file=f)
            print(seq,file=f)

from Bio import SeqIO
import argparse 


parser = argparse.ArgumentParser(description='Make small changes to the fasta header')

parser.add_argument('genome',type=str,help='Enter the fasta genome')

parser.add_argument('name',type=str,help='Enter the new name')

args = parser.parse_args()


newFasta = open(args.name+'.fasta','w')

for seq in SeqIO.parse(args.genome,'fasta'):

    seq.id = seq.id + '-' + args.name

    SeqIO.write(seq,newFasta,'fasta')




def genomeSizeSelect(genome,size,outfile):

    f = open(outfile,'w')

    for seq in SeqIO.parse(genome,'fasta'):
        if len(seq) > size:
            SeqIO.write(seq,f,'fasta')


import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='This size selects out contigs from a genome')

parser.add_argument('genome',type=str,help='Path to genome fasta file')

parser.add_argument('size',type=int,help='Size threshold')

parser.add_argument('outFile',type=str,help='Out file name')

args = parser.parse_args()

if __name__ == '__main__':

    genomeSizeSelect(args.genome,args.size,args.outFile)

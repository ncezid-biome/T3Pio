from Bio import SeqIO
import glob
import argparse
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser(description='This checks emboss or really any single fasta file for variability in the sequence and returns min/max lengths')

parser.add_argument('fastaDir',type=str,help='Path to directiory containing files')

args = parser.parse_args()

if __name__ == '__main__':

    lengthList = []

    fastas = sorted(glob.glob(args.fastaDir+'*'))

    for fasta in fastas:
        for record in SeqIO.parse(fasta,'fasta'):
            for letter in record.seq:
                if letter not in ('ATGC'):
                    lengthList.append(len(record))
                    break

    print('The number of files with var is: ' + str(len(lengthList)))
    print('The min length is: ' + str(min(lengthList)))
    print('The max length is: ' + str(max(lengthList)))

    plt.hist(lengthList,bins=100)
    plt.show()

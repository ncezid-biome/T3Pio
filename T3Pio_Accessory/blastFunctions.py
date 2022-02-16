from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.Seq import Seq
import os


def isSeqRecord(seq):

    ans = type(seq) == Bio.SeqRecord.SeqRecord

    return ans

def isFile(seq):

    ans = os.path.isfile(seq)

    return ans

def makeSeqRecord(seq):

    if isFile(seq):
        
        try:

        fasta_record = SeqIO.read(seq,'fasta')

        except ValueError:

            continue
    
    else:
        
        try:
        
            fasta_record = Seq(seq)

        except TypeError:

            continue

    return fasta_record


def makeSeqTranslation(fasta_record):

    try:
        
        fasta_record = fasta_record.translate()

    except ValueError:

        continue



def callNCBIWWWqblast(seq,addArgs=None):

    

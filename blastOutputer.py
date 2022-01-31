class BlastOut:

    def __init__(self,subject,identity,length,start,end,evalue):

        self.subject = subject
        self.identity = identity
        self.length = length
        self.start = start
        self.end = end
        self.evalue = evalue



def blastnRunner(blastDB,query):

    r1 = subprocess.Popen(('blastn','-db',blastDB,'-query',query,'-outfmt','7'),stdout=subprocess.PIPE)

    r2 = (r1.stdout.read().decode('ascii'))

    
    blastClassList = blastOutputParser(r2)

    return(blastClassList)


def blastOutputParser(blastOut):

    blastClassList = []

    blastResults = blastOut.split('\n')

    for result in blastResults:

        try:
        
            result = result.split('\t')

            r = BlastOut(result[1],result[2],result[3],result[8],result[9],result[10])

            blastClassList.append(r)

        except IndexError:

    return(blastClassList)


def blastClassListReducer(blastClassList,threshold):

    reducedBlastClassList = []

    for blastClass in blastClassList:

        if blastClass.evalue <= threshold:

            reducedBlastClassList.append(blastClass)

    return(reducedBlastClassList)


import argparse
import subprocess

parser = argparse.ArgumentParser(description='Script designed to run command line blast on a local database and return hits below a user input e-value')

parser.add_argument('blastDB',type=str,help='Enter the blast database name to use')

parser.add_argument('query',type=str,help='Enter the fasta file to search the DB with')

parser.add_argument('threshold',type=float,help='Enter the e-value limit')

parser.add_argument('outfile',type=str,help='Enter the name of the outfile')

args = parser.parse_args()

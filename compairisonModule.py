
class PairWise:

    """Calculated pairwise difference counts between two isolate multifasta files

    Attributes:
        isolate1: String representing the name of the first isolate in the comparison
        isolate2: String representing the name of the second isolate in the compairison
        diffSites: Integer representing the number of different sequences between the two isolates
        snpDiff: Integer representing the number of difference SNPs between the two isolates
        orthoDiff: Interger representing the number of different orthogroup sites. Calculated 
                   differently than diffSites as multiple sequences could come from the same 
                   orthogroup.
    """

    def __init__(self,isolate1,isolate2,diffSites,orthoDiff):

        """Returns a new Pairwise object"""
        
        self.isolate1 = isolate1
        self.isolate2 = isolate2
        self.diffSites = diffSites
        self.orthoDiff = orthoDiff

class DictionaryRecord:

    """SeqIO dictionary and name of the isolate the dictionary was created from
    
    Atrributes:
        name: String representing the name of the isolate 
        dictionary: SeqIO dictionary created from a multifasta file
    """

    def __init__(self,name,dictionary):
        
        """Returns a new DictionaryRecord object"""

        self.name = name
        self.dictionary = dictionary


#Shelved for later consideration
'''
class Amplicons:

    def __init__(self):

    ampliconDict = {}
        
    def AmpliconAccumulator(self,key,seq,se2):

        try:

            ampliconDict[key].append(seq)
            
            ampliconDict[key].append(seq2)

        except KeyError:

            ampliconDict.setdefault(key,[])

            ampliconDict[key].append(seq)
            
            ampliconDict[key].append(seq2)
'''


def FastaToDicts(fastaFiles):

    """Takes multifasta files and turns them into SeqIO dictionary records

    Params
    ------
    fastaFiles: List of Strings
         File locations of the multifasta files


    Returns
    -------
    dictClassList: List 
         List of DictionaryRecord objects 
    """

    dictClassList = []

    for fasta in fastaFiles:

        fastaName = fasta.split('/')[-1].split('.')[0]

        fastaDict = SeqIO.to_dict(SeqIO.parse(fasta,'fasta'))

        dictRecord = DictionaryRecord(fastaName,fastaDict)
    
        dictClassList.append(dictRecord)

    return(dictClassList)



#Retired
#Plays no role in code now. Kept for posterity
'''
def SnpDifferenceCount(seq1,seq2):

    """Takes two strings and counts the differences between them

    Params
    ------
    seq1: String
         DNA sequence
    seq2: String
         DNA sequence

    Returns
    -------
    count: Integer
         Total number of differences between two strings
    """

    count = sum(1 for a,b in zip(seq1,seq2) if a != b) + abs(len(seq1)-len(seq2))

    return(count)
'''


def OrthogroupDiffCount(orths,fastaDict,fastaDict2):

    """Takes a list of partial keys and determines if there are differences between two dictionaries

    Params
    ------
    orths: List
         List of strings that are partial keys for fastaDict
    fastaDict: SeqIO dictionary
         SeqIO dictionary of a multifasta file
    fastaDict2: SeqIO dictionary
         SeqIO dictionary of a multifasta file

    
    Searches both fastaDict and fastaDict2 using partial key matching to gather a list of values 
    consisting of DNA sequences. Then compares the two lists to determine if they are the same or
    not.

    Returns
    -------
    orthDiffs: Integer
          Total number of differences found between compared lists
    """

    orthDiffs = 0

    for orth in orths:

        boolList = []

        dictList = [str(val.seq) for key,val in fastaDict.items() if orth in key]

        dictList2 = [str(val.seq) for key,val in fastaDict2.items() if orth in key]

        if dictList2:
            
            if len(dictList) != len(dictList2):
            
                orthDiffs += 1

            else:

                count = 0

                while count < len(dictList):

                    reverseCheck = ReverseComplementCheck(dictList[count],dictList2[count])

                    boolList.append(reverseCheck)

                    count += 1

                if all(boolList) == False:

                    orthDiffs += 1

    return(orthDiffs)


def ReverseComplementCheck(seq1,seq2):

    """Takes two sequences and compares them for similarity. If different, checks to see if seq2 is a reverse_complement of seq1

    Params
    ------
    seq1: String
         Nucleotide sequence
    seq2: String 
         Nucleotide sequence



    """

    compare = False

    if seq1 == seq2:

        compare = True
    
    else:

       seq2 = Seq(seq2)

       seq2 = str(seq2.reverse_complement())

       if seq1 == seq2:
           
           compare = True


    return(compare)



def FastaCompare(fastaDict,fastaDict2,fastaName,fastaName2):

    """Takes dictionaries and compares values for differences

    Params
    ------
    fastaDict: Dictionary
         SeqIO dictionary of a multifasta file
    fastaDict2: Dictionary
         SeqIO dictionary of a multifasta file
    fastaName: String
         Name of isolate from which fastaDict was constructed
    fastaName2: String 
         Name of isolate from which fastaDict2 was constructed

    
    Iterates through fastaDict using it's keys to get the associated values from fastaDict2.
    Convertes both fastaDict and fastaDict2 values to strings from SeqIO dictionary values objects.
    Sends stings to be compared by SnpDifferenceCount function and stores partial key values from 
    fastaDict in orths List. Converts orthos List to a Set and back to List to ensure unique values
    and sends orths List, fastaDict, and fastaDict2 to OrthogroupDiffCount.

    Returns
    -------
    pairwise: PairWise object
         Creates PairWise object using integer values generated during comparison
    """

    diffSites = 0

    orths = []

    for k,v in fastaDict.items():

        orth = k.split('-')[0]

        orths.append(orth)

        try:
                
            seq1=str(v.seq)

            seq2 = str(fastaDict2[k].seq)
                
            siteCheck = ReverseComplementCheck(seq1,seq2)

            if siteCheck == False:

                diffSites += 1

        except KeyError:
            continue


    orths = list(set(orths))

    orthoDiff = OrthogroupDiffCount(orths,fastaDict,fastaDict2)

    pairwise = PairWise(str(fastaName),str(fastaName2),int(diffSites),int(orthoDiff))

    return(pairwise)



def FastaDictCheck(isolatePair,dictionaryRecords,q):

    """Checks that the SeqIO dictionaries created from fasta files are the correct size

    Params
    ------
    isolatePair: Tuple
         File paths to two fasta files
    dictionaryRecords: Tuple
         Two SeqIO dictionary records corrosponding to the two fasta files in isolatePair
    q: queue
         queue for logging events

    Checks the size of the fasta files and compairs them to the size of the SeqIO dictionaries to make
    sure the dictionaries were created correctly.

    Returns
    -------
    ok: boolean
         True/False value indicating whether the dictionaries were the correct size

    """

    ok = True

    r1 = subprocess.Popen(('wc','-l',isolatePair[0]),stdout=subprocess.PIPE)
    r1Len = int((r1.stdout.read().decode('ascii')).split(' ')[0])

    r2 = subprocess.Popen(('wc','-l',isolatePair[1]),stdout=subprocess.PIPE)
    r2Len = int((r2.stdout.read().decode('ascii')).split(' ')[0])

    if r1Len != int(len(dictionaryRecords[0].dictionary)*2):
        
        q.put('ERROR '+isolatePair[0]+' may not be a fasta file. \n It did not form an appropriately sized SeqIO dictionary object')

        ok = False

    if r2Len != int(len(dictionaryRecords[1].dictionary)*2):
      
        q.put('ERROR '+isolatePair[1]+' may not be a fasta file. \n It did not form an appropriately sized SeqIO dictionary object')
        
        ok = False

    return(ok)



def listener_process(q):

    """Sets the logging file for the listener process

    Params
    ------
    q: queue 
         queue for loggeing events

    Returns
    -------
    None

    """

    with open(logFile,'w') as f:

        while 1:
            m = q.get()

            if m == 'kill':
                
                break

            f.write(str(m)+'\n')
            f.flush()




def ParallelFunctions(isolatePair,q):

    """Sets up functions to be run during multiprocessing

    Params
    ------
    isolatePair: Tuple
         A set of two file paths to multifasta files

    Returns
    -------
    pairwiseInfo: PairWise object
         PairWise object generated from comparison
    """
    
    q.put(isolatePair[0]+'\t'+isolatePair[1])

    dictionaryRecords = FastaToDicts(isolatePair)

    fastaCheck = FastaDictCheck(isolatePair,dictionaryRecords,q)

    if fastaCheck == False:
        
        q.put('ERROR: An error has occured with either ' +isolate[0]+' or '+isolate[1]+' please check the log')
        

    pairwiseInfo = FastaCompare(dictionaryRecords[0].dictionary,dictionaryRecords[1].dictionary,dictionaryRecords[0].name,dictionaryRecords[1].name)

    return(pairwiseInfo)



def ParallelRunner(isolatePairs):

    """Sets the multiprocessing pool and process

    Params
    ------
    isolatePairs: List
         List of pairs of file paths to multifasta files


    Creates a manager queue for a listener process that will be used for logging across all spawned processes. Only the listener process will have write acess to the log file. Queue ensures only one write at a time. 

    Returns
    -------
    result_list: List
         List of PairWise objects
    """

    manager = multiprocessing.Manager()

    q = manager.Queue()
        
    pool = multiprocessing.Pool(processes=args.NumberCores)

    watcher = pool.apply_async(listener_process,(q,))

    parallelStatic = partial(ParallelFunctions,q=q)

    result_list = pool.map(parallelStatic,isolatePairs)

    q.put('kill')

    pool.close()

    pool.join()

    return(result_list)


#Shelved for now
'''
def OutfileWriter(results,outfileName,dataType):

    f = open(outfileName+dataType+'.tsv','w')

    for pairwise in results:
        print(pairwise.isolate1+'\t'+pairwise.isolate2+'\t'+str(pairwise.dataType),file=f)

    f.close()
'''

from Bio import SeqIO
from Bio.Seq import Seq
import glob
import argparse
from itertools import combinations
import multiprocessing
import subprocess
import logging
import logging.handlers
from functools import partial

parser = argparse.ArgumentParser(description='Comparse multifasta files and returns pairwise distance files of both number of different sites and total number of snp differences')

parser.add_argument('multiFastaFolder',type=str,help='Enter the path to a folder containing multi fasta files')

parser.add_argument('outfileName',type=str,help='Enter the name of the outfiles')

parser.add_argument('NumberCores',type=int,help='Enter number of cores to run program with')

args = parser.parse_args()


fastaFiles = sorted(glob.glob(args.multiFastaFolder+'*.fasta'))


isolatePairings = list(combinations(fastaFiles,2))

logFile = args.outfileName+'.log'

if __name__ == '__main__':

    results = ParallelRunner(isolatePairings)



site = open(args.outfileName+'sites.tsv','w')
ortho = open(args.outfileName+'ortho.tsv','w')

for pairwise in results:

    print(pairwise.isolate1+'\t'+pairwise.isolate2+'\t'+str(pairwise.diffSites),file=site)
    print(pairwise.isolate1+'\t'+pairwise.isolate2+'\t'+str(pairwise.orthoDiff),file=ortho)

site.close()
ortho.close()

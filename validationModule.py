class Primers:

    """Relevent information regarding primers used in the analysis

    Attributes:
        orthogroupInfo: String representing which orthogroup the primer belongs to
        number: String representing which primer it is for a respective orthogroup
        leftSeq: String representing the forward primer sequence 5'-3'
        rightSeq: String representing the reverse primer sequence 5'-3'
        leftHit: Integer representing where the first base of the forward primer lands in a genome
        rightHit: Integer representing where the last base of the reverse primer lands in a genome
        leftLen: Integer representing the length of the forward primer
        rightLen: Integer representing the length of the reverse primer

    """
    
    def __init__(self,orthogroupInfo,number,leftSeq,rightSeq,leftHit,rightHit,leftLen,rightLen):
        
        """Returns a new Primers object"""

        self.orthogroupInfo = orthogroupInfo
        self.number = number
        self.leftSeq = leftSeq
        self.rightSeq = rightSeq
        self.leftHit = leftHit
        self.rightHit = rightHit
        self.leftLen = leftLen
        self.rightLen = rightLen

class PrimerSearchResults:

    """Relevent information from in silico PCR program Emboss/primersearch

    Atributes:
       primerInfo: String representing the primer used in the in silico analysis
       sequenceName: String representing the contig the primers hit
       ampLen: Integer representing the length of the predicted amplicon
       leftHit: Integer representing where the first base of the forward primer lands
       rightHit: Integer representing where the last base of the reverse primer lands
       sequence: String representing the predicted amplicon

    """
    
    def __init__(self,primerInfo,sequenceName,ampLen,leftHit,rightHit,sequence):

        """Returns a new PrimerSearchResults object"""

        self.primerInfo = primerInfo
        self.sequenceName = sequenceName
        self.ampLen = ampLen
        self.leftHit = leftHit
        self.rightHit = rightHit
        self.sequence = sequence
  

#########################################################
#('Def') OutputLoggingFileSetter() comments
#Takes as input ('String') outfile path 
##### ('String') primerFile name
#Sets Logger file name 
#Sets input logging level to Debug
#Sets logger file format to include Time,logging level, message
#Set logger file handler
#Returns nothing
#####Just sets the logging file
def OutputLoggingFileSetter(outfile,primerFile):

    logger = logging.getLogger(outfile+primerFile)

    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter('%(asctime)s:%(levelname)s:%(message)s')

    file_handler = logging.FileHandler(outfile+primerFile+'.log','w')

    file_handler.setFormatter(formatter)

    logger.addHandler(file_handler)    



#########################################################
#('Def') PrimmersearchRunner() comments
#Takes as input ('String') trimalFile path
##### ('List') of ('Primers') ('Objects')
#Creates ('IO') primer file for primersearch 
##### Pulls primer number/forward/reverse sequences from ('Primers') ('Object')
#####Stores information in ('List') primerInfoList
#Runs primersearch using subprocess
#Returns ('String') primersearchFile.ps file path
def PrimersearchRunner(assemblyFile, primerFile,outfile,logger):

    primersearchFile = primerFile.split('/')[-1].split('.')[0]

    primersearchFile = outfile+primersearchFile
    
    try:

        return_code = subprocess.check_output(['primersearch','-seqall',assemblyFile,'-infile',primerFile,'-mismatchpercent','6','-outfile',primersearchFile+'.ps'],stderr=subprocess.STDOUT)

    except subprocess.CalledProcessError as error:
        
        errorMsg = ('Status : Fail',error.returncode,error.output.strip('\n'))
        
        logger.error(str(errorMsg))

    return(primersearchFile+'.ps')



#########################################################
#('Def') PrimersearchValidator() comments:
#Takes as input ('List') ampliconInfo
##### ('Int') numberIsolates
#Uses list comprehension to pull sequence names from ('List') ampliconInfo
#####Creates ('Set') sequencesSet of ('List') sequences
#####Checks for equal length to ensure no duplicate isolates pulled
#Checks correct number sequences present against ('Int') numberIsolates
#Returns ('True') if proper number of isolates and no duplicates present
#Returns ('False') if improper number of isolates and/or duplicates present
def PrimersearchValidator(ampliconInfo,numberIsolates,logger):

    sequenceStrings = [ seq for seq in ampliconInfo if 'Sequence' in seq]
    
    sequences = []

    for seqs in sequenceStrings:
        
        seqs = seqs.strip('\t')
        seqs = seqs.strip('\n')
        seqs = seqs.split(':')[1]
        seqs = seqs.strip(' ')
        sequences.append(seqs.split('.')[0])

    sequencesSet = set(sequences)

    if len(sequences) != len(sequencesSet):
        
        return(False)

    if len(sequences) != numberIsolates:

        return(False)
    else:
        
        return(True)


#########################################################
#('Def') PrimersearchComber() comments:
#Takes as input ('List') ampliconInfo
#####('Primer') ('Object') primer
#####('Dictionary') sequenceRecordDict
#Parses ('List') ampliconInfo to find forwardHit/reverseHit/ampliconLen/sequenceName
#####for each isolate record present in ('List') ampliconInfo
#Pulls isolate sequence from ('Dictionary') sequenceRecordDict and slices ('String') for amplicon
#####sequence 
#Stores information in ('Object') ('PrimerSearchResults') for each isolate
#Returns ('List') of ('PrimerSearchResults') ('Objects')
def PrimersearchComber(ampliconInfo,primer,sequenceRecordDict,logger):
    
    primersearchObjects = []

    for line in ampliconInfo:
        
        if line.startswith('Amplimer'):

            try:
            
                forwardHit = int(ampliconInfo[ampliconInfo.index(line)+3].split(' ')[5])
                reverseHit = int((ampliconInfo[ampliconInfo.index(line)+4].split(' ')[5]).replace('[','').replace(']',''))
                ampliconLen = int(ampliconInfo[ampliconInfo.index(line)+5].split(' ')[2])
                sequenceName = str(ampliconInfo[ampliconInfo.index(line)+1].split(' ')[1].strip(' ').strip('\n'))
                
                try:

                    sequence = str(sequenceRecordDict[sequenceName].seq)

                except KeyError:
                    
                    logger.error('The contig headers for this file differ from the headers primersearch captured. Parsing the sequenceRecordDict for a partial match to the key using the sequenceName.')

                    seqName = [key for key in sequenceRecordDict.keys() if sequenceName in key]

                    if len(seqName) == 1:

                        logger.error('Single partial match was found. Proceeding as intended')

                        sequence = str(sequenceRecordDict[seqName[0]].seq)

                    if len(seqName) == 0:

                        logger.error('Contig header issue could not be corrected. No matches were found corresponding to primersearch file')

                        logger.error(sequenceName+' Could not be found in sequenceRecordDict')
                        
                        break

                    if len(seqName) > 1:

                        logger.error('Contig header issue could not be corrected. Too many matches to found.')

                        logger.error(sequenceName+' Is not unique in sequenceRecordDict')
                    
                        break


                sequence = (sequence.replace('-',''))

                sequenceLen = len(sequence)
                
                sequence = sequence[forwardHit+primer.leftLen-1:sequenceLen-reverseHit-primer.rightLen+1]
                
                primersearchObject = PrimerSearchResults(primer,sequenceName,ampliconLen,forwardHit,reverseHit,sequence)
                
                primersearchObjects.append(primersearchObject)

                logger.debug(str(forwardHit)+'\s'+str(primer.leftLen)+'\s'+':'+'\s'+str(sequenceLen)+'\s'+str(reverseHit)+'\s'+str(primer.rightLen))

                logger.debug(str(primersearchObject.primerInfo)+'\t'+str(primersearchObject.sequenceName))
                logger.debug(str(primersearchObject.primerInfo.orthogroupInfo)+'\t'+str(primersearchObject.primerInfo.number))

                logger.debug(primer.leftSeq+'\t'+primer.rightSeq)

                logger.debug(str(primersearchObject.leftHit)+'\t'+str(primersearchObject.rightHit))

                logger.debug(sequence+'\n'+primersearchObject.sequence)

            except IndexError:
                
                primersearchObject = PrimerSearchResults(primer,None,None,None,None,None)

                primersearchObjects.append(primersearchObject)

                logger.debug(primersearchObject.primerInfo.number)
                logger.debug(ampliconInfo)
                logger.debug(primer.leftSeq+'\t'+primer.rightSeq)
                logger.warning(str(primersearchObject.primerInfo.number)+' has no amplcions for this assembly')

    return(primersearchObjects)


#########################################################
#('Def') PrimersearchParser() comments:
#Takes as input ('String') primersearchFile file path
#####('Int') numberIsolates
#####('List') of ('Primer') ('Objects')
#####('String') trimalFile file path
#('IO') reads in primersearchFile to ('List') primersearchInfo
#Stores trimalFile into ('SeqIO') ('Dictionary') ('Object') sequenceRecordDict
#Indexes ('List') primersearchInfo on 'Primer name' + ('Primer') ('Object') primer.numer
#Slices ('List') primersearchInfo using index and ('Int') numerIsolates * 6 into ('List') ampliconInfo
#Sends ('List') ampliconInfo and ('Int') numberIsolates to ('Def') PrimersearchValidator()
#Sends ('List') ampliconInfo and ('Primer') ('Object') to ('Def') PrimersearchComber
#Returns ('List') of ('PrimerSearchResults') ('Objects')
def PrimersearchParser(primersearchFile,numberIsolates,primerObjectList,sequenceRecordDict,logger):

    try:

        with open(primersearchFile,'r') as f:
            primersearchInfo = f.readlines()
        f.close()

        primersearchObjectLists = []
    
        primersearchObjectList = []


        for primer in primerObjectList:
        
            ampliconInfoStartIndex = primersearchInfo.index('Primer name '+str(primer.number)+'\n')

            ampliconInfo = primersearchInfo[ampliconInfoStartIndex:ampliconInfoStartIndex+(numberIsolates*6)+1]

            validationMark = PrimersearchValidator(ampliconInfo,numberIsolates,logger)

            if validationMark == True:
            
                parsedInfo = PrimersearchComber(ampliconInfo,primer,sequenceRecordDict,logger)
        
                primersearchObjectLists.append(parsedInfo)

            else:

                break

        for primersearchObjects in primersearchObjectLists:

            primersearchObjectList = primersearchObjectList + primersearchObjects

    except FileNotFoundError:

        logger.error('No primersearch file created')


    return(primersearchObjectList)

#########################################################
#('Def') ValidationPrimerIntake() comments
#Takes as input ('String') primerFile primer file path
#####Logger logging file
#Creates ('Class') Primers object 
#Returns primerObjects
def ValidationPrimerIntake(primerFile,logger):

    primerObjects = []
        
    try:

        with open(primerFile,'r') as f:
            primers = f.readlines()
            f.close()

    except FileNotFoundError:

        logger.error('No primer file found')

    try:

        for primer in primers:
            
            primerInfo = primer.split('\t')
            
            primerObj = Primers(primerInfo[0][:9],primerInfo[0],primerInfo[1],primerInfo[2].strip('\n'),None,None,len(primerInfo[1]),len(primerInfo[2].strip('\n')))
            
            primerObjects.append(primerObj)

    except IndexError:

        logger.error('Primer file not formed as expected')

    return(primerObjects)



#########################################################
#('Def') AssemblyDictionary() comments
#Takes as input ('String') assemblyFile file path
#####Logger logging file
#Turns assembly file into ('Dict') sequenceRecordDict SeqIO dictionary object
#Returns ('Dict') sequenceRecordDict
def AssemblyDictionary(assemblyFile,logger):

       sequenceRecordDict =SeqIO.to_dict(SeqIO.parse(assemblyFile,'fasta'))

       logger.debug(assemblyFile)

       return(sequenceRecordDict)


#########################################################


#########################################################
#('Def') ParallelFunctions() comments
#Takes as input('String') primerFile file path
#####('String') assemblyFile file path
#####('String') outfile file path
#Sets functions to run in parallel
#Returns ('List') of ('Objects') PrimerSearchResults
def ParallelFunctions(primerFile,assemblyFile,outfile):

    primerFileName = primerFile.split('/')[-1]
    
    OutputLoggingFileSetter(outfile,primerFileName)

    logger = logging.getLogger(outfile+primerFileName)

    logger.debug(str(primerFile + ' is in the parallel universe'))

    primerObjects = ValidationPrimerIntake(primerFile,logger)

    if len(primerObjects) == 0:
        
        logger.error('No primers were found in primer file')

    sequenceRecordDict = AssemblyDictionary(assemblyFile,logger)

#Put check here to make sure the sequence record exists

    primersearchFile = PrimersearchRunner(assemblyFile,primerFile,outfile,logger)

    if path.exists(primersearchFile) == False:

        logger.error('Primersearch file does not exist')

    primersearchObjectList = PrimersearchParser(primersearchFile,1,primerObjects,sequenceRecordDict,logger)

    subprocess.check_output(['rm',primersearchFile])

    return(primersearchObjectList)


#########################################################
def ParallelRunner(primers):

    pool = multiprocessing.Pool(processes=args.NumberCores)

    parallelStatic = partial(ParallelFunctions,assemblyFile=assembly,outfile=outfile)

    result_list = pool.map(parallelStatic,primers)

    pool.close()

    return(result_list)

########################################################
def MultifastaCreator(primersearchObjectList,fastaFile):

    for primersearchObject in primersearchObjectList:

        if primersearchObject.sequenceName != None:

            print('>'+primersearchObject.primerInfo.orthogroupInfo +'-'+str(primersearchObject.primerInfo.number),file=fastaFile)
            print(primersearchObject.sequence,file=fastaFile)

    return(None)
    
#########################################################
#Main Body


import argparse
import glob
import multiprocessing 
from functools import partial
import subprocess
import os
import sys
from Bio import SeqIO
import logging
import os.path
from os import path



#DEBUG: Detailed information, typically of interest only when diagnosing problems.

#INFO: Confirmation that things are working as expected.

#WARNING: And indication that something unexpected happened, or indicative of some problem in the near
#future (e.g. 'disk space low'). The software is still working as expected

#ERROR: Due to more serious problems, the software has not been able to perform some function

#CRITICAL: A serious error, indicating that the program itself may be unable to continue running



parser = argparse.ArgumentParser(description='Validation script for T3Pio. \n Generates amplcions from primers generated through T3Pio on assemblies not used in the T3Pio script')

parser.add_argument('assemblyDir',type=str,help='Input the directory containing assembly files')

parser.add_argument('fileType',type=str,help='Input the file extention to be used. eg. fna,fasta,')

parser.add_argument('primerDir',type=str,help='Input the directory containing primer files')

parser.add_argument('NumberCores',type=int,help='Input number of cores to use')

parser.add_argument('outDir',type=str,help='Input the out directory')

args = parser.parse_args()

assemblies = sorted(glob.glob(args.assemblyDir+'*.'+args.fileType))

primers = sorted(glob.glob(args.primerDir+'*'))


#########################################################
#Execution script
#Takes ('List') of assemblies and executes the above functions
#Runs on assembly at a time
#Returns nothing
#Outputs multifasta files

for assembly in assemblies:

    assemblyName = assembly.split('/')[-1].split('.')[0]

    if not os.path.exists(args.outDir+'/'+assemblyName+'/'):

        os.makedirs(args.outDir+'/'+assemblyName+'/',exist_ok=True)

    outfile = args.outDir+'/'+assemblyName+'/'

    if __name__== '__main__':

        results = ParallelRunner(primers)

    fastaFile = open(outfile+assemblyName+'.fasta','w')

    for primersearchObjectList in results:

        MultifastaCreator(primersearchObjectList,fastaFile)

    fastaFile.close()

    r1 = subprocess.Popen(('grep','-r','ERROR',outfile),stdout=subprocess.PIPE)
    
    r2 =(r1.stdout.read().decode('ascii'))
    
    r2 = r2.split('\n')

    logFiles = glob.glob(outfile+'*.log')
    
    for logs in logFiles:

        subprocess.check_output(['rm',logs])

    errorLog = open(outfile+assemblyName+'.log','w')

    for errors in r2:

        print(errors,file=errorLog)

    errorLog.close()

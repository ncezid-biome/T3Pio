


class Orthogroup:
    """"Parsed orthogroup info

    Attributes:
        orthogroup: String representing the orthogroup
        sequence: String representing the sequence
    """
    def __init__(self,orthogroup,sequence):

        #Returns an orthogroup object
        self.orthogroup = orthogroup
        self.sequence = sequence


class Primers:

    """
    Attributes:
        orthogroupInfo: String representing which orthogroup
        number: Int representing primer number
        leftSeq: String representing thge left sequence
        rightSeq: String representing the right sequence
    """
    
    def __init__(self,orthogroupInfo,number,leftSeq,rightSeq,leftHit,rightHit,leftLen,rightLen):
        self.orthogroupInfo = orthogroupInfo
        self.number = number
        self.leftSeq = leftSeq
        self.rightSeq = rightSeq
        self.leftHit = leftHit
        self.rightHit = rightHit
        self.leftLen = leftLen
        self.rightLen = rightLen



class OrthofinderResults:

    def __init__(self,orthoNumber,seqNames):
        self.orthoNumber = orthoNumber
        self.seqNames = seqNames

class PrimerSearchResults:
    
    def __init__(self,primerInfo,sequenceName,ampLen,leftHit,rightHit,sequence):
        self.primerInfo = primerInfo
        self.sequenceName = sequenceName
        self.ampLen = ampLen
        self.leftHit = leftHit
        self.rightHit = rightHit
        self.sequence = sequence

#########################################################
#('Def') orthoFinderParser() comments
#Takes as input the Orthogroup.txt file that is generated from OrthoFinder
#Returns orthDict dictionary 
#Keys are orthogroup names and values are lists of insolates in the orthogroup
#Very simple file parsing
#Removes trailing '\n' from each line and splits on spaces ' '
#Removes ':' from orthogroup name

def orthoFinderParser(OrthofinderTxt):

    with open(OrthofinderTxt) as f:
        orthogroupInfo = f.readlines()
    f.close()
    
    orthogroupResultsList = []

    for line in orthogroupInfo:
        line = line.strip('\n').split(' ')
        orthogroup = OrthofinderResults(line[0].strip(':'),line[1:])
        orthogroupResultsList.append(orthogroup)
        
    return(orthogroupResultsList)

#########################################################
def OrthofinderClassCheck(orthogroupResultsList):

    for orthogroupObject in orthogroupResultsList:

        if not isinstance(orthogroupObject.orthoNumber,str):

            sys.exit('One of your orthogroup objects is malformed')

        if not isinstance(orthogroupObject.seqNames,list):
            
            sys.exit('One of your orthogroup objects is malformed')

#########################################################
#('Def') dictionaryCleaner() comments
#Takes as input a ('Dictionary')
#Tracks through the input ('Dictionary') looking for keys with empty values
#Deletes keys that have empty values
#Returns the ('Dictionary') 

def dictionaryCleaner(dictionary):
    
    for k in list(dictionary.keys()):
        if dictionary[k] == []:
            del dictionary[k]
    
    return(dictionary)


#########################################################
#('Def') primer3Parser() comments
#Takes as input a ('Dictionary') containing the STDOUT from primer3
#Initiates the ('List') primerPairObjectList which will hold the ('Class Objects') Primers
#Initiates the ('List') primerInfoList to hold the ('Strings') from the input ('Dictionary')
#Uses indexing on the ('List') primerInfoList to pull out ('Int') totalPrimersReturned
#('Int') totalPrimersReturned is used for the while loop to iterate through ('List') primerInfoList to pull out:
#
#('String') leftSequence
#('String') rightSequence
#('Int') leftHit
#('Int') rightHit
#('Int') leftLength
#('Int') rightLength
#Pulled information is stored in ('Class Object') Primers 
#Returns ('List') of ('Class Objects') Primers

def primer3Parser(primer3Dict):
    
    primerPairObjectList = []

    for k in primer3Dict.keys():
        
        primerInfoList = primer3Dict[k]
        try:

            totalPrimersReturned = int(primerInfoList[primerInfoList.index('PRIMER_PAIR_NUM_RETURNED') + 1])

            if totalPrimersReturned != 0:
                
                orthogroupObj = Orthogroup(primerInfoList[0][:9],primerInfoList[2])

                count = 0
                while count < totalPrimersReturned:
                    leftSequence = (primerInfoList[primerInfoList.index('PRIMER_LEFT_'+str(count)+'_SEQUENCE') + 1])
                    rightSequence = (primerInfoList[primerInfoList.index('PRIMER_RIGHT_'+str(count)+'_SEQUENCE') + 1])
                    leftHit = int((primerInfoList[primerInfoList.index('PRIMER_LEFT_'+str(count)) + 1]).split(',')[0])
                    rightHit = int((primerInfoList[primerInfoList.index('PRIMER_RIGHT_'+str(count)) + 1]).split(',')[0])
                    leftLength = int((primerInfoList[primerInfoList.index('PRIMER_LEFT_'+str(count)) + 1]).split(',')[1])
                    rightLength = int((primerInfoList[primerInfoList.index('PRIMER_RIGHT_'+str(count)) + 1]).split(',')[1])
                    
                    primerObj = Primers(orthogroupObj,count,leftSequence,rightSequence,leftHit,rightHit,leftLength,rightLength)
                    primerPairObjectList.append(primerObj)
                    count += 1

        except ValueError as e:
            print (f" ******************* caught error:  {str(e)}")
            continue


    return(primerPairObjectList)




#########################################################################################
##################### GBK FILE PROCESSING BLOCK #########################################
#########################################################################################
#Code block designed to extract nucleotide and protein sequences from GBK files
#Uses bioPython and the SeqIO.parse function to accomplish this from the CDSs in the GBK files
#Accepts as input a list of absolute file paths to GBK files
#the dictionary nuceotideDictionary is returned and contains the sampleName, and prokka locus information as a key and the nucleotide sequence as the value for all the CDSs


def gbkParser(gbkFiles,outFileName):

    nucleotideDictionary = {}    

    if not os.path.exists(outFileName+'/faaFiles/'):            
        os.makedirs(outFileName+'/faaFiles/',exist_ok=True)

    totalFeatures = 0

    for i in gbkFiles:            
        sampleName = i.split('/')     
        sampleName = sampleName[len(sampleName)-1]    
        sampleName = sampleName.split('.')[0]         
        inFile = open(i,'r')                          
        faaOut = open(outFileName+'/faaFiles/'+sampleName+'.faa','w')     
        
        featureCount = 0

        for prot in SeqIO.parse(inFile,'genbank'):
    
            for feature in prot.features:            
 
                if feature.type == 'CDS':

                    featureCount += 1
                    
                    totalFeatures += 1
                    
                    assert len(feature.qualifiers['translation']) == 1
 
                    faaOut.write('>'+sampleName+'.'+'%s.%s\n%s\n' %(prot.name,feature.qualifiers['locus_tag'][0],feature.qualifiers['translation'][0])) 

                    nucleotideSeq = feature.location.extract(prot).seq       

                    nucleotideDictionary['>'+sampleName+'.'+prot.name+'.'+feature.qualifiers['locus_tag'][0]] = str(nucleotideSeq) 
        
        if len(nucleotideDictionary) != totalFeatures:
            sys.exit('Number of stored features does not equal number of found features')

        if featureCount == 0:
            sys.exit('Something is wrong with file '+i+' no CDS\'s were found.')

    return(nucleotideDictionary,outFileName+'/faaFiles/')            
    
    
#########################################################
#('Def') OrthofinderRunner() comments
#Takes as input ('Int') number of cores to use ('numCores')
#####('String') file path to the folder containing amino acid fasta files ('faaFiles')
#####('String') user supplied date ('date')
#Runs orthofinder program
#Returns ('String') file path to the orthofinder output file ('orthoText')
def OrthofinderRunner(numberOfCores,fastaAminoAcidFolder):

    date = datetime.datetime.now()
    # date = date.strftime('%b%d%H%M')
    date = date.strftime('%b%d')
    
    try:

        # return_code = subprocess.check_output(
        #     ['orthofinder','-f',fastaAminoAcidFolder,'-t',numberOfCores,'-og', '-n', date],
        #     stderr=subprocess.STDOUT)
        return_code = subprocess.check_output(['orthofinder','-f',fastaAminoAcidFolder,'-t',numberOfCores,'-og'],stderr=subprocess.STDOUT)
    
    except subprocess.CalledProcessError as error:
        
        print('Status : Fail',error.returncode,error.output.strip('\n'))

    #orthofinderText = fastaAminoAcidFolder+'/OrthoFinder/Results_'+date+'/Orthogroups/Orthogroups.txt'
    orthofinderText = fastaAminoAcidFolder+'/Results_'+date+'/Orthogroups.txt'
    
    return(orthofinderText)

#########################################################
#('Def') OrthofinderParing() comments
#Takes as input ('Int') number of isolates originally submitted ('numberOfIsolates')
#####('Float') percent inclusion to be considered an orthogroup ('percentInclusion')
#####('List') of ('OrthofinderObjectList')('Objects') ('obList')
#Calculates the number is isolates in an orthogroup needed in order to pass filter ('numberNeeded')
#('For') loops through ('List') of ('OrthofinderResults')('Objects') and checks to see if an orthogroup has a sufficient number
#####First checks for multiple of same isolate in an orthogroup
#####appends seqNames into a ('List')('SeqNameList') after split on ('.') to get isolate name
#####creates a ('Set')('seqNameSet') and compares ('Set') and ('List') lengths 
#####if the same then checks to see if orthogroup has appropriate number of isolates in it
#returns ('List') ('paredDownList') of ('OrthofinderResults')('Objects') 
def OrthofinderParing(numberOfIsolates,percentInclusion,orthofinderObjectList):

    paredDownList = []

    numberNeeded = numberOfIsolates*percentInclusion
    
    for orthogroups in orthofinderObjectList:
        seqNameList = []

        for sequences in orthogroups.seqNames:
            seqNameList.append(sequences.split('.')[0])
        seqNameSet = set(seqNameList)

        if len(seqNameList) == len(seqNameSet):
        
            if len(orthogroups.seqNames) >= numberNeeded:
            
                paredDownList.append(orthogroups)

    return(paredDownList)

#########################################################
#('Def') MultifastaGenerator() comments
#Takes as input single ('OrthofinderResults') ('Object') ('ofrObject')
#####('Dictionary') from ('Def') gbkParser() (dictionary[SequenceName]:'Sequence') ('nucDict')
#####('String') user supplies out directory ('outName')
#opens multifasta file 
#Uses ('Dictionary') search to print the Sequence name and sequence in fasta format
#Returns ('String') file path to the multifasta file 
def MultifastaGenerator(orthofinderObject,nucleicAcidDict,outFile):
    
    f = open(outFile+'/'+orthofinderObject.orthoNumber+'.fasta','w')

    for sequence in orthofinderObject.seqNames:
        
        print('>'+sequence.split('.')[0],file=f)
        print(nucleicAcidDict['>'+sequence],file=f)

    f.close()

    numberOrthogroupIsolates = len(orthofinderObject.seqNames)

    return(outFile+'/'+orthofinderObject.orthoNumber+'.fasta',numberOrthogroupIsolates)

#########################################################
#('Def') MuscleFileGenerator() comments
#Takes as input ('String') file path to a multifasta file ('multifastaFile')
#Runs muscle
#Returns ('String') file path to the muscle file
def MuscleFileGenerator(multifastaFile):
    
    muscleFile = multifastaFile.split('.')[0]
    
    # return_code = subprocess.check_output(['muscle','-align',multifastaFile,'-output',muscleFile+'.muscle', '-threads', '2', '-quiet'])
    # return_code = subprocess.check_output(['muscle','-align',multifastaFile,'-output',muscleFile+'.muscle'])
    return_code = subprocess.check_output(['muscle','-in',multifastaFile,'-out',muscleFile+'.muscle','-seqtype','dna', '-quiet'])

    return(muscleFile+'.muscle')

#########################################################
#('Def') TrimAlFileGenerator() comments
#Takes as input ('String') file path to a muscle file ('muscleFile')
#Runs trimAl 
#Returns ('String') file path to the trimAl file
def TrimAlFileGenerator(muscleFile):

    trimAlFile = muscleFile.split('.')[0]

    return_code = subprocess.check_output(['trimal','-in',muscleFile,'-out',trimAlFile+'.trimAl','-fasta','-automated1'])

    return(trimAlFile+'.trimAl')

#########################################################
#('Def') ConsambigFileGenerator() comments
#Light string manipulation to gather consambig inputs
#####('fileName') is the name of the output file
#####('directoryPath') is the path to the directory where the output will be 
#Runs EMBOSS consambig 
#Returns ('String') file path to the consambig file
def ConsambigFileGenerator(trimAlFile):

    consambigFile = trimAlFile.split('.')[0]

    fileName = consambigFile.split('/')[-1]

    directoryPath = '/'.join(map(str,consambigFile.split('/')[:-1]))

    return_code = subprocess.check_output(['consambig','-sequence',trimAlFile,'-outseq',consambigFile+'.fa','-name',fileName,'-snucleotide','-sformat1','fasta','-osformat2','fasta','-osdirectory2', directoryPath])

    return(consambigFile+'.fa')

#########################################################
#('Def') RunPrimer3() comments
#Takes as input ('String') file path to a consambig file ('consambigFile')
#####('List') of ('Strings') containing all the global primer3 values used in a boulder file
#reads in ('consambigFile') and creates the Sequence information for primer3
#generates the boulder file
#Runs primer3
#####Captures stdout from primer3 and generates a ('Dictionary') ('primer3Dict')
#Returns ('Dictionary') (dictionary[Orthogroup]:['Primer3',...,'Output']
def RunPrimer3(consambigFile,primer3List):
    
    primer3BoulderFile = consambigFile.split('.')[0]

    consensusInfo = list(SeqIO.parse(consambigFile,'fasta'))

    f = open(primer3BoulderFile,'w')

    orthogroup = ('SEQEUNCE_ID='+ consensusInfo[0].id)
    
    consensusSeq = ('SEQUENCE_TEMPLATE='+ re.sub(r'[a-z]','N',str(consensusInfo[0].seq).strip('\n')))
    
    print(orthogroup +'\n'+consensusSeq,file=f)
    for i in primer3List:
        print(i,file=f)
    f.close()
    
    primer3Dict = {}
    
    r1 = subprocess.Popen(('primer3_core',primer3BoulderFile),stdout=subprocess.PIPE)
    r2 = (r1.stdout.read().decode('ascii'))
    
    primer3Dict[(r2.split('\n')[0]).split('=')[1]]=r2.replace('=','\n').split('\n')[1:]

    return(primer3Dict)

#########################################################
#('Def') PrimerFlankingRegionCheck() comments
#Takes as input ('List') of ('Primers') ('Objects')
#Slices out the primer flanked region of the sequence
#Checks the sequence for inclusion of SNPs
#Returns ('List') of ('Primers') ('Objects') ('clearedPrimers')
def PrimerFlankingRegionCheck(primerObjectList):

    clearedPrimers = []
    
    for primer in primerObjectList:
        
        sequence = primer.orthogroupInfo.sequence[(primer.leftHit+primer.leftLen)-1:(primer.rightHit)-1]

        for letter in sequence:
            
            if letter not in ('ATGCN'):

                clearedPrimers.append(primer)
                break
    
    return(clearedPrimers)
                

#########################################################
#('Def') PrimmersearchRunner() comments
#Takes as input ('String') trimalFile path
##### ('List') of ('Primers') ('Objects')
#Creates ('IO') primer file for primersearch 
##### Pulls primer number/forward/reverse sequences from ('Primers') ('Object')
#####Stores information in ('List') primerInfoList
#Runs primersearch using subprocess
#Returns ('String') primersearchFile.ps file path
def PrimersearchRunner(trimalFile, primerObjectList):

    primersearchFile = trimalFile.split('.')[0]

    primerFileList = []

    for primers in primerObjectList:
        
        primerInfoList = []
        
        primerInfoList.append(str(primers.number)+'\t'+primers.leftSeq+'\t'+primers.rightSeq)

        primerFileList.append(primerInfoList)

    
    f = open(primersearchFile+'Primers','w')

    for primerInfo in primerFileList:
        for primers in primerInfo:
            print(primers,file=f)

    f.close()
        
    command = ['primersearch','-seqall',trimalFile,'-infile',primersearchFile+'Primers','-mismatchpercent','6','-outfile',primersearchFile+'.ps']
    process = subprocess.run(command, capture_output=True, text=True)
    # primers file might be empty, which would cause primersearch throw out error
    if process.returncode == 0:
        return(primersearchFile+'.ps')
    # return_code = subprocess.check_output(['primersearch','-seqall',trimalFile,'-infile',primersearchFile+'Primers','-mismatchpercent','6','-outfile',primersearchFile+'.ps'])

    # return(primersearchFile+'.ps')


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
def PrimersearchValidator(ampliconInfo,numberIsolates):

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
def PrimersearchComber(ampliconInfo,primer,sequenceRecordDict):

    primersearchObjects = []
    
    for line in ampliconInfo:
        
        if line.startswith('Amplimer'):
            try: 
                forwardHit = int(ampliconInfo[ampliconInfo.index(line)+3].split(' ')[5])
                reverseHit = int((ampliconInfo[ampliconInfo.index(line)+4].split(' ')[5]).replace('[','').replace(']',''))
                print(ampliconInfo[ampliconInfo.index(line)+1])
                ampliconLen = int(ampliconInfo[ampliconInfo.index(line)+5].split(' ')[2])
                sequenceName = str(ampliconInfo[ampliconInfo.index(line)+1].split(' ')[1].strip(' ').strip('\n'))

                sequence = str(sequenceRecordDict[sequenceName].seq)
                sequence = (sequence.replace('-',''))

                sequence = sequence[forwardHit+primer.leftLen:ampliconLen-reverseHit-primer.rightLen]

                primersearchObject = PrimerSearchResults(primer,sequenceName,ampliconLen,forwardHit,reverseHit,sequence)

                primersearchObjects.append(primersearchObject)
            except IndexError as ie:
                print (str(ie))
                continue


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
def PrimersearchParser(primersearchFile,numberIsolates,primerObjectList,trimalFile):

    with open(primersearchFile,'r') as f:
        primersearchInfo = f.readlines()
    f.close()

    primersearchObjectLists = []
    
    primersearchObjectList = []


    sequenceRecordDict =SeqIO.to_dict(SeqIO.parse(trimalFile,'fasta'))


    for primer in primerObjectList:
        
        ampliconInfoStartIndex = primersearchInfo.index('Primer name '+str(primer.number)+'\n')

        ampliconInfo = primersearchInfo[ampliconInfoStartIndex:ampliconInfoStartIndex+(numberIsolates*6)+1]

        validationMark = PrimersearchValidator(ampliconInfo,numberIsolates)

        if validationMark == True:
            # print (f"ampliconInfo is: {ampliconInfo}")
            # print (f"primer is: {primer}")
            # print (f"sequenceRecordDict is: {sequenceRecordDict}")
            parsedInfo = PrimersearchComber(ampliconInfo,primer,sequenceRecordDict)
        
            primersearchObjectLists.append(parsedInfo)

        else:

            break

    for primersearchObjects in primersearchObjectLists:

        primersearchObjectList = primersearchObjectList + primersearchObjects

    return(primersearchObjectList)


#########################################################
def ParallelFunctions(orthofinderObject,nucleicAcidDict,outfile,primerDesignInfo):

    print(str(orthofinderObject) +'is in the parallel universe')
    multifasta = MultifastaGenerator(orthofinderObject,nucleicAcidDict,outfile)

    muscleFile = MuscleFileGenerator(multifasta[0])

    trimalFile = TrimAlFileGenerator(muscleFile)

    consambigFile = ConsambigFileGenerator(trimalFile)
    primer3Info = RunPrimer3(consambigFile,primerDesignInfo)
    primerObjectList = primer3Parser(primer3Info)
    checkedPrimerObjectList = PrimerFlankingRegionCheck(primerObjectList)
    primersearchFile = PrimersearchRunner(trimalFile,checkedPrimerObjectList)

    primersearchObjectList = []
    if primersearchFile: ### PrimersearchRunner might return None if there the primer is not 'valid', containing no sequences
        primersearchObjectList = PrimersearchParser(primersearchFile,multifasta[1],checkedPrimerObjectList,trimalFile)

    return(primersearchObjectList)


#########################################################
#('Def') BoulderIOParser() comments
#Takes as input ('String') boulderFile file path
#Parses ('List') designInfo
#Returns ('List') primerDesignInfo
def BoulderIOParser(boulderFile):

    with open(boulderFile, 'r') as f:
        
        designInfo = f.readlines()

    f.close()
    
    primerDesignInfo = []

    for info in designInfo:
        
        primerDesignInfo.append(info.strip('\n'))

    return(primerDesignInfo)


#########################################################
#('Def') IsolateNumberVerification() comments
#Takes as input ('List') gbkFiles
#####('Int') userSuppliedNumber
#Confirms number of files taken in is equal to number of files expected
#Returns ('True') is numbers match
#Returns ('False') if number do not match
def IsolateNumberVerification(gbkFiles,userSuppliedNumber):

    correctNumber = False
    
    if len(gbkFiles) == userSuppliedNumber:

        correctNumber = True

    return(correctNumber)


#########################################################
#('Def') parallelRunner() comments
#Takes as input ('List') ('OrthofinderResults') ('Objects')
#Sets pool to number of processes to be used
#Uses partial for ('Def') ParallelFunctions() to set static arguments
#####Allows for looping through ('List') ('OrthofinderResults') ('Objects')
#Returns ('lists') of ('PrimersearchResults') ('Objects')
def parallelRunner(orthofinderResultObjects):

    pool = multiprocessing.Pool(processes=args.NumberCores)

    parallelStatic = partial(ParallelFunctions,nucleicAcidDict=parsedGbk[0],outfile=args.OutDirectory,primerDesignInfo=boulderio)

    result_list = pool.map(parallelStatic,orthofinderResultObjects)
    
    print(result_list)

#########################################################
import argparse
import re
import os
import glob
from Bio import SeqIO
import subprocess
import datetime
import multiprocessing
from functools import partial
import sys


parser = argparse.ArgumentParser(description='T3Pio has been designed to generate primers for a partial core genome MLST scheme')

parser.add_argument('NumberIsolates',type=int,help='Input the number of isolates that will be analyed')

parser.add_argument('PercentInclusion',type=float,help='Input the percent inclusion you would to have for your partial core MLST scheme design')

parser.add_argument('BoulderIOFile',type=str,help='Provide the path to a partially filled out boulder.io file for primer3. Leave the SEQUENCE_ID and SEQUENCE_TEMPLATE lines out of the file')

parser.add_argument('GbkFolder',type=str,help='Provide the path to the folder containing the GBK files to analyzed')

parser.add_argument('OutDirectory',type=str,help='Provide the directory where you wish intermediate and final files to be stored')

parser.add_argument('NumberCores',type=int,help='Input the number of cores you would like dedicate to orthofinder and analysis')

args = parser.parse_args()



gbkFiles = sorted(glob.glob(args.GbkFolder+'*.gbk'))


numberOfFilesCheck = IsolateNumberVerification(gbkFiles,args.NumberIsolates)

if numberOfFilesCheck == False:
    
    sys.exit('The number of files gathered and the number of files entered are not the same.  \n Please correct and restart.')

    

print('Starting GBK File Parsing')
parsedGbk = gbkParser(gbkFiles,args.OutDirectory)

print('Starting orthofinder run')
orthofinderRun = OrthofinderRunner(str(args.NumberCores),parsedGbk[1])

print('Starting orthofinder file parsing')
parsedOrthofinder = orthoFinderParser(orthofinderRun)

print('Starting orthofinder output paring')
paredOrthofinderResults = OrthofinderParing(args.NumberIsolates,args.PercentInclusion,parsedOrthofinder)

boulderio = BoulderIOParser(args.BoulderIOFile)



if __name__ == '__main__':

    parallelRunner(paredOrthofinderResults)
    
    

import argparse
import os
import subprocess
import glob
import re
import numpy
import csv
from Bio import SeqIO


#Setting up the parser to parse input
parser = argparse.ArgumentParser(description='This script is for the generation of primer candidates')

#Adding the parser for the InputGbkFolder
parser.add_argument('InputGbkFolder',type=str,help='Please input a folder containing prokka annotated gbk files')

#Closing the parser
args = parser.parse_args()

#Creates array of file locations
GbkFiles = glob.glob(args.InputGbkFolder + '/*.gbk')

numGBKFiles= len(GbkFiles)

#Create the fna and faa directories if they do not already exist
if not os.path.exists('fnaFiles/'):
    os.makedirs('fnaFiles/',exist_ok=True)
if not os.path.exists('faaFiles/'):
    os.makedirs('faaFiles/',exist_ok=True)

print('Parsing GBK Files:')

#Creat while loop to go through each gbk file and extract the nucleotide sequence and amino acid sequence in fasta format. Use basic SeqIO.convert to get fna and then loops and features funtion to pull out faa
i=0
while i < len(GbkFiles):
    sample = GbkFiles[i].split('/')[16]
    sample = sample.split('.')[0]
    GbkIn = open(GbkFiles[i],'r')
    faaOut = open('faaFiles/'+sample+'.faa','w')
    fnaOut =open('fnaFiles/'+sample+'.fna','w')
#Still convinced is magic
    for prot in SeqIO.parse(GbkIn,'genbank'):
        for feature in prot.features:
            if feature.type=='CDS':
                assert len(feature.qualifiers['translation'])==1
                faaOut.write('>'+sample+'.'+'%s.%s\n%s\n' %(prot.name,feature.qualifiers['locus_tag'][0],feature.qualifiers['translation'][0]))
                l2 = feature.location.extract(prot).seq
                fnaOut.write('>'+sample+'.'+'%s.%s\n%s\n' %(prot.name,feature.qualifiers['locus_tag'][0],str(l2)))
    GbkIn.close()
    faaOut.close()
    fnaOut.close()
    i += 1

print('Starting OrthoFinder')
#Call subprocess to run orthofinder. Deposits Results directory into the faaFiles folder

return_code = subprocess.check_output(['orthofinder', '-f', 'faaFiles','-t','21','-og'])


#Get path to single copy orthogroups text file to be used to pull out for core orthogroups

with open('faaFiles/Results_Jan10/Orthogroups.csv') as f:
    reader = csv.reader(f)
    orths = [r for r in reader]

#Move through list of lists orths[i][0]

numIso = len(orths[0][0].split('\t')[1:])


testog = []
finalOrthGroupDict = {}

i=1
while i < len(orths):
    orthGroupDict = {}
    og = orths[i][0].split('\t')[0]
    isoList = orths[i][0].split('\t')[1:]
    isoList = list(filter(lambda a: a != '',isoList))
    orthGroupDict[og] = isoList
    if (len(orthGroupDict[og])/numIso) > 0.95 and len(orthGroupDict[og]) <= numIso:
        setList = []
        for iso in isoList:
            baseName = iso.split('.')[0]
            setList.append(baseName)
        setList = set(setList)
        if len(setList) == len(isoList):
            testog.append(og)
            finalOrthGroupDict[og] = orthGroupDict[og]
    i += 1

if not os.path.exists('/scratch/NewHMAS/ogFastaFiles/'):
    os.makedirs('/scratch/NewHMAS/ogFastaFiles/',exist_ok=True)

for k in finalOrthGroupDict.keys():
    ogFasta = open('/scratch/NewHMAS/ogFastaFiles/'+k+'.fasta','a')
    print(k)
    for iso in finalOrthGroupDict[k]:
        baseName = iso.split('.')[0]
        with open('fnaFiles/'+baseName+'.fna','r') as line:
            sequence = line.readlines()
        line.close()
        seq = '>'+iso+'\n'
        seqIndex = sequence.index(seq)
        ogFasta.write(seq+sequence[seqIndex+1])

#On to MSA using muscle

if not os.path.exists('muscleAlignmentFiles/'):
    os.makedirs('muscleAlignmentFiles/',exist_ok=True)

fastaFiles = sorted(glob.glob('/scratch/NewHMAS/ogFastaFiles/*.fasta'))


#Set up script to run Muscle Aignment
#Basic file name handling and subprocess input
#Loops through all multifasta files
r=0
while r < len(fastaFiles):
    inFileName = fastaFiles[r]
    outFileName=fastaFiles[r].split('/')[4]
    outFileName=('muscleAlignment'+outFileName)
    return_code=subprocess.check_output(['muscle','-in', inFileName ,'-out', 'muscleAlignmentFiles/'+outFileName,'-seqtype', 'dna'])
    r += 1


#Trumper Section

if not os.path.exists('trimalTrimmedFastaFiles/'):
    os.makedirs('trimalTrimmedFastaFiles/',exist_ok=True)

msaFastaFiles = sorted(glob.glob('muscleAlignmentFiles/*.fasta'))

t=0
while t < len(msaFastaFiles):
    inFileName = msaFastaFiles[t]
    outFileName = msaFastaFiles[t].split('/')[1]
    outFileName = (outFileName)[15:]
    return_code = subprocess.check_output(['trimal','-in',inFileName,'-out','trimalTrimmedFastaFiles/trimmed'+outFileName,'-fasta','-automated1'])
    t += 1

#Emboss consensus section

if not os.path.exists('embossConsensusAlignments/'):
    os.makedirs('embossConsensusAlignments/',exist_ok=True)

trimmedMsaFiles = sorted(glob.glob('trimalTrimmedFastaFiles/*.fasta'))

e=0
while e < len(trimmedMsaFiles):
    inFile = trimmedMsaFiles[e]
    outFileName = trimmedMsaFiles[e].split('/')[1]
    outFileName = (outFileName)[7:]
    outFileName = outFileName.split('.')[0]
    return_code = subprocess.check_output(['consambig','-sequence',inFile,'-outseq','embossConsensusAlignments/'+outFileName+'Consensus.fasta','-name',outFileName+'Consensus.fasta','-snucleotide1','-sformat1', 'fasta','-osformat2','fasta','-osdirectory2','embossConsensusAlignments/'])
    e += 1
    print('Working on file number: ' + str(e) + ' Which is file name: ' + outFileName )

if not os.path.exists('primer3FilesTest/'):
    os.makedirs('primer3FilesTest/',exist_ok=True)

consensusFiles = sorted(glob.glob('embossConsensusAlignments/*.fasta'))


if not os.path.exists('primer3Files/'):
    os.makedirs('primer3Files/',exist_ok=True)

consensusFiles = sorted(glob.glob('embossConsensusAlignments/*.fasta'))

#Primer3 block. Takes in the emboss consensus files and goes through one at a time. Extracts the sequence name and sequence and stores them in boulder IO
#format. Then prints them to a file called boulderFile with all of the other options that are wanted for primer3. Once the file is created primer3 is run
#through subprocess and the output, which is normally to stdout is captured and saved to a file that includes the OG# in the title. This process repeats,
#with the boulder file being overwritten each time until all consensus sequences have been processed. Also replaces lowercase letters with N's in the seq

p=0
while p < len(consensusFiles):
    primer3Out=((consensusFiles[p])[25:]).split('.')[0]
    print(primer3Out)
    with open(consensusFiles[p], 'r') as seq:
        seqTmp='SEQUENCE_TEMPLATE='
        for line in seq:
            if '>' in line:
                seqID =('SEQUENCE_ID='+line.strip('>')).split('.')[0]
            else:
                seqTmp =(seqTmp + re.sub(r'[a-z]','N',line).strip('\n'))
    f = open('boulderFile','w')
    print(seqID, file=f)
    print(seqTmp, file=f)
    print('PRIMER_TASK=generic', file=f)
    print('PRIMER_NUM_RETURN=10', file=f)
    print('PRIMER_PICK_LEFT_PRIMER=1', file=f)
    print('PRIMER_PICK_RIGHT_PRIMER=1', file=f)
    print('PRIMER_OPT_SIZE=21', file=f)
    print('PRIMER_MIN_SIZE=18', file=f)
    print('PRIMER_MAX_SIZE=27', file=f)
    print('PRIMER_MAX_NS_ACCEPTED=3', file=f)
    print('PRIMER_PRODUCT_SIZE_RANGE=180-250', file=f)
    print('P3_FILE_FLAG=0', file=f)
    print('PRIMER_EXPLAIN_FLAG=0', file=f)
    print('PRIMER_FIRST_BASE_INDEX=1', file=f)
    print('PRIMER_MIN_TM=57.0', file=f)
    print('PRIMER_MAX_TM=63.0', file=f)
    print('PRIMER_OPT_TM=60.0', file=f)
    print('PRIMER_PAIR_MAX_DIFF_TM=1.0', file=f)
    print('PRIMER_LIBERAL_BASE=3', file=f)
    print('PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/apps/x86_64/primer3/primer3-2.3.4/src/primer3_config/', file=f)
    print('=', file=f)
    f.close()
    return_code=subprocess.call(['primer3_core','boulderFile'],stdout=open('primer3Files/'+primer3Out+'primers.txt','w'))
    p += 1


#Amp Filter Process

if not os.path.exists('filteredPrimers/'):
    os.makedirs('filteredPrimers/',exist_ok=True)

rawPrimer3Files = sorted(glob.glob('primer3Files/*.txt'))

#Want to inital sort to see if file has any actual primers or sequence_template was too short
#Code takes the input primer3 folder and iterates through until done. To start, a file is read in and then checked to see if it is the correct size. A length of 22 indicates that the sequence_template was too short, a length of 26 indicates that no primers were returned. If the file is appropriately large, the sequence template is saved for later use. Then the number of primer pairs is found to make sure all primer pairs are assesed. Left stop and right start positions are found in the while loop which looks for the each primer pairs left primer and then using line and indexing the right primer is found. The coordinates are calculated and then used to slice the sequence into the amplicon region. This region is then checked for and bases that are not ATGCN and if there are variable regions the Orthogroup and primer pair are saved to a file for primersearch later. This is done for all primer pairs in all files that make it through the filters

q=0
while q < len(rawPrimer3Files):
    OG = (rawPrimer3Files[q]).split('/')[1]
    OG = (OG)[:9]
    print(OG)
    with open(rawPrimer3Files[q],'r') as pI:
        primerInfo = pI.readlines()
        if len(primerInfo) != 22 and len(primerInfo) != 26:
            sequence = primerInfo[1]
            sequence = sequence.split('=')[1]
            numIndex = [primerInfo.index(i) for i  in primerInfo if 'PRIMER_PAIR_NUM_RETURNED=' in i]
            numPairReturned = primerInfo[numIndex[0]]
            numPairReturned = numPairReturned.split('=')[1]
            count = 0
            while count < int(numPairReturned):
                primNum = ('PRIMER_LEFT_'+str(count)+'=')
                for line in primerInfo:
                    if primNum in line:
                        leftPrimer = primerInfo[primerInfo.index(line)-2].split('=')[1].strip('\n')
                        if 'N' in leftPrimer:
                            break
                        rightPrimer = primerInfo[primerInfo.index(line)-1].split('=')[1].strip('\n')
                        if 'N' in rightPrimer:
                            break
                        leftStop = line
                        leftStop = leftStop.split('=')[1]
                        leftStop = leftStop.split(',')
                        leftStop = int(leftStop[0])+int(leftStop[1])
                        rightStart = (primerInfo[primerInfo.index(line)+1])
                        rightStart = rightStart.split('=')[1]
                        rightStart = rightStart.split(',')
                        rightStart = int(rightStart[0])-int(rightStart[1])
                        amplicon = sequence[leftStop:rightStart]
                        for letter in amplicon:
                            if letter not in ('ATGC'):
                                primersOut = open('filteredPrimers/'+OG+'filteredPrimers','a')
                                primersOut.write(OG+'primerGroup'+ str(count) + '\t' + ((primerInfo[primerInfo.index(line)-2]).split('=')[1]).strip('\n') + '\t' + ((primerInfo[primerInfo.index(line)-1]).split('=')[1]).strip('\n')+'\n')
                                primersOut.close()
                                break
                count += 1
    q+=1


#Primer Specificity

if not os.path.exists('specificityPrimerSearch/'):
    os.makedirs('specificityPrimerSearch/',exist_ok=True) #create directory to hold files

specPrimerFiles = sorted(glob.glob('filteredPrimers/*')) #grab filtered primer pairs ( amplicons contain SNPs and no N bases in the primer pairs) 
#Grab metaGenome Assemblies
#Loop through both, meta being outer loop

krakenFiles = sorted(glob.glob('krakenSequences/*')) #grabs kraken files for dicovery of unclassified contigs
metaFiles = sorted(glob.glob('metaGenomeAssemblies/*')) #grabs the metagenome fasta files
print(krakenFiles)


#code block to run primerSearch on the metagenomes with the filtered primers
#uses subprocess to call primerSearch (EMBOSS/6.4.0)

mf = 0
while mf < len(metaFiles):
    spf = 0
    while spf < len(specPrimerFiles):
        metaName = metaFiles[mf].split('/')[1]
        metaName = metaName.split('.')[0]
        ogName = ((specPrimerFiles[spf])).split('/')[1]
        ogName = ogName[:9]
        print(ogName)
        return_code= subprocess.check_output(['primersearch','-seqall',metaFiles[mf],'-infile',specPrimerFiles[spf],'-mismatchpercent','6','-outfile','specificityPrimerSearch/'+ogName+metaName])
        spf += 1
    mf += 1


#Final Specificity Part

transFiles = sorted(glob.glob('translatedKrakenSequences/*')) #grabs the translated kraken files 
unclassFiles = sorted(glob.glob('unclassifiedSequences/*')) #grabs the unclassified contig files

specDict={} #specification dictionary to tag contigs that are either unclassified or corresponding to an organism that isn't salmonella/bug of interest

tf = 0 
while tf < len(transFiles):
    tranFileName = transFiles[tf].split('/')[1]
    tranFileName = tranFileName.split('.')[0]
    contigList = []
    with open(transFiles[tf],'r') as tran:
        tranList = tran.readlines()
    tran.close()
    with open(unclassFiles[tf],'r') as unc:
        unclassList = unc.readlines()
    unc.close()
    for line in tranList:
        contig = line.split('\t')[0]
        species = line.split(';')
        if 'Salmonella' not in species[len(species)-1]:
            contigList.append(contig)
    for line in unclassList:
        contigList.append(line.strip('\n'))
    specDict[tranFileName] = contigList #fills the dictionary with list of undesirable contigs keyed by the metagenome name
    tf += 1

specPrimerSearchFiles = sorted(glob.glob('specificityPrimerSearch/*')) #grabs the primerSearch result files

#test = specPrimerSearchFiles[2].split('/')[1]
#test = test[9:len(test)-7]
badPrimerPairs = [] #list to hold any primer pairs that are found to be non-specific. each element is formated as OG#######primerGroup# so it can be tracked and used properly further down the line

sps = 0 
while sps < len(specPrimerSearchFiles):
    metaName = specPrimerSearchFiles[sps].split('/')[1]
    metaName = metaName[9:]                     
    with open(specPrimerSearchFiles[sps],'r') as spsf:
        searchLines = spsf.readlines()
    spsf.close()
    count = 0
    dictContigs = specDict[metaName]
    while count < len(searchLines):
        if searchLines[count].startswith('Amplimer'):
            checkContig = searchLines[count +1].split(':')[1]
            checkContig = checkContig.strip(' ')
            checkContig = checkContig.strip('\n')
            if (checkContig in dictContigs):
                badPrimerPairs.append(searchLines[count - 1])
                print(badPrimerPairs) 
            elif searchLines[count].startswith('Amplimer 2'):
                badPrimerPairs.append(searchLines[count - 7])
                print(badPrimerPairs)
            count += 1
                                      
        else:
            count += 1
    sps += 1

badPrimerPairs = list(set(badPrimerPairs))

outSpec = open('specificityFilteredPrimerPairList.txt','w')

for line in badPrimerPairs:
    outSpec.write(line)

outSpec.close()


#Primer scoring block. Weights the score based on position of variable site(3' end bad) and length of primer so longer primers with a variable base in the same 5' position as a shorted primer will be be prefered. 
#Both left and right primer scores are combined for the total score for comparison
def primerScore(primer):
    spot = 0
    position = 1
    for letter in primer:
        if letter not in ('ATGC'):
            spot = position
        position += 1
    score = spot/(len(primer))
    return score

rawPrimer3Files = sorted(glob.glob('primer3Files/*')) #grabs the primer3 generated files

if not os.path.exists('overlapFilteredPrimers/'): #creates directory to store the files generated from the script
    os.makedirs('overlapFilteredPrimers/',exist_ok=True)

badDict = {} #dictionary to hold the orthogroup and associated bad primers as decided from specificity testing
with open('specificityFilteredPrimerPairList.txt','r') as b: #file that contains the bad primer pair information
    badPrimer = b.readlines() #b is for bad
b.close()
for line in badPrimer:
    name = line.split(' ')[2]
    orth = name[:9] #important to key by orthogroup since the primer3 files are named from just the orthogroup
    pNum = name[9:].strip('\n') #the bad primer pair for an orthogroup
    badDict.setdefault(orth, []) #allows dictionary value to be an appendable list
    badDict[orth].append(pNum) #fills that list 

q=0
while q < len(rawPrimer3Files):
    OG = (rawPrimer3Files[q]).split('/')[1]
    OG = (OG)[:9] #used for naming convention and associating with bad primers
    print(OG)
    badPrimerList = [] #makes a list of bad primers for this particular orthogroup 
    try: 
        badPrimerList = badDict[OG] #get correct list for orthogroup
    except KeyError:
        print('not in list') #some(most) orthogroups do not have bad primer pairs so expception is needed
    snpDict = {} #contain the snps for each primer pair found for an orthogroup
    primerDict = {} #dictionary to hold the left/right primer for each primer pair
#Pretty much the same filtering script body as from the first pass with additional logic to remove primer pairs that were tagged as bad from specificity
    with open(rawPrimer3Files[q],'r') as pI:
        primerInfo = pI.readlines()
        if len(primerInfo) != 22 and len(primerInfo) != 26: #length of 22 indicates the target sequence was to short and length of 26 indicates no primer pairs were found
            sequence = primerInfo[1]
            sequence = sequence.split('=')[1]
            numIndex = [primerInfo.index(i) for i  in primerInfo if 'PRIMER_PAIR_NUM_RETURNED=' in i] #finds the line with number of primer pairs found
            numPairReturned = primerInfo[numIndex[0]] 
            numPairReturned = numPairReturned.split('=')[1] #gets integer value of how many primer pairs returned
            count = 0
            while count < int(numPairReturned):
                primNum = ('PRIMER_LEFT_'+str(count)+'=')
                for line in primerInfo:
                    if primNum in line:
                        groupNum = primNum[len(primNum)-2:len(primNum)-1]
                        ogPrimeNum = OG + 'primerGroup' + groupNum #Formating name
                        checkName = 'primerGroup'+groupNum #used to see if this is bad primer pair
#Need redundant N filtering because both the initial filter for specificity pull from the same place, primer3 files and therefore there is no master currated list to pull from
                        if checkName in badPrimerList: #logic to skip a bad primer pair
                            break
                        snpList = [] #list of snps for a primer pair
                        primerList = [] #list of left and right primer 
                        a = 0
                        leftPrimer = primerInfo[primerInfo.index(line)-2].split('=')[1].strip('\n')
                        if 'N' in leftPrimer:
                            break
                        rightPrimer = primerInfo[primerInfo.index(line)-1].split('=')[1].strip('\n') #logic to skip primers with any N's, left or right
                        if 'N' in rightPrimer:
                            break
                        primerList.append(leftPrimer)
                        primerList.append(rightPrimer) #puts primers into list
                        primerDict[ogPrimeNum] = primerList #fills dictionary with primer list keyed to primer number
                        leftStop = line
                        leftStop = leftStop.split('=')[1]
                        leftStop = leftStop.split(',')
                        leftStop = int(leftStop[0])+int(leftStop[1])
                        rightStart = (primerInfo[primerInfo.index(line)+1])
                        rightStart = rightStart.split('=')[1]
                        rightStart = rightStart.split(',')
                        rightStart = int(rightStart[0])-int(rightStart[1]) #finds the stop and start for each primer pair
                        amplicon = sequence[leftStop:rightStart] #extracts the amplicon
                        while a < len(amplicon):
                            site = amplicon[a]
                            if site not in ('ATGC'):
                                snpSite = leftStop + a 
                                snpList.append(snpSite) #logic for finding snp positions in the amplicon
                            a += 1
                        snpDict[ogPrimeNum] = snpList #adds snp list to snp dictionary
                count += 1 #loop through all primer pairs found
#Necessary because of the problem stated above, that the inital filtering and therefore primer sets is just to get a list of bad primers from specificity testing.
#This needs to be here to account for empty primer sets that might arrise from the 'redundant' filtering above
    for k in list(snpDict.keys()):
        if len(snpDict.get(k)) == 0:
            del snpDict[k]
#This is for filtering primer pairs that have the same captured SNPs
    olf=0
    while olf < len(snpDict):
        print(len(snpDict)) 
        keyRemoval = [] #track which primer pairs should be removed
        keyTrac = [] #track which primers have already been used to check 
        for k1 in list(snpDict.keys()):
            keyTrac.append(k1) #keep track of used primers
            left0 = primerDict.get(k1)[0]
            right0 = primerDict.get(k1)[1]
            left0score = primerScore(left0) #logic to get the primers for the first primer group to search with and generate primer score
            right0score = primerScore(right0)
            total0 = left0score + right0score
            for k2 in list(snpDict.keys()):
                if k1 != k2:
                    if snpDict.get(k1) == snpDict.get(k2):
                        #print(k1)
                        #print(k2)
                        left1 = primerDict.get(k2)[0]
                        right1 = primerDict.get(k2)[1]
                        left1score = primerScore(left1)
                        right1score = primerScore(right1) #same as above but for the primers being searched against
                        total1 = left1score + right1score
                        if total0 < total1: 
                            keyRemoval.append(k2) #remove the searched against primer is the search primer is better
                        elif total0 == total1:
                            if k2 not in keyTrac:
                                keyRemoval.append(k2) #remove searched agaisnt primer if they are equal
        keyRemoval = list(set(keyRemoval)) #remove any duplicate values that might arrise
        print(keyRemoval)
        for key in keyRemoval:
            del snpDict[key] #delete actual primers
        olf += 1
#This is for filtering primer pairs that are subgroups of another primer pair
    sub=0
    while sub < len(snpDict):
        print(len(snpDict))
        for k1 in list(snpDict.keys()):
            left0 = primerDict.get(k1)[0]
            right0 = primerDict.get(k1)[1]
            left0score = primerScore(left0)
            right0score = primerScore(right0)
            total0 = left0score + right0score
            for k2 in list(snpDict.keys()):
                if k1 != k2:
                    list1 = snpDict.get(k1) 
                    list2 = snpDict.get(k2)
                    if list1 != None:
                        if set(list1).issubset(set(list2)):
                            left1 = primerDict.get(k2)[0]
                            right1 = primerDict.get(k2)[1]
                            left1score = primerScore(left1)
                            right1score = primerScore(right1)
                            total1 = left1score + right1score
                            if total0 > total1: #Similar to above but looking are primer pairs that are subsets of another and keeping them if they have a better score
                                del snpDict[k1]
        sub += 1

    for k,v in snpDict.items():
        print(k)
        print(v)
    if len(snpDict) != 0:
        for k in list(snpDict.keys()):
            lprimer = primerDict.get(k)[0]
            rprimer = primerDict.get(k)[1]
            primerGroup = k.strip('=')[9:]
            fileName = k[:9]
            primersOut = open('overlapFilteredPrimers/'+fileName+'filteredPrimers','a')
            primersOut.write(fileName+primerGroup + '\t' + lprimer + '\t' + rprimer + '\n') #print out remaining primer pairs in primersearch format 
    q+=1


#PRIMER SEARCH! ALMOST DONE!But not really :'(

if not os.path.exists('primerSearchFiles/'):
    os.makedirs('primerSearchFiles/',exist_ok=True)

filteredPrimerFiles = sorted(glob.glob('overlapFilteredPrimers/*'))

f=0
while f < len(filteredPrimerFiles):
    ogName =((filteredPrimerFiles[f])).split('/')[1]
    ogName = ogName[:9]
    print(ogName)
    return_code= subprocess.check_output(['primersearch','-seqall','trimalTrimmedFastaFiles/trimmed'+ogName+'.fasta','-infile',filteredPrimerFiles[f],'-mismatchpercent','6','-outfile','primerSearchFiles/'+ogName+'primerSearched'])
    f += 1


#PrimerSearchFilter(1)

if not os.path.exists('filteredPrimerSearch/'):
    os.makedirs('filteredPrimerSearch/',exist_ok=True)

searchFiles = sorted(glob.glob('primerSearchFiles/*'))

#GbkFiles = sorted(glob.glob('GbkFolder/*'))
#numGBKFiles = len(GbkFiles)


fi = 0
while fi < len(searchFiles):
    numPrimerNames = []
#Read in primersearch files, one at a time and save them into a list
    with open(searchFiles[fi],'r') as sF:
        parseList = sF.readlines()
    sF.close()
#Initiate a variable to keep track of where in parseList the program is in
#Initiate a dictionary that will keep track of the orthogroup and later hold the primerpair(with relevent info as values) as its values 
    count = 0
    countDict = {}
    while count < len(parseList):
#step through the file, depending on what the line starts with
        if parseList[count].startswith('\n'):
            checkDict={}
            count += 1
        elif parseList[count].startswith('Primer name'):
            outFileName = parseList[count].split(' ')[2]
            count += 1
#Most information is contained under Amplimer. Once it sees that as the start of a line it will grab the information wanted and clean the information and store in a list and add it to the primerpair dictionary. 
        elif parseList[count].startswith('Amplimer'):
            tmpList = []
            sequence = parseList[count+1]
            sequence = (sequence.split(' ')[1]).strip('\n')
            size = parseList[count+2]
            size = (size.split(' ')[0]).strip('\t')
            forward = parseList[count+3]
            forward = (forward.split(' ')[5])
            reverse = parseList[count+4]
            reverse = (reverse.split(' ')[5]).replace('[','').replace(']','')
            tmpList = [size,forward,reverse]
            checkDict[sequence] = tmpList
            countDict[outFileName]= checkDict
            count += 6
        else:
            continue
#Worst part is here. Itterates over key value pairs and extracts the amplicons from the trimal trimmed msa files and stores the amplicons into a new msa fasta file as long as they meet the criteria. Such as not to many primerpair hits for a single primerpair
    for k,v in countDict.items():
        outFile = open('filteredPrimerSearch/'+k.strip('\n')+'.fasta','a')
        orthoName = k
        orthoName = orthoName[:9]
        with open('trimalTrimmedFastaFiles/trimmed'+orthoName+'.fasta','r') as ff:
            msa = ff.readlines()
        ff.close()
        msa = [s.strip('\n') for s in msa]
        msa= ''.join(msa)
#Most annoying part. Stupid trimal and its damned \n EVERYWHERE        
        delimiters = 'bp','>'
        regex = '|'.join(map(re.escape, delimiters))
        msa = re.split(regex, msa)
        msalen = int((len(msa)-1)/2)
        for k2,v2 in countDict[k].items():
            if len(v) == msalen:
                if len(v2) == 3:
                    k2 = k2+' '+v2[0]+' '
                    isolate = msa[msa.index(k2)]
                    sequence = msa[msa.index(k2)+1]
                    outFile.write('>'+isolate+'\n'+sequence[int(v2[1])-1:int(v2[0])-int(v2[2])-1]+'\n')
                    print(orthoName)
    
    fi += 1


#Creating the two directories where the bulk of the files will be stored

if not os.path.exists('alleleCallFiles/'):
    os.makedirs('alleleCallFiles/',exist_ok=True)
if not os.path.exists('orthogroupPrimerPairAlleles/'):
    os.makedirs('orthogroupPrimerPairAlleles/',exist_ok=True)
#Note that the sorted aspect allows for order when grabbing (a,b,c/1,2,3) and so on
primerParis = sorted(glob.glob('overlapFilteredPrimers/*'))#Primer pairs to attach to the allele matrix for each OGprimerPair
amplicons = sorted(glob.glob('filteredPrimerSearch/*')) #Fasta files contaning sequences
sequenceNames = sorted(glob.glob('GbkFolder/*')) # Need for naming convention
#seqNameList = [] 
#for name in sequenceNames:
#    realName = name.split('/')[1]
#    realName = realName.split('.')[0]
#    seqNameList.append(realName)
amp = 0 #Setting variable to loop through all files
OGPPNumAlleles = {} #dictionary to track how many alleles there are for all of the primer pairs
OGPPComparitorDict = {} #dictionary to track which alleles each isolate has for each primer pair
while amp < len(amplicons):
    ampDict = {} #dictionary to track which amplicon belongs to which isolate
    orderList = [] #list used to move through the dictionaries in an orderly fashion
    ampName = amplicons[amp].split('/')[1]
    ampName = ampName.split('.')[0]
    OG = ampName[:9] + 'filteredPrimers' #Name used to get primer pair file for each orthogroup
    primerFile = 'overlapFilteredPrimers/'+OG #Used to open filtered primer file to get primer pairs
    primerGroup = ampName[9:]
    print(ampName)
    
    with open(amplicons[amp],'r') as a: #Reading in amplicon file
        ampList = a.readlines()
    a.close()
    with open(primerFile,'r') as p: #Reading in associated primer pairs file
        plist = p.readlines()
    alleleList = [] #list to be filled with allele group calls
    cnt = 0
    while cnt < len(ampList):
        ampDict[ampList[cnt].split('.')[0]] = ampList[cnt+1] #filling in the amplicon dictionary with isolate and associated amplicon
#important to split on . here due to erratic(read here as schizoid) naming convention
        orderList.append(ampList[cnt].split('.')[0]) #filling in the order list
        alleleList.append(ampList[cnt+1]) #filling in the allele list with amplicons
        cnt += 2

    alleleSet = set(alleleList) #remove duplicate amplicons leaving a set of unique amplicons
    alleleDict = {} #dictionary to keep track of which amplicon is associated with which allele group
    alcnt = 1
    for allele in alleleSet:
        alleleDict['Allele'+str(alcnt)] = allele #Filling the allele dictionary with allele groups and associated amplicon sequence
        alcnt += 1
    allOut = open('orthogroupPrimerPairAlleles/'+ampName+'Alleles','a') #create file that contains the allele group and amplicon sequence for each primer pair
    for k,v in alleleDict.items():
        print(k+'\t'+v,file=allOut)
    allOut.close()

    OGPPNumAlleles[ampName]= len(alleleDict) #fills dictionary with how many alleles there are for a given primer pair
    alleleMatrix = numpy.zeros(((len(ampDict)+1),(len(ampDict)+1)),dtype=object) #creates the matrix for pair wise comparison
    for line in plist:
        if primerGroup in line:
            lprimer = line.split('\t')[1] #pulls out the associated primers to add into the matrix
            rprimer = line.split('\t')[2]
            rprimer = rprimer.strip('\n')
    alleleMatrix[0][0] = lprimer + ' ' + rprimer #fills the first box with the primer pair
    z = 1
    w = 0
    while w < len(orderList):
        alleleMatrix[0][z] = orderList[w].strip('\n') #fills the 1 column/row with the name of the isolates
        alleleMatrix[z][0] = orderList[w].strip('\n')
        z += 1
        w += 1
    x = 0
    while x < len(ampDict):
        y=0
        while y < len(ampDict):
            if ampDict[orderList[x]] == ampDict[orderList[y]]: #compares the amplicons and if the are equal a 1 is placed in the associated matrix space
                alleleMatrix[x+1][y+1] = 1
            y += 1
        x += 1

    bugAlleleNumDict = {} #dictionary to track which allele is associated with each isolate
    for k,v in ampDict.items():
        bugAlleleNumDict[k] = list(alleleDict.keys())[list(alleleDict.values()).index(v)] #using indexing to grab the correct allele group for the amplicon that the isolate has
        OGPPComparitorDict.setdefault(k,[]) #enables a list to be created and appended in the dictionary that keeps track of all the alleles for all of the isolates
        OGPPComparitorDict[k].append(bugAlleleNumDict[k]) #adds the allele to the list

    bugAllOut = open('orthogroupPrimerPairAlleles/'+ampName+'SampleTags','a') #open file to hold the isolate and allele group info for each prim for k,v in bugAlleleNumDict.items():
    print(k.strip('\n') + '\t' + v,file=bugAllOut) #add to the file
    bugAllOut.close()

    #print(alleleMatrix[2][2])
    amp += 1
    numpy.savetxt('alleleCallFiles/'+ampName,alleleMatrix, delimiter='\t',fmt='%s') #saves the pair wise matrix in tab deliminated format

numAlleleOut = open('OrthogroupPrimerPairAlleleCount.txt','a') #generate file that has the number of alleles for each primer pair
for k,v in OGPPNumAlleles.items():
    print(k+'\t'+str(v),file=numAlleleOut)


w = csv.writer(open('Seans_Awesome_Alleles.csv','w'))
for k,v in OGPPComparitorDict.items():
    w.writerow([k,v])

sameAlleles = []
sampleToSampleList = []

for k1 in list(OGPPComparitorDict.keys()):
    for k2 in list(OGPPComparitorDict.keys()):
        if k1 != k2:
            if OGPPComparitorDict.get(k1) == OGPPComparitorDict.get(k2): #logic to go through and find identical isolates and tag them
                if k2 not in sameAlleles:
                    sameAlleles.append(k2)
                    sampleToSampleList.append(k1 + 'is the same as ' + k2)

sameOut = open('SameGenomeList.txt','a') #prints file containing all of the isolates that have at least one other identical isolate
for line in sameAlleles:
    print(line,file=sameOut)

sameSampOut = open('SampleToSampleSame.txt','a') #prints file that says which isolate is identical to which other isolate
for line in sampleToSampleList:
    print(line,file=sameSampOut)

print(len(sameAlleles))


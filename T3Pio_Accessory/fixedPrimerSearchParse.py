import re
import os
import glob
import subprocess
import argparse



#PRIMER SEARCH! ALMOST DONE!But not really :'(

if not os.path.exists('validationPrimerSearchFiles/'):
    os.makedirs('validationPrimerSearchFiles/',exist_ok=True)

filteredPrimerFiles = sorted(glob.glob('overlapFilteredPrimers/*'))


for i in filteredPrimerFiles:
    for k in assemblies:
        ogName =((filteredPrimerFiles[f])).split('/')[1]
        ogName = ogName[:9]
        print(ogName)
        return_code= subprocess.check_output(['primersearch','-seqall',k+'.fasta','-infile',i,'-mismatchpercent','6','-outfile','validationPrimerSearchFiles/'+ogName+'primerSearched'])
    


#PrimerSearchFilter(1)

if not os.path.exists('validationFilteredPrimerSearch/'):
    os.makedirs('validationFilteredPrimerSearch/',exist_ok=True)

searchFiles = sorted(glob.glob('validationPrimerSearchFiles/*'))

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
            leftAmp = parseList[count+3].split(' ')[0].strip('\t')
            leftSize = re.sub('([\(\[]).*?([\)\]])','\g<1>\g<2>', leftAmp)
            leftSize = leftSize.replace('[','N').replace(']','')
            leftSize = len(leftSize)
            rightAmp = parseList[count+4].split(' ')[0].strip('\t')
            rightSize = re.sub('([\(\[]).*?([\)\]])','\g<1>\g<2>', rightAmp)
            rightSize = rightSize.replace('[','N').replace(']','')
            rightSize = len(rightSize)
            tmpList = [size,forward,reverse,leftSize,rightSize]
            checkDict[sequence] = tmpList
            countDict[outFileName]= checkDict
            count += 6
        else:
            continue
#Worst part is here. Itterates over key value pairs and extracts the amplicons from the trimal trimmed msa files and stores the amplicons into a new msa fasta file as long as they meet the criteria. Such as not to many primerpair hits for a single primerpair
    for k,v in countDict.items():
        outFile = open('validationFilteredPrimerSearch/'+k.strip('\n')+'.fasta','a')
        orthoName = k
        orthoName = orthoName[:9]
        with open('../Final_Set_Pipeline_Run/trimalTrimmedFastaFiles/trimmed'+orthoName+'.fasta','r') as ff:
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
                if len(v2) == 5:
                    k2 = k2+' '+v2[0]+' '
                    isolate = msa[msa.index(k2)]
                    sequence = msa[msa.index(k2)+1]
                    sequence = sequence.replace('-','')
                    outFile.write('>'+isolate+'\n'+sequence[int(v2[1])+int(v2[3]):int(v2[0])-int(v2[2])-int(v2[4])]+'\n')
                    print(orthoName)
    
    fi += 1


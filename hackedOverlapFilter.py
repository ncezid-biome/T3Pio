import glob
import os
import pdb

def primerScore(primer):
    spot = 0
    position = 1
    for letter in primer:
        if letter not in ('ATGC'):
            spot = position
        position += 1
    score = spot/(len(primer))
    return score

rawPrimer3Files = sorted(glob.glob('../../../../primer3Files/*')) #grabs the primer3 generated files

if not os.path.exists('overlapFilteredPrimers/'): #creates directory to store the files generated from the script
    os.makedirs('overlapFilteredPrimers/',exist_ok=True)

badDict = {} #dictionary to hold the orthogroup and associated bad primers as decided from specificity testing
with open('offTargetPrimerHits','r') as b: #file that contains the bad primer pair information
    badPrimer = b.readlines() #b is for bad
b.close()
for line in badPrimer:
    name = line
    orth = name[:9] #important to key by orthogroup since the primer3 files are named from just the orthogroup
    pNum = name[9:].strip('\n') #the bad primer pair for an orthogroup
    badDict.setdefault(orth, []) #allows dictionary value to be an appendable list
    badDict[orth].append(pNum) #fills that list 

q=0

while q < len(rawPrimer3Files):
    OG = (rawPrimer3Files[q]).split('/')[-1]
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

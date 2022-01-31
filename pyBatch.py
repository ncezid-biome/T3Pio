import glob
import argparse
import subprocess
import os
import re

parser = argparse.ArgumentParser(description = 'Runs primersearch on input assembly list')

parser.add_argument('assembly',type=str,help='Assembly name')

#parser.add_argument('fileName',type=str,help='File name to find assembly')

args = parser.parse_args()

primerPairs = sorted(glob.glob('../overlapFilteredPrimers/*'))

if not os.path.exists('/scratch/yim0/primersearchFiles/'):
    os.makedirs('/scratch/yim0/primersearchFiles/',exist_ok=True)


#assemName = args.assembly.replace('-','_')
#fileName = args.fileName.split('.')[0]

assembly = args.assembly
assemblyName = assembly.split('/')[-1]
assemblyName = assemblyName.split('.')[0]

'''
assemName = ''.join(f for f in args.assembly if f.isalnum())
print(assemName)
for i in assemblies:
    a = i.split('/')
    b = a[len(a)-1]
    a = ''.join(f for f in b if f.isalnum())
    print(a)
    if a == assemName:
        assemfilepath=i
        assemfile = b
        break
    else:
        continue
#make searchFiles empty list before primersearch then append with outfile name. No need to search directory with all the damn files in it
'''

searchFiles = []
for i in primerPairs:
    OG = (i.split('/')[2])[:9]
    return_code = subprocess.check_output(['primersearch','-seqall',assembly,'-infile',i,'-mismatchpercent','6','-outfile','/scratch/yim0/primersearchFiles/'+OG+'--'+assemblyName])
    searchFiles.append('/scratch/yim0/primersearchFiles/'+OG+'--'+assemblyName)


#psAssemFile = assemfilepath+'/'+assemfile+'.fasta'

if not os.path.exists('filteredPrimerSearch/'):
    os.makedirs('filteredPrimerSearch/',exist_ok=True)

sf = 0
while sf < len(searchFiles):
    searchList = searchFiles[sf].split('--')
    if len(searchList) > 2:
        isolate = searchList[1]+'--'+searchList[2]
    else:
        isolate = searchList[1]
    with open(searchFiles[sf],'r') as f:
        parseList = f.readlines()
    f.close()
    with open(assembly,'r') as a:
        seq = a.readlines()
    a.close()
    print(searchFiles[sf])
    print(isolate)
    count = 0
    while count < len(parseList):
#Need way to account for primerSearch files with no information and multiple primers
#Probably add extra conditional and can make use of existing try/except block to handle out of bounds issues
#Not sure if anything else should be accounted for
#Plan is to see if line below 'Primer name' is new line character: adding and parseList[count+1] != '\n'
        try:
            if parseList[count].startswith('Primer name') and parseList[count+1] != '\n':
                try:
                    primerPair = (parseList[count].split(' ')[2]).strip('\n')
                    outfile = open('filteredPrimerSearch/'+ isolate + '.fasta','a')
                    contig = (parseList[count+2].split(' ')[1])
                    contig = '>'+contig
                    for line in seq:
                        if line.startswith(contig):
                            sequence = seq[seq.index(line)+1]
                            size =len(sequence)
                            forward = parseList[count+4]
                            forward = forward.split(' ')[5]
                            reverse = parseList[count+5]
                            reverse = (reverse.split(' ')[5]).replace('[','').replace(']','')
                            leftAmp = parseList[count+3].split(' ')[0].strip('\t')
                            leftSize = re.sub('([\(\[]).*?([\)\]])','\g<1>\g<2>', leftAmp)
                            leftSize = leftSize.replace('[','N').replace(']','')
                            leftSize = len(leftSize)
                            rightAmp = parseList[count+4].split(' ')[0].strip('\t')
                            rightSize = re.sub('([\(\[]).*?([\)\]])','\g<1>\g<2>', rightAmp)
                            rightSize = rightSize.replace('[','N').replace(']','')
                            rightSize = len(rightSize)
                            outfile.write(primerPair+'.'+'>'+isolate+ '\n' + sequence[int(forward)+leftSize:(size-int(reverse)-rightSize)] + '\n')
                            count += 7
                except IndexError:                
                            count += 1
            else:
                count += 1
        except IndexError:
            count += 1
    sf += 1

for i in searchFiles:
    return_code = subprocess.check_output(['rm',i])

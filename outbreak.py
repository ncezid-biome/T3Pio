import glob
import argparse
import subprocess
import os

parser = argparse.ArgumentParser(description = 'Runs primersearch on input assembly list')

parser.add_argument('assemblies',type=str,help='File containing the names of assemblies')

args = parser.parse_args()

qFolder = args.assemblies.split('.')[0]
with open(args.assemblies,'r') as a:
    assem = a.readlines()

assem = list(map(lambda s: s.strip('\n'), assem))
assemPath = sorted(glob.glob(qFolder+'/*'))
trueAssemblyName = []
for j in assemPath:
    j = j.split('/')[1]
    k = j.replace('-','_')
    if k in assem:
        trueAssemblyName.append(j)

primerPairs = sorted(glob.glob('overlapFilteredPrimers/*'))
if not os.path.exists('primersearchFiles/'):
    os.makedirs('primersearchFiles/',exist_ok=True)

for i in trueAssemblyName:
    pp = 0
    while pp <  len(primerPairs):
        OG = (primerPairs[pp].split('/')[1])[:9]
        return_code = subprocess.check_output(['primersearch','-seqall',qFolder+'/'+i+'/'+i+'.fasta','-infile',primerPairs[pp],'-mismatchpercent','6','-outfile','primersearchFiles/'+OG+'--'+i])
        pp += 1        


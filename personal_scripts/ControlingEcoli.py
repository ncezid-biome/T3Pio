import subprocess
import time
import argparse

parser = argparse.ArgumentParser(description='This script is meant to control the submission of a set of SRR and SRS numbers to NCBI for read downloads without making the whole damn world angry')

parser.add_argument('InputFile',type=str,help='The file')

parser.add_argument('qstatLen',type=int,help='The normal length of your qstat return when split on \n')

parser.add_argument('sleepLen',type=int,help='How long do you want your script to sleep between checking the queue')

args = parser.parse_args()

with open(args.InputFile,'r') as f:
    Accessions = f.readlines()
f.close()


while len(Accessions) >=20:
    queue = True
    submission = Accessions[0:20]
    del(Accessions[0:20])
    f = open('HPCSubFile','w')
    for i in submission:
        i = i.strip('\n')
        print(i,file=f)
    f.close()
    return_code = subprocess.check_output(['./2.0batch.sh','20','HPCSubFile'])
    while queue == True:
        time.sleep(args.sleepLen)
        r1 = subprocess.Popen(('qstat'),stdout=subprocess.PIPE)
        qLen = (r1.stdout.read().decode('ascii'))
        qLen = qLen.split('\n')
        if len(qLen) == args.qstatLen:
            queue = False
        else:
            continue

for i in Accessions:
    f = open('HPCSubFile','w')
    for i in submission:
        i = i.strip('\n')
        print(i,file=f)
    f.close()
    return_code = subprocess.check_output(['./2.0batch.sh',len(Accessions),'HPCSubFile'])

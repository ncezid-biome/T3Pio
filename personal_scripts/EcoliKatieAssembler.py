import argparse
import subprocess
import os
import time
import urllib.request

parser = argparse.ArgumentParser(description='This script is to be used in conjunction with bash script to mass gather and assemble reads from SRS and SRR numbers on the HPC')

parser.add_argument('accession',type=str,help='The accession number')

parser.add_argument('sleepTime',type=int,help='Pass a variable letting script know how long to sleep')

args = parser.parse_args()

accession = args.accession
accession = accession.strip('\n')

print(args.sleepTime)
a=True

while a == True:
    time.sleep(args.sleepTime)
    try:
        if accession[:3] == 'SRS':
            print(accession)
            r1 = subprocess.Popen(('esearch','-db','sra','-query',accession,'--api_key','ba1d246a5eff7d139247f581024ac7a3ca08'),stdout=subprocess.PIPE)
            r2 = subprocess.Popen(('efetch','--format','runinfo','--api_key','ba1d246a5eff7d139247f581024ac7a3ca08'),stdin=r1.stdout,stdout=subprocess.PIPE)
            r3 = subprocess.Popen(('cut','-d',',','-f','1'),stdin=r2.stdout,stdout=subprocess.PIPE)
            r4 = subprocess.Popen(('grep','SRR'),stdin=r3.stdout,stdout=subprocess.PIPE)
            SRR=(r4.stdout.read().decode('ascii'))
            SRR=SRR.strip('\n')
            os.makedirs('Ecoli_Reads/'+accession)
        
            return_code = subprocess.check_output(['prefetch','-v',SRR])
            return_code = subprocess.check_output(['fastq-dump','--accession',SRR,'--outdir','Ecoli_Reads/'+accession,'--defline-seq','@$ac.$si/$ri','--defline-qual','+','--split-files','--skip-technical','--dumpbase','--clip','--gzip'])

            fastqArgs = '-o Ecoli_Reads/'+accession+' --noextract'
        
            return_code = subprocess.check_output(['trim_galore','--fastqc_args',fastqArgs,'--clip_R1','10','--clip_R2','10','--three_prime_clip_R1','2','--three_prime_clip_R2','2','--paired','Ecoli_Reads/'+accession+'/'+SRR+'_1.fastq.gz','Ecoli_Reads/'+accession+'/'+SRR+'_2.fastq.gz','-o','Ecoli_Reads/'+accession])
            return_code = subprocess.check_output(['./skesaRunner.sh','Ecoli_Reads/'+accession+'/'+SRR+'_1_val_1.fq.gz','Ecoli_Reads/'+accession+'/'+SRR+'_2_val_2.fq.gz','Ecoli_Reads/'+accession+'/'+accession+'.fasta'])
            
            a = False
        else:
            os.makedirs('Ecoli_Reads/'+accession)
    
            return_code = subprocess.check_output(['prefetch','-v',accession])
            return_code = subprocess.check_output(['fastq-dump','--accession',accession,'--outdir','Ecoli_Reads/'+accession,'--defline-seq','@$ac.$si/$ri','--defline-qual','+','--split-files','--skip-technical','--dumpbase','--clip','--gzip'])


            fastqArgs = '-o Ecoli_Reads/'+accession+' --noextract'


            return_code = subprocess.check_output(['trim_galore','--fastqc_args',fastqArgs,'--clip_R1','10','--clip_R2','10','--three_prime_clip_R1','2','--three_prime_clip_R2','2','--paired','Ecoli_Reads/'+accession+'/'+accession+'_1.fastq.gz','Ecoli_Reads/'+accession+'/'+accession+'_2.fastq.gz','-o','Ecoli_Reads/'+accession])

            return_code = subprocess.check_output(['./skesaRunner.sh','Ecoli_Reads/'+accession+'/'+accession+'_1_val_1.fq.gz','Ecoli_Reads/'+accession+'/'+accession+'_2_val_2.fq.gz','Ecoli_Reads/'+accession+'/'+accession+'.fasta'])
            
            a = False
    except urllib.error.HTTPError:
        time.sleep(10)
        continue
















#return_code = subprocess.Popen(('skesa','--cores','','--fastq','Ecoli_Reads/'+accession+'_1_val_1.fq.gz','Ecoli_Reads/'+accession+'_2_val_2.fq.gz','--use_paired_ends'),stdout=subprocess.PIPE)
#return_assem = subprocess.Popen(('>','Ecoli_Reads/'+accession+'.fasta'),stdin=retrun_code.stdout)

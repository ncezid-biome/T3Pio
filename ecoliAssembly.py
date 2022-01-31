import glob
import subprocess
import os

rawReads = sorted(glob.glob('Raw_Reads/*'))

i = 0

while i < len(rawReads):
    isoName = rawReads[i].split('/')[1]
    isoName = isoName.split('_')[0]
    os.makedirs(isoName, exist_ok=True)
    print(rawReads[i]+' '+ rawReads[i+1])
    return_code = subprocess.check_output(['trim_galore','--fastqc','--fastqc_args','"--outdir trimedFastqc"','--clip_R1','15','--clip_R2','15','--three_prime_clip_R1','15','--three_prime_clip_R2','15','--paired',rawReads[i],rawReads[i+1]])


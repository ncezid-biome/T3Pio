#!/bin/bash 

source /etc/profile.d/modules.sh
module purge
module load Python/3.7
module load trim_galore
module load cutadapt
module load fastqc
module load Entrez
module load sratoolkit/2.9.1

FILE=~/CIMS/Katie_Ecoli/$1;
ASSEMBLY=$(sed "${SGE_TASK_ID}q;d" $FILE);
NUM=$SGE_TASK_ID

echo $ASSEMBLY
#echo $FILE
echo $NUM

python ~/CIMS/Katie_Ecoli/EcoliKatieAssembler.py $ASSEMBLY $NUM



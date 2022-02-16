#!/bin/bash 

source /etc/profile.d/modules.sh

FILE=~/CIMS/Katie_Ecoli/$1;
ASSEMBLY=$(sed "${SGE_TASK_ID}q;d" $FILE);

echo $ASSEMBLY

~/CIMS/Katie_Ecoli/localAssembleBashBatch.sh $ASSEMBLY

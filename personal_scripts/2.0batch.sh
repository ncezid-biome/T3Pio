#!/bin/bash 

source /etc/profile.d/modules.sh
module purge
module load Python/3.7

TMP=$(mktemp --tmpdir='.' --directory commensalEcoli.XXXXXXXX)
mkdir -p $TMP/log

qsub -q edlb.q -o $TMP/log -j y -cwd -pe smp 1 -V -N commensalEcoli -t 1:$1 ~/CIMS/Katie_Ecoli/batch.sh $2





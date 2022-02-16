#!/bin/bash

source /etc/profile.d/modules.sh
module load Skesa

skesa --cores 4 --fastq $1,$2 --use_paired_ends > $3

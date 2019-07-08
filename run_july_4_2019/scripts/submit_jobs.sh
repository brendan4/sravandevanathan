#!/bin/bash

scripts=`ls 12*sh`

for s in $scripts;
do
   sbatch $s
done

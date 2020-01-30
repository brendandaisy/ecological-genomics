#!/bin/bash

cd ~/ecological-genomics/mydata

date=$(date +%F)
dir=fastqc-${date}
mkdir ${dir} # makes dir if not already exists

for f in ${1}*PRK*
do
    fastqc ${f} -o ${dir}
done


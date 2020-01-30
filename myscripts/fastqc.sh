#!/bin/bash

cd ~/ecological-genomics

mkdir fastqc # makes dir if not already exists

for f in /data/project_data/RS_ExomeSeq/fastq/edge_fastq/PRK*
do
    fastqc ${f} -o fastqc/
done


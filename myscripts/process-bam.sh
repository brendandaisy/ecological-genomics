#!/bin/bash

# set the output folder and other globals
myrepo="/users/b/c/bcase/ecological-genomics"
outfolder="${myrepo}/myresults/BWA"

mkdir ${outfolder}

mypop="PRK"
ref="/data/project_data/RS_ExomeSeq/ReferenceGenomes/Pabies1.0-genome_reduced.fa"

# move to the folder containing trimmed reads
cd /data/project_data/RS_ExomeSeq/fastq/edge_fastq

for R1 in ${mypop}*R1_fastq.gz
do
    R2=${R1/_R1_fastq.gz/_R2_fastq.gz}
    f=${R1/_R1_fastq.gz/}
    name=`basename ${f}` # sets output file name corresp. to this individual
    bwa mem -t 1 -M -a ${ref} ${R1} ${R2} > ${outfolder}/${name}.sam
    sambamba-0.7.1-linux-static view -S --format=bam ${outfolder}/${name}.sam -o ${outfolder}/${name}.bam # convert to binary
    sambamba-0.7.1-linux-static markdup -r -t 1 ${outfolder}/${name}.bam ${outfolder}/${name}.rmdup.bam # remove PCR duplicates
    samtools sort ${outfolder}/${name}.rmdup.bam -o ${outfolder}/${name}.sorted.rmdup.bam
done

